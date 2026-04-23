/******************************************************************************
 * canonical_decomposition.h
 * *
 * Decomposes a hypergraph in its canonical decomposition. For this, we first
 * compute the prime decomposition of the hypergraph with the help of a split
 * oracle and we then re-assemble the prime decomposition into the canonical
 * decomposition, which can be used to generate the hypercactus.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_CANONICAL_DECOMPOSITION_H
#define SMHM_CANONICAL_DECOMPOSITION_H

#include <iostream>
#include <cassert>
#include <vector>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/orderer/orderer.h"
#include "lib/conversion/hypergraph_converter.h"
#include "lib/coarsening/strongly_connected_components.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

// The result of the split oracle query
struct SplitOracleResult
{
    bool hasFoundSplit;
    NodeID sourceID;
    NodeID sinkID;
    std::vector<bool> splitSide;
};

// The marker edge that connects two decomposition elements in the decomposition tree
// NB: isFromGoodSplit indicates whether the marker edge comes from a good split (i.e. is present in the canonical decomposition)
struct MarkerEdge
{
    NodeID markerID;
    NodeID representativeNodeID;
    ClusterIndex targetDecompositionIndex;
    EdgeIndex markerEdgeIndexAtTarget;
    bool isFromGoodSplit;
};

// The decomposition element of a hypergraph, which contains the included nodes and the marker edges
// NB: isDeleted indicates whether the decomposition element is deleted, i.e. it is not part of the canonical decomposition
// NB2: isSolidPolygon indicates whether the decomposition element is a solid polygon
// NB3: hasGraphCycle indicates whether the decomposition element has a graph cycle (only used for solid polygons)
struct DecompositionElement
{
    std::vector<NodeID> includedNodes;
    std::vector<MarkerEdge> markerEdges;
    bool isDeleted = false;
    bool isSolidPolygon = false;
    bool hasGraphCycle = false;
};

class CanonicalDecomposition
{
private:
    // Store the minimum edge cut of the original hypergraph
    CutValue minEdgeCut;
    // Store the number of nodes in the kernelized hypergraph
    NodeIndex kernelNumNodes;
    // Store the mapping of the kernelization, i.e. which node in the original hypergraph corresponds to which node in the kernelized hypergraph
    std::vector<ClusterID> kernelizationMapping;
    // Store the nodes of the kernelized hypergraph for which we should check if isolating them leads to a split (due to the previous kernelization)
    std::vector<NodeID> nodesToCheckForIsolatedSplit;
    // Computes the tight ordering of the hypergraph passed to the split oracle
    Orderer<StaticHypergraph, EdgeWeight> tightOrderer;
    // Store the tight ordering of the hypergraph passed to the split oracle
    std::vector<NodeID> nodeOrdering;
    // Store the representive node id for each node in the kernelized hypergraph (only used for computing the prime decomposition of the kernelized hypergraph)
    std::vector<NodeID> representativeNodeID;
    // Store the tree of decomposition elements computed during the recursive calls
    std::vector<DecompositionElement> decompositionTree;
    // Store the number of decomposition elements computed so far
    ClusterIndex numDecompositionElements = 0;
    // Store the number of markers
    EdgeIndex numMarkers = 0;
    // Store the initial number of markers (used if we want to use the IDs of the markers as indices)
    EdgeIndex initialNumMarkers = 0;
    // Whether we should print verbose output
    bool verbose = false;

    // Query the split oracle to either find a split or return a pair {s, t} of nodes for which there does not exist a (s-t)-split
    SplitOracleResult query_split_oracle(StaticHypergraph &hypergraph);

    // Compute max flow on the tight graph between the last two nodes in the MA-ordering (Arikati & Melhorn)
    // NB: The trick is that we can directly use the tight ordering, as MA-ordering and tight ordering are equivalent for graphs
    //     and the tight ordering of the hypergraph is identical to the tight ordering of the tight graph
    // NB2: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
    bool compute_max_flow_on_tight_graph(StaticGraph &tightGraph, const NodeID sourceID, const NodeID sinkID, std::vector<EdgeWeight> &flow, std::vector<EdgeID> &directedFlowEdge);

    // Perform BFS to find the nodes that are contracted into the marker node in the decomposition tree
    // We then assign the same cluster ID to all nodes that are contracted into the marker node and return this cluster ID
    // NB: The marker edge is directed, i.e. we find the nodes of one side of the split that are contracted into the marker node
    NodeID find_nodes_contracted_to_marker_node_in_decomposition_tree(mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                      const ClusterIndex sourceDecompositionIndex,
                                                                      const MarkerEdge markerEdge);

    // Check if the hypergraph is a solid polygon
    // A solid polygon is a hypergraph that consists of a GRAPH cycle (including all nodes) where each edge has the same weight A and a hyperedge containing all nodes with weight B.
    // NB: If A = 0 (i.e. there is no graph cycle) OR if B = 0 (i.e. there is no hyperedge containing all nodes), then the hypergraph is still a solid polygon.
    bool is_solid_polygon(StaticHypergraph &hypergraph) const;

    // Return for a marker edge the new node ID of the marker node in the hypercactus
    NodeID get_new_node_id_of_marker_in_hypercactus(MarkerEdge &markerEdge, std::vector<bool> &hasSeenMarkerNode, NodeIndex &numHypercactusNodes);

public:
    CanonicalDecomposition(const NodeIndex kernelNumNodes,
                           const EdgeIndex kernelNumEdges,
                           std::vector<ClusterID> kernelizationMapping,
                           const CutValue minEdgeCut,
                           const bool verbose = false);

    // Recursively compute the prime decomposition of the kernelized hypergraph
    void compute_prime_decomposition_of_kernel(StaticHypergraph hypergraph, NodeIndex decompositionIndex = 0);

    // Convert the prime decomposition of the kernelized hypergraph into the prime decomposition of the original hypergraph
    void convert_to_prime_decomposition_of_original_hypergraph();

    // Compute the canonical decomposition from the prime decomposition by trying to merge the decompositions elements by their markers
    void compute_canonical_decomposition_from_prime_decomposition(StaticHypergraph &hypergraph);

    // Create the hypercactus from the canonical decomposition
    // NB: The weight of a node in the hypercactus indicates how many nodes from the original hypergraph it represents
    mt_kahypar_hypergraph_t create_hypercactus(StaticHypergraph &hypergraph, std::vector<ClusterID> &hypercactusMapping, bool allowOutputFloatEdgeWeights = false);

    // Print the decomposition tree
    void print_decomposition_tree(const bool isCanonical = false) const;
};

// Query the split oracle to either find a split or return a pair {s, t} of nodes for which there does not exist a (s-t)-split
inline SplitOracleResult CanonicalDecomposition::query_split_oracle(StaticHypergraph &hypergraph)
{
    // Get the number of nodes in the hypergraph
    const NodeID hypergraphNumNodes = hypergraph.initialNumNodes();
    // Get the number of edges in the hypergraph
    const EdgeID hypergraphNumEdges = hypergraph.initialNumEdges();

    // Compute the tight ordering of the nodes in the hypergraph
    tightOrderer.compute_ordering(hypergraph, hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes(), &nodeOrdering, nullptr, nullptr);

    // Create the (static) tight graph
    // NB1: The tight graph can have parallel edges, since the same pair of nodes can be the last two pins in the node ordering for different hyperedges
    // NB2: The uniqueID of each edge in the tight graph corresponds to ID of the edge in the hypergraph
    // NB3: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
    // NB4: When computing the max flow on the tight graph, we need that the incident edges of each nodes are sorted, which is achieved by setting the stable_construction_of_incident_edges flag to true
    mt_kahypar_hypergraph_t tightGraphWrapper = HypergraphConverter::hypergraph_to_tight_graph<StaticHypergraph>(hypergraph, nodeOrdering, DETERMINISTIC, true);
    // Make sure that the created graph is a static graph
    assert(tightGraphWrapper.type == STATIC_GRAPH);
    // Cast the static graph wrapper to a static graph
    StaticGraph &tightGraph = mt_kahypar::utils::cast<StaticGraph>(tightGraphWrapper);

    // Print the structure of the tight graph
    if (verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "################### TIGHT GRAPH ###################" << std::endl;
        std::cout << "===================================================" << std::endl;
        for (NodeID nodeID : nodeOrdering)
            std::cout << nodeID + 1 << " ";
        std::cout << std::endl;
        std::cout << "tight_graph_num_nodes \t\t\t" << tightGraph.initialNumNodes() << std::endl;
        std::cout << "tight_graph_num_edges \t\t\t" << tightGraph.initialNumEdges() / 2 << std::endl;
        for (EdgeID edgeID : tightGraph.edges())
            if (tightGraph.edgeSource(edgeID) < tightGraph.edgeTarget(edgeID))
                std::cout << tightGraph.edgeSource(edgeID) + 1 << " <-> " << tightGraph.edgeTarget(edgeID) + 1 << " (weight: " << tightGraph.edgeWeight(edgeID) << " " << "unique_id: " << tightGraph.uniqueEdgeID(edgeID) << ")" << std::endl;
    }

    // Compute max flow on the tight graph between the last two nodes in the MA-ordering (Arikati & Melhorn)
    // NB: The trick is that we can directly use the tight ordering, as MA-ordering and tight ordering are equivalent for graphs
    //     and the tight ordering of the hypergraph is identical to the tight ordering of the tight graph
    // NB2: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
    const NodeID sourceID = hypergraphNumNodes - 2;
    const NodeID sinkID = hypergraphNumNodes - 1;
    std::vector<EdgeWeight> flow;
    // Store the mapping from the unique edge ID to the directed edge ID in the tight graph
    // NB: We only use this if the flow is non-zero, since the flow grows
    std::vector<EdgeID> directedFlowEdge;
    bool canFindSplit = compute_max_flow_on_tight_graph(tightGraph, sourceID, sinkID, flow, directedFlowEdge);
    // Check if we can find a split
    if (!canFindSplit)
    {
        mt_kahypar_free_hypergraph(tightGraphWrapper);
        return {false, nodeOrdering[sourceID], nodeOrdering[sinkID], {}};
    }

    // Output the results of the max flow computation
    if (verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "##################### MAX FLOW ####################" << std::endl;
        std::cout << "===================================================" << std::endl;
        for (EdgeID edgeID : tightGraph.edges())
        {
            EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
            EdgeID directedEdgeID = directedFlowEdge[uniqueEdgeID];
            if (flow[uniqueEdgeID] > 0 && edgeID == directedEdgeID)
                std::cout << tightGraph.edgeSource(directedEdgeID) + 1 << " -> " << tightGraph.edgeTarget(directedEdgeID) + 1 << " (flow: " << flow[uniqueEdgeID] << ")" << std::endl;
        }
    }

    // Create the (static) digraph
    // NB: The digraph representation is usually directed, but the graph datastructure only supports undirected graphs (i.e. forward and backward edges)
    //     Therefore, on must make sure to not use any of the backward edges, which can be done by checking the IDs of the nodes in the digraph representation.
    // NB2: The IDs of the nodes and edges are mapped in a specific way from the hypergraph to the digraph (see HypergraphConverter::hypergraph_to_digraph for details)
    mt_kahypar_hypergraph_t digraphWrapper = HypergraphConverter::hypergraph_to_digraph<StaticHypergraph>(hypergraph, DETERMINISTIC);
    // Make sure that the created graph is a static graph
    assert(digraphWrapper.type == STATIC_GRAPH);
    // Cast the static graph wrapper to a static graph
    StaticGraph &digraph = mt_kahypar::utils::cast<StaticGraph>(digraphWrapper);

    // Compute the residual digraph
    // NB: We use the specific mapping of the IDs as described in HypergraphConverter::hypergraph_to_digraph
    for (EdgeID edgeID : digraph.edges())
    {
        NodeID sourceID = digraph.edgeSource(edgeID);
        NodeID targetID = digraph.edgeTarget(edgeID);
        // Check if we have an edge of the form (v -> e+)
        if ((sourceID < hypergraphNumNodes && targetID >= hypergraphNumNodes + hypergraphNumEdges))
        {
            EdgeID uniqueEdgeIDInTightGraph = targetID - hypergraphNumNodes - hypergraphNumEdges;
            EdgeID directedFlowEdgeIDInTightGraph = directedFlowEdge[uniqueEdgeIDInTightGraph];
            // Check if v receives some flow in the tight graph through the unique edge
            // If so, set the edge weight of the backward edge (v -> e+) to the received flow, otherwise disable it by setting its weight to zero
            EdgeWeight edgeWeight = (flow[uniqueEdgeIDInTightGraph] > 0 && nodeOrdering[tightGraph.edgeTarget(directedFlowEdgeIDInTightGraph)] == sourceID) ? flow[uniqueEdgeIDInTightGraph] : 0;
            digraph.setEdgeWeight(edgeID, edgeWeight);
        }
        // Check if we have an edge of the form (e- -> v)
        else if (sourceID >= hypergraphNumNodes && sourceID < hypergraphNumNodes + hypergraphNumEdges && targetID < sourceID)
        {
            EdgeID uniqueEdgeIDInTightGraph = sourceID - hypergraphNumNodes;
            EdgeID directedFlowEdgeIDInTightGraph = directedFlowEdge[uniqueEdgeIDInTightGraph];
            // Check if v gives some flow in the tight graph through the unique edge
            // If so, set the edge weight of the backward edge (e- -> v) to the given flow, otherwise disable it by setting its weight to zero
            EdgeWeight edgeWeight = (flow[uniqueEdgeIDInTightGraph] > 0 && nodeOrdering[tightGraph.edgeSource(directedFlowEdgeIDInTightGraph)] == targetID) ? flow[uniqueEdgeIDInTightGraph] : 0;
            digraph.setEdgeWeight(edgeID, edgeWeight);
        }
        // Check if we have an edge of the form (e+ -> e-)
        else if (targetID >= hypergraphNumNodes && targetID < hypergraphNumNodes + hypergraphNumEdges && targetID < sourceID)
        {
            EdgeID uniqueEdgeIDInTightGraph = targetID - hypergraphNumNodes;
            // Set the edge weight of the edge (e+ -> e-) to the flow of the unique edge in the tight graph
            digraph.setEdgeWeight(edgeID, flow[uniqueEdgeIDInTightGraph]);
        }
        // Check if we have an edge of the form (e- -> e+)
        else if (sourceID >= hypergraphNumNodes && sourceID < hypergraphNumNodes + hypergraphNumEdges && sourceID < targetID)
        {
            EdgeID uniqueEdgeIDInTightGraph = sourceID - hypergraphNumNodes;
            // Decrease the edge weight of the edge (e- -> e+) by the flow of the unique edge in the tight graph
            digraph.setEdgeWeight(edgeID, digraph.edgeWeight(edgeID) - flow[uniqueEdgeIDInTightGraph]);
        }
    }

    // Print the structure of the residual digraph
    if (verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "################# RESIDUAL DIGRAPH ################" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "residual_digraph_num_nodes \t\t" << digraph.initialNumNodes() << std::endl;
        EdgeIndex residualDigraphNumEdges = 0;
        for (EdgeID edgeID : digraph.edges())
            if (digraph.edgeWeight(edgeID) > 0)
                residualDigraphNumEdges++;
        std::cout << "residual_digraph_num_edges \t\t" << residualDigraphNumEdges << std::endl;
        for (EdgeID edgeID : digraph.edges())
            if (digraph.edgeWeight(edgeID) > 0)
                std::cout << digraph.edgeSource(edgeID) + 1 << " -> " << digraph.edgeTarget(edgeID) + 1 << " (weight: " << digraph.edgeWeight(edgeID) << " " << "unique_id: " << digraph.uniqueEdgeID(edgeID) << ")" << std::endl;
    }

    // Store the SCC ID for each node in the graph
    // NB: The SCC IDs also correspond to the finishing times of the SCCs in the DFS
    std::vector<ClusterID> componentID;
    // Computes the strongly connected components (SCCs) of a directed graph
    // NB: We consider only edges with a positive weight, i.e. we ignore edges with a zero or negative weight
    // NB2: The component IDs are stored in the passed vector and the number of components is returned
    StronglyConnectedComponentsFinder sccFinder(digraph.initialNumNodes(), digraph.initialNumEdges());
    ClusterIndex componentCounter = sccFinder.find_sccs(digraph, componentID);

    // Since the sink of the flow is the last node in the tight ordering of the orignal hypergraph,
    // we observe that the SCC of the sink has no incoming edges in the residual digraph.
    // Therefore, we can assign to the SCC of the sink the highest ID (i.e. the last finishing time in the DFS).
    // All SCCs with a higher ID than the SCC of the sink will decrement their ID by one.
    // NB: We use the nodeOrdering mapping because sourceID and sinkID are the IDs in the tight graph, which is equal to their position in the tight ordering
    // NB2: We also directly count here how many nodes (v) of the original hypergraph are in the SCC of the source and in the SCC of the source.
    NodeIndex numOriginalNodesInSourceComponent = 0;
    NodeIndex numOriginalNodesInSinkComponent = 0;
    ClusterID oldSourceComponentID = componentID[nodeOrdering[sourceID]];
    ClusterID oldSinkComponentID = componentID[nodeOrdering[sinkID]];

    for (NodeID nodeID : digraph.nodes())
        if (componentID[nodeID] == oldSinkComponentID)
        {
            if (nodeID < hypergraphNumNodes)
                numOriginalNodesInSinkComponent++;
            // Assign the highest ID to the SCC of the sink
            componentID[nodeID] = componentCounter - 1;
        }
        else
        {
            if (nodeID < hypergraphNumNodes && componentID[nodeID] == oldSourceComponentID)
                numOriginalNodesInSourceComponent++;
            // Decrease the ID of the SCCs that were higher than the SCC of the sink
            if (componentID[nodeID] > oldSinkComponentID)
                componentID[nodeID]--;
        }

    // Since the sink is the last node in the tight ordering, {sinkID} is a trivial mincut between the source and the sink, i.e. the sink cannot reach any other node of the original hypergraph.
    // This means that the SCC of the sink contains only one the sink itself (as nodes in the original hypergraph).
    assert(numOriginalNodesInSinkComponent == 1);
    // Make sure that the SCC of the sink has the highest ID
    assert(componentID[nodeOrdering[sinkID]] == componentCounter - 1);

    // Print the SCC assignment of the residual digraph
    if (verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "####### SCC ASSIGNMENT OF RESIDUAL DIGRAPH ########" << std::endl;
        std::cout << "===================================================" << std::endl;
        for (NodeID nodeID : digraph.nodes())
            std::cout << componentID[nodeID] + 1 << " ";
        std::cout << std::endl;
        std::cout << "num_sccs \t\t\t\t" << componentCounter << std::endl;
        std::cout << "num_original_nodes_in_source_scc \t" << numOriginalNodesInSourceComponent << std::endl;
        std::cout << "id_of_source_scc \t\t\t" << componentID[nodeOrdering[sourceID]] + 1 << std::endl;
        std::cout << "num_original_nodes_in_sink_scc \t\t" << numOriginalNodesInSinkComponent << std::endl;
        std::cout << "id_of_sink_scc \t\t\t\t" << componentID[nodeOrdering[sinkID]] + 1 << std::endl;
    }

    // Now we try to find a non-trivial mincut between the source and the sink in the residual digraph
    std::vector<bool> splitSide(hypergraphNumNodes, 0);
    // NB: We do not explicitly contract the SCCs in the digraph, but we implicitly use the SCCs as the elements in the enumeration algorithm of Schrage & Baker
    // NB2: We must flip the idea of the enumeration algorithm, so here a set is feasible if for every element in the set, all SUCCESSORS of the element are also in the set (i.e. we have a closed node set)
    //      Therefore, we order the SCCs according to their finishing times (= IDs of SCCs) in the DFS (i.e. all successors of a SCCs come before the SCC itself)
    // NB3: Every SCC contains at least one node (v) of the original hypergraph.
    //      We observe that the SCC cannot be empty, so if it does not contain v, it must contain either at least one e- or at least one e+ node.
    //      Case 1: If there exists an SCC including any e-, then the edge (e-, e+) either pushes non-zero flow or not. If it does, this means e- receives non-zero flow from a node v, i.e. both
    //              edges (v, e-) and (e-, v) exist in the residual digraph, which means that the SCC contains v. If it does not, then the edge (e-, e+) exists in the residual digraph, which means
    //              we have a cycle back to e- over e+ and some v, i.e. the SCC contains v.
    //      Case 2: If there exists an SCC including any e+, then the edge (e-, e+) either pushes non-zero flow or not. If it does, this means e+ givess non-zero flow to a node v, i.e. both
    //              edges (e+, v) and (v, e-) exist in the residual digraph, which means that the SCC contains v. If it does not, then the edge (e-, e+) exists in the residual digraph, which means
    //              we have a cycle back to e+ over some v and e-, i.e. the SCC contains v.
    // NB4: As a reminder, the SCC of the sink contains only the sink itself (as nodes in the original hypergraph).
    // First we make quick checks to find non-trivial mincuts. If one of these quick checks succeeds, we do not need to run the enumeration algorithm of Schrage & Baker
    // We keep track where we need to cut in the order to get the non-trivial mincut of a successful quick check. Every component ID greater or equal to quickNonTrivialCut will be assigned to the right split side (= 1).
    // This means that quickNonTrivialCut = 0 means that we do not have a found a non-trivial mincut with the quick checks.
    ClusterIndex quickNonTrivialCut = 0;
    // Quick check 1: If the SCC of the source contains at least two nodes (v) of the original hypergraph AND there lies a SCC between the SCC of the source and the SCC of the sink
    //                in the order, then we have a non-trivial mincut (cut after SCC of source)
    // Quick check 2: If there lies a SCC before the SCC of the source in the order AND there lies a SCC between the SCC of the source and the SCC of the sink in the order, then we have a
    //                non-trivial mincut (cut after SCC of source)
    if ((numOriginalNodesInSourceComponent >= 2 || componentID[nodeOrdering[sourceID]] > 0) &&
        componentID[nodeOrdering[sourceID]] < componentID[nodeOrdering[sinkID]] - 1)
        quickNonTrivialCut = componentID[nodeOrdering[sourceID]] + 1;
    // Quick check 3: If there are at least two SCCs between the SCC of the source and the SCC of the sink in the order, then we have a non-trivial mincut (cut after SCC that follows the SCC of the source)
    // NB: We need to make sure that there are at least 3 SCCs in total, otherwise we might get an overflow when subtracting 2 from the component ID of the sink
    else if (componentCounter >= 3 && componentID[nodeOrdering[sourceID]] < componentID[nodeOrdering[sinkID]] - 2)
        quickNonTrivialCut = componentID[nodeOrdering[sourceID]] + 2;

    // Check if we have found a non-trivial mincut with the quick checks
    if (quickNonTrivialCut > 0)
    {
        for (NodeIndex i = 0; i < hypergraphNumNodes; ++i)
            splitSide[i] = (componentID[i] >= quickNonTrivialCut);
        if (verbose)
            std::cout << "Found split: Quick check was successful" << std::endl;
        mt_kahypar_free_hypergraph(tightGraphWrapper);
        mt_kahypar_free_hypergraph(digraphWrapper);
        return {true, nodeOrdering[sourceID], nodeOrdering[sinkID], splitSide};
    }
    // If we have less than 4 SCCs, then there exist the following possible forms: [SCC(source), SCC(sink)] or [SCC(source), SCC(other), SCC(sink)] or [SCC(other), SCC(source), SCC(sink)]
    // In the first case, we know that SCC(sink) contains only the sink node (invariant). Thus there is no non-trivial mincut separating the source and the sink.
    // In the second case, we know that both SCC(source) and SCC(sink) contain only the source and sink node, respectively (invariant + quick check 1). Thus we cannot find a non-trivial mincut separating the source and the sink.
    // In the third case, we know that SCC(sink) contains only the sink node (invariant). Thus we cannot find a non-trivial mincut separating the source and the sink.
    // Therefore, we can directly stop if there are less than 4 SCCs
    if (componentCounter < 4)
    {
        if (verbose)
            std::cout << "No split: There are less than 4 SCCs and all quick checks failed." << std::endl;
        mt_kahypar_free_hypergraph(tightGraphWrapper);
        mt_kahypar_free_hypergraph(digraphWrapper);
        return {false, nodeOrdering[sourceID], nodeOrdering[sinkID], {}};
    }

    // Enumerate the minimum cuts between the source and the sink using the dynamic programming algorithm of Schrage & Baker
    // NB: Since we have at least 4 SCCs, we know that the SCC of the source comes directly before the SCC of the sink in the order (quick checks 2 & 3)
    assert(componentID[nodeOrdering[sourceID]] == componentID[nodeOrdering[sinkID]] - 1);
    // NB2: We are only interested in the mincut between the source and the sink, and we know that the SCC of the source comes directly before the SCC of the sink and the SCC of the sink
    //     is the last SCC in the order. The first mincut between the source and the sink that the algorithm considers is when every SCC is included in the set from left to tight in the order
    //     up until the SCC of the source included (i.e. only exclude SCC of sink). The algorithm then checks if we can remove any of the previous SCCs from the set, which is the case if the
    //     SCC has no PREDECESSOR in the set. If we can find such a SCC, we have found a non-trivial mincut between the source and the sink. The order in which we scan is not important.
    std::vector<bool> hasPredecessorInSet(componentCounter, false);
    for (NodeID nodeID : digraph.nodes())
        for (EdgeID edgeID : digraph.incidentEdges(nodeID))
        {
            NodeID targetID = digraph.edgeTarget(edgeID);
            // Check if the target belongs to a different SCC than the node
            // NB: If the node comes from the same SCC as the sink, then skip (since the SCC of the sink is not in the set)
            if (digraph.edgeWeight(edgeID) > 0 && componentID[nodeID] != componentID[targetID] && componentID[nodeID] != componentID[nodeOrdering[sinkID]])
                hasPredecessorInSet[componentID[targetID]] = true;
        }

    ClusterID componentIDWithNoPredecessorInSet = std::numeric_limits<ClusterID>::max();
    // Check if we can find a SCC with no predecessor in the set (excluding the SCCs of the source and of the sink)
    for (ClusterIndex i = 0; i < componentCounter - 2; ++i)
        if (!hasPredecessorInSet[i])
        {
            componentIDWithNoPredecessorInSet = i;
            break;
        }

    // If we have found such a component, we have found a non-trivial mincut between the source and the sink
    if (componentIDWithNoPredecessorInSet == std::numeric_limits<ClusterID>::max())
    {
        if (verbose)
            std::cout << "No split: Cannot find componennt with no predecessor in set." << std::endl;
        mt_kahypar_free_hypergraph(tightGraphWrapper);
        mt_kahypar_free_hypergraph(digraphWrapper);
        return {false, nodeOrdering[sourceID], nodeOrdering[sinkID], {}};
    }

    // Assign the split side based on the found component ID with no predecessor in the set
    for (NodeIndex i = 0; i < hypergraphNumNodes; ++i)
        splitSide[i] = (componentID[i] == componentIDWithNoPredecessorInSet || componentID[i] == componentID[nodeOrdering[sinkID]]);

    mt_kahypar_free_hypergraph(tightGraphWrapper);
    mt_kahypar_free_hypergraph(digraphWrapper);
    return {true, nodeOrdering[sourceID], nodeOrdering[sinkID], splitSide};
}

// Compute max flow on the tight graph between the last two nodes in the MA-ordering (Arikati & Melhorn)
// NB: The trick is that we can directly use the tight ordering, as MA-ordering and tight ordering are equivalent for graphs
//     and the tight ordering of the hypergraph is identical to the tight ordering of the tight graph
// NB2: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
inline bool CanonicalDecomposition::compute_max_flow_on_tight_graph(StaticGraph &tightGraph, const NodeID sourceID, const NodeID sinkID, std::vector<EdgeWeight> &flow, std::vector<EdgeID> &directedFlowEdge)
{
    // Store whether an edge contains flow from the source node (true) or flow to the sink node (false)
    // NB: We only use this flag if the flow is non-zero, since the flow grows from both the source and the sink node and we need to distinguish between the two directions
    std::vector<bool> containsFlowFromSource(tightGraph.maxUniqueID(), false);
    // Store the max flow value and the remaining flow to push over the edges
    FlowValue maxFlowValue = 0;
    FlowValue remainingFlowToPush = 0;
    // Store the flow for each unique edge in the tight graph
    flow.resize(tightGraph.maxUniqueID());
    std::fill(flow.begin(), flow.end(), 0);
    // Store the mapping from the unique edge ID to the directed edge ID in the tight graph
    // NB: We only use this if the flow is non-zero, since the flow grows
    directedFlowEdge.resize(tightGraph.maxUniqueID());

    // Saturate all incident edges to the sink node and compute the remaining flow to push over the edges
    for (EdgeID edgeID : tightGraph.incidentEdges(sinkID))
    {
        NodeID targetID = tightGraph.edgeTarget(edgeID);
        EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
        // Saturate all incident edges to the sink
        flow[uniqueEdgeID] = tightGraph.edgeWeight(edgeID);
        containsFlowFromSource[uniqueEdgeID] = false;
        // Increase the max flow value
        maxFlowValue += flow[uniqueEdgeID];
        // Only increase the remaining flow to push if the edges does not connect the sink with the source
        if (targetID != sourceID)
            remainingFlowToPush += flow[uniqueEdgeID];
    }

    // If the max flow value between the source and the sink is larger than the mincut value, there is no split between the source and the sink
    if (maxFlowValue > minEdgeCut)
    {
        if (verbose)
            std::cout << "No split: The computed max flow value (" << maxFlowValue << ") is larger than the mincut value (" << minEdgeCut << ")." << std::endl;
        return false;
    }

    if (remainingFlowToPush > 0)
    {
        // Push the entire remaining flow from the source node over the edges by prioritizing the edges with the highest target position in the tight ordering (excluding edges that connect the source with the sink)
        // NB: By construction of the tight graph, the incident edges of each node are already sorted in ascending order by the target ID (= position in the tight ordering), so we use the reversed order
        FlowValue flowToPushFromSource = remainingFlowToPush;
        auto edgeIterator = tightGraph.incidentEdges(sourceID).end();
        while (edgeIterator != tightGraph.incidentEdges(sourceID).begin() && flowToPushFromSource > 0)
        {
            // NB: We need to first decrement the iterator since .end() points at the element AFTER the vector
            EdgeID edgeID = *(--edgeIterator);
            EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
            NodeID targetID = tightGraph.edgeTarget(edgeID);

            // Ignore edges that connect the source with the sink
            if (targetID == sinkID)
                continue;

            // Compute the flow to push over the edge
            FlowValue flowToPush = std::min((FlowValue)tightGraph.edgeWeight(edgeID), flowToPushFromSource);
            // Push the flow over the edge
            flow[uniqueEdgeID] += flowToPush;
            containsFlowFromSource[uniqueEdgeID] = true;
            // Decrease the remaining flow to push from the source
            flowToPushFromSource -= flowToPush;
        }

        // Make sure that we pushed all the remaining flow from the source node
        assert(flowToPushFromSource == 0);
        // Make sure that we have at least 3 nodes if we have still remaining flow to push
        assert(tightGraph.initialNumNodes() > 2);

        // Route the remaining flow through the graph from the source to the sink
        NodeIndex currentNodeID = tightGraph.initialNumNodes() - 3;
        while (true)
        {
            // Store how much flow is already coming from the source or going to the sink
            // NB: For this, we consider only the incident edges with a target that has a HIGHER position in the tight ordering than the current node
            FlowValue flowToSink = 0;
            FlowValue flowFromSource = 0;

            for (EdgeID edgeID : tightGraph.incidentEdges(currentNodeID))
            {
                EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
                // Get the ID of the target (= position of the target in the tight ordering)
                NodeID targetID = tightGraph.edgeTarget(edgeID);

                // Ignore edges where the target position is LOWER than current node position
                if (targetID < currentNodeID)
                    continue;

                // Check if the edge contains non-zero flow
                if (flow[uniqueEdgeID] > 0)
                {
                    // If so, increase the respective flow counter
                    if (containsFlowFromSource[uniqueEdgeID])
                        flowFromSource += flow[uniqueEdgeID];
                    else
                        flowToSink += flow[uniqueEdgeID];
                }
            }

            // Update the remaining flow to push
            remainingFlowToPush -= std::min(flowToSink, flowFromSource);
            // Stop if there is no remaining flow to push
            if (remainingFlowToPush == 0)
                break;

            // Check whether the current node gets more flow from the source than it sends to the sink
            bool hasMoreFlowFromSourceThanToSink = (flowFromSource > flowToSink);
            // Compute the flow to push over the edges with a target that has a LOWER position in the tight ordering than the current node
            FlowValue flowToPushToNeighbours = hasMoreFlowFromSourceThanToSink ? (flowFromSource - flowToSink) : (flowToSink - flowFromSource);

            // Push the flow to the neighbours of the current node by prioritizing the edges with the highest target position in the tight ordering (but lower than the current node position)
            // NB: By construction of the tight graph, the incident edges of each node are already sorted in ascending order by the target ID (= position in the tight ordering), so we use the reversed order
            auto edgeIterator = tightGraph.incidentEdges(currentNodeID).end();
            while (edgeIterator != tightGraph.incidentEdges(currentNodeID).begin() && flowToPushToNeighbours > 0)
            {
                // NB: We need to first decrement the iterator since .end() points at the element AFTER the vector
                EdgeID edgeID = *(--edgeIterator);
                EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
                // Get the ID of the target (= position of the target in the tight ordering)
                NodeID targetID = tightGraph.edgeTarget(edgeID);

                // Ignore edges where the target position is HIGHER than current node position
                if (targetID > currentNodeID)
                    continue;

                // Compute the flow to push over the edge
                FlowValue flowToPush = std::min((FlowValue)tightGraph.edgeWeight(edgeID), flowToPushToNeighbours);
                // Push the flow over the edge
                flow[uniqueEdgeID] += flowToPush;
                containsFlowFromSource[uniqueEdgeID] = hasMoreFlowFromSourceThanToSink;
                // Decrease the remaining flow to push
                flowToPushToNeighbours -= flowToPush;
            }

            // Make sure that we pushed all the remaining flow to the neighbours of the current node
            assert(flowToPushToNeighbours == 0);

            // Stop if we have handled all nodes in the tight ordering
            if (currentNodeID == 0)
                break;
            else
                // Otherwise, decrease the node id (= position in the tight ordering)
                currentNodeID--;
        }
    }

    // Store the mapping from the unique edge ID to the directed edge ID in the tight graph
    // NB: We only use this if the flow is non-zero, since the flow grows
    for (EdgeID edgeID : tightGraph.edges())
    {
        // Get the ID of the source and the target (= position of the sourc and target in the tight ordering)
        NodeID sourceID = tightGraph.edgeSource(edgeID);
        NodeID targetID = tightGraph.edgeTarget(edgeID);
        EdgeID uniqueEdgeID = tightGraph.uniqueEdgeID(edgeID);
        // If the flow comes from the source, it flows downwards in the tight ordering, otherwise it flows upwards
        if (flow[uniqueEdgeID] > 0 &&
            ((containsFlowFromSource[uniqueEdgeID] && sourceID > targetID) ||
             (!containsFlowFromSource[uniqueEdgeID] && sourceID < targetID)))
            directedFlowEdge[uniqueEdgeID] = edgeID;
    }

    return true;
}

#endif // end of SMHM_CANONICAL_DECOMPOSITION_H