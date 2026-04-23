/******************************************************************************
 * hypergraph_converter.h
 * *
 * Converts a hypergraph into different formats (e.g. static, dynamic) or into
 * graph representations (e.g. tight graph).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_HYPERGRAPH_CONVERTER_H
#define SMHM_HYPERGRAPH_CONVERTER_H

// Own headers
#include "lib/utils/definitions.h"
#include "lib/io/mt_kahypar_io.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"

class HypergraphConverter
{
private:
public:
    // Converts a static hypergraph into a dynamic hypergraph and return whether the hypergraph has weighted edges
    static mt_kahypar_hypergraph_t static_hypergraph_to_dynamic_hypergraph(StaticHypergraph &hypergraph, bool &hasWeightedEdges);

    // Converts a (dynamic or static) hypergraph into a (dynamic or static) tight graph representation
    // NB: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
    template <typename Hypergraph>
    static mt_kahypar_hypergraph_t hypergraph_to_tight_graph(Hypergraph &hypergraph,
                                                             std::vector<NodeID> &nodeOrdering,
                                                             const mt_kahypar_preset_type_t tightGraphPresetType = DETERMINISTIC,
                                                             const bool stable_construction_of_incident_edges = false);

    // Converts a (dynamic or static) hypergraph into its (dynamic or static) digraph representation
    // NB: The digraph representation is usually directed, but the graph datastructure only supports undirected graphs (i.e. forward and backward edges)
    //     Therefore, on must make sure to not use any of the backward edges, which can be done by checking the IDs of the nodes in the digraph representation.
    template <typename Hypergraph>
    static mt_kahypar_hypergraph_t hypergraph_to_digraph(Hypergraph &hypergraph, const mt_kahypar_preset_type_t digraphPresetType = DETERMINISTIC);
};

// Converts a (dynamic or static) hypergraph into a (dynamic or static) tight graph representation
// NB: The ID of each node in the tight graph representation corresponds to the position of the node in the node ordering
template <typename Hypergraph>
mt_kahypar_hypergraph_t HypergraphConverter::hypergraph_to_tight_graph(Hypergraph &hypergraph, std::vector<NodeID> &nodeOrdering, const mt_kahypar_preset_type_t tightGraphPresetType, const bool stable_construction_of_incident_edges)
{
    // Get the number of nodes and edges in the tight graph
    const NodeIndex numNodes = hypergraph.initialNumNodes();
    const EdgeIndex numEdges = hypergraph.initialNumEdges();
    // Initialize the edge list list of the tight graph and the list of the weight of the edges and of the nodes
    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edgeList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(2 * numEdges);
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edgesWeights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(numEdges);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> nodesWeights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(numNodes);

    // Store for each node the position in the node ordering, so that we can directly access the position of each hyperedge pin
    std::vector<NodeIndex> nodeOrderingPosition(numNodes);
    forRangeSequentialOrParallel(0, numNodes, i, NodeIndex)
    {
        nodeOrderingPosition[nodeOrdering[i]] = i;
    }
    endFor;

    // Fill the edge list
    // NB: While doing so, we can also directly fill the edge and nod weights weights
    // NB2: We do not ignore disabled edges here, since there usually should not be any and we must build a valid consecutive adjacency list
    forRangeSequentialOrParallel(0, numEdges, i, EdgeIndex)
    {
        // Store the weight of the edge
        edgesWeights[i] = hypergraph.edgeWeight(i);

        // Store the positions of the last two pins of the hyperedge in the ordering
        NodeIndex positionOfLastPinInOrdering = 0;
        NodeIndex positionOfSecondLastPinInOrdering = 0;

        // Loop over the pins of the edge
        for (NodeID pinID : hypergraph.pins(i))
        {
            // Get the position of the pin in the ordering
            NodeIndex position = nodeOrderingPosition[pinID];

            // Store the weight of the pin
            nodesWeights[position] = hypergraph.nodeWeight(pinID);

            // Check if we found a pin with a higher/later position
            if (position > positionOfLastPinInOrdering)
            {
                positionOfSecondLastPinInOrdering = positionOfLastPinInOrdering;
                positionOfLastPinInOrdering = position;
            }
            else if (position > positionOfSecondLastPinInOrdering)
            {
                positionOfSecondLastPinInOrdering = position;
            }
        }
        // Store an edge connecting the last two pins of the hyperedge in the ordering
        edgeList[2 * i] = positionOfSecondLastPinInOrdering;
        edgeList[2 * i + 1] = positionOfLastPinInOrdering;
    }
    endFor;

    return MtKaHyParIO::mt_kahypar_create_graph(tightGraphPresetType,
                                                numNodes,
                                                numEdges,
                                                edgeList.get(),
                                                edgesWeights.get(),
                                                nodesWeights.get(),
                                                stable_construction_of_incident_edges);
}

// Converts a (dynamic or static) hypergraph into its (dynamic or static) digraph representation
// NB: The digraph representation is usually directed, but the graph datastructure only supports undirected graphs (i.e. forward and backward edges)
//     Therefore, on must make sure to not use any of the backward edges, which can be done by checking the IDs of the nodes in the digraph representation.
template <typename Hypergraph>
mt_kahypar_hypergraph_t HypergraphConverter::hypergraph_to_digraph(Hypergraph &hypergraph, const mt_kahypar_preset_type_t digraphPresetType)
{
    // Extract the number of nodes, edges and pins from the hypergraph
    const NodeIndex hypergraphNumNodes = hypergraph.initialNumNodes();
    const EdgeIndex hypergraphNumEdges = hypergraph.initialNumEdges();
    const NodeIndex hypergraphNumPins = hypergraph.initialNumPins();
    // Get the number of nodes and edges in the residual digraph
    const NodeIndex numNodes = hypergraphNumNodes + 2 * hypergraphNumEdges;
    const EdgeIndex numEdges = hypergraphNumEdges + 2 * hypergraphNumPins;

    // Initialize the edge list list of the residual digraph and the list of the weight of the edges and of the nodes
    // We denote n := hypergraphNumNodes, m := hypergraphNumEdges, p := hypergraphNumPins
    // NB: The node weights list is build as follows: | 0 to (n-1) | 0 to (m-1) | 0 to (m-1) |
    //     So first we list the nodes (v) of the hypergraph, then the source nodes (e-) of the hyperedges, and finally the target nodes (e+) of the hyperedges.
    // NB2: The edge weights list is build as follows: | 0 to (m-1) | 0 to (p-1) | 0 to (p-1) |
    //     So first we list the edges (e- -> e+) of the hypergraph, then the outgoing edges (v -> e-) of the nodes, and finally the incoming edges (e+ -> v) of the nodes.
    // NB3: The edge list follows the same format as the edge weights, we just double the space (i.e. 2 * i and 2 * i + 1 in edge list map both to i in edge weights list).
    std::unique_ptr<mt_kahypar_hypernode_id_t[]> edgeList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(2 * numEdges);
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edgesWeights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(numEdges);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> nodesWeights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(numNodes);

    // Compute the prefix sum of the node degrees
    // NB: While doing so, we can also directly fill (some) node weights
    // NB2: We do not ignore disabled nodes here, since there usually should not be any and we must build a valid consecutive adjacency list
    std::unique_ptr<size_t[]> degreePrefixSum = std::make_unique<size_t[]>(numNodes + 1);
    forRangeSequentialOrParallel(0, hypergraphNumNodes, i, NodeIndex)
    {
#ifdef SMHM_PARALLEL
        // In the parallel case, we need to compute the prefix sum later
        degreePrefixSum[i + 1] = hypergraph.nodeDegree(i);
#else
        // In the sequential case, we can directly compute the prefix sum
        degreePrefixSum[i + 1] = hypergraph.nodeDegree(i) + (i > 0) * degreePrefixSum[i];
#endif
        // Store the weight of the node
        nodesWeights[i] = hypergraph.nodeWeight(i);
    }
    endFor;

#ifdef SMHM_PARALLEL
    // Compute the prefix sum of the number of pins for the parallel case
    mt_kahypar::parallel_prefix_sum(degreePrefixSum.get(), degreePrefixSum.get() + numNodes + 1, degreePrefixSum.get(), std::plus<>(), UL(0));
#endif

    // Add the outgoing edges (v -> e-) and the incoming edges (e+ -> v) of the nodes
    // NB: We do not ignore disabled nodes or edges here, since there usually should not be any and we must build a valid consecutive adjacency list
    forAllNodesSequentialOrParallel(hypergraph, nodeID)
    {
        NodeIndex numHandledIncidentEdges = 0;
        // Loop over the pins of the edge
        for (NodeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            EdgeID outgoingEdgeID = hypergraphNumEdges + degreePrefixSum[nodeID] + numHandledIncidentEdges;
            EdgeID incomingEdgeID = hypergraphNumEdges + hypergraphNumPins + degreePrefixSum[nodeID] + numHandledIncidentEdges;

            // Store the edge connecting the node with the source of the hyperedge in the digraph
            edgeList[2 * outgoingEdgeID] = nodeID;
            edgeList[2 * outgoingEdgeID + 1] = hypergraphNumNodes + edgeID;
            edgesWeights[outgoingEdgeID] = std::numeric_limits<mt_kahypar_hyperedge_weight_t>::max();
            // Store the edge connecting the target of the hyperedge with the node in the digraph
            edgeList[2 * incomingEdgeID] = hypergraphNumNodes + hypergraphNumEdges + edgeID;
            edgeList[2 * incomingEdgeID + 1] = nodeID;
            edgesWeights[incomingEdgeID] = std::numeric_limits<mt_kahypar_hyperedge_weight_t>::max();
            // Increment the number of handled incident edges
            numHandledIncidentEdges++;
        }
    }
    endFor;

    // Add the edges between the source and target of the hyperedges (e- -> e+)
    // NB: We do not ignore disabled edges here, since there usually should not be any and we must build a valid consecutive adjacency list
    forRangeSequentialOrParallel(0, hypergraphNumEdges, i, EdgeIndex)
    {
        // Store the weight of the edge
        edgesWeights[i] = hypergraph.edgeWeight(i);

        NodeID edgeSourceID = hypergraphNumNodes + i;
        NodeID edgeTargetID = hypergraphNumNodes + hypergraphNumEdges + i;

        // Store an edge connecting the source and target of the hyperedge in the digraph
        edgeList[2 * i] = edgeSourceID;
        edgeList[2 * i + 1] = edgeTargetID;

        // Store the weight of the source and target nodes of the hyperedge
        // NB: We simply use the edge weight as the weight of the source and target nodes
        nodesWeights[edgeSourceID] = hypergraph.edgeWeight(i);
        nodesWeights[edgeTargetID] = hypergraph.edgeWeight(i);
    }
    endFor;

    return MtKaHyParIO::mt_kahypar_create_graph(digraphPresetType,
                                                numNodes,
                                                numEdges,
                                                edgeList.get(),
                                                edgesWeights.get(),
                                                nodesWeights.get());
}

#endif // end of SMHM_HYPERGRAPH_CONVERTER_H