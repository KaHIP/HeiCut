/******************************************************************************
 * hypergraph_converter.cpp
 * *
 * Converts a hypergraph into different formats (e.g. static, dynamic) or into
 * graph representations (e.g. tight graph).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "hypergraph_converter.h"

// Converts a static hypergraph into a dynamic hypergraph and return whether the hypergraph has weighted edges
mt_kahypar_hypergraph_t HypergraphConverter::static_hypergraph_to_dynamic_hypergraph(StaticHypergraph &hypergraph, bool &hasWeightedEdges)
{
    // Assume that the hypergraph is unweighted by default
    hasWeightedEdges = false;
    // Get the number of nodes, edges and pins in the dynamic hypergraph
    const NodeIndex numNodes = hypergraph.initialNumNodes();
    const EdgeIndex numEdges = hypergraph.initialNumEdges();
    const NodeIndex numPins = hypergraph.initialNumPins();
    // Initialize the adjacency list of the dynamic hypergraph and the list of the weight of the edges and of the nodes
    std::unique_ptr<size_t[]> edgesAdjIndices = std::make_unique<size_t[]>(numEdges + 1);
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> edgesAdjList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(numPins);
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edgesWeights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(numEdges);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> nodesWeights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(numNodes);

    // Fill the adjacency indices by computing the prefix sum of the number of pins
    // NB: While doing so, we can also directly fill the edge weights and determine whether the hypergraph has weighted edges
    // NB2: We do not ignore disabled edges here, since there usually should not be any and we must build a valid consecutive adjacency list
    forRangeSequentialOrParallel(0, numEdges, i, EdgeIndex)
    {
#ifdef SMHM_PARALLEL
        // In the parallel case, we need to compute the prefix sum later
        edgesAdjIndices[i + 1] = hypergraph.edgeSize(i);
#else
        // In the sequential case, we can directly compute the prefix sum
        edgesAdjIndices[i + 1] = hypergraph.edgeSize(i) + (i > 0) * edgesAdjIndices[i];
#endif
        // Store the weight of the edge
        edgesWeights[i] = hypergraph.edgeWeight(i);

        // Check if the hypergraph has weighted edges
        if (hypergraph.edgeWeight(i) != 1 && !hasWeightedEdges)
            hasWeightedEdges = true;
    }
    endFor;

#ifdef SMHM_PARALLEL
    // Compute the prefix sum of the number of pins for the parallel case
    mt_kahypar::parallel_prefix_sum(edgesAdjIndices.get(), edgesAdjIndices.get() + numEdges + 1, edgesAdjIndices.get(), std::plus<>(), UL(0));
#endif

    // Build the adjacency list (will be used to build the dynamc hypergraph)
    // NB: We do not ignore disabled edges or pins here, since there usually should not be any and we must build a valid consecutive adjacency list
    forAllEdgesSequentialOrParallel(hypergraph, edgeID)
    {
        NodeIndex numPinsOfEdgeAdded = 0;
        // Loop over the pins of the edge
        for (NodeID pinID : hypergraph.pins(edgeID))
        {
            // Sdd the pin to the adjacency list and store its weight
            edgesAdjList[edgesAdjIndices[edgeID] + (numPinsOfEdgeAdded++)] = pinID;
            nodesWeights[pinID] = hypergraph.nodeWeight(pinID);
        }
    }
    endFor;

    // Create the dynamic hypergraph
    return mt_kahypar_create_hypergraph(HIGHEST_QUALITY,
                                        numNodes,
                                        numEdges,
                                        edgesAdjIndices.get(),
                                        edgesAdjList.get(),
                                        edgesWeights.get(),
                                        nodesWeights.get());
}