/******************************************************************************
 * pruner.h
 * *
 * Applies the pruning rules to the hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_PRUNER_H
#define SMHM_PRUNER_H

#include <cassert>
#include <vector>
#include <deque>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
#include "lib/data_structure/union_find/union_find_sequential.h"
#include "lib/data_structure/union_find/union_find_parallel.h"
#include "lib/data_structure/union_find/union_find_base.h"
// Google headers
#include <sparsehash/dense_hash_map>
// Mt-KaHyPar headers
#include "mt-kahypar/parallel/stl/scalable_vector.h"
// oneTBB headers
#ifdef SMHM_PARALLEL
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/blocked_range.h>
#endif

class Pruner
{
private:
    // Factory function to create the appropriate union find data structure object based on isParallel
    std::unique_ptr<UnionFindBase> create_union_find(const NodeIndex initialNumNodes);

    // Disable a hyperedge and increase the number of removed hyperedges
    // NB: For the parallel version, we need to use this because removeEdge is not thread-save, but even
    // for the sequential case, this is much faster than removeEdge
    void disable_hyperedge(StaticHypergraph &hypergraph, const EdgeID edgeID);

    // Get the weighted node degrees of the hypergraph
    std::vector<EdgeWeight> get_weighted_node_degrees(const StaticHypergraph &hypergraph);

public:
    // Remove hyperedges of size one or with weight zero
    // NB: The passed hypergraph is modified in place
    void remove_hyperedges_of_size_one_or_weight_zero(StaticHypergraph &hypergraph);

    // Compute a naive estimate of the minimum cut using the minimum weighted node degree
    CutValue compute_naive_mincut_estimate(const StaticHypergraph &hypergraph);

    // Contract hyperedges that are not lighter than the given estimate (VieCut pruning rule 1)
    StaticHypergraph contract_hyperedges_not_lighter_than_estimate(StaticHypergraph &hypergraph,
                                                                   mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                   CutValue mincutEstimate,
                                                                   const bool findAllMincuts = false);

    // Contract overlaps that are not lighter than the given estimate
    StaticHypergraph contract_overlaps_not_lighter_than_estimate(StaticHypergraph &hypergraph,
                                                                 mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                 CutValue mincutEstimate,
                                                                 const bool findAllMincuts = false);

    // Contract strictly nested isolated substructures of parent hyperedges
    // NB: Nested isolated substructures are never part of any mincut, thus we can also contract them if we want to find all mincuts without any modification
    StaticHypergraph contract_strictly_nested_isolated_substructures(StaticHypergraph &hypergraph,
                                                                     mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID);

    // Contract hyperedges of size two that can be shifted to one side of the cut (VieCut pruning rule 2)
    StaticHypergraph contract_shiftable_hyperedges_of_size_two(StaticHypergraph &hypergraph,
                                                               mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                               const CutValue mincutEstimate,
                                                               const bool findAllMincuts = false);

    // Contract hyperedges of size two that meet the triangle conditions (VieCut pruning rules 3 and 4)
    StaticHypergraph contract_triangle_hyperedges_of_size_two(StaticHypergraph &hypergraph,
                                                              mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                              const CutValue mincutEstimate,
                                                              const bool findAllMincuts = false);
};

// Create the appropriate union find datastructure object based on the number of threads used
inline std::unique_ptr<UnionFindBase> Pruner::create_union_find(const NodeIndex initialNumNodes)
{
#ifdef SMHM_PARALLEL
    return std::make_unique<UnionFindParallel>(initialNumNodes);
#else
    return std::make_unique<UnionFindSequential>(initialNumNodes);
#endif
}

// Disable a hyperedge and increase the number of removed hyperedges
// NB: For the parallel version, we need to use this because removeEdge is not thread-save, but even
// for the sequential case, this is much faster than removeEdge
inline void Pruner::disable_hyperedge(StaticHypergraph &hypergraph, const EdgeID edgeID)
{
    hypergraph.disableHyperedge(edgeID);
    hypergraph.setNumRemovedHyperedges(hypergraph.numRemovedHyperedges() + 1);
}

// Remove hyperedges of size one or with weight zero
// NB: The passed hypergraph is modified in place
inline void Pruner::remove_hyperedges_of_size_one_or_weight_zero(StaticHypergraph &hypergraph)
{
    // Remove hyperedges of size one or with weight zero
    forAllEdgesSequentialOrParallel(hypergraph, edgeID)
    {
        if (hypergraph.edgeIsEnabled(edgeID) && (hypergraph.edgeSize(edgeID) == 1 || hypergraph.edgeWeight(edgeID) == 0))
            // Disable the hyperedge and increase the number of removed hyperedges
            disable_hyperedge(hypergraph, edgeID);
    }
    endFor;
}

// Compute a naive estimate of the minimum cut using the minimum weighted node degree
inline CutValue Pruner::compute_naive_mincut_estimate(const StaticHypergraph &hypergraph)
{
#ifdef SMHM_PARALLEL
    CutValue minCutEstimate = tbb::parallel_reduce(
        tbb::blocked_range<NodeID>(0, hypergraph.initialNumNodes()),
        std::numeric_limits<CutValue>::max(),
        [&](const tbb::blocked_range<NodeID> &range, CutValue localMin) -> CutValue
        {
            for (NodeID nodeID = range.begin(); nodeID != range.end(); ++nodeID)
            {
                // Ignore disabled nodes
                if (!hypergraph.nodeIsEnabled(nodeID))
                    continue;
                // Compute the weighted node degree of the node
                CutValue weightedNodeDegree = 0;
                for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
                    if (hypergraph.edgeIsEnabled(edgeID))
                        weightedNodeDegree += hypergraph.edgeWeight(edgeID);
                // Update thelocal  minimum cut estimate (if necessary)
                localMin = std::min(localMin, weightedNodeDegree);
            }
            return localMin;
        },
        // Reduce local minimum cut estimates to a global minimum cut estimate
        [](CutValue a, CutValue b) -> CutValue
        {
            return std::min(a, b);
        });
#else
    CutValue minCutEstimate = std::numeric_limits<CutValue>::max();
    for (NodeID nodeID : hypergraph.nodes())
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Compute the weighted node degree of the node
        CutValue weightedNodeDegree = 0;
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
            if (hypergraph.edgeIsEnabled(edgeID))
                weightedNodeDegree += hypergraph.edgeWeight(edgeID);
        // Update the minimum cut estimate (if necessary)
        minCutEstimate = std::min(minCutEstimate, weightedNodeDegree);
    }
#endif

    return minCutEstimate;
}

// Get the weighted node degrees of the hypergraph
inline std::vector<EdgeWeight> Pruner::get_weighted_node_degrees(const StaticHypergraph &hypergraph)
{
    std::vector<EdgeWeight> weightedNodeDegrees(hypergraph.initialNumNodes(), 0);

    forAllNodesSequentialOrParallel(hypergraph, nodeID)
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continueFor;

        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
            if (hypergraph.edgeIsEnabled(edgeID))
                weightedNodeDegrees[nodeID] += hypergraph.edgeWeight(edgeID);
    }
    endFor;

    return weightedNodeDegrees;
}

// Contract hyperedges that are not lighter than the given estimate (VieCut pruning rule 1)
inline StaticHypergraph Pruner::contract_hyperedges_not_lighter_than_estimate(StaticHypergraph &hypergraph,
                                                                              mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                              CutValue mincutEstimate,
                                                                              const bool findAllMincuts)
{
    // If we want to find all mincuts, the inequality must be strict, which is the same as increasing the mincut estimate by 1
    if (findAllMincuts)
        mincutEstimate += 1;

    // Get the initial number of nodes
    // NB: For static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();

    // Make sure that clusterID has the correct size
    assert(clusterID.size() == initialNumNodes);

    // Initialize the union-find data structure, which keeps track of the nodes that can be contracted together
    auto unionFindPtr = create_union_find(initialNumNodes);

    // Iterate over all hyperedges
    forAllEdgesSequentialOrParallel(hypergraph, edgeID)
    {
        // Ignore disabled hyperedges or hyperedges that are lighter than the mincut estimate or hyperedges with only one pin
        if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeWeight(edgeID) < mincutEstimate || hypergraph.edgeSize(edgeID) <= 1)
            continueFor;

        // Store the first enabled pin of the hyperedge
        NodeID firstPinID = std::numeric_limits<NodeID>::max();
        // Iterate over the pins of the hyperedge
        for (NodeID pinID : hypergraph.pins(edgeID))
        {
            // Ignore disabled pins
            if (!hypergraph.nodeIsEnabled(pinID))
                continue;

            // Store the first enabled pin of the hyperedge
            if (firstPinID == std::numeric_limits<NodeID>::max())
                firstPinID = pinID;
            // Otherwise, merge the pins
            else
                unionFindPtr->Union(firstPinID, pinID);
        }
        // Disable the hyperedge and increase the number of removed hyperedges
        disable_hyperedge(hypergraph, edgeID);
    }
    endFor;

    // Get for each node the cluster id by finding the representative in the union-find data structure
    forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
    {
        clusterID[i] = unionFindPtr->Find(i);
    }
    endFor;

    return hypergraph.contract(clusterID);
}

// Contract overlaps that are not lighter than the given estimate
inline StaticHypergraph Pruner::contract_overlaps_not_lighter_than_estimate(StaticHypergraph &hypergraph,
                                                                            mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                            CutValue mincutEstimate,
                                                                            const bool findAllMincuts)
{
    // If we want to find all mincuts, the inequality must be strict, which is the same as increasing the mincut estimate by 1
    if (findAllMincuts)
        mincutEstimate += 1;

    // Get the weighted node degrees of the hypergraph
    std::vector<EdgeWeight> weightedNodeDegrees = get_weighted_node_degrees(hypergraph);

    // Get the initial number of nodes
    // NB: For static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();

    // Make sure that clusterID has the correct size
    assert(clusterID.size() == initialNumNodes);

    // Initialize the union-find data structure, which keeps track of the nodes that can be contracted together
    auto unionFindPtr = create_union_find(initialNumNodes);

    // Intitialize the hash vector to store the total weight of shared hyperedges between the current node and any other node.
    // NB: The first element of the pair is the total weight to the other node (whose id is equal to index), the second element
    //     is the current node id for which the total weight is stored. This allows us to use the same vector for all nodes, since
    //     we exactly know when to to reset the values in the vector.
    std::pair<EdgeWeight, NodeID> hashVectorInitialValue = {0, std::numeric_limits<NodeID>::max()};
#ifdef SMHM_PARALLEL
    tbb::enumerable_thread_specific<std::vector<std::pair<EdgeWeight, NodeID>>> hashVectorPerThread(
        [&]
        { return std::vector<std::pair<EdgeWeight, NodeID>>(initialNumNodes, hashVectorInitialValue); });
#else
    std::vector<std::pair<EdgeWeight, NodeID>> hashVector(initialNumNodes, hashVectorInitialValue);
#endif

    // Iterate over all hypernodes
    forAllNodesSequentialOrParallel(hypergraph, nodeID)
    {
        // Ignore disabled nodes or nodes with a weighted degree smaller than the mincut estimate
        if (!hypergraph.nodeIsEnabled(nodeID) || weightedNodeDegrees[nodeID] < mincutEstimate)
            continueFor;

#ifdef SMHM_PARALLEL
        // Get the hash vector for the current thread
        std::vector<std::pair<EdgeWeight, NodeID>> &hashVector = hashVectorPerThread.local();
#endif

        // Iterate over the incident hyperedges of the node
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            // Ignore disabled hyperedges
            if (!hypergraph.edgeIsEnabled(edgeID))
                continue;

            // Iterate over the neighbors of the node
            for (NodeID pinID : hypergraph.pins(edgeID))
            {
                // Ignore disabled pins, pins that are equal to the node or pins with a weighted degree smaller than the mincut estimate
                if (!hypergraph.nodeIsEnabled(pinID) || pinID == nodeID || weightedNodeDegrees[pinID] < mincutEstimate)
                    continue;

                // Reset the total weight if we look for the first time
                // at the connection beween nodeID and pinID
                if (hashVector[pinID].second != nodeID)
                {
                    hashVector[pinID].first = 0;
                    hashVector[pinID].second = nodeID;
                }

                // Do not merge the nodes if they were already merged
                if (hashVector[pinID].first >= mincutEstimate)
                    continue;

                // Increase the total weight of shared hyperedges between nodeID and pinID
                hashVector[pinID].first += hypergraph.edgeWeight(edgeID);

                // Merge the nodes of the pins if the total weight is now larger than the mincut estimate
                if (hashVector[pinID].first >= mincutEstimate)
                    unionFindPtr->Union(nodeID, pinID);
            }
        }
    }
    endFor;

    // Get for each node the cluster id by finding the representative in the union-find data structure
    forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
    {
        clusterID[i] = unionFindPtr->Find(i);
    }
    endFor;

    return hypergraph.contract(clusterID);
}

// Contract strictly nested isolated substructures of parent hyperedges
// NB: Nested isolated substructures are never part of any mincut, thus we can also contract them if we want to find all mincuts without any modification
inline StaticHypergraph Pruner::contract_strictly_nested_isolated_substructures(StaticHypergraph &hypergraph,
                                                                                mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID)
{
    // Get the initial number of nodes
    // NB: For static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();

    // Make sure that clusterID has the correct size
    assert(clusterID.size() == initialNumNodes);

    // Initialize the union-find data structure, which keeps track of the nodes that can be contracted together
    auto unionFindPtr = create_union_find(initialNumNodes);

    // Intitialize the hash vector to store the total number of pins that the current hyperedge shares with any other hyperedge.
    // NB: The first element of the pair is the total number of pins shared with the other hyperedge (whose id is equal to index),
    //     the second element is the current hyperedge id for which the total number is stored. This allows us to use the same vector
    //     for all hyperedges, since we exactly know when to to reset the values in the vector.
    //
    // Initialize also the mark vector, which marks the pins that can leave the parent hyperedge when starting at the pin and using a
    // path that does not include a super parent hyperedge
    // NB: We only use one vector (per thread) since we go once over the hyperedges and flag for each hyperedge its pins
    //     by setting at the corresponding index (= pin id) the value to the id of the hyperedge. This means
    //     that we simply overwrite old values and do not have to reset the vector for each hyperedge.
    std::pair<NodeIndex, EdgeID> hashVectorInitialValue = {0, std::numeric_limits<EdgeID>::max()};
    NodeID markedPinsInitialValue = std::numeric_limits<NodeID>::max();

#ifdef SMHM_PARALLEL
    tbb::enumerable_thread_specific<std::vector<std::pair<NodeIndex, EdgeID>>> hashVectorPerThread(
        [&]
        { return std::vector<std::pair<NodeIndex, EdgeID>>(hypergraph.initialNumEdges(), hashVectorInitialValue); });

    tbb::enumerable_thread_specific<std::vector<NodeID>> markedPinsPerThread(
        [&]
        { return std::vector<NodeID>(initialNumNodes, markedPinsInitialValue); });
#else
    std::vector<std::pair<NodeIndex, EdgeID>> hashVector(hypergraph.initialNumEdges(), hashVectorInitialValue);
    std::vector<NodeID> markedPins(initialNumNodes, markedPinsInitialValue);
#endif

    // Iterate over all hyperedges, which will be considered as the parent hyperedges
    forAllEdgesSequentialOrParallel(hypergraph, edgeID)
    {
        // Ignore disabled hyperedges or hyperedges with only two pins (since they cannot have strictly nested isolated hyperedges with at least two pins)
        if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeSize(edgeID) <= 2)
            continueFor;

#ifdef SMHM_PARALLEL
        // Get the hash vector for the current thread
        std::vector<std::pair<NodeIndex, EdgeID>> &hashVector = hashVectorPerThread.local();
        // Get the marked pins vector for the current thread
        std::vector<NodeID> &markedPins = markedPinsPerThread.local();
#endif
        // Stores the nested hyperedges of the current parent hyperedge
        std::deque<EdgeID> nestedHyperedges;

        // Go a first time over the pins of the hyperedge to compute the total number of shared pins with other hyperedges
        for (NodeID pinID : hypergraph.pins(edgeID))
        {
            // Ignore disabled pins
            if (!hypergraph.nodeIsEnabled(pinID))
                continue;

            // Iterate over the incident hyperedges of the pin
            for (EdgeID incidentEdgeID : hypergraph.incidentEdges(pinID))
            {
                // Ignore disabled hyperedges or hyperedges that are the current parent hyperedge
                if (!hypergraph.edgeIsEnabled(incidentEdgeID) || incidentEdgeID == edgeID)
                    continue;

                // Reset the total sum if we look for the first time
                // at the shared pins beween edgeID and incidentEdgeID
                if (hashVector[incidentEdgeID].second != edgeID)
                {
                    hashVector[incidentEdgeID].first = 0;
                    hashVector[incidentEdgeID].second = edgeID;
                }

                // Increase the total sum of shared pins between edgeID and incidentEdgeID
                hashVector[incidentEdgeID].first += 1;

                // Check if we found a strictly nested hyperedge and if so, store it (if not already done)
                if (hypergraph.edgeSize(incidentEdgeID) < hypergraph.edgeSize(edgeID) && hashVector[incidentEdgeID].first == hypergraph.edgeSize(incidentEdgeID))
                    nestedHyperedges.push_back(incidentEdgeID);
            }
        }

        // Stop early if the current parent hyperedge has no nested hyperedge
        if (nestedHyperedges.empty())
            continueFor;

        // Stores the pins that can leave the current parent hyperedge when starting at the pin and using a path
        // that does not include the parent hyperedge or a super parent hyperedge (i.e. a hyperedge that entirely contains the current parent hyperedge)
        // NB: Every node is added at most once to the queue.
        std::deque<EdgeID> pinsThatCanLeaveParentEdgeByNonSuperParentEdge;

        // Go a second time over the pins of the hyperedge to find the pins that can DIRECTLY leave the current parent hyperedge by using a non (super) parent hyperedge
        for (NodeID pinID : hypergraph.pins(edgeID))
        {

            // Ignore disabled pins
            if (!hypergraph.nodeIsEnabled(pinID))
                continue;

            // Iterate over the incident hyperedges of the pin
            for (EdgeID incidentEdgeID : hypergraph.incidentEdges(pinID))
            {
                // Ignore hyperedges that were not considered in the first pass
                // NB: This also inplicity ignores the parent hyperedge (i.e. incidentEdgeID != edgeID)
                if (hashVector[incidentEdgeID].second != edgeID)
                    continue;

                // The pin can directly leave the current parent hyperedge via the incident hyperedge if the shared number of pins is
                // smaller than both the size of the incident hyperedge and the size of the current parent hyperedge. The latter
                // assures that the incident hyperedge is not a super parent hyperedge.
                if (markedPins[pinID] != edgeID &&
                    hashVector[incidentEdgeID].first < hypergraph.edgeSize(incidentEdgeID) &&
                    hashVector[incidentEdgeID].first < hypergraph.edgeSize(edgeID))
                {
                    markedPins[pinID] = edgeID;
                    pinsThatCanLeaveParentEdgeByNonSuperParentEdge.push_back(pinID);
                }
            }
        }

        // We first assume that all nested hyperedges are valid
        EdgeID numValidNestedHyperedges = nestedHyperedges.size();

        // Iterate over the queue until it is empty
        while (!pinsThatCanLeaveParentEdgeByNonSuperParentEdge.empty())
        {
            // Get the next pin in the queue
            NodeID pinID = pinsThatCanLeaveParentEdgeByNonSuperParentEdge.front();
            pinsThatCanLeaveParentEdgeByNonSuperParentEdge.pop_front();
            // Iterate over the incident hyperedges of the pin
            for (EdgeID incidentEdgeID : hypergraph.incidentEdges(pinID))
            {
                // Ignore hyperedges that were not considered in the first pass
                // NB: This also inplicity ignores the parent hyperedge (i.e. incidentEdgeID != edgeID)
                if (hashVector[incidentEdgeID].second != edgeID)
                    continue;

                // Check if we found a nested hyperedge and if so, invalidate it (since its pins can leave the current parent hyperedge)
                if (hypergraph.edgeSize(incidentEdgeID) < hypergraph.edgeSize(edgeID) && hashVector[incidentEdgeID].first == hypergraph.edgeSize(incidentEdgeID))
                {
                    // Invalidate the nested hyperedge
                    hashVector[incidentEdgeID].first = 0;
                    numValidNestedHyperedges--;
                    // Add the pins of the nested hyperedge to the queue (if not already done)
                    for (NodeID nestedPinID : hypergraph.pins(incidentEdgeID))
                        if (hypergraph.nodeIsEnabled(nestedPinID) && markedPins[nestedPinID] != edgeID)
                        {
                            markedPins[nestedPinID] = edgeID;
                            pinsThatCanLeaveParentEdgeByNonSuperParentEdge.push_back(nestedPinID);
                        }
                }
            }
        }

        // Stop early if the current parent hyperedge has no VALID nested hyperedge
        if (numValidNestedHyperedges == 0)
            continueFor;

        // Iterate over the nested hyperedges and contract them if they are not invalid
        for (EdgeID nestedEdgeID : nestedHyperedges)
        {
            // Ingore invalid nested hyperedges
            if (hashVector[nestedEdgeID].first != hypergraph.edgeSize(nestedEdgeID))
                continue;

            // Store the first enabled pin of the hyperedge
            NodeID firstPinID = std::numeric_limits<NodeID>::max();
            // Iterate over the pins of the hyperedge
            for (NodeID pinID : hypergraph.pins(nestedEdgeID))
            {
                // Ignore disabled pins
                if (!hypergraph.nodeIsEnabled(pinID))
                    continue;

                // Store the first enabled pin of the hyperedge
                if (firstPinID == std::numeric_limits<NodeID>::max())
                    firstPinID = pinID;
                // Otherwise, merge the pins
                else
                    unionFindPtr->Union(firstPinID, pinID);
            }
            // Disable the hyperedge and increase the number of removed hyperedges
            disable_hyperedge(hypergraph, edgeID);
        }
    }
    endFor;

    // Get for each node the cluster id by finding the representative in the union-find data structure
    forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
    {
        clusterID[i] = unionFindPtr->Find(i);
    }
    endFor;

    return hypergraph.contract(clusterID);
}

// Contract hyperedges of size two that can be shifted to one side of the cut (VieCut pruning rule 2)
inline StaticHypergraph Pruner::contract_shiftable_hyperedges_of_size_two(StaticHypergraph &hypergraph,
                                                                          mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                          const CutValue mincutEstimate,
                                                                          const bool findAllMincuts)
{
    // Get the weighted node degrees of the hypergraph
    // NB:  In the parallel version of VieCut, only uncontracted nodes are allowed to be contracted for PR 2, because they argue that
    //      the weighted node degree changes when nodes are contracted together. However, this is not done in the sequential version.
    //      Theoretically, the change of the weighted node degree should not break the rule, even if it is applied to adjacent edges,
    //      thus we omit this restriction (i.e. we follow the sequential version of VieCut even for our parallel version).
    //      It is only important that the inequality is strict, as explained below, which is the case for VieCut in the code but not in the paper.
    std::vector<EdgeWeight> weightedNodeDegrees = get_weighted_node_degrees(hypergraph);

    // Get the initial number of nodes
    // NB: For static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();

    // Make sure that clusterID has the correct size
    assert(clusterID.size() == initialNumNodes);

    // Initialize the union-find data structure, which keeps track of the nodes that can be contracted together
    auto unionFindPtr = create_union_find(initialNumNodes);

    forAllEdgesSequentialOrParallel(hypergraph, edgeID)
    {
        // Ignore disabled hyperedges and hyperedges that are not of size two
        if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeSize(edgeID) != 2)
            continueFor;

        // Get the weight of the hyperedge
        EdgeWeight edgeWeight = hypergraph.edgeWeight(edgeID);
        // Store whether we can shift the hyperedge
        bool canShiftHyperedge = false;
        // Store the first enabled pin of the hyperedge
        NodeID firstPinID = std::numeric_limits<NodeID>::max();

        // Iterate over the pins of the hyperedge
        for (NodeID pinID : hypergraph.pins(edgeID))
        {
            // Ignore disabled pins
            if (!hypergraph.nodeIsEnabled(pinID))
                continue;

            // Check if we can shift the hyperedge by moving the current pin to the other side of the cut
            // IMPORTANT: We would need to verify that the current pin has still other incident hyperedges apart from the current hyperedge,
            //            because otherwise the shifting idea of the rule does not work and the reduced hypergraph might have a larger mincut value.
            //            Example: An unweighted 2-uniform dumbbell hypergraph but where one side is a single node instead of a clique has obviously still
            //                     a mincut of 1, but the shifting rule is applicable to the "bridge" hyperedge, giving us the clique as a reduced hypergraph.
            //            However, this additional check is only necessary if we wan to find all minimum cuts. This is not necessary if we only compute the mincut value,
            //            because we update the mincut value before and after every pruning rule using the naive estimate, i.e. we still capture the correct mincut
            //            value and the increased mincut in the reduced hypergraph does not affect the result as we update the mincut only if we found a new lower value.
            // IMPORTANT 2: The inequality (2 * edgeWeight > weightedNodeDegrees[pinID]) MUST be strict because otherwise it can happen that two neighboring 2-uniform
            //              hyperedges both assume that they can be shifted but their common pin cannot be shifted to both sides at the same time.
            //              Example: A cycle of 7 nodes where we have the following 5 hyperedges with the respective weights in parenthesis:
            //                       {1, 2} (5), {2, 3} (5), {3, 4, 5} (10), {5, 6} (1), {6, 7, 1} (10). The minimum cut is 6 by crossing either {1, 2} & {5, 6} or {2, 3} & {5, 6}.
            //                       If the inequality is not strict, then {1, 2} and {2, 3} would be contracted and the minimum cut of the contracted hypergraph is 10 (isolating node 4 or 7).
            //                       Note that even updating the mincut estimate after each pruning rule would not help here, because before and after contraction the min weighted node degree is larger than 6.
            //              One could make the inequality non-strict if one would allow pins to take only part of a single contraction (similar to pruning rule 3 of VieCut)
            if (!canShiftHyperedge && 2 * edgeWeight > weightedNodeDegrees[pinID] && (!findAllMincuts || weightedNodeDegrees[pinID] > mincutEstimate))
                canShiftHyperedge = true;

            // Store the first enabled pin of the hyperedge
            if (firstPinID == std::numeric_limits<NodeID>::max())
                firstPinID = pinID;
            // Otherwise, merge the pins if we can shift the hyperedge
            else if (canShiftHyperedge)
                unionFindPtr->Union(firstPinID, pinID);
        }
        // Disable the hyperedge and increase the number of removed hyperedges (if we can shift the hyperedge)
        if (canShiftHyperedge)
            disable_hyperedge(hypergraph, edgeID);
    }
    endFor;

    // Get for each node the cluster id by finding the representative in the union-find data structure
    forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
    {
        clusterID[i] = unionFindPtr->Find(i);
    }
    endFor;

    return hypergraph.contract(clusterID);
}

// Contract hyperedges of size two that meet the triangle conditions (VieCut pruning rules 3 and 4)
inline StaticHypergraph Pruner::contract_triangle_hyperedges_of_size_two(StaticHypergraph &hypergraph,
                                                                         mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                         const CutValue mincutEstimate,
                                                                         const bool findAllMincuts)
{
    // Get the weighted node degrees of the hypergraph
    std::vector<EdgeWeight> weightedNodeDegrees = get_weighted_node_degrees(hypergraph);

    // Get the initial number of nodes
    // NB: For static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();

    // Make sure that clusterID has the correct size
    assert(clusterID.size() == initialNumNodes);

    // Initialize the union-find data structure, which keeps track of the nodes that can be contracted together
    auto unionFindPtr = create_union_find(initialNumNodes);

    // Intitialize the hash vector to store the id of the hyperedge of size two connecting the current node with any other node.
    // NB: The first element of the pair is the id of the hyperedge of size two to the other node (whose id is equal to index),
    //     the second element is the current node id for which the hyperedge id is stored. This allows us to use the same vector
    //     for all nodes, since we exactly know when to to reset the values in the vector.
    //
    // Also keep track of the already finished and the already contracted nodes

    std::pair<EdgeID, NodeID> hashVectorInitialValue = {std::numeric_limits<EdgeID>::max(), std::numeric_limits<NodeID>::max()};
#ifdef SMHM_PARALLEL
    tbb::enumerable_thread_specific<std::vector<std::pair<EdgeID, NodeID>>> hashVectorPerThread(
        [&]
        { return std::vector<std::pair<EdgeID, NodeID>>(initialNumNodes, hashVectorInitialValue); });

    std::vector<uint8_t> finished(initialNumNodes, false);
    std::vector<uint8_t> contracted(initialNumNodes, false);
#else
    std::vector<std::pair<EdgeID, NodeID>> hashVector(initialNumNodes, hashVectorInitialValue);
    std::vector<bool> finished(initialNumNodes, false);
    std::vector<bool> contracted(initialNumNodes, false);
#endif

    // Iterate over all hypernodes
    forAllNodesSequentialOrParallel(hypergraph, nodeID)
    {
        // Ignore disabled nodes or already finished nodes
        if (!hypergraph.nodeIsEnabled(nodeID) || finished[nodeID])
            continueFor;

#ifdef SMHM_PARALLEL
        // Get the hash vector for the current thread
        std::vector<std::pair<NodeIndex, EdgeID>> &hashVector = hashVectorPerThread.local();
#endif

        // Set the node as finished
        finished[nodeID] = true;

        // Iterate over the incident hyperedges of the node to mark the edges and the neighboring pins
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            // Ignore disabled hyperedges and hyperedges that are not of size two
            if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeSize(edgeID) != 2)
                continue;

            // Get the iterator range of the pins of the edge
            auto pinsIterator = hypergraph.pins(edgeID).begin();
            // Get the target of the hyperedge of size two
            NodeID targetID = 0;
            do
            {
                targetID = *(pinsIterator++);
            } while ((targetID == nodeID || !hypergraph.nodeIsEnabled(targetID)) && pinsIterator != hypergraph.pins(edgeID).end());

#ifndef SMHM_PARALLEL
            // Make sure that the implication (targetID < nodeID) => finished[targetID] == true holds (only for the sequential version, since
            // for the parallel version the processing order is not guaranteed to be the same as the order of the nodes in the hypergraph)
            // NB:  In the sequential version, it can happen that (targetID < nodeID) if the target was set to finished by another node with
            //      an even smaller ID, thus the target did not perform its own iteration. Therefore it can happen that target is finished
            //      and the current node is not finished, even when they are connected by a hyperedge of size 2.
            assert(targetID > nodeID || finished[targetID] == true);
#endif

            // Ignore those targets that have a smaller ID than the ID of the current node
            // NB: This guarantees that each hyperedge of size 2 is only put twice in the hashVector. HOWEVER this alone does not guarantee that
            //     each hyperedge of size 2 is only looked at twice, since we can reach the same target from different nodes and each time the target
            //     goes over all its incident edges . To ensure that each hyperedge of size 2 is only looked at twice, the finished array is used.
            //     Note that for the parallel version, the weak constraint targetID > nodeID is sufficient. We only impose for the sequential
            //     version that each hyperedge of size 2 is only looked at twice using the finished array (cf. below).
            // NB2: For the sequential version, targetID < nodeID would mean that finished[targetID] == true. But even if targetID > nodeID,
            //     it can be that finished[targetID] == true but then the edge (nodeID, targetID) was only considered once yet
            if (targetID > nodeID)
                // If so, mark the target by storing the hyperedge of size two connecting it to the node
                hashVector[targetID] = std::make_pair(edgeID, nodeID);
        }

        // Iterate over the incident hyperedges of the node to check the triangle conditions
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            // Ignore disabled hyperedges and hyperedges that are not of size two
            if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeSize(edgeID) != 2)
                continue;

            // Get the iterator range of the pins of the edge
            auto pinsIterator = hypergraph.pins(edgeID).begin();
            // Get the target of the hyperedge of size two
            NodeID targetID = 0;
            do
            {
                targetID = *(pinsIterator++);
            } while ((targetID == nodeID || !hypergraph.nodeIsEnabled(targetID)) && pinsIterator != hypergraph.pins(edgeID).end());

#ifdef SMHM_PARALLEL
            // For the parallel version, we omit the constraint that each edge should only be considered twice (and therefore we also allow
            // SOME targets to be finished, more specifically those with targetID > nodeID), because due to the parallel execution a single iteration
            // can be more expensive to find more potential triangles. Therefore for the parallel version, we only make sure that targetID > nodeID,
            // which is a weaker condition.
            if (targetID < nodeID)
                continue;
#else
            // For the sequential version, Ignore the edge if the target is already finished (only for the sequential version)
            // NB1: Finished targets already scanned their incident edges and thus we cannot do it again, since each edge should only be
            //      considered twice, i.e. once for each direction. The reason why we do not directly exclude all finished targets when
            //      building hashVector[targetID] above is that the targetOfTarget (i.e. the third node closing the triangle) is allowed to be
            //      finished but not the target (i.e. the second node which also scans its incident edges).
            // NB2: This includes the case that targetID < nodeID, since then the target was also already processed and marked as finished
            if (finished[targetID])
            {
                // Unmark the target
                hashVector[targetID] = std::make_pair(std::numeric_limits<EdgeID>::max(), std::numeric_limits<NodeID>::max());
                continue;
            }
            // Make sure that the implication (targetID < nodeID) => finished[targetID] == true holds (only for the sequential version)
            // Since we know here that finished[targetID] == 0 (because we did not reach the continue), it must be that targetID > nodeID
            assert(targetID > nodeID);
#endif

            // Set the target as finished
            finished[targetID] = true;

            // Store the weight of the main edge in the 2-uniform triangle
            EdgeWeight triangleMainEdgeWeight = hypergraph.edgeWeight(edgeID);
            // Store the lowest cut leaving the triangles where the current edge is part of
            EdgeWeight lowestCutLeavingTriangles = triangleMainEdgeWeight;
            // Store if we found any triangle (otherwise we do not have to check for VieCut pruning rule 4)
            bool hasFoundAnyTriangle = false;

            // Go over the incident edges of the target to find the triangles
            for (EdgeID backwardEdgeID : hypergraph.incidentEdges(targetID))
            {
                // Ignore disabled hyperedges and hyperedges that are not of size two
                if (!hypergraph.edgeIsEnabled(backwardEdgeID) || hypergraph.edgeSize(backwardEdgeID) != 2)
                    continue;

                // Get the iterator range of the pins of the edge
                auto targetPinsIterator = hypergraph.pins(backwardEdgeID).begin();
                // Get the target of the hyperedge of size two
                EdgeID targetOfTargetID = 0;
                do
                {
                    targetOfTargetID = *(targetPinsIterator++);
                } while ((targetOfTargetID == targetID || !hypergraph.nodeIsEnabled(targetOfTargetID)) && targetPinsIterator != hypergraph.pins(backwardEdgeID).end());

                // Ignore the edge if it was not marked by the node
                if (hashVector[targetOfTargetID].first == std::numeric_limits<EdgeID>::max() || hashVector[targetOfTargetID].second != nodeID)
                    continue;

                // Store the weights of the other two edges of the 2-uniform triangle
                EdgeWeight triangleForwardEdgeWeight = hypergraph.edgeWeight(hashVector[targetOfTargetID].first);
                EdgeWeight triangleBackwardEdgeWeight = hypergraph.edgeWeight(backwardEdgeID);

                // Increment the lowest cut leaving the triangles
                lowestCutLeavingTriangles += std::min(triangleForwardEdgeWeight, triangleBackwardEdgeWeight);
                hasFoundAnyTriangle = true;

                // Check VieCut pruning rule 3
                // IMPORTANT: We would need to verify that the nodeID (resp. targetID) has still other incident hyperedges apart from the current hyperedge,
                //            because otherwise the shifting idea of the rule does not work and the reduced hypergraph might have a larger mincut value.
                //            However, this additional check is only necessary if we wan to find all minimum cuts. This is not necessary if we only compute the mincut value,
                //            because we update the mincut value before and after every pruning rule using the naive estimate, i.e. we still capture the correct mincut
                //            value and the increased mincut in the reduced hypergraph does not affect the result as we update the mincut only if we found a new lower value.
                // IMPORTANT 2: Contrary to the pruning rule 2 of VieCut (cf. above), the inequalities do not have to be strict here, as we only allow nodes to be part of at most one contraction.
                //              Nonetheless it remains open if enforcing strict inequalities allows nodes to be part of more than one contraction without breaking the rule,
                //              as done for the adaptation of pruning rule 2 of VieCut (cf. above).
                //              HOWEVER: If we want to find all minimum cuts, we need to enforce strict inequalities, because otherwise it can happen that we destroy with
                //              the contraction a mincut in favour of another mincut, i.e. we do not keep all minimum cuts.
                EdgeWeight twiceForwardTriangleCut = 2 * (triangleMainEdgeWeight + triangleForwardEdgeWeight);
                EdgeWeight twiceBackwardTriangleCut = 2 * (triangleMainEdgeWeight + triangleBackwardEdgeWeight);
                if ((!findAllMincuts && twiceForwardTriangleCut >= weightedNodeDegrees[nodeID] && twiceBackwardTriangleCut >= weightedNodeDegrees[targetID]) ||
                    (findAllMincuts && twiceForwardTriangleCut > weightedNodeDegrees[nodeID] && twiceBackwardTriangleCut > weightedNodeDegrees[targetID] &&
                     weightedNodeDegrees[nodeID] > mincutEstimate && weightedNodeDegrees[targetID] > mincutEstimate))
                {
#ifdef SMHM_PARALLEL
                    // For the parallel version, we use a CAS operation to make sure that the nodes are not already contracted
                    // NB: We use two sepearate if statements to avoid the second CAS operation if the first one already failed
                    if (__sync_bool_compare_and_swap(&contracted[nodeID], false, true))
                        if (__sync_bool_compare_and_swap(&contracted[targetID], false, true))
                        {
                            unionFindPtr->Union(nodeID, targetID);
                            break;
                        }
#else
                    // For the sequential version, we can simply make sure that the nodes are not already contracted
                    if (!contracted[nodeID] && !contracted[targetID])
                    {
                        contracted[nodeID] = true;
                        contracted[targetID] = true;
                        unionFindPtr->Union(nodeID, targetID);
                        break;
                    }
#endif
                }
            }
            // Check VieCut pruning rule 4
            if (hasFoundAnyTriangle && ((!findAllMincuts && lowestCutLeavingTriangles >= mincutEstimate) || (findAllMincuts && lowestCutLeavingTriangles > mincutEstimate)))
            {
                contracted[nodeID] = true;
                contracted[targetID] = true;
                unionFindPtr->Union(nodeID, targetID);
            }

            // Unmark the target so that the other edge of the triangle is not considered
            hashVector[targetID] = std::make_pair(std::numeric_limits<EdgeID>::max(), std::numeric_limits<NodeID>::max());
        }
    }
    endFor;

    // Get for each node the cluster id by finding the representative in the union-find data structure
    forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
    {
        clusterID[i] = unionFindPtr->Find(i);
    }
    endFor;

    return hypergraph.contract(clusterID);
}

#endif // end of SMHM_PRUNER_H