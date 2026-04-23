/******************************************************************************
 * orderer.h
 * *
 * Computes an ordering of the nodes of a hypergraph. Note that one instance of
 * the Orderer class can only be used for one hypergraph (and its contracted versions).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_ORDERER_H
#define SMHM_ORDERER_H

#include <vector>
#include <type_traits>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
#include "lib/utils/random.h"
#include "lib/data_structure/bucket_max_queue/bucket_max_queue.h"
// Mt-KaHyPar headers
#include "kahypar-resources/datastructure/binary_heap.h"

// Wrapper class for the orderer so that we can define the KeyType of the orderer at runtime based on the config (only needed for the parallel version)
template <typename Hypergraph>
class OrdererWrapper
{
public:
    virtual ~OrdererWrapper() = default;
    virtual void compute_ordering(const Hypergraph &hypergraph,
                                  const NodeIndex numNodes,
                                  std::vector<NodeID> *nodeOrdering,
                                  std::vector<EdgeID> *edgeHeadOrdering,
                                  std::vector<NodeID> *edgeHead) = 0;
};

template <typename Hypergraph, typename KeyType>
class Orderer : public OrdererWrapper<Hypergraph>
{
private:
    // The weight used in the ordering function
    KeyType orderingWeight;
    // Whether the binary heap should be used
    const bool useBinaryHeap;
    // Whether the orderer works on a static hypergraph
    bool worksOnStaticHypergraph = false;
    // Initial number of edges in the hypergraph
    const EdgeIndex initialNumEdges;
    // Initial number of nodes in the hypergraph
    const NodeIndex initialNumNodes;
    // Count the number of nodes added to the ordering
    NodeIndex numNodesAddedToOrdering = 0;
    // Count the number of edges added to the ordering
    EdgeIndex numEdgesAddedToOrdering = 0;
    // Type of the ordering
    const OrderingType orderingType;
    // Store for each edge the number of pins added to the ordering
    std::vector<NodeIndex> edgeNumPinsAdded;
    // Store the permutation of the nodes when initializing the priority queue
    std::vector<NodeID> permutation;
    // Store the binary max heap (only if useBinaryHeap == true)
    kahypar::ds::BinaryMaxHeap<NodeID, KeyType> binaryMaxHeap;
    // Store the bucket max priority queue (only if useBinaryHeap == false)
    BucketMaxQueue<NodeID, NodeIndex, EdgeIndex> bucketMaxQueue;

    // Get the number of levels to elevate the non-ordered pins of an edge in the priority queue
    // NB: A value of 0 means that the non-ordered pins should not be elevated
    KeyType get_num_levels_to_elevate_non_ordered_pins_of_edge(const NodeIndex edgeNumOfPinsAdded,
                                                               const NodeIndex edgeSize,
                                                               const EdgeWeight edgeWeight);

    // Elevate the non-ordered pins of an edge in the bucket priority queue
    void elevate_non_ordered_pins_of_edge(const Hypergraph &hypergraph,
                                          const EdgeID edgeID,
                                          const KeyType numLevelsToElevate);

    // Reset the priority queue
    NodeID reset_priority_queue(const Hypergraph &hypergraph, const NodeIndex numNodes);

public:
    Orderer(const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges, const OrderingType orderingType, const bool hasWeightedEdges, MersenneTwister randEngine, const KeyType orderingWeight = 0);

    // Compute the ordering of the nodes of a hypergraph
    // The result is stored in the vectors nodeOrdering, edgeHeadOrdering, and edgeHead
    void compute_ordering(const Hypergraph &hypergraph,
                          const NodeIndex numNodes,
                          std::vector<NodeID> *nodeOrdering,
                          std::vector<EdgeID> *edgeHeadOrdering,
                          std::vector<NodeID> *edgeHead);
};

template <typename Hypergraph, typename KeyType>
Orderer<Hypergraph, KeyType>::Orderer(const NodeIndex initialNumNodes,
                                      const EdgeIndex initialNumEdges,
                                      const OrderingType orderingType,
                                      const bool hasWeightedEdges,
                                      MersenneTwister randEngine,
                                      const KeyType orderingWeight)
    : initialNumNodes(initialNumNodes),
      initialNumEdges(initialNumEdges),
      orderingType(orderingType),
      edgeNumPinsAdded(initialNumEdges, 0),
      orderingWeight(orderingWeight),
      worksOnStaticHypergraph(std::is_same<Hypergraph, StaticHypergraph>::value),
      useBinaryHeap(hasWeightedEdges || std::is_floating_point<KeyType>::value),
      // NB: It is important to use the initial number of nodes for the allocation here,
      //     because we use the node ID as the index of the vector, which is might be larger
      //     than the number of nodes (e.g. if some nodes were contracted)
      bucketMaxQueue(!useBinaryHeap ? initialNumNodes : 0),
      // NB: It is important to use the initial number of nodes for the allocation here,
      //     because we want to reuse the binary max heap also for the contracted versions
      //     of the hypergraph, which might have a different number of nodes
      binaryMaxHeap(useBinaryHeap ? initialNumNodes : 0),
      // NB: It is important to use the initial number of nodes for the allocation here,
      //     because we permutate the node IDs. Removed/disabled nodes are filtered out
      //     when initializing the priority queue
      permutation(initialNumNodes)
{
    // Initialize the random permutation of the nodes
    RandomFunctions::permutate_vector_good(&permutation, true, randEngine);
};

// Compute the ordering of the nodes of a hypergraph
// The result is stored in the vectors nodeOrdering, edgeHeadOrdering, and edgeHead
template <typename Hypergraph, typename KeyType>
inline void Orderer<Hypergraph, KeyType>::compute_ordering(const Hypergraph &hypergraph,
                                                           const NodeIndex numNodes,
                                                           std::vector<NodeID> *nodeOrdering,
                                                           std::vector<EdgeID> *edgeHeadOrdering,
                                                           std::vector<NodeID> *edgeHead)
{
    // Reset the data structures
    NodeID startNodeID = reset_priority_queue(hypergraph, numNodes);
    std::fill(edgeNumPinsAdded.begin(), edgeNumPinsAdded.end(), 0);
    numNodesAddedToOrdering = 0;
    numEdgesAddedToOrdering = 0;

    // Stop if all nodes are ordered
    while (numNodesAddedToOrdering < numNodes)
    {
        // Get the next node with the highest connectivity to the so-far ordered nodes (depending on the ordering type)
        // NB: When no node has been added to the ordering yet, we take the start node
        NodeID nodeID;
        if (numNodesAddedToOrdering == 0)
            nodeID = startNodeID;
        else if (useBinaryHeap)
        {
            nodeID = binaryMaxHeap.top();
            // Remove the node from the binary max heap
            binaryMaxHeap.pop();
        }
        else
            // We purposefully do not decrease the highest non-empty bucket index here,
            // because we do it at the end of the while iteration, which is after the non-ordered
            // pins of the edge have been elevated. Since the elevation itself might change the
            // highest non-empty bucket index, it is sometimes unnecessary to decrease it here.
            nodeID = bucketMaxQueue.topAndPop(false);

        // Add the node to the ordering (if it is not a nullptr)
        if (nodeOrdering != nullptr)
            (*nodeOrdering)[numNodesAddedToOrdering] = nodeID;
        // Increase the number of nodes added to the ordering
        numNodesAddedToOrdering++;

        // Go over all incident edges of the node
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            // Ignore disabled edges
            if (!hypergraph.edgeIsEnabled(edgeID))
                continue;

            // Increase for the edge the number of pins added to the ordering
            edgeNumPinsAdded[edgeID]++;
            // Check if this is the first pin of the edge added to the ordering
            if (edgeNumPinsAdded[edgeID] == 1)
            {
                // Add the edge to the head ordering (if it is not a nullptr)
                if (edgeHeadOrdering != nullptr)
                    (*edgeHeadOrdering)[numEdgesAddedToOrdering] = edgeID;
                // Increase the number of edges added to the ordering
                numEdgesAddedToOrdering++;
                // Set the head of the edge (if it is not a nullptr)
                if (edgeHead != nullptr)
                    (*edgeHead)[edgeID] = nodeID;
            }

            // Get the number of levels to elevate the non-ordered pins of the edge
            // NB: A value of 0 means that the non-ordered pins should not be elevated
            KeyType numLevelsToElevate = get_num_levels_to_elevate_non_ordered_pins_of_edge(edgeNumPinsAdded[edgeID],
                                                                                            hypergraph.edgeSize(edgeID),
                                                                                            hypergraph.edgeWeight(edgeID));

            if (numLevelsToElevate > 0)
                // Elevate the non-ordered pins of an edge in the priority queue
                elevate_non_ordered_pins_of_edge(hypergraph, edgeID, numLevelsToElevate);
        }

        // Decrease the highest non-empty bucket index if the corresponding bucket is empty (only if useBinaryHeap == false)
        if (!useBinaryHeap && numNodesAddedToOrdering < numNodes)
            bucketMaxQueue.decreaseHighestNonEmptyBucketIndex();
    }
}

// Get the number of levels to elevate the non-ordered pins of an edge in the priority queue
// NB: A value of 0 means that the non-ordered pins should not be elevated
template <typename Hypergraph, typename KeyType>
inline KeyType Orderer<Hypergraph, KeyType>::get_num_levels_to_elevate_non_ordered_pins_of_edge(const NodeIndex edgeNumOfPinsAdded, const NodeIndex edgeSize, const EdgeWeight edgeWeight)
{
    switch (orderingType)
    {
    case OrderingType::MA:
        return (edgeNumOfPinsAdded == 1) * edgeWeight;
    case OrderingType::TIGHT:
        return (edgeNumOfPinsAdded == edgeSize - 1) * edgeWeight;
    case OrderingType::QUEYRANNE:
        // NB: It can happen that the non-ordered pins of an edge must be directly elevated (2 * edgeWeight) levels,
        //     i.e. when the hyperedge has a size of 2 and the first pin is added to the ordering. In this case,
        //     the second pin must be elevated directly (2 * edgeWeight) levels. In other words, we simply omit the 1/2 ordering
        //     weight, which comes with the advantage that we have always integer levels.
        return ((edgeNumOfPinsAdded == 1) + (edgeNumOfPinsAdded == edgeSize - 1)) * edgeWeight;
#ifdef SMHM_PARALLEL
    case OrderingType::MIX_UNIFORM:
        return (orderingWeight * (edgeNumOfPinsAdded == 1) + (1 - orderingWeight) * (edgeNumOfPinsAdded == edgeSize - 1)) * edgeWeight;
#endif
    default:
        return 0;
    }
}

// Elevate the non-ordered pins of an edge in the bucket priority queue
template <typename Hypergraph, typename KeyType>
inline void Orderer<Hypergraph, KeyType>::elevate_non_ordered_pins_of_edge(const Hypergraph &hypergraph,
                                                                           const EdgeID edgeID,
                                                                           const KeyType numLevelsToElevate)
{
    // Go over all pins of the edge
    for (NodeID pinID : hypergraph.pins(edgeID))
    {
        // Ignore disabled pins
        if (!hypergraph.nodeIsEnabled(pinID))
            continue;

        // Increase the key of the pin in the priority queue
        if (useBinaryHeap)
        {
            if (binaryMaxHeap.contains(pinID))
                binaryMaxHeap.increaseKeyBy(pinID, numLevelsToElevate);
        }
        else
        {
            if (bucketMaxQueue.contains(pinID))
                bucketMaxQueue.increaseByKey(pinID, numLevelsToElevate);
        }
    }
}

// Reset the priority queue
template <typename Hypergraph, typename KeyType>
inline NodeID Orderer<Hypergraph, KeyType>::reset_priority_queue(const Hypergraph &hypergraph, const NodeIndex numNodes)
{
    if (useBinaryHeap)
        // Clear the binary max heap
        binaryMaxHeap.clear();
    else
    {
        // Get the maximum node degree
        EdgeIndex maxDegree = 0;
        for (NodeID nodeID : hypergraph.nodes())
            if (hypergraph.nodeIsEnabled(nodeID))
                maxDegree = std::max(maxDegree, hypergraph.nodeDegree(nodeID));

        // Reset/resize the bucket priority queue and reserve the space of the first bucket to hold all elements
        // NB1: Since the maxDegree may be different for contracted versions of the hypergraph, we need to resize the bucket priority queue entirely each time
        // NB2: The start node does not need to be inserted, as it will be handled at the very beginning
        bucketMaxQueue.reset(orderingType == OrderingType::QUEYRANNE ? (2 * maxDegree + 1) : (maxDegree + 1), numNodes - 1);
    }

    NodeIndex indexInFirstBucket = 0;
    NodeID startNodeID = std::numeric_limits<NodeID>::max();
    for (NodeIndex i = 0; i < initialNumNodes; i++)
    {
        // Get the node id in the permutation
        NodeID nodeID = permutation[i];

        // Ignore disabled nodes
        // NB: numNodes is exactly the number of non-disabled nodes
        // NB2: For static hypergraphs, the node IDs are remapped to the range [0, numNodes - 1], so we need to use a different check
        if ((worksOnStaticHypergraph && nodeID >= numNodes) || !hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Set the start node if it is not set yet (i.e. the start node is the first non-disabled node in the permutation)
        // NB: The start node is not added to the priority queue
        if (startNodeID == std::numeric_limits<NodeID>::max())
        {
            startNodeID = nodeID;
            if (!useBinaryHeap)
                // Mark the start node as being not in a bucket
                bucketMaxQueue.markElementAsNotInBucket(startNodeID);
            continue;
        }

        if (useBinaryHeap)
            binaryMaxHeap.push(nodeID, 0);
        else
        {
            bucketMaxQueue.insertInFirstBucket(nodeID, indexInFirstBucket);
            indexInFirstBucket++;
        }
    }
    return startNodeID;
}

#endif // end of SMHM_ORDERER_H