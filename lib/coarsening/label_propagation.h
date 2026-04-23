/******************************************************************************
 * label_propagation.h
 * *
 * Label propagation algorithm to coarsen a hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_LABELPROPAGATION_H
#define SMHM_LABELPROPAGATION_H

#include <cassert>
#include <vector>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
#include "lib/utils/random.h"
// Mt-KaHyPar headers
#include "mt-kahypar/parallel/stl/scalable_vector.h"

class LabelPropagation
{
private:
    // Number of iterations to perform the label propagation
    const IterationIndex numIterations = DEFAULT_LP_NUM_ITERATIONS;

    // Mode of the label propagation algorithm
    const LabelPropagationMode mode = DEFAULT_LP_MODE;

    // Number of pins to sample (only used if mode is probabilistic)
    const NodeIndex numPinsToSample = DEFAULT_LP_NUM_PINS_TO_SAMPLE;

public:
    LabelPropagation(const IterationIndex numIterations, const LabelPropagationMode mode, const NodeIndex numPinsToSample);

    // Perform the label propagation and the contraction of the label clusters
    StaticHypergraph propagate_and_contract_labels(StaticHypergraph &hypergraph, mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID);

    // Get the score of an edge (which contributes to the score between a node and a cluster)
    ScoreValue get_score_of_edge(const StaticHypergraph &hypergraph, const EdgeID edgeID);
};

// Get the score of an edge (which contributes to the score between a node and a cluster)
inline ScoreValue LabelPropagation::get_score_of_edge(const StaticHypergraph &hypergraph, const EdgeID edgeID)
{
    ScoreValue score = static_cast<ScoreValue>(hypergraph.edgeWeight(edgeID)) / (hypergraph.edgeSize(edgeID) - 1);

    // Correct the score for the sampling of the pins (only used if mode is probabilistic)
    if (mode == LabelPropagationMode::PROBABILISTIC)
        score = (score * hypergraph.edgeSize(edgeID)) / numPinsToSample;

    return score;
}

// Perform the label propagation and the contraction of the label clusters
inline StaticHypergraph LabelPropagation::propagate_and_contract_labels(StaticHypergraph &hypergraph, mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID)
{
    // Get the initial number of nodes and the initial number of edges of the hypergraph, since we use the node/edge id indices.
    // Besides, for static hypergraphs the length of the cluster id vector must be equal to the initial number of nodes.
    NodeIndex initialNumNodes = hypergraph.initialNumNodes();
    NodeIndex initialNumEdges = hypergraph.initialNumEdges();

    // Initialize the cluster ids (= labels) of the nodes to contain 0, 1, 2, ..., n-1
    // NB: We purposefully ignore race conditions here in the parallel version
    assert(clusterID.size() == initialNumNodes);
    std::iota(clusterID.begin(), clusterID.end(), 0);

    // Initialize the random permutation of the nodes
    std::vector<NodeID> permutation(initialNumNodes);
    RandomFunctions::permutate_vector_local(&permutation, true);

    // Initialize the incidence array to store the pins of the edges (only used if mode is probabilistic)
    std::vector<std::vector<NodeID>> incidenceArray(mode == LabelPropagationMode::PROBABILISTIC ? initialNumEdges : 0);

    if (mode == LabelPropagationMode::PROBABILISTIC)
    {
        // Fill the incidence array with the pins of the edges
        forRangeSequentialOrParallel(0, initialNumEdges, i, EdgeIndex)
        {
            // Ignore disabled hyperedges
            if (!hypergraph.edgeIsEnabled(i))
                continueFor;

            incidenceArray[i].reserve(hypergraph.edgeSize(i));
            for (NodeID pinID : hypergraph.pins(i))
                // NB: It it important to add here all pins of the edge (also if they are disabled) because
                //     we will use the size of the edge to iterate, which also includes the disabled pins.
                incidenceArray[i].push_back(pinID);
        }
        endFor;
    }

    // Intitialize the hash vector to store the scores (= connection strengths of the current node to the clusters).
    // NB: The first element of the pair is the score of the cluster (i.e. its id is equal to index), the second element
    //     is the current node id for which the score is stored. This allows to use the same vector for all nodes, since
    //     we exactly know when to to reset the values in the vector.
    //
    // For the parallel version also initialize the random number generator per thread (since MersenneTwister is not thread safe)
    std::pair<ScoreValue, NodeID> hashVectorInitialValue = {0, 0};
#ifdef SMHM_PARALLEL
    tbb::enumerable_thread_specific<std::vector<std::pair<ScoreValue, NodeID>>> hashVectorPerThread(
        [&]
        { return std::vector<std::pair<ScoreValue, NodeID>>(initialNumNodes, hashVectorInitialValue); });

    std::atomic<int> counter{0};
    tbb::enumerable_thread_specific<MersenneTwister> randEnginePerThread(
        [&]()
        {
            return MersenneTwister(RandomFunctions::get_seed() + counter++);
        });
#else
    std::vector<std::pair<ScoreValue, NodeID>> hashVector(initialNumNodes, hashVectorInitialValue);
#endif

    for (IterationIndex it = 0; it < numIterations; it++)
    {
        // Sample a constant number of pins of every edges (only used if mode is probabilistic)
        // The sampled pins of an edge are moved to the front of the respective range in the incidence array
        // NB: We do not have to reset the incidence array for the next iteration, since we will again move the sampled pins simply to the front
        if (mode == LabelPropagationMode::PROBABILISTIC)
        {
            forRangeSequentialOrParallel(0, initialNumEdges, i, EdgeIndex)
            {
                // Do nothing if the edge is disabled or if the size of the edge is smaller than the number of pins to sample
                if (!hypergraph.edgeIsEnabled(i) || hypergraph.edgeSize(i) <= numPinsToSample)
                    continueFor;

#ifdef SMHM_PARALLEL
                // Get the random number generator for the current thread
                MersenneTwister &randEngine = randEnginePerThread.local();
#else
                // Get the random number generator
                MersenneTwister &randEngine = RandomFunctions::get_random_engine();
#endif

                // Sample numPinsToSample pins of the edge by moving them to the front in the incidence array
                for (NodeIndex j = 0; j < numPinsToSample; j++)
                {
                    // Get a random node index and swap it to the the position j
                    NodeIndex randomNodeIndex = RandomFunctions::get_uniform_random_int_in_bounds<NodeIndex>(j, hypergraph.edgeSize(i) - 1, randEngine);
                    if (randomNodeIndex != j)
                        std::swap(incidenceArray[i][j], incidenceArray[i][randomNodeIndex]);
                }
            }
            endFor;
        }

        // Iterate over all nodes in the fixed random order using the permutation vector
        forRangeSequentialOrParallel(0, initialNumNodes, i, NodeIndex)
        {

#ifdef SMHM_PARALLEL
            // Get the hash vector for the current thread
            std::vector<std::pair<ScoreValue, NodeID>> &hashVector = hashVectorPerThread.local();
            // Get the random number generator for the current thread
            MersenneTwister &randEngine = randEnginePerThread.local();
#else
            // Get the random number generator
            MersenneTwister &randEngine = RandomFunctions::get_random_engine();
#endif

            // Get the node id in the permutation
            NodeID nodeID = permutation[i];

            // Ignore the node if it is disabled
            if (!hypergraph.nodeIsEnabled(nodeID))
                continueFor;

            // Initialize the cluster id with the highest score (and its respective score)
            ClusterID maxClusterID = clusterID[nodeID];
            ScoreValue maxClusterScore = 0;
            // Go over the edges of the node
            for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
            {
                // Ignore disabled hyperedges or hyperedges with only one pin
                if (!hypergraph.edgeIsEnabled(edgeID) || hypergraph.edgeSize(edgeID) <= 1)
                    continue;

                // Compute the score of the edge
                ScoreValue scoreOfEdge = get_score_of_edge(hypergraph, edgeID);

                // Get the iterator range of the pins of the edge (only used if mode is clique expanded)
                auto pinsIterator = hypergraph.pins(edgeID).begin();

                // Go over the pins of the edge
                for (NodeIndex i = 0; i < hypergraph.edgeSize(edgeID); i++)
                {
                    // Stop if we have already looked at all sampled pins of the edge (only used if mode is probabilistic)
                    if (mode == LabelPropagationMode::PROBABILISTIC && i >= numPinsToSample)
                        break;

                    // Get the pin id of the edge based on the mode
                    NodeID pinID = (mode == LabelPropagationMode::PROBABILISTIC) ? incidenceArray[edgeID][i] : *(pinsIterator++);

                    // Ignore the pin if it is disabled or if it is the node itself
                    if (!hypergraph.nodeIsEnabled(pinID) || pinID == nodeID)
                        continue;

                    // Get the cluster id that contains the pin
                    ClusterID pinClusterID = clusterID[pinID];
                    // Reset the score if we look for the first time
                    // at the connection beween nodeID and pinClusterID
                    if (hashVector[pinClusterID].second != nodeID)
                    {
                        hashVector[pinClusterID].first = 0;
                        hashVector[pinClusterID].second = nodeID;
                    }

                    // Increase the score between nodeID and pinClusterID
                    hashVector[pinClusterID].first += scoreOfEdge;
                    // Check if the score is the highest computed so far
                    if (hashVector[pinClusterID].first > maxClusterScore ||
                        (hashVector[pinClusterID].first == maxClusterScore &&
                         RandomFunctions::get_uniform_random_int(randEngine) % 2))
                    {
                        maxClusterScore = hashVector[pinClusterID].first;
                        maxClusterID = pinClusterID;
                    }
                }
            }
            // Update the cluster of the node
            clusterID[nodeID] = maxClusterID;
        }
        endFor;
    }

    // Contract the hypergraph based on the cluster ids
    // NB: We do not have to remap the cluster ids since this is already done during the contraction
    // NB2: We do not have to free the hypergraph since we optain directly the hypergraph and not a wrapper
    // NB3: Single-pin or parallel hyperedges are removed from the contracted hypergraph automatically
    return hypergraph.contract(clusterID);
}

#endif // end of SMHM_LABELPROPAGATION_H