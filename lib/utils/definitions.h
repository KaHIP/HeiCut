/******************************************************************************
 * definitions.h
 * *
 * General definitions for the algorithms.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_DEFINITIONS_H
#define SMHM_DEFINITIONS_H

// Mt-KaHyPar header
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/dynamic_graph.h"

// Makros
#ifdef SMHM_PARALLEL
#define forAllEdgesSequentialOrParallel(hypergraph, edgeID) { hypergraph.doParallelForAllEdges([&](EdgeID edgeID)
#define forAllNodesSequentialOrParallel(hypergraph, nodeID) { hypergraph.doParallelForAllNodes([&](NodeID nodeID)
#define forRangeSequentialOrParallel(startValue, endValue, i, Index) { tbb::parallel_for(static_cast<Index>(startValue), endValue, [&](Index i)
#define forBlockedRangeSequentialOrParallel(startValue, endValue, i, Index) { tbb::parallel_for(tbb::blocked_range<Index>(static_cast<Index>(startValue), endValue), \
                                                                              [&](const tbb::blocked_range<Index> &range) { for (Index i = range.begin(); i != range.end(); ++i)
#define continueFor return
#define endFor );}
#define endForBlocked });}
#else
#define forAllEdgesSequentialOrParallel(hypergraph, edgeID) { for (EdgeID edgeID : hypergraph.edges())
#define forAllNodesSequentialOrParallel(hypergraph, nodeID) { for (NodeID nodeID : hypergraph.nodes())
#define forRangeSequentialOrParallel(startValue, endValue, i, Index) { for (Index i = startValue; i < endValue; i++)
#define forBlockedRangeSequentialOrParallel(startValue, endValue, i, Index) { for (Index i = startValue; i < endValue; i++)
#define continueFor continue
#define endFor }
#define endForBlocked }
#endif

// Type of the ordering
enum class OrderingType : uint8_t
{
#ifdef SMHM_PARALLEL
    MIX_DISCRETE, // Mix orderings in a discrete way (i.e. different threads get different orderings)
    MIX_UNIFORM,  // Mix orderings in a uniform way (i.e. all threads get different, unformly spaced ordering weights)
#endif
    MA,       // Klimmek and Wagner
    TIGHT,    // Mak and Wong
    QUEYRANNE // Queyranne
};

// Mode of the contraction performed with the ordering
enum class OrderingMode : uint8_t
{
    SINGLE,
    MULTI
};

// Mode of the label propagation algorithm
enum class LabelPropagationMode : uint8_t
{
    CLIQUE_EXPANDED,
    PROBABILISTIC
};

// Underlying exact solver of the kernelizer
enum class BaseSolver : uint8_t
{
    ILP,
    SUBMODULAR
};

// Mode of the pruning rules
enum class PruningMode : uint8_t
{
    BEST,
    ALL,
};

// Mode of the ILP solver
enum class ILPMode : uint8_t
{
    BIP,
    MILP
};

// Type definitions
typedef uint32_t NodeID;
typedef uint32_t NodeIndex;
typedef uint32_t EdgeID;
typedef uint32_t EdgeIndex;
typedef uint64_t NodeWeight;
typedef uint64_t EdgeWeight;
typedef uint32_t ClusterID;
typedef uint32_t ClusterIndex;
typedef uint32_t CutValue;
typedef uint32_t FlowValue;
typedef uint32_t TrimmerValue;
typedef uint32_t IterationIndex;
typedef size_t Fingerprint;
typedef double ScoreValue;
typedef mt_kahypar::ds::StaticHypergraph StaticHypergraph;
typedef mt_kahypar::ds::DynamicHypergraph DynamicHypergraph;
typedef mt_kahypar::ds::StaticGraph StaticGraph;
typedef mt_kahypar::ds::DynamicGraph DynamicGraph;

#endif // end of SMHM_DEFINITIONS_H