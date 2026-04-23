/******************************************************************************
 * submodular.cpp
 * *
 * Solver that finds the minimum cut of a hypergraph via submodular optimization.
 * More specifically, it implements the algorithms of Klimmek and Wagner, Mak
 * and Wong and Queyranne. Note that one instance of the SubmodularMincut class
 * can only be used for one hypergraph (and its contracted versions).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <vector>
// Own headers
#include "submodular.h"
#include "lib/utils/definitions.h"
#include "lib/orderer/orderer.h"
#include "lib/utils/random.h"
// Mt-KaHyPar headers
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"

SubmodularMincut::SubmodularMincut(const NodeIndex initialNumNodes,
                                   const EdgeIndex initialNumEdges,
                                   const OrderingType orderingType,
                                   const OrderingMode orderingMode,
                                   const bool hasWeightedEdges,
                                   const size_t numThreads)
#ifdef SMHM_PARALLEL
    : numThreads(numThreads),
      orderingMode(orderingMode),
      nodeOrderingPerThread(numThreads, std::vector<NodeID>(initialNumNodes)),
      markedEdgePinsPerThread(numThreads, std::vector<MarkedEdgePins>(initialNumEdges, {0, std::numeric_limits<NodeID>::max(), false, false})),
      ordererPerThread(numThreads),
      orderingTypePerThread(numThreads)
{
    // Create a list of all possible ordering types (only if the ordering type is ALL)
    constexpr OrderingType discreteOrderingTypes[] = {OrderingType::TIGHT, OrderingType::QUEYRANNE, OrderingType::MA};

    auto initializeThread = [&](size_t threadIndex)
    {
        MersenneTwister threadRandEngine(RandomFunctions::get_seed() + threadIndex);
        OrderingType threadOrderingType = (orderingType == OrderingType::MIX_DISCRETE) ? discreteOrderingTypes[threadIndex % 3] : orderingType;
        // IMPORTANT: Even if we perform pairwise contractions on the dynamic hypergraph, we can still guarantee that if the original hypergraph is unweighted,
        //            every contracted version is also unweighted. This is because contrary to the static hypergraph, the dynamic hypergraph does not
        //            merge parallel edges or delete single-pin edges during contractions. Thus, the unweighted node degree is identical to the weighted node degree
        //            and there exists as many parallel hyperedges (with weight 1) as the weight of the merged hyperedge would be.
        if (threadOrderingType == OrderingType::MIX_UNIFORM)
            ordererPerThread[threadIndex] = std::make_unique<Orderer<DynamicHypergraph, double>>(initialNumNodes, initialNumEdges, threadOrderingType, hasWeightedEdges, threadRandEngine, (double)threadIndex / (numThreads - 1));
        else
            ordererPerThread[threadIndex] = std::make_unique<Orderer<DynamicHypergraph, EdgeWeight>>(initialNumNodes, initialNumEdges, threadOrderingType, hasWeightedEdges, threadRandEngine);
        orderingTypePerThread[threadIndex] = threadOrderingType;
    };

    // Create a vector to hold the threads that will initialize the orderers
    std::vector<std::thread> initializingThreads;

    // Initialize the orderers for each thread in parallel
    for (size_t threadIndex = 1; threadIndex < numThreads; ++threadIndex)
        initializingThreads.emplace_back([&, threadIndex]()
                                         { initializeThread(threadIndex); });

    // Also use the main thread
    initializeThread(0);

    // Wait for all threads to finish
    for (auto &t : initializingThreads)
        t.join();
};
#else
    : orderingMode(orderingMode),
      nodeOrdering(initialNumNodes),
      markedEdgePins(initialNumEdges, {0, std::numeric_limits<NodeID>::max(), false, false}),
      // IMPORTANT: Even if we perform pairwise contractions on the dynamic hypergraph, we can still guarantee that if the original hypergraph is unweighted,
      //            every contracted version is also unweighted. This is because contrary to the static hypergraph, the dynamic hypergraph does not
      //            merge parallel edges or delete single-pin edges during contractions. Thus, the unweighted node degree is identical to the weighted node degree
      //            and there exists as many parallel hyperedges (with weight 1) as the weight of the merged hyperedge would be.
      orderer(initialNumNodes, initialNumEdges, orderingType, hasWeightedEdges, RandomFunctions::get_random_engine()),
      orderingType(orderingType) {};
#endif