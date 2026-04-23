/******************************************************************************
 * hypercactus.cpp
 * *
 * Computes the hypercactus of a hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <cassert>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/utils/definitions.h"
#include "lib/io/mt_kahypar_io.h"
#include "lib/solvers/kernelizer.h"
#include "lib/utils/random.h"
#include "lib/utils/output.h"
#include "lib/decomposition/canonical_decomposition.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
// KaHIP headers
#include "kahip/timer.h"

// Output the hypercactus and the mapping from nodes in the original hypergraph to nodes in the hypercactus
double output_hypercactus(StaticHypergraph &hypercactus,
                          std::vector<ClusterID> &hypercactusMapping,
                          const char *outputFileName,
                          timer &t,
                          CutValue minEdgeCut,
                          bool allowOutputFloatEdgeWeights)
{
    NodeIndex hypercactusNumNodes = hypercactus.initialNumNodes() - hypercactus.numRemovedHypernodes();
    EdgeIndex hypercactusNumEdges = hypercactus.initialNumEdges() - hypercactus.numRemovedHyperedges();

    double hypercactusTime = t.elapsed();

    std::cout << "===================================================" << std::endl;
    std::cout << "################### HYPERCACTUS ###################" << std::endl;
    std::cout << "===================================================" << std::endl;
    // Print the number of edges in the hypercactus
    std::cout << "hypercactus_num_edges \t\t\t" << hypercactusNumEdges << std::endl;
    // Print the number of nodes in the hypercactus
    std::cout << "hypercactus_num_nodes \t\t\t" << hypercactusNumNodes << std::endl;
    // Print the time it took to create the hypercactus
    std::cout << "hypercactus_time \t\t\t" << hypercactusTime << std::endl;

    // Create a string stream buffer for the hypercactus
    std::ostringstream hypercactusBuffer;
    std::string hypercactusFileName = outputFileName;
    hypercactusFileName += ".hypercactus";
    // Open output file of hypercactus
    std::ofstream fHypercactus(hypercactusFileName);

    // Iterate over all hyperedges
    for (EdgeID edgeID : hypercactus.edges())
    {
        // Ignore disabled hyperedges
        if (!hypercactus.edgeIsEnabled(edgeID))
            continue;

        EdgeWeight edgeWeight = hypercactus.edgeWeight(edgeID);
        // If we allow floating point edge weights in the output, we need to correct the potentially wrong integer rounding that
        // has happened during the construction of the hypercactus.
        if (allowOutputFloatEdgeWeights && minEdgeCut % 2 == 1 && edgeWeight < minEdgeCut)
            hypercactusBuffer << minEdgeCut / 2.0 << " ";
        else
            hypercactusBuffer << edgeWeight << " ";

        // Iterate over the pins of the hyperedge
        for (NodeID pinID : hypercactus.pins(edgeID))
        {
            // Ignore disabled pins
            if (!hypercactus.nodeIsEnabled(pinID))
                continue;

            hypercactusBuffer << pinID + 1 << " ";
        }
        hypercactusBuffer << std::endl;
    }

    // Iterate over all hypernodes
    for (NodeID nodeID : hypercactus.nodes())
    {
        // Ignore disabled nodes
        if (!hypercactus.nodeIsEnabled(nodeID))
            continue;

        hypercactusBuffer << hypercactus.nodeWeight(nodeID) << std::endl;
    }

    std::cout << "hypercactus_output_file \t\t" << hypercactusFileName << std::endl;

    // Write the header to the file
    fHypercactus << hypercactusNumEdges << " " << hypercactusNumNodes << " 11" << std::endl;
    // Write the hyperedges to the file
    fHypercactus << hypercactusBuffer.str();
    // Close the file
    fHypercactus.close();

    // Create a string stream buffer for the hypercactus mapping
    std::ostringstream hypercactusMappingBuffer;
    std::string hypercactusMappingFileName = outputFileName;
    hypercactusMappingFileName += ".hypercactus.mapping";
    // Open output file of hypercactus mapping
    std::ofstream fHypercactusMapping(hypercactusMappingFileName);

    // Iterate over the mapping
    for (NodeID nodeID : hypercactusMapping)
        hypercactusMappingBuffer << nodeID + 1 << std::endl;

    std::cout << "hypercactus_mapping_output_file \t" << hypercactusMappingFileName << std::endl;

    // Write the mapping to the file
    fHypercactusMapping << hypercactusMappingBuffer.str();
    // Close the file
    fHypercactusMapping.close();

    return hypercactusTime;
}

int main(int argc, char *argv[])
{
#ifdef NDEBUG
    std::cout << "build_mode \t\t\t\tRELEASE" << std::endl;
#else
    std::cout << "build_mode \t\t\t\tDEBUG" << std::endl;
#endif

    // Stores the config
    HypercactusGeneratorConfig config;

    int parseReturnCode = ParseParams::parse_parameters_hypercactus(argc, argv, config);

    if (parseReturnCode)
        return 1;

    // Initialize the timer
    timer t;

    // Set the seed for the random number generator
    RandomFunctions::set_seed(config.seed);

    // Read the hypergraph from the given file
    mt_kahypar_hypergraph_t hypergraphWrapper = MtKaHyParIO::read_hypergraph_from_file(config);

    // Cast the hypergraph wrapper to a static hypergraph
    StaticHypergraph &hypergraph = mt_kahypar::utils::cast<StaticHypergraph>(hypergraphWrapper);

    // Extract the number of nodes and edges
    EdgeIndex inputNumEdges = hypergraph.initialNumEdges() - hypergraph.numRemovedHyperedges();
    NodeIndex inputNumNodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();

    if (config.unweighted)
    {
        for (EdgeID edgeID : hypergraph.edges())
            if (hypergraph.edgeIsEnabled(edgeID))
                hypergraph.setEdgeWeight(edgeID, 1);
    }

    std::cout << "seed \t\t\t\t\t" << config.seed << std::endl;
    std::cout << "unweighted \t\t\t\t" << config.unweighted << std::endl;
    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;
    // Initialize the total time
    double totalComputingTime = 0;

    if (config.verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "############### MINCUT COMPUTATION ################" << std::endl;
    }

    // Initialize the kernelizer config
    KernelizerConfig kernelizerConfig;
    kernelizerConfig.seed = config.seed;
    kernelizerConfig.numThreads = config.numThreads;
    kernelizerConfig.baseSolver = DEFAULT_BASE_SOLVER;
    kernelizerConfig.ilpTimeout = DEFAULT_TIMEOUT_ILP;
    kernelizerConfig.ilpMode = DEFAULT_ILP_MODE;
    kernelizerConfig.orderingType = DEFAULT_ORDERING_TYPE;
    kernelizerConfig.orderingMode = DEFAULT_ORDERING_MODE;
    kernelizerConfig.pruningMode = DEFAULT_PRUNING_MODE;
    // Deactivate LP to get an exact mincut value
    kernelizerConfig.LPNumIterations = 0;
    kernelizerConfig.LPMode = DEFAULT_LP_MODE;
    kernelizerConfig.LPNumPinsToSample = DEFAULT_LP_NUM_PINS_TO_SAMPLE;
    // Set the verbose flag
    kernelizerConfig.verbose = config.verbose;

    // Initialize the kernelizer
    Kernelizer kernelizer(kernelizerConfig);

    // Compute the minimum cut of the hypergraph using kernelization
    // NB: We need to copy the hypergraph, since the kernelizer modifies/contracts it
    StaticHypergraph mincutHypergraph = hypergraph.copy();
    KernelizerResult result = kernelizer.compute_mincut(mincutHypergraph, inputNumNodes, inputNumEdges);

    std::cout << "===================================================" << std::endl;
    // Print the final minimum edge cut value
    std::cout << "mincut_value \t\t\t\t" << result.minEdgeCut << std::endl;
    // Print the mincut computing time
    std::cout << "mincut_computing_time \t\t\t" << result.time << std::endl;
    totalComputingTime += result.time;

    if (config.verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "################## KERNELIZATION ##################" << std::endl;
        std::cout << "===================================================" << std::endl;
    }

    double kernelizationTime = 0;
    // Initialize the kernelization mapping to contain 0, 1, 2, ..., inputNumNodes-1
    std::vector<ClusterID> kernelizationMapping(inputNumNodes);
    std::iota(kernelizationMapping.begin(), kernelizationMapping.end(), 0);
    // NB: We need to copy the hypergraph, since the kernelizer modifies/contracts it
    StaticHypergraph kernelizedHypergraph = hypergraph.copy();
    bool didStopEarly = kernelizer.apply_kernelization(kernelizedHypergraph, inputNumNodes, inputNumEdges, kernelizationTime, true, &kernelizationMapping);
    // Extract the number of nodes and edges of the kernelized hypergraph
    EdgeIndex kernelNumEdges = kernelizedHypergraph.initialNumEdges();
    NodeIndex kernelNumNodes = kernelizedHypergraph.initialNumNodes();

    // Print the results of the kernelization
    if (!config.verbose || didStopEarly)
        std::cout << "===================================================" << std::endl;
    // Print the number of edges after the kernelization
    std::cout << "kernel_num_edges  \t\t\t" << kernelNumEdges << std::endl;
    // Print the number of nides after the kernelization
    std::cout << "kernel_num_nodes  \t\t\t" << kernelNumNodes << std::endl;
    // Print the kernelization time
    std::cout << "kernelization_time \t\t\t" << kernelizationTime << std::endl;
    totalComputingTime += kernelizationTime;

    // Print the mapping
    if (config.verbose)
    {
        std::cout << "kernelization_mapping \t\t\t";
        for (ClusterID i = 0; i < kernelizationMapping.size(); ++i)
            std::cout << kernelizationMapping[i] << " ";
        std::cout << std::endl;
    }

    // Restart timer
    t.restart();

    // Make sure that if the hypergraph is disconnected (i.e. minimum edge cut is 0), the kernelized hypergraph has no edges
    assert(result.minEdgeCut > 0 || kernelNumEdges == 0);

    // Stop early if the minimum edge cut is 0, since then the hypercactus is simply the kernelized hypergraph
    if (result.minEdgeCut == 0)
    {
        // We need to set the weight of the hypercactus nodes to the number of nodes they represent in the original hypergraph
        for (NodeID nodeID : kernelizedHypergraph.nodes())
            kernelizedHypergraph.setNodeWeight(nodeID, 0);
        for (NodeID nodeID : kernelizationMapping)
            kernelizedHypergraph.setNodeWeight(nodeID, kernelizedHypergraph.nodeWeight(nodeID) + 1);
        std::cout << "===================================================" << std::endl;
        std::cout << "Minimum edge cut is 0, hypercactus is simply the kernelized hypergraph." << std::endl;
        totalComputingTime += output_hypercactus(kernelizedHypergraph,
                                                 kernelizationMapping,
                                                 config.outputFileName,
                                                 t,
                                                 result.minEdgeCut,
                                                 config.allowOutputFloatEdgeWeights);
        std::cout << "===================================================" << std::endl;
        std::cout << "total_computing_time \t\t\t" << totalComputingTime << std::endl;
        return 0;
    }

    // Initialize the canonical decomposition of the hypergraph
    CanonicalDecomposition canonicalDecomposition(kernelNumNodes, kernelNumEdges, kernelizationMapping, result.minEdgeCut, config.verbose);
    // First compute the prime decomposition of the kernelized hypergraph
    canonicalDecomposition.compute_prime_decomposition_of_kernel(std::move(kernelizedHypergraph));
    if (config.verbose)
        canonicalDecomposition.print_decomposition_tree();
    // Convert the prime decomposition of the kernelized hypergraph into the prime decomposition of the original hypergraph
    canonicalDecomposition.convert_to_prime_decomposition_of_original_hypergraph();
    if (config.verbose)
        canonicalDecomposition.print_decomposition_tree();
    // Compute the canonical decomposition from the prime decomposition of the original hypergraph
    canonicalDecomposition.compute_canonical_decomposition_from_prime_decomposition(hypergraph);
    canonicalDecomposition.print_decomposition_tree(true);

    // Create the hypercactus from the canonical decomposition
    // NB: The weight of a node in the hypercactus indicates how many nodes from the original hypergraph it represents
    std::vector<ClusterID> hypercactusMapping(inputNumNodes);
    mt_kahypar_hypergraph_t hypercactusWrapper = canonicalDecomposition.create_hypercactus(hypergraph, hypercactusMapping, config.allowOutputFloatEdgeWeights);

    StaticHypergraph &hypercactus = mt_kahypar::utils::cast<StaticHypergraph>(hypercactusWrapper);

    // Output the hypercactus and the mapping from nodes in the original hypergraph to nodes in the hypercactus
    totalComputingTime += output_hypercactus(hypercactus,
                                             hypercactusMapping,
                                             config.outputFileName,
                                             t,
                                             result.minEdgeCut,
                                             config.allowOutputFloatEdgeWeights);
    std::cout << "===================================================" << std::endl;
    std::cout << "total_computing_time \t\t\t" << totalComputingTime << std::endl;

    // Make sure that the minimum edge cut of the hypercactus is equal to the minimum edge cut of the original hypergraph
    // NB: If we allow floating point edge weights in the output, we cannot guarantee that the mincut of the hypercactus is equal to the mincut of the original hypergraph (due inter edge weights of Mt-KaHyPar)
    assert([&]()
           {
               kernelizerConfig.verbose = false;
               Kernelizer silentKernelizer(kernelizerConfig);
               return config.allowOutputFloatEdgeWeights || silentKernelizer.compute_mincut(hypercactus, hypercactus.initialNumNodes(), hypercactus.initialNumEdges()).minEdgeCut == result.minEdgeCut; }());

    // Free the hypercactus wrapper
    mt_kahypar_free_hypergraph(hypercactusWrapper);
    // Free the hypergraph wrapper
    mt_kahypar_free_hypergraph(hypergraphWrapper);

    return 0;
}
