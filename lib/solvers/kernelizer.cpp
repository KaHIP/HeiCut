/******************************************************************************
 * kernelizer.cpp
 * *
 * Computes the (estimated) mincut of a hypergraph. For this, the kernel of the
 * hypergraph is built by coarsening the hypergraph using the label propagation
 * and pruning rules. Afterwards, the mincut of the coarsened hypergraph is
 * computed using a submodular solver or an ILP (with floating point relaxation).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
#include <functional>
#include <vector>
#include <cassert>
// Own headers
#include "kernelizer.h"
#include "lib/utils/definitions.h"
#include "lib/coarsening/label_propagation.h"
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/coarsening/pruner.h"
#ifdef USE_GUROBI
#include "lib/solvers/ilp.h"
#endif
#include "lib/solvers/submodular.h"
#include "lib/conversion/hypergraph_converter.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/utils/cast.h"
#ifdef USE_GUROBI
// Gurobi headers
#include "gurobi_c++.h"
#endif
// KaHIP headers
#include "kahip/timer.h"

struct ContractionMethod
{
    std::function<void()> execute;
    std::string suffix;
    std::string timeDescription;
};

Kernelizer::Kernelizer(KernelizerConfig config) : config(config) {}

// Apply kernelization to the hypergraph and return whether we can stop early (i.e. mincut can already be determined)
bool Kernelizer::apply_kernelization(StaticHypergraph &hypergraph,
                                     const NodeIndex initialNumNodes,
                                     const EdgeIndex initialNumEdges,
                                     double &totalComputingTime,
                                     const bool findAllMincuts,
                                     std::vector<ClusterID> *mapping)
{
    // Remove hyperedges of size one or with weight zero
    t.restart();
    pruner.remove_hyperedges_of_size_one_or_weight_zero(hypergraph);
    totalComputingTime += print_stats_and_return_used_time(hypergraph, initialNumNodes, initialNumEdges, "_no_trivial_edges \t\t", "trivial_edges_pruning_time \t\t");

    // Initialize the label propagater
    LabelPropagation labelPropagater(config.LPNumIterations, config.LPMode, config.LPNumPinsToSample);

    // Initialize the round counter
    IterationIndex roundCounter = 0;
    // Store the number of nodes before the pruning rules in each round (to check if the pruning rules have reduced the hypergraph)
    NodeIndex numNodesBeforePruningRules = 0;
    // Define the pruning rules applied during the while loop
    std::vector<ContractionMethod> contractionMethods;

    // Add the best pruning rules to the list of contraction methods
    contractionMethods.insert(contractionMethods.begin(),
                              {// Contract hyperedges that are not lighter than the given estimate (VieCut pruning rule 1)
                               {[&]()
                                {
                                    // Store the number of nodes before the pruning rules in each round (to check if the pruning rules have reduced the hypergraph)
                                    numNodesBeforePruningRules = get_current_num_nodes(hypergraph);
                                    hypergraph = pruner.contract_hyperedges_not_lighter_than_estimate(hypergraph, clusterID, minEdgeCut, findAllMincuts);
                                },
                                "_no_heavy_edges \t\t", "heavy_edges_pruning_time \t\t"},
                               // Contract overlaps of hyperedges that are not lighter than the given estimate
                               {[&]()
                                { hypergraph = pruner.contract_overlaps_not_lighter_than_estimate(hypergraph, clusterID, minEdgeCut, findAllMincuts); },
                                "_no_heavy_overlaps \t", "heavy_overlaps_pruning_time \t\t"},
                               // Contract hyperedges of size two that can be shifted to one side of the cut (VieCut pruning rule 2)
                               {[&]()
                                { hypergraph = pruner.contract_shiftable_hyperedges_of_size_two(hypergraph, clusterID, minEdgeCut, findAllMincuts); },
                                "_no_shiftable_2-edges \t", "shiftable_2-edges_pruning_time \t\t"},
                               // Contract hyperedges of size two that meet the triangle conditions (VieCut pruning rules 3 and 4)
                               {[&]()
                                { hypergraph = pruner.contract_triangle_hyperedges_of_size_two(hypergraph, clusterID, minEdgeCut, findAllMincuts); },
                                "_no_triangle_2-edges \t", "triangle_2-edges_pruning_time \t\t"}});

    // Add the remaining pruning rules to the list of contraction methods if all pruning rules should be applied
    if (config.pruningMode == PruningMode::ALL)
        contractionMethods.insert(contractionMethods.end() - 2,
                                  // Contract strictly nested isolated substructures of parent hyperedges
                                  {[&]()
                                   { hypergraph = pruner.contract_strictly_nested_isolated_substructures(hypergraph, clusterID); },
                                   "_no_nested_subgraphs \t", "nested_subgraphs_pruning_time \t\t"});

    // Perform the label propagation (if necessary)
    if (config.LPNumIterations > 0)
        contractionMethods.insert(contractionMethods.begin(),
                                  {[&]()
                                   { hypergraph = labelPropagater.propagate_and_contract_labels(hypergraph, clusterID); },
                                   "_coarsened \t\t", "coarsening_time \t\t\t"});

    // Apply the pruning rules (and optionally the label propagation) until the number of nodes does not change anymore
    // NB: For the first round we do not check the number of nodes, since we want to perform at least one round
    while (roundCounter == 0 || get_current_num_nodes(hypergraph) < numNodesBeforePruningRules)
    {
        // Print the current round
        if (config.verbose)
        {
            std::cout << "<<<<<<<<<<<<<<<<<<<<< Round " << roundCounter << " >>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cout << "===================================================" << std::endl;
        }
        // Increment the round counter
        ++roundCounter;

        // Go over all contraction methods and apply them to the hypergraph
        for (const ContractionMethod &method : contractionMethods)
        {
            // Restart timer for measuring the time of the current contraction method
            t.restart();
            // Resize the clusterID vector to the initial number of nodes of the hypergraph
            clusterID.resize(hypergraph.initialNumNodes());
            // Perform the contraction method
            method.execute();
            // Update the mapping (if necessary)
            update_mapping(mapping, clusterID);
            // Update the mincut value
            update_min_cut_value(pruner, hypergraph);
            // Print the statistics of the current hypergraph
            totalComputingTime += print_stats_and_return_used_time(hypergraph, initialNumNodes, initialNumEdges, method.suffix, method.timeDescription);
            // Check if we can stop early
            if (can_stop_early(hypergraph, findAllMincuts))
                return true;
        }
    }
    return false;
}

// Compute the minimum cut of the hypergraph using kernelization
KernelizerResult Kernelizer::compute_mincut(StaticHypergraph &hypergraph, const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges)
{
    // Initialize the total computing time
    double totalComputingTime = 0;

    // Get the naive estimate of the minimum cut and the weighted node degrees
    t.restart();
    CutValue naiveEstimate = pruner.compute_naive_mincut_estimate(hypergraph);
    // Initialize the minimum edge cut value as the naive estimate
    minEdgeCut = naiveEstimate;
    double naiveEstimateTime = t.elapsed();
    totalComputingTime += naiveEstimateTime;
    if (config.verbose)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<< Preprocessing >>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "naive_mincut_value \t\t\t" << minEdgeCut << std::endl;
        std::cout << "naive_mincut_time \t\t\t" << naiveEstimateTime << std::endl;
        std::cout << "===================================================" << std::endl;
    }

    // Stop early if the naive estimate is already 0
    if (naiveEstimate == 0)
    {
        if (config.verbose)
            std::cout << "Naive mincut is already 0, no base solver necessary." << std::endl;
        return {0, 0, totalComputingTime};
    }

    // Apply the kernelization to the hypergraph
    bool canStopEarly = apply_kernelization(hypergraph, initialNumNodes, initialNumEdges, totalComputingTime);

    // Check if we can stop early
    if (canStopEarly)
        return {naiveEstimate, minEdgeCut, totalComputingTime};

    if (config.verbose)
    {
        std::cout << "<<<<<<<<<<<<<<<<<< Exact Solver >>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "===================================================" << std::endl;
    }

    // Restart timer for measuring the time of the exact solver
    t.restart();

    // NB: If we reach this point, we could not stop early and the hypergraph needs to be solved using the exact solver
    switch (config.baseSolver)
    {
    case BaseSolver::ILP:
    {
#ifdef USE_GUROBI
        try
        {
            // Create an environment
            GRBEnv env = GRBEnv(true);
            env.start();
            // Set the log file
            env.set("LogFile", "gurobi.log");

            // Initialize the ILP solver
            ILPMincut mincutILPSolver(get_current_num_nodes(hypergraph),
                                      get_current_num_edges(hypergraph),
                                      env,
                                      config.seed,
                                      config.ilpTimeout,
                                      config.ilpMode,
                                      config.numThreads);
            // Add the node variables and constraints to the ILP
            mincutILPSolver.add_node_variables_and_constraints(hypergraph);
            // Add the edge variables and constraints to the ILP
            mincutILPSolver.add_edge_variables_and_constraints(hypergraph);
            // Set the objective function of the ILP
            mincutILPSolver.set_objective();
            // Optimize the ILP
            mincutILPSolver.optimize();
            // Solve the minimum cut problem via ILP
            CutValue ilpMinEdgeCut = mincutILPSolver.get_result(hypergraph);
            minEdgeCut = std::min(minEdgeCut, ilpMinEdgeCut);
            double ilpTime = t.elapsed();
            totalComputingTime += ilpTime;
            // Print the results of the ILP
            if (config.verbose)
            {
                std::cout << "base_solver_mincut_value \t\t" << ilpMinEdgeCut << std::endl;
                std::cout << "base_solver_time \t\t\t" << ilpTime << std::endl;
            }
        }
        catch (GRBException e)
        {
            std::cerr << "Error code = " << e.getErrorCode() << std::endl;
            std::cerr << e.getMessage() << std::endl;
        }
#else
        std::cerr << "Error: ILP base solver requires Gurobi. Build with -DUSE_GUROBI=ON or use --base_solver=submodular." << std::endl;
#endif
        break;
    }
    case BaseSolver::SUBMODULAR:
    {
        // Before we can use the submodular solver, we need to convert the static hypergraph to a dynamic hypergraph
        // Stores whether the hypergraph has weighted edges
        bool hasWeightedEdges;
        // Create the dynamic hypergraph
        mt_kahypar_hypergraph_t dynamicHypergraphWrapper = HypergraphConverter::static_hypergraph_to_dynamic_hypergraph(hypergraph, hasWeightedEdges);
        // Make sure that the created hypergraph is a dynamic hypergraph
        assert(dynamicHypergraphWrapper.type == DYNAMIC_HYPERGRAPH);
        // Cast the dynamic hypergraph wrapper to a dynamic hypergraph
        DynamicHypergraph &dynamicHypergraph = mt_kahypar::utils::cast<DynamicHypergraph>(dynamicHypergraphWrapper);

        double conversionTime = t.elapsed();
        totalComputingTime += conversionTime;
        // Print the time needed for the conversion of the static hypergraph to the dynamic hypergraph
        if (config.verbose)
            std::cout << "static_to_dynamic_conversion_time \t" << conversionTime << std::endl;

        // Restart timer
        t.restart();

        // Initialize the submodular mincut solver
        SubmodularMincut mincutSubmodularSolver(dynamicHypergraph.initialNumNodes(), dynamicHypergraph.initialNumEdges(), config.orderingType, config.orderingMode, hasWeightedEdges, config.numThreads);

        // Solve the minimum cut problem via submodular optimization
        SubmodularMincutResult result = mincutSubmodularSolver.solve(dynamicHypergraph);
        minEdgeCut = std::min(minEdgeCut, result.minEdgeCut);
        double submodularTime = t.elapsed();
        totalComputingTime += submodularTime;
        // Print the results of the submodular solver
        if (config.verbose)
        {
            std::cout << "base_solver_mincut_value \t\t" << result.minEdgeCut << std::endl;
            std::cout << "base_solver_num_iterations \t\t" << result.numIterations << std::endl;
            std::cout << "base_solver_contractions_per_it \t" << result.meanContractionsPerIteration << std::endl;
            std::cout << "base_solver_time \t\t\t" << submodularTime << std::endl;
        }

        // Free the dynamic hypergraph wrapper
        mt_kahypar_free_hypergraph(dynamicHypergraphWrapper);
        break;
    }
    default:
        std::cerr << "Error: Unknown exact solver." << std::endl;
        break;
    }

    return {naiveEstimate, minEdgeCut, totalComputingTime};
}
