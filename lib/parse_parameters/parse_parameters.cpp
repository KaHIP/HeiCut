/******************************************************************************
 * parse_parameters.cpp
 * *
 * Parses the parameters provided from the command line.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
#include <string.h>
#include <sstream>
#include <regex.h>
// Own headers
#include "parse_parameters.h"
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
// Argtable headers
#include "argtable3.h"

// Define the base argument table
struct BaseArgtable
{
    bool isStaticHypergraph = true;

    struct arg_lit *help = arg_lit0("h", "help", "Print help.");
    struct arg_str *hypergraphFileName = arg_strn(NULL, NULL, "PATH_TO_HYPERGRAPH", 1, 1, "Path to hypergraph file.");
    struct arg_rex *hypergraphFileFormat = arg_rex0(NULL, "file_format", "^(HMETIS|METIS)$", "fileFormat", REG_EXTENDED, "Format [HMETIS|METIS] of the hypergraph file. Default: HMETIS.");
    struct arg_int *seed = arg_int0(NULL, "seed", NULL, "Seed to use for PRNG. Default: 0.");
    struct arg_rex *presetType = arg_rex0(NULL, "preset_type", "^(deterministic|default)$", "fileFormat", REG_EXTENDED, "Preset type [deterministic|default] when reading the hypergraph. Default: deterministic.");
    struct arg_int *numThreads = arg_int0("t", "num_threads", NULL, "Number of threads used. 0 means all available cores are used. Default: 0.");
    struct arg_lit *unweighted = arg_lit0(NULL, "unweighted", "Parse hypergraph as unweighted. Default: false.");

    BaseArgtable(bool isStaticHypergraph = true) : isStaticHypergraph(isStaticHypergraph) {};

    std::vector<void *> to_vector()
    {
        std::vector<void *> argtable = {help, hypergraphFileName, hypergraphFileFormat, seed, unweighted};
        if (isStaticHypergraph)
            argtable.push_back(presetType);
#ifdef SMHM_PARALLEL
        argtable.push_back(numThreads);
#endif
        return argtable;
    };

    void extract_arguments(BaseConfig &config)
    {
        config.hypergraphFileName = hypergraphFileName->sval[0]; // hypergraphFileName is a required argument
        config.seed = (seed->count > 0) ? seed->ival[0] : DEFAULT_SEED;
        config.unweighted = (unweighted->count > 0) ? true : DEFAULT_UNWEIGHTED;

#ifdef SMHM_PARALLEL
        config.numThreads = (numThreads->count > 0) ? numThreads->ival[0] : DEFAULT_PARALLEL_NUM_THREADS;
#else
        config.numThreads = 1;
#endif

        if (isStaticHypergraph)
            if (presetType->count > 0)
                if (strcmp("deterministic", presetType->sval[0]) == 0)
                    config.presetType = mt_kahypar_preset_type_t::DETERMINISTIC;
                else if (strcmp("default", presetType->sval[0]) == 0)
                    config.presetType = mt_kahypar_preset_type_t::DEFAULT;
                else
                    config.presetType = DEFAULT_PRESET_TYPE;
            else
                config.presetType = DEFAULT_PRESET_TYPE;
        else
            config.presetType = mt_kahypar_preset_type_t::HIGHEST_QUALITY;

        if (hypergraphFileFormat->count > 0)
            if (strcmp("HMETIS", hypergraphFileFormat->sval[0]) == 0)
                config.hypergraphFileFormat = mt_kahypar_file_format_type_t::HMETIS;
            else if (strcmp("METIS", hypergraphFileFormat->sval[0]) == 0)
                config.hypergraphFileFormat = mt_kahypar_file_format_type_t::METIS;
            else
                config.hypergraphFileFormat = DEFAULT_FILE_FORMAT;
        else
            config.hypergraphFileFormat = DEFAULT_FILE_FORMAT;
    }
};

// Print the help content
void ParseParams::print_help_content(const char *nameOfProgram, std::vector<void *> &argtable)
{
    printf("Usage: %s", nameOfProgram);
    arg_print_syntax(stdout, argtable.data(), "\n");
    arg_print_glossary(stdout, argtable.data(), "  %-40s %s\n");
    arg_freetable(argtable.data(), argtable.size());
}

// Print the error
void ParseParams::print_error(const char *nameOfProgram, std::vector<void *> &argtable, struct arg_end *end)
{
    arg_print_errors(stderr, end, nameOfProgram);
    printf("Try '%s --help' for more information.\n", nameOfProgram);
    arg_freetable(argtable.data(), argtable.size());
}

// Check if the help was requested or if there were any errors
int ParseParams::check_for_help_or_errors(int argc, char **argv, const char *nameOfProgram, std::vector<void *> &argtable, struct arg_lit *help, struct arg_end *end)
{
    // Parse the arguments
    int nmbErrors = arg_parse(argc, argv, argtable.data());

    // Check if the help was requested
    if (help->count > 0)
    {
        print_help_content(nameOfProgram, argtable);
        return 1;
    }

    // Check if there were any errors
    if (nmbErrors > 0)
    {
        print_error(nameOfProgram, argtable, end);
        return 1;
    }
    return 0;
}

// Parse the parameters for the hypercactus algorithm
int ParseParams::parse_parameters_hypercactus(int argc, char **argv, HypercactusGeneratorConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable;
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
    struct arg_str *outputFileName = arg_strn(NULL, NULL, "PATH_TO_OUTPUT", 1, 1, "Path to output file to store the generated hypercactus in HMETIS format.");
    struct arg_lit *allowOutputFloatEdgeWeights = arg_lit0("f", "output_float_edge_weights", "Allow floating point edge weights in the output. Default: false.");
    struct arg_lit *verbose = arg_lit0("v", "verbose", "Output additional information. Default: false.");
    struct arg_end *end = arg_end(100);

    // Add the rest of the arguments to the argument table
    argtable.push_back(outputFileName);
    argtable.push_back(allowOutputFloatEdgeWeights);
    argtable.push_back(verbose);
    argtable.push_back(end);

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);
    config.outputFileName = outputFileName->sval[0]; // outputFileName is a required argument
    config.allowOutputFloatEdgeWeights = (allowOutputFloatEdgeWeights->count > 0) ? true : DEFAULT_ALLOW_OUTPUT_FLOAT_EDGE_WEIGHTS;
    config.verbose = (verbose->count > 0) ? true : DEFAULT_VERBOSE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the label propagation
int ParseParams::parse_parameters_kernelizer(int argc, char **argv, KernelizerConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable;
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
    struct arg_rex *baseSolver = arg_rex0(NULL, "base_solver", "^(submodular|ilp)$", "BaseSolver", REG_EXTENDED, "Underlying base solver [submodular|ilp] after the pruning rules. Default: submodular.");
    struct arg_int *ilpTimeout = arg_int0(NULL, "ilp_timeout", NULL, "Timeout of ILP in seconds. Only used if base_solver is ilp. Default: 7200s (= 2h).");
    struct arg_rex *ilpMode = arg_rex0(NULL, "ilp_mode", "^(BIP|MILP)$", "ILPMode", REG_EXTENDED, "Mode [BIP|MILP] for the ILP. Only used if base_solver is ilp. Default: BIP.");
#ifdef SMHM_PARALLEL
    struct arg_rex *orderingType = arg_rex0(NULL, "ordering_type", "^(mix_discrete|mix_uniform|max_adjacency|tight|queyranne)$", "OrderingType", REG_EXTENDED, "Ordering type [mix_discrete|mix_uniform|max_adjacency|tight|queyranne] for the submodular solver. Only used if base_solver is submodular. Default: tight.");
#else
    struct arg_rex *orderingType = arg_rex0(NULL, "ordering_type", "^(max_adjacency|tight|queyranne)$", "OrderingType", REG_EXTENDED, "Ordering type [max_adjacency|tight|queyranne] for the submodular solver. Only used if base_solver is submodular. Default: tight.");
#endif
    struct arg_rex *orderingMode = arg_rex0(NULL, "ordering_mode", "^(single|multi)$", "OrderingMode", REG_EXTENDED, "Contraction mode [single|multi] for the ordering. Only used if base_solver is submodular. Default: single.");
    struct arg_rex *pruningMode = arg_rex0(NULL, "pruning_mode", "^(best|all)$", "PruningMode", REG_EXTENDED, "Mode [best|all] for applying the pruning rules. Default: best.");
    struct arg_int *LPNumIterations = arg_int0(NULL, "lp_num_iterations", NULL, "Number of iterations for the label propagation. Default: 0.");
    struct arg_rex *LPMode = arg_rex0(NULL, "lp_mode", "^(clique_expanded|probabilistic)$", "LabelPropagationMode", REG_EXTENDED, "Mode [clique_expanded|probabilistic] of the label propagation algorithm. Default: clique_expanded.");
    struct arg_int *LPNumPinsToSample = arg_int0(NULL, "lp_num_pins_to_sample", NULL, "Number of pins to sample from each edge. Only used if lp_mode is probabilistic. Default: 25.");
    struct arg_lit *verbose = arg_lit0("v", "verbose", "Output additional information. Default: false.");
    struct arg_end *end = arg_end(100);

    // Add the rest of the arguments to the argument table
    argtable.push_back(baseSolver);
    argtable.push_back(ilpTimeout);
    argtable.push_back(ilpMode);
    argtable.push_back(orderingType);
    argtable.push_back(orderingMode);
    argtable.push_back(pruningMode);
    argtable.push_back(LPNumIterations);
    argtable.push_back(LPMode);
    argtable.push_back(LPNumPinsToSample);
    argtable.push_back(verbose);
    argtable.push_back(end);

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);

    // Extract the arguments
    config.LPNumIterations = (LPNumIterations->count > 0) ? LPNumIterations->ival[0] : DEFAULT_LP_NUM_ITERATIONS;
    config.LPNumPinsToSample = (LPNumPinsToSample->count > 0) ? LPNumPinsToSample->ival[0] : DEFAULT_LP_NUM_PINS_TO_SAMPLE;
    config.ilpTimeout = (ilpTimeout->count > 0) ? ilpTimeout->ival[0] : DEFAULT_TIMEOUT_ILP;
    config.verbose = (verbose->count > 0) ? true : DEFAULT_VERBOSE;

    if (LPMode->count > 0)
        if (strcmp("clique_expanded", LPMode->sval[0]) == 0)
            config.LPMode = LabelPropagationMode::CLIQUE_EXPANDED;
        else if (strcmp("probabilistic", LPMode->sval[0]) == 0)
            config.LPMode = LabelPropagationMode::PROBABILISTIC;
        else
            config.LPMode = DEFAULT_LP_MODE;
    else
        config.LPMode = DEFAULT_LP_MODE;

    if (baseSolver->count > 0)
        if (strcmp("ilp", baseSolver->sval[0]) == 0)
            config.baseSolver = BaseSolver::ILP;
        else if (strcmp("submodular", baseSolver->sval[0]) == 0)
            config.baseSolver = BaseSolver::SUBMODULAR;
        else
            config.baseSolver = DEFAULT_BASE_SOLVER;
    else
        config.baseSolver = DEFAULT_BASE_SOLVER;

    if (orderingType->count > 0)
        if (strcmp("max_adjacency", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MA;
        else if (strcmp("tight", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::TIGHT;
        else if (strcmp("queyranne", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::QUEYRANNE;
#ifdef SMHM_PARALLEL
        else if (strcmp("mix_discrete", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MIX_DISCRETE;
        else if (strcmp("mix_uniform", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MIX_UNIFORM;
#endif
        else
            config.orderingType = DEFAULT_ORDERING_TYPE;
    else
        config.orderingType = DEFAULT_ORDERING_TYPE;

    if (orderingMode->count > 0)
        if (strcmp("single", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::SINGLE;
        else if (strcmp("multi", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::MULTI;
        else
            config.orderingMode = DEFAULT_ORDERING_MODE;
    else
        config.orderingMode = DEFAULT_ORDERING_MODE;

    if (pruningMode->count > 0)
        if (strcmp("all", pruningMode->sval[0]) == 0)
            config.pruningMode = PruningMode::ALL;
        else if (strcmp("best", pruningMode->sval[0]) == 0)
            config.pruningMode = PruningMode::BEST;
        else
            config.pruningMode = DEFAULT_PRUNING_MODE;
    else
        config.pruningMode = DEFAULT_PRUNING_MODE;

    if (ilpMode->count > 0)
        if (strcmp("BIP", ilpMode->sval[0]) == 0)
            config.ilpMode = ILPMode::BIP;
        else if (strcmp("MILP", ilpMode->sval[0]) == 0)
            config.ilpMode = ILPMode::MILP;
        else
            config.ilpMode = DEFAULT_ILP_MODE;
    else
        config.ilpMode = DEFAULT_ILP_MODE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the mincut using the ILP
int ParseParams::parse_parameters_mincut_ilp(int argc, char **argv, ILPConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable;
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
    struct arg_int *ilpTimeout = arg_int0(NULL, "ilp_timeout", NULL, "Timeout of ILP in seconds. Default: 7200s (= 2h).");
    struct arg_rex *ilpMode = arg_rex0(NULL, "ilp_mode", "^(BIP|MILP)$", "ILPMode", REG_EXTENDED, "Mode [BIP|MILP] for the ILP. Default: BIP.");
    struct arg_end *end = arg_end(100);

    // Add the rest of the arguments to the argument table
    argtable.push_back(ilpTimeout);
    argtable.push_back(ilpMode);
    argtable.push_back(end);

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);

    // Extract the arguments
    config.ilpTimeout = (ilpTimeout->count > 0) ? ilpTimeout->ival[0] : DEFAULT_TIMEOUT_ILP;

    if (ilpMode->count > 0)
        if (strcmp("BIP", ilpMode->sval[0]) == 0)
            config.ilpMode = ILPMode::BIP;
        else if (strcmp("MILP", ilpMode->sval[0]) == 0)
            config.ilpMode = ILPMode::MILP;
        else
            config.ilpMode = DEFAULT_ILP_MODE;
    else
        config.ilpMode = DEFAULT_ILP_MODE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the mincut using the k-trimmer
int ParseParams::parse_parameters_mincut_trimmer(int argc, char **argv, TrimmerConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable;
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
    struct arg_rex *orderingType = arg_rex0(NULL, "ordering_type", "^(max_adjacency|tight|queyranne)$", "OrderingType", REG_EXTENDED, "Ordering type [max_adjacency|tight|queyranne] for the submodular solver. Default: tight.");
    struct arg_rex *orderingMode = arg_rex0(NULL, "ordering_mode", "^(single|multi)$", "OrderingMode", REG_EXTENDED, "Contraction mode [single|multi] for the ordering. Default: single.");
    struct arg_end *end = arg_end(100);
    // Add the rest of the arguments to the argument table
    argtable.insert(argtable.end(), {orderingType, orderingMode, end});

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);

    if (orderingType->count > 0)
        if (strcmp("max_adjacency", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MA;
        else if (strcmp("tight", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::TIGHT;
        else if (strcmp("queyranne", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::QUEYRANNE;
        else
            config.orderingType = DEFAULT_ORDERING_TYPE;
    else
        config.orderingType = DEFAULT_ORDERING_TYPE;

    if (orderingMode->count > 0)
        if (strcmp("single", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::SINGLE;
        else if (strcmp("multi", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::MULTI;
        else
            config.orderingMode = DEFAULT_ORDERING_MODE;
    else
        config.orderingMode = DEFAULT_ORDERING_MODE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the mincut using the submodular solver
int ParseParams::parse_parameters_mincut_submodular(int argc, char **argv, SubmodularConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable(false);
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
#ifdef SMHM_PARALLEL
    struct arg_rex *orderingType = arg_rex0(NULL, "ordering_type", "^(mix_discrete|mix_uniform|max_adjacency|tight|queyranne)$", "OrderingType", REG_EXTENDED, "Ordering type [mix_discrete|mix_uniform|max_adjacency|tight|queyranne] for the submodular solver. Only used if base_solver is submodular. Default: tight.");
#else
    struct arg_rex *orderingType = arg_rex0(NULL, "ordering_type", "^(max_adjacency|tight|queyranne)$", "OrderingType", REG_EXTENDED, "Ordering type [max_adjacency|tight|queyranne] for the submodular solver. Default: tight.");
#endif
    struct arg_rex *orderingMode = arg_rex0(NULL, "ordering_mode", "^(single|multi)$", "OrderingMode", REG_EXTENDED, "Contraction mode [single|multi] for the ordering. Default: single.");
    struct arg_end *end = arg_end(100);
    // Add the rest of the arguments to the argument table
    argtable.insert(argtable.end(), {orderingType, orderingMode, end});

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);

    if (orderingType->count > 0)
        if (strcmp("max_adjacency", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MA;
        else if (strcmp("tight", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::TIGHT;
        else if (strcmp("queyranne", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::QUEYRANNE;
#ifdef SMHM_PARALLEL
        else if (strcmp("mix_discrete", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MIX_DISCRETE;
        else if (strcmp("mix_uniform", orderingType->sval[0]) == 0)
            config.orderingType = OrderingType::MIX_UNIFORM;
#endif
        else
            config.orderingType = DEFAULT_ORDERING_TYPE;
    else
        config.orderingType = DEFAULT_ORDERING_TYPE;

    if (orderingMode->count > 0)
        if (strcmp("single", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::SINGLE;
        else if (strcmp("multi", orderingMode->sval[0]) == 0)
            config.orderingMode = OrderingMode::MULTI;
        else
            config.orderingMode = DEFAULT_ORDERING_MODE;
    else
        config.orderingMode = DEFAULT_ORDERING_MODE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the dumbbell generator
int ParseParams::parse_parameters_dumbbell_generator(int argc, char **argv, DumbbellGeneratorConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];

    // Define the arguments
    struct arg_lit *help = arg_lit0("h", "help", "Print help.");
    struct arg_str *outputFileName = arg_strn(NULL, NULL, "PATH_TO_OUTPUT", 1, 1, "Path to output file to store the generated dumbbell hypergraph in HMETIS format.");
    struct arg_int *seed = arg_int0(NULL, "seed", NULL, "Seed to use for PRNG. Default: 0.");
    struct arg_int *minDegree = arg_int0(NULL, "min_degree", NULL, "Minimum degree per node. Passing 0 will connect the node to every other node of the same side (complete hypergraph). Default: 0.");
    struct arg_int *numNodesPerSide = arg_int0(NULL, "num_nodes_per_side", NULL, "Number of nodes per side of the dumbbell hypergraph. Default: 100.");
    struct arg_lit *weighted = arg_lit0("w", "weighted", "Generate a weighted hypergraph. Default: false.");
    struct arg_int *edgeSize = arg_int0(NULL, "k", NULL, "Fixed hyperedge size. Default: 3.");
    struct arg_end *end = arg_end(100);

    // Add the arguments to the argument table
    std::vector<void *> argtable = {help, outputFileName, seed, minDegree, numNodesPerSide, weighted, edgeSize, end};

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, help, end))
        return 1;

    // Extract the arguments
    config.outputFileName = outputFileName->sval[0]; // outputFileName is a required argument
    config.seed = (seed->count > 0) ? seed->ival[0] : DEFAULT_SEED;
    config.minDegree = (minDegree->count > 0) ? minDegree->ival[0] : DEFAULT_DUMBBELL_MIN_DEGREE;
    config.numNodesPerSide = (numNodesPerSide->count > 0) ? numNodesPerSide->ival[0] : DEFAULT_DUMBBELL_NUM_NODES_PER_SIDE;
    config.weighted = (weighted->count > 0) ? true : DEFAULT_DUMBBELL_WEIGHTED;
    config.edgeSize = (edgeSize->count > 0) ? edgeSize->ival[0] : DEFAULT_DUMBBELL_EDGE_SIZE;

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}

// Parse the parameters for the k-core generator
int ParseParams::parse_parameters_kcore_generator(int argc, char **argv, KCoreGeneratorConfig &config)
{
    // Get the name of the program
    const char *nameOfProgram = argv[0];
    // Get the base argument table
    BaseArgtable baseArgtable;
    std::vector<void *> argtable = baseArgtable.to_vector();

    // Define the rest of the arguments
    struct arg_str *outputFileName = arg_strn(NULL, NULL, "PATH_TO_OUTPUT", 1, 1, "Path to output file to store the generated dumbbell hypergraph in HMETIS format.");
    struct arg_end *end = arg_end(100);

    // Add the rest of the arguments to the argument table
    argtable.insert(argtable.end(), {outputFileName, end});

    // Check if the help was requested or if there were any errors
    if (check_for_help_or_errors(argc, argv, nameOfProgram, argtable, baseArgtable.help, end))
        return 1;

    // Extract the arguments
    baseArgtable.extract_arguments(config);

    // Extract the arguments
    config.outputFileName = outputFileName->sval[0]; // outputFileName is a required argument

    arg_freetable(argtable.data(), argtable.size());
    return 0;
}
