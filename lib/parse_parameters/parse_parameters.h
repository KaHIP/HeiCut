/******************************************************************************
 * parse_parameters.h
 * *
 * Parses the parameters provided from the command line.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_PARSE_PARAMS_H
#define SMHM_PARSE_PARAMS_H

// Own headers
#include "lib/utils/definitions.h"

struct BaseConfig
{
    const char *hypergraphFileName;
    mt_kahypar_file_format_type_t hypergraphFileFormat;
    int seed;
    mt_kahypar_preset_type_t presetType;
    size_t numThreads;
    bool unweighted;
};

struct KernelizerConfig : BaseConfig
{
    BaseSolver baseSolver;
    double ilpTimeout;
    ILPMode ilpMode;
    OrderingType orderingType;
    OrderingMode orderingMode;
    PruningMode pruningMode;
    IterationIndex LPNumIterations;
    LabelPropagationMode LPMode;
    NodeIndex LPNumPinsToSample;
    bool verbose;
};

struct ILPConfig : BaseConfig
{
    double ilpTimeout;
    ILPMode ilpMode;
};

struct TrimmerConfig : BaseConfig
{
    OrderingType orderingType;
    OrderingMode orderingMode;
};

struct SubmodularConfig : BaseConfig
{
    OrderingType orderingType;
    OrderingMode orderingMode;
};

struct DumbbellGeneratorConfig
{
    int seed;
    NodeIndex numNodesPerSide;
    NodeIndex edgeSize;
    EdgeIndex minDegree;
    bool weighted;
    const char *outputFileName;
};

struct HypercactusGeneratorConfig : BaseConfig
{
    bool verbose;
    bool allowOutputFloatEdgeWeights;
    const char *outputFileName;
};

struct KCoreGeneratorConfig : BaseConfig
{
    const char *outputFileName;
};

class ParseParams
{
public:
    // Print the help content
    static void print_help_content(const char *nameOfProgram, std::vector<void *> &argtable);

    // Print the error
    static void print_error(const char *nameOfProgram, std::vector<void *> &argtable, struct arg_end *end);

    // Check if the help was requested or if there were any errors
    static int check_for_help_or_errors(int argc, char **argv, const char *nameOfProgram, std::vector<void *> &argtable, struct arg_lit *help, struct arg_end *end);

    // Parse the parameters for the hypercactus algorithm
    static int parse_parameters_hypercactus(int argc, char **argv, HypercactusGeneratorConfig &config);

    // Parse the parameters for the kernelizer of the mininum cut value
    static int parse_parameters_kernelizer(int argc, char **argv, KernelizerConfig &config);

    // Parse the parameters for the mincut using the ILP
    static int parse_parameters_mincut_ilp(int argc, char **argv, ILPConfig &config);

    // Parse the parameters for the mincut using the k-trimmer
    static int parse_parameters_mincut_trimmer(int argc, char **argv, TrimmerConfig &config);

    // Parse the parameters for the mincut using the submodular solver
    static int parse_parameters_mincut_submodular(int argc, char **argv, SubmodularConfig &config);

    // Parse the parameters for the dumbbell generator
    static int parse_parameters_dumbbell_generator(int argc, char **argv, DumbbellGeneratorConfig &config);

    // Parse the parameters for the k-core generator
    static int parse_parameters_kcore_generator(int argc, char **argv, KCoreGeneratorConfig &config);
};

#endif // end of SMHM_PARSE_PARAMS_H
