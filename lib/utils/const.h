/******************************************************************************
 * const.h
 * *
 * General constants for the algorithms.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_CONSTANTS_H
#define SMHM_CONSTANTS_H

// Own headers
#include "lib/utils/definitions.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"

// Default seed
const int DEFAULT_SEED = 0;

// Default timeout for the ILP in seconds
const double DEFAULT_TIMEOUT_ILP = 7200;

// Default preset type when reading the hypergraph
const mt_kahypar_preset_type_t DEFAULT_PRESET_TYPE = mt_kahypar_preset_type_t::DETERMINISTIC;

// Default option for the base solver of the kernelizer
const BaseSolver DEFAULT_BASE_SOLVER = BaseSolver::SUBMODULAR;

// Default option for the verbose mode
const bool DEFAULT_VERBOSE = false;

// Default number of threads (0 means all available cores)
const size_t DEFAULT_PARALLEL_NUM_THREADS = 0;

// Default ordering type
const OrderingType DEFAULT_ORDERING_TYPE = OrderingType::TIGHT;

// Default ordering contraction mode
const OrderingMode DEFAULT_ORDERING_MODE = OrderingMode::SINGLE;

// Default pruning mode
const PruningMode DEFAULT_PRUNING_MODE = PruningMode::BEST;

// Default ilp mode
const ILPMode DEFAULT_ILP_MODE = ILPMode::BIP;

// Default number of iterations for the label propagation algorithm
const IterationIndex DEFAULT_LP_NUM_ITERATIONS = 0;

// Default number of pins to sample for the probabilistic label propagation algorithm
const NodeIndex DEFAULT_LP_NUM_PINS_TO_SAMPLE = 25;

// Default mode for the label propagation algorithm
const LabelPropagationMode DEFAULT_LP_MODE = LabelPropagationMode::CLIQUE_EXPANDED;

// Default file format
const mt_kahypar_file_format_type_t DEFAULT_FILE_FORMAT = mt_kahypar_file_format_type_t::HMETIS;

// Default edge size for the dumbbell generator
const NodeIndex DEFAULT_DUMBBELL_EDGE_SIZE = 3;

// Default number of nodes per side for the dumbbell generator
const NodeIndex DEFAULT_DUMBBELL_NUM_NODES_PER_SIDE = 100;

// Default minimum degree for the dumbbell generator
// NB: 0 means that each side of the dumbbell is a complete hypergraph
const EdgeIndex DEFAULT_DUMBBELL_MIN_DEGREE = 0;

// Default option on whether the dumbbell hypergraph should be weighted
const bool DEFAULT_DUMBBELL_WEIGHTED = false;

// Default option on whether the hypergraph should be processed as unweighted
const bool DEFAULT_UNWEIGHTED = false;

// Default option on whether the hypercactus generator should allow floating point edge weights in the output
const bool DEFAULT_ALLOW_OUTPUT_FLOAT_EDGE_WEIGHTS = false;

#endif // end of SMHM_CONSTANTS_H
