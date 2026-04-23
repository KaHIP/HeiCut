/******************************************************************************
 * output.h
 * *
 * Operators for the output of results.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_OUTPUT_H
#define SMHM_OUTPUT_H

#include <iostream>
// Own headers
#include "lib/utils/definitions.h"

// Operator for the output of the ordering type
inline std::ostream &operator<<(std::ostream &os, OrderingType orderingType)
{
    switch (orderingType)
    {
#ifdef SMHM_PARALLEL
    case OrderingType::MIX_DISCRETE:
        return os << "mix_discrete";
    case OrderingType::MIX_UNIFORM:
        return os << "mix_uniform";
#endif
    case OrderingType::MA:
        return os << "max_adjacency";
    case OrderingType::TIGHT:
        return os << "tight";
    case OrderingType::QUEYRANNE:
        return os << "queyranne";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the ordering mode
inline std::ostream &operator<<(std::ostream &os, OrderingMode orderingMode)
{
    switch (orderingMode)
    {
    case OrderingMode::SINGLE:
        return os << "single";
    case OrderingMode::MULTI:
        return os << "multi";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the label propagation mode
inline std::ostream &operator<<(std::ostream &os, LabelPropagationMode labelPropagationMode)
{
    switch (labelPropagationMode)
    {
    case LabelPropagationMode::CLIQUE_EXPANDED:
        return os << "clique_expanded";
    case LabelPropagationMode::PROBABILISTIC:
        return os << "probabilistic";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the pruning mode
inline std::ostream &operator<<(std::ostream &os, PruningMode pruningMode)
{
    switch (pruningMode)
    {
    case PruningMode::ALL:
        return os << "all";
    case PruningMode::BEST:
        return os << "best";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the ILP mode
inline std::ostream &operator<<(std::ostream &os, ILPMode ilpMode)
{
    switch (ilpMode)
    {
    case ILPMode::BIP:
        return os << "BIP";
    case ILPMode::MILP:
        return os << "MILP";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the base solver
inline std::ostream &operator<<(std::ostream &os, BaseSolver baseSolver)
{
    switch (baseSolver)
    {
    case BaseSolver::ILP:
        return os << "ilp";
    case BaseSolver::SUBMODULAR:
        return os << "submodular";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the file format
inline std::ostream &operator<<(std::ostream &os, mt_kahypar_file_format_type_t fileFormat)
{
    switch (fileFormat)
    {
    case mt_kahypar_file_format_type_t::HMETIS:
        return os << "HMETIS";
    case mt_kahypar_file_format_type_t::METIS:
        return os << "METIS";
    default:
        return os << "unknown";
    }
}

// Operator for the output of the preset type
inline std::ostream &operator<<(std::ostream &os, mt_kahypar_preset_type_t presetType)
{
    switch (presetType)
    {
    case mt_kahypar_preset_type_t::DETERMINISTIC:
        return os << "deterministic";
    case mt_kahypar_preset_type_t::LARGE_K:
        return os << "large_k";
    case mt_kahypar_preset_type_t::DEFAULT:
        return os << "default";
    case mt_kahypar_preset_type_t::QUALITY:
        return os << "quality";
    case mt_kahypar_preset_type_t::HIGHEST_QUALITY:
        return os << "highest_quality";
    default:
        return os << "unknown";
    }
}

#endif // SMHM_OUTPUT_H