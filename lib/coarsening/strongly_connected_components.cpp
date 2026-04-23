/******************************************************************************
 * strongly_connected_components.cpp
 * *
 * Computes the strongly connected components (SCCs) of a directed graph. The
 * SCCs are computed using a depth-first search (DFS) based algorithm by Dijkstra.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "strongly_connected_components.h"
#include "lib/utils/definitions.h"

StronglyConnectedComponentsFinder::StronglyConnectedComponentsFinder(const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges)
    : dfsNum(initialNumNodes, std::numeric_limits<NodeIndex>::max()),
      traversedEdge(initialNumEdges, false) {};