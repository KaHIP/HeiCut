/******************************************************************************
 * strongly_connected_components.h
 * *
 * Computes the strongly connected components (SCCs) of a directed graph. The
 * SCCs are computed using a depth-first search (DFS) based algorithm by Dijkstra.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_SCCS_H
#define SMHM_SCCS_H

#include <iostream>
#include <stack>
#include <vector>
// Own headers
#include "lib/utils/definitions.h"

class StronglyConnectedComponentsFinder
{
private:
    // Store the DFS number for each node in the graph
    std::vector<NodeIndex> dfsNum;
    // Store the visited nodes counter (will be used to assign the DFS number)
    NodeIndex visitedNodesCounter = 0;
    // Store the component counter (will be used to assign the SCC ID)
    ClusterIndex componentCounter = 0;
    // Store the open representatives ordered by increasing DFS number
    std::stack<NodeID> openRepresentatives;
    // Store the open nodes ordered by increasing DFS number
    std::stack<NodeID> openNodes;
    // Store the iteration stack for the DFS
    std::stack<EdgeID> dfsStack;
    // Store for each edge if we have traversed it (will be used to detect when to backtrack)
    std::vector<bool> traversedEdge;

    // Set the node as a root
    void root(NodeID nodeID);

    // Traverse the tree edge 
    void traverseTreeEdge(NodeID sourceID, NodeID targetID);

    // Traverse the non-tree edge
    void traverseNonTreeEdge(NodeID sourceID, NodeID targetID, std::vector<ClusterID> &componentID);

    // Backtrack the edge
    void backtrack(NodeID sourceID, NodeID targetID, std::vector<ClusterID> &componentID);

public:
    StronglyConnectedComponentsFinder(const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges);

    // Computes the strongly connected components (SCCs) of a directed graph
    // NB: We consider only edges with a positive weight, i.e. we ignore edges with a zero or negative weight
    ClusterIndex find_sccs(StaticGraph &graph, std::vector<ClusterID> &componentID);
};

// Computes the strongly connected components (SCCs) of a directed graph
// NB: We consider only edges with a positive weight, i.e. we ignore edges with a zero or negative weight
// NB2: The component IDs are stored in the passed vector and the number of components is returned
inline ClusterIndex StronglyConnectedComponentsFinder::find_sccs(StaticGraph &graph, std::vector<ClusterID> &componentID)
{
    // Store the SCC ID for each node in the graph
    // NB: The SCC IDs also correspond to the finishing times of the SCCs in the DFS
    componentID.resize(graph.initialNumNodes());
    std::fill(componentID.begin(), componentID.end(), std::numeric_limits<ClusterID>::max());

    // Perform DFS to find the SCCs in the graph
    for (NodeID nodeID : graph.nodes())
    {
        // Start DFS from every non-visited node
        if (dfsNum[nodeID] == std::numeric_limits<NodeIndex>::max())
        {
            // Set the node as a root
            root(nodeID);

            // Push all incident edges (i.e. with a positive weight) of the node onto the stack
            for (EdgeID edgeID : graph.incidentEdges(nodeID))
                if (graph.edgeWeight(edgeID) > 0)
                    dfsStack.push(edgeID);

            // Handle all edges in the stack
            while (!dfsStack.empty())
            {
                // Get the top edge from the stack
                EdgeID edgeID = dfsStack.top();
                // Get the source of the edge
                NodeID sourceID = graph.edgeSource(edgeID);
                // Get the target of the edge
                NodeID targetID = graph.edgeTarget(edgeID);

                // Check if we have already traversed the edge
                if (traversedEdge[edgeID])
                {
                    // If so, we backtrack the edge
                    backtrack(sourceID, targetID, componentID);
                    // Remove the edge from the stack
                    dfsStack.pop();
                }
                else
                {
                    // Otherwise, mark the edge as traversed
                    traversedEdge[edgeID] = true;
                    // Check if we see the target of the edge for the first time
                    if (dfsNum[targetID] == std::numeric_limits<NodeIndex>::max())
                    {
                        // If so, traverse the tree edge
                        traverseTreeEdge(sourceID, targetID);
                        // Push all incident edges (i.e. with a positive weight) of the  target onto the stack
                        for (EdgeID targetEdgeID : graph.incidentEdges(targetID))
                            if (graph.edgeWeight(targetEdgeID) > 0)
                                dfsStack.push(targetEdgeID);
                    }
                    else
                    {
                        // Otherwise, traverse the non-tree edge
                        traverseNonTreeEdge(sourceID, targetID, componentID);
                        // Remove the edge from the stack
                        dfsStack.pop();
                    }
                }
            }

            // Backtrack the imaginary edge (nodeID, nodeID)
            backtrack(nodeID, nodeID, componentID);
        }
    }

    assert(openRepresentatives.empty());
    assert(openNodes.empty());
    assert(visitedNodesCounter == graph.initialNumNodes());
    assert(componentCounter > 0);

    return componentCounter;
}

// Set the node as a root
inline void StronglyConnectedComponentsFinder::root(NodeID nodeID)
{
    dfsNum[nodeID] = visitedNodesCounter++;
    openRepresentatives.push(nodeID);
    openNodes.push(nodeID);
}

// Traverse the tree edge
inline void StronglyConnectedComponentsFinder::traverseTreeEdge(NodeID sourceID, NodeID targetID)
{
    dfsNum[targetID] = visitedNodesCounter++;
    openRepresentatives.push(targetID);
    openNodes.push(targetID);
}

// Traverse the non-tree edge
inline void StronglyConnectedComponentsFinder::traverseNonTreeEdge(NodeID sourceID, NodeID targetID, std::vector<ClusterID> &componentID)
{
    // Ensure that the target is an open node
    if (componentID[targetID] != std::numeric_limits<ClusterID>::max())
        return;

    while (!openRepresentatives.empty() && dfsNum[openRepresentatives.top()] > dfsNum[targetID])
        openRepresentatives.pop();
}

// Backtrack the edge
inline void StronglyConnectedComponentsFinder::backtrack(NodeID sourceID, NodeID targetID, std::vector<ClusterID> &componentID)
{
    // Ensure that the target is the top of the open representatives stack
    if (targetID != openRepresentatives.top())
        return;

    openRepresentatives.pop();
    NodeID topOpenNodeID;
    do
    {
        topOpenNodeID = openNodes.top();
        openNodes.pop();
        componentID[topOpenNodeID] = componentCounter;
    } while (topOpenNodeID != targetID);
    componentCounter++;
}

#endif // end of SMHM_SCCS_H