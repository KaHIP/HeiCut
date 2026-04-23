/******************************************************************************
 * canonical_decomposition.cpp
 * *
 * Decomposes a hypergraph in its canonical decomposition. For this, we first
 * compute the prime decomposition of the hypergraph with the help of a split
 * oracle and we then re-assemble the prime decomposition into the canonical
 * decomposition, which can be used to generate the hypercactus.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <queue>
// Own headers
#include "canonical_decomposition.h"
#include "lib/utils/random.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

CanonicalDecomposition::CanonicalDecomposition(const NodeIndex kernelNumNodes,
                                               const EdgeIndex kernelNumEdges,
                                               std::vector<ClusterID> kernelizationMapping,
                                               const CutValue minEdgeCut,
                                               const bool verbose)
    // NB: Since the tight orderer also computes the tight ordering of the contracted versions of the static hypergraph and since
    //     parallel hyperedges are merged during the contraction of the static hypergraph, the contracted versions of the static hypergraph
    //     are likely to have weighted hyperedges even if the original static hypergraph has unweighted hyperedges. Therefore, we always pass true for hasWeightedEdges.
    : tightOrderer(kernelNumNodes, kernelNumEdges, OrderingType::TIGHT, true, RandomFunctions::get_random_engine()),
      nodeOrdering(kernelNumNodes),
      representativeNodeID(kernelNumNodes),
      kernelizationMapping(kernelizationMapping),
      kernelNumNodes(kernelNumNodes),
      minEdgeCut(minEdgeCut),
      verbose(verbose)
{
  // Initialize the decomposition tree with the root hypergraph
  DecompositionElement decompositionElement;
  decompositionElement.includedNodes.resize(kernelNumNodes);
  // Initialize the included nodes vector to contain 0, 1, 2, ..., kernelNumNodes-1
  std::iota(decompositionElement.includedNodes.begin(), decompositionElement.includedNodes.end(), 0);
  decompositionTree.push_back(decompositionElement);
  numDecompositionElements = 1;
  // Initialize the representative node ID vector to contain 0, 1, 2, ..., kernelNumNodes-1
  std::iota(representativeNodeID.begin(), representativeNodeID.end(), 0);

  // Determine which nodes in the kernelized hypergraph are contracted nodes (i.e. incorporate more than 1 original node)
  // We need to check for these contracted nodes if isolating them leads to a split in the hypergraph.
  // Observation: It is NOT possible that a split runs through a contracted nodes, due to the way the reduction rules are defined.
  std::vector<NodeIndex> kernelNodeSizes(kernelNumNodes, 0);
  for (NodeID nodeID : kernelizationMapping)
  {
    kernelNodeSizes[nodeID]++;
    // Push the node ID into the contracted nodes vector only once (i.e. when we know it incorporates (at least) 2 original nodes)
    if (kernelNodeSizes[nodeID] == 2)
      nodesToCheckForIsolatedSplit.push_back(nodeID);
  }
};

// Recursively compute the prime decomposition of the kernelized hypergraph
void CanonicalDecomposition::compute_prime_decomposition_of_kernel(StaticHypergraph hypergraph, NodeIndex decompositionIndex)
{
  // Get the number of nodes in the hypergraph
  const NodeID numNodes = hypergraph.initialNumNodes();
  // Get the number of edges in the hypergraph
  const EdgeID numEdges = hypergraph.initialNumEdges();
  // Store if we have found a split in the hypergraph
  bool hasFoundSplit = false;
  // Store whether the split is good or not (only relevant if hasFoundSplit is true)
  bool isGoodSplit = false;
  // Store the split (if found)
  std::vector<bool> splitSide(numNodes, 0);
  // Store the source ID and sink ID of the split oracle (if necessary)
  NodeID splitOracleSourceID = std::numeric_limits<NodeID>::max();
  NodeID splitOracleSinkID = std::numeric_limits<NodeID>::max();

  // Make sure that the number of decomposition elements computed so far is equal to the size of the decomposition tree
  assert(numDecompositionElements == decompositionTree.size());

  // Due to the previous kernelization, we need to first check if isolating any of the nodes in nodesToCheckForIsolatedSplit leads to a split
  // If so, we need to put them into their own decomposition element.
  // NB: We only do this check for the first decomposition element (i.e. decompositionIndex == 0), since we do not want to check for isolated splits in
  //     the other decomposition elements. Only if nodesToCheckForIsolatedSplit is empty, we search for a split in the hypergraph using the split oracle.
  while (decompositionIndex == 0 && !hasFoundSplit && !nodesToCheckForIsolatedSplit.empty())
  {
    NodeID nodeID = nodesToCheckForIsolatedSplit.back();
    nodesToCheckForIsolatedSplit.pop_back();

    // If there are only 2 nodes in the hypergraph, then isolating the current node only leads to a split if the other node is also a contracted node (i.e. nodesToCheckForIsolatedSplit is not empty)
    if (numNodes <= 2 && nodesToCheckForIsolatedSplit.empty())
      continue;

    // Compute the cut when isolating the node
    CutValue isolationCut = 0;
    for (EdgeID edgeID : hypergraph.incidentEdges(representativeNodeID[nodeID]))
      if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
        isolationCut += hypergraph.edgeWeight(edgeID);

    // Check if {representativeNodeID[nodeID]} is a split
    if (isolationCut == minEdgeCut)
    {
      // If so, a split has been found
      // NB: The split is good, because if there would be another split that crosses the isolating split, then it would need to go through the isolating node, which is not possible
      //     This means that we can skip checking whether the split is good or not when building the canonical decomposition from the prime decomposition
      hasFoundSplit = true;
      isGoodSplit = true;
      splitSide[representativeNodeID[nodeID]] = 1;
    }
  }

  // Any hypergraph with less than 4 nodes is prime, since it contains no split (i.e. any cut is a trivial cut)
  // NB: We need to exclude the case where we have found a split by isolating a node from the the previous kernelization
  if (numNodes < 4 && !hasFoundSplit)
    return;

  // If we have not found a split yet, query the split oracle to either find a split or return a pair {s, t} of nodes for which there does not exist a (s-t)-split
  if (!hasFoundSplit)
  {
    SplitOracleResult splitOracleResult = query_split_oracle(hypergraph);
    assert(!splitOracleResult.hasFoundSplit || splitOracleResult.splitSide.size() == numNodes);
    if (splitOracleResult.hasFoundSplit)
      splitSide = std::move(splitOracleResult.splitSide);
    splitOracleSourceID = splitOracleResult.sourceID;
    splitOracleSinkID = splitOracleResult.sinkID;
    hasFoundSplit = splitOracleResult.hasFoundSplit;
    // NB: Even if we have found a split, we don't know if the split is good or not, so isGoodSplit remains false, i.e. it will be later checked when building the canonical decomposition
  }

  // If we have not found a split yet, check if {sinkID, sourceID} is a split.
  // NB: We first need to ensure that {sinkID, sourceID} is is not a split before contracting them
  // NB2: Contrary to the algorithm described by Chekuri and Xu, we make this check when going down the recursion stack (i.e. before contracting the nodes)
  //      and not when going up the recursion stack, because this is conceptionally much easier and also faster, since we do not need to uncontract the hypergraph
  if (!hasFoundSplit)
  {
    std::vector<bool> markedHyperedges(numEdges, false);
    CutValue sourceSinkCut = 0;

    // For the source, mark all (enabled) incident hyperedges that have more than one pin and increase the cut value
    for (EdgeID edgeID : hypergraph.incidentEdges(splitOracleSourceID))
      if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
      {
        markedHyperedges[edgeID] = true;
        sourceSinkCut += hypergraph.edgeWeight(edgeID);
      }

    // For the sink, increase the cut value only for all (enabled) UNMARKED incident hyperedges that have more than one pin
    // NB: Marked hyperedges with size 2 are not going over the cut {sinkID, sourceID}, thus we need to remove them from the cut value
    for (EdgeID edgeID : hypergraph.incidentEdges(splitOracleSinkID))
      if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
      {
        if (!markedHyperedges[edgeID])
          sourceSinkCut += hypergraph.edgeWeight(edgeID);
        else if (hypergraph.edgeSize(edgeID) == 2)
          sourceSinkCut -= hypergraph.edgeWeight(edgeID);
      }

    // Check if {sinkID, sourceID} is a split
    if (sourceSinkCut == minEdgeCut)
    {
      // If so, a split has been found
      // NB: The split is good, because if there would be another split that crosses this split, then it would seperate the source and the sink, which is not possible
      //     This means that we can skip checking whether the split is good or not when building the canonical decomposition from the prime decomposition
      hasFoundSplit = true;
      isGoodSplit = true;
      splitSide[splitOracleSourceID] = 1;
      splitSide[splitOracleSinkID] = 1;
    }
  }

  // Print the result of the split oracle
  if (verbose)
  {
    std::cout << "===================================================" << std::endl;
    std::cout << "################### SPLIT RESULT ##################" << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "has_found_split \t\t\t" << (hasFoundSplit ? "true" : "false") << std::endl;
    std::cout << "is_good_split \t\t\t\t" << (isGoodSplit ? "true" : "false") << std::endl;
    if (splitOracleSourceID != std::numeric_limits<NodeID>::max() && splitOracleSinkID != std::numeric_limits<NodeID>::max())
    {
      std::cout << "source_id \t\t\t\t" << splitOracleSourceID + 1 << std::endl;
      std::cout << "sink_id \t\t\t\t" << splitOracleSinkID + 1 << std::endl;
    }
    else
    {
      assert(hasFoundSplit);
      std::cout << "No split oracle called, i.e. split was found by isolating a node from the previous kernelization." << std::endl;
    }

    if (hasFoundSplit)
    {
      for (NodeIndex i = 0; i < numNodes; ++i)
        std::cout << splitSide[i] << " ";
      std::cout << std::endl;
    }
  }

  // Check if a split has been found (either by isolating a node from the previous kernelization or by the split oracle or by isolating the sink and source nodes)
  if (hasFoundSplit)
  {
    // Increase the number of decomposition elements computed so far
    numDecompositionElements++;

    // Store the cluster IDs of the left and right partition when using the split
    // NB: splitSide indicates whether a node is in the left partition (= 0) or right partition (= 1 ) after the split. This means that for the left partition we contract
    // the nodes with split side of 1 and for the right partition we contract the nodes with split side 0
    // NB2: We determine the first node ID that is contracted into the marker node in the left and right partition and then use this node ID for the respective cluster ID
    mt_kahypar::parallel::scalable_vector<ClusterID> leftClusterID(numNodes);
    mt_kahypar::parallel::scalable_vector<ClusterID> rightClusterID(numNodes);
    NodeID firstNodeIDContractedIntoLeftMarker = std::numeric_limits<NodeID>::max();
    NodeID firstNodeIDContractedIntoRightMarker = std::numeric_limits<NodeID>::max();
    for (NodeIndex i = 0; i < numNodes; ++i)
    {
      bool isAssignedToRightSide = (splitSide[i] == 1);
      if (isAssignedToRightSide)
      {
        if (firstNodeIDContractedIntoLeftMarker == std::numeric_limits<NodeID>::max())
          firstNodeIDContractedIntoLeftMarker = i;
      }
      else if (firstNodeIDContractedIntoRightMarker == std::numeric_limits<NodeID>::max())
        firstNodeIDContractedIntoRightMarker = i;
      leftClusterID[i] = isAssignedToRightSide ? firstNodeIDContractedIntoLeftMarker : i;
      rightClusterID[i] = isAssignedToRightSide ? i : firstNodeIDContractedIntoRightMarker;
    }

    // Contract the left and right partition of the hypergraph
    StaticHypergraph leftHypergraph = hypergraph.contract(leftClusterID);
    StaticHypergraph rightHypergraph = hypergraph.contract(rightClusterID);

    // Iterate over all included nodes of old the decomposition element and split them into left and right included nodes
    // NB: We also update the representativeNodeID of the included nodes to point to the representative node in the left or right hypergraph
    std::vector<NodeID> leftIncludedNodes;
    std::vector<NodeID> rightIncludedNodes;
    for (NodeID &nodeID : decompositionTree[decompositionIndex].includedNodes)
      if (splitSide[representativeNodeID[nodeID]] == 0)
      {
        leftIncludedNodes.push_back(nodeID);
        representativeNodeID[nodeID] = leftClusterID[representativeNodeID[nodeID]];
      }
      else
      {
        rightIncludedNodes.push_back(nodeID);
        representativeNodeID[nodeID] = rightClusterID[representativeNodeID[nodeID]];
      }

    // Iterate over all marker edges of the old decomposition element and split them into left and right marker edges
    // NB: We also update the representativeNodeID of the marker edges to point to the representative node in the left or right hypergraph
    std::vector<MarkerEdge> leftMarkerEdges;
    std::vector<MarkerEdge> rightMarkerEdges;
    for (MarkerEdge &markerEdge : decompositionTree[decompositionIndex].markerEdges)
      if (splitSide[markerEdge.representativeNodeID] == 0)
      {
        markerEdge.representativeNodeID = leftClusterID[markerEdge.representativeNodeID];
        decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].targetDecompositionIndex = decompositionIndex;
        // NB: We use leftMarkerEdges.size() as the index since the marker edge was not yet added to the leftMarkerEdges vector
        decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].markerEdgeIndexAtTarget = leftMarkerEdges.size();
        leftMarkerEdges.push_back(markerEdge);
      }
      else
      {
        markerEdge.representativeNodeID = rightClusterID[markerEdge.representativeNodeID];
        decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].targetDecompositionIndex = numDecompositionElements - 1;
        // NB: We use rightMarkerEdges.size() as the index since the marker edge was not yet added to the rightMarkerEdges vector
        decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].markerEdgeIndexAtTarget = rightMarkerEdges.size();
        rightMarkerEdges.push_back(markerEdge);
      }

    // Create two new marker edges that connect the left and right decomposition elements
    MarkerEdge newLeftMarkerEdge = {numMarkers,
                                    leftClusterID[firstNodeIDContractedIntoLeftMarker],
                                    numDecompositionElements - 1,
                                    static_cast<EdgeIndex>(rightMarkerEdges.size()),
                                    isGoodSplit};
    MarkerEdge newRightMarkerEdge = {numMarkers,
                                     rightClusterID[firstNodeIDContractedIntoRightMarker],
                                     decompositionIndex,
                                     static_cast<EdgeIndex>(leftMarkerEdges.size()),
                                     isGoodSplit};
    numMarkers++;
    initialNumMarkers++;
    leftMarkerEdges.push_back(newLeftMarkerEdge);
    rightMarkerEdges.push_back(newRightMarkerEdge);

    // Replace the old decomposition element with the left decomposition element
    decompositionTree[decompositionIndex].includedNodes = std::move(leftIncludedNodes);
    decompositionTree[decompositionIndex].markerEdges = std::move(leftMarkerEdges);
    // Add the right decomposition element as a new decomposition element in the decomposition tree
    DecompositionElement rightDecompositionElement;
    rightDecompositionElement.includedNodes = std::move(rightIncludedNodes);
    rightDecompositionElement.markerEdges = std::move(rightMarkerEdges);
    decompositionTree.push_back(std::move(rightDecompositionElement));

    // Recursively compute the prime decomposition of the left and right decomposition
    compute_prime_decomposition_of_kernel(std::move(leftHypergraph), decompositionIndex);
    compute_prime_decomposition_of_kernel(std::move(rightHypergraph), numDecompositionElements - 1);
    return;
  }

  // If we have not found a split, we need to contract the source and sink nodes
  // Initialize the clusterID to contain 0, 1, 2, ..., numNodes-1
  mt_kahypar::parallel::scalable_vector<ClusterID> clusterID(numNodes);
  std::iota(clusterID.begin(), clusterID.end(), 0);
  // Contract the source and the sink by setting the cluster ID of the source to the sink ID
  clusterID[splitOracleSourceID] = splitOracleSinkID;
  StaticHypergraph simplifiedHypergraph = hypergraph.contract(clusterID);
  // Update the representative node ID of the included nodes in the decomposition element
  std::vector<NodeID> &includedNodes = decompositionTree[decompositionIndex].includedNodes;
  for (NodeID &nodeID : decompositionTree[decompositionIndex].includedNodes)
    representativeNodeID[nodeID] = clusterID[representativeNodeID[nodeID]];
  // Update the representative node ID of the marker edges in the decomposition element
  for (MarkerEdge &markerEdge : decompositionTree[decompositionIndex].markerEdges)
    markerEdge.representativeNodeID = clusterID[markerEdge.representativeNodeID];
  // Recurse on the simplified hypergraph
  compute_prime_decomposition_of_kernel(std::move(simplifiedHypergraph), decompositionIndex);
}

// Convert the prime decomposition of the kernelized hypergraph into the prime decomposition of the original hypergraph
void CanonicalDecomposition::convert_to_prime_decomposition_of_original_hypergraph()
{
  // Store for each kernel node ID the original node IDs that are contracted into it
  std::vector<std::vector<NodeID>> originalNodeIDsPerKernelNodeID(kernelNumNodes);
  for (NodeID originalNodeID = 0; originalNodeID < kernelizationMapping.size(); ++originalNodeID)
    originalNodeIDsPerKernelNodeID[kernelizationMapping[originalNodeID]].push_back(originalNodeID);

  // Loop over all decomposition elements in the decomposition tree
  for (ClusterIndex i = 0; i < numDecompositionElements; ++i)
  {
    // Construct the original included nodes of the decomposition element
    std::vector<NodeID> originalIncludedNodes;
    for (NodeID kernelNodeID : decompositionTree[i].includedNodes)
    {
      // Get the original node IDs that are contracted into the kernel node ID
      const std::vector<NodeID> &originalNodeIDs = originalNodeIDsPerKernelNodeID[kernelNodeID];
      // Add the original node IDs to the original included nodes of the decomposition element
      originalIncludedNodes.insert(originalIncludedNodes.end(), originalNodeIDs.begin(), originalNodeIDs.end());
    }
    // Update the included nodes of the decomposition element
    decompositionTree[i].includedNodes = std::move(originalIncludedNodes);
  }
}

// Compute the canonical decomposition from the prime decomposition by trying to merge the decompositions elements by their markers
void CanonicalDecomposition::compute_canonical_decomposition_from_prime_decomposition(StaticHypergraph &hypergraph)
{
  // Make sure that we have at least one decomposition element (i.e. there exists a prime decomposition)
  assert(numDecompositionElements > 0 && numDecompositionElements == decompositionTree.size());

  // Get the initial number of nodes in the hypergraph
  const NodeID initialNumNodes = hypergraph.initialNumNodes();

  // Initialize the cluster ID vector
  mt_kahypar::parallel::scalable_vector<ClusterID> clusterID(initialNumNodes);

  // Loop over all decomposition elements
  for (ClusterIndex i = 0; i < numDecompositionElements; ++i)
  {
    // Skip the decomposition element if it is deleted or if it has no marker edges
    if (decompositionTree[i].isDeleted || decompositionTree[i].markerEdges.empty())
      continue;

    // Loop over all marker edges of the current decomposition element
    EdgeIndex numMarkerEdges = decompositionTree[i].markerEdges.size();
    EdgeIndex markerEdgeIndex = 0;
    while (markerEdgeIndex < numMarkerEdges)
    {
      // Get the marker edge
      MarkerEdge &markerEdge = decompositionTree[i].markerEdges[markerEdgeIndex];
      // Skip the marker edge if it is already certain that it from a good split
      if (markerEdge.isFromGoodSplit)
      {
        markerEdgeIndex++;
        continue;
      }

      // Make sure that the target decomposition element is not deleted
      assert(!decompositionTree[markerEdge.targetDecompositionIndex].isDeleted);
      // Make sure that the reversed marker edge is not marked as a good split
      assert(!decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].isFromGoodSplit);

      // Reset the cluster ID vector to contain 0, 1, 2, ..., initialNumNodes-1
      std::iota(clusterID.begin(), clusterID.end(), 0);

      // Loop over the other marker edges of both decomposition elements to find which original nodes are contracted into the marker node via BFS
      // Observation: Two decomposition elements cannot share more than one marker node, hence we can perform two independent searches (i.e. we get no cycles)
      for (const MarkerEdge &otherMarkerEdge : decompositionTree[i].markerEdges)
        if (otherMarkerEdge.markerID != markerEdge.markerID)
          find_nodes_contracted_to_marker_node_in_decomposition_tree(clusterID, i, otherMarkerEdge);
      for (const MarkerEdge &otherMarkerEdge : decompositionTree[markerEdge.targetDecompositionIndex].markerEdges)
        if (otherMarkerEdge.markerID != markerEdge.markerID)
          find_nodes_contracted_to_marker_node_in_decomposition_tree(clusterID, markerEdge.targetDecompositionIndex, otherMarkerEdge);

      // Contract the hypergraph with the cluster ID vector to obtain the merged hypergraph from the two decomposition elements
      StaticHypergraph mergedHypergraph = hypergraph.contract(clusterID);

      // Check if the merged hypergraph is a solid polygon
      // A solid polygon is a hypergraph that consists of a GRAPH cycle (including all nodes) where each edge has the same weight A and a hyperedge containing all nodes with weight B.
      // NB: If A = 0 (i.e. there is no graph cycle) OR if B = 0 (i.e. there is no hyperedge containing all nodes), then the hypergraph is still a solid polygon.
      bool isSolidPolygon = is_solid_polygon(mergedHypergraph);

      if (verbose)
      {
        std::cout << "===================================================" << std::endl;
        std::cout << "################ MERGED HYPERGRAPH ################" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "num_nodes \t\t\t\t" << mergedHypergraph.initialNumNodes() << std::endl;
        std::cout << "num_edges \t\t\t\t" << mergedHypergraph.initialNumEdges() << std::endl;
        // Print the cluster ID vector
        std::cout << "cluster_id \t\t\t\t";
        for (NodeID i = 0; i < initialNumNodes; ++i)
          std::cout << clusterID[i] + 1 << " ";
        std::cout << std::endl;
        std::cout << "is_solid_polygon \t\t\t" << (isSolidPolygon ? "true" : "false") << std::endl;
      }

      // If the merged hypergraph is a solid polygon, then the split is not good, i.e. merge the two decomposition elements and remove the marker edges
      if (isSolidPolygon)
      {
        // IMPORTANT: We need to copy the marker edge to avoid using an invalid reference after the marker edge is removed from the decomposition tree
        MarkerEdge oldMarkerEdge = markerEdge;
        // Merge the included nodes of the target decomposition element into the current decomposition element
        std::vector<NodeID> &includedNodes = decompositionTree[i].includedNodes;
        std::vector<NodeID> &targetIncludedNodes = decompositionTree[oldMarkerEdge.targetDecompositionIndex].includedNodes;
        includedNodes.insert(includedNodes.end(), std::make_move_iterator(targetIncludedNodes.begin()), std::make_move_iterator(targetIncludedNodes.end()));
        targetIncludedNodes.clear();
        // Delete the marker edge in the current decomposition element
        std::swap(decompositionTree[i].markerEdges[markerEdgeIndex], decompositionTree[i].markerEdges.back());
        decompositionTree[i].markerEdges.pop_back();
        numMarkerEdges--;
        // Delete the marker edge in the target decomposition element
        std::swap(decompositionTree[oldMarkerEdge.targetDecompositionIndex].markerEdges[oldMarkerEdge.markerEdgeIndexAtTarget], decompositionTree[oldMarkerEdge.targetDecompositionIndex].markerEdges.back());
        decompositionTree[oldMarkerEdge.targetDecompositionIndex].markerEdges.pop_back();
        // Merge the remaining marker edges of the target decomposition element into the current decomposition element
        std::vector<MarkerEdge> &markerEdges = decompositionTree[i].markerEdges;
        std::vector<MarkerEdge> &targetMarkerEdges = decompositionTree[oldMarkerEdge.targetDecompositionIndex].markerEdges;
        numMarkerEdges += targetMarkerEdges.size();
        markerEdges.insert(markerEdges.end(), std::make_move_iterator(targetMarkerEdges.begin()), std::make_move_iterator(targetMarkerEdges.end()));
        targetMarkerEdges.clear();
        // Update the reversed marker edges for all marker edges in the range [markerEdgeIndex, numMarkerEdges - 1]
        for (EdgeIndex j = markerEdgeIndex; j < numMarkerEdges; ++j)
        {
          MarkerEdge &otherMarkerEdge = decompositionTree[i].markerEdges[j];
          decompositionTree[otherMarkerEdge.targetDecompositionIndex].markerEdges[otherMarkerEdge.markerEdgeIndexAtTarget].targetDecompositionIndex = i;
          decompositionTree[otherMarkerEdge.targetDecompositionIndex].markerEdges[otherMarkerEdge.markerEdgeIndexAtTarget].markerEdgeIndexAtTarget = j;
        }
        // Decrease the TOTAL number of markers
        numMarkers--;
        // Store that the current decomposition element is a polygon
        decompositionTree[i].isSolidPolygon = true;
        // Store whether the current decomposition element has a graph cycle
        decompositionTree[i].hasGraphCycle = mergedHypergraph.initialNumEdges() > 1;
        // Set the target decomposition element as deleted
        decompositionTree[oldMarkerEdge.targetDecompositionIndex].isDeleted = true;
        // NB: We do not decrease the number of decomposition elements, since we only mark them as deleted and do not remove them from the decomposition tree
        // NB2: There is no need to update the markerEdgeIndex, since we have swapped the last marker edge to the current position and we will look at it in the next iteration
      }
      // Otherwise if the merged hypergraph is NOT a solid polygon, then the split is good
      else
      {
        // Mark the marker edges as coming from a good split
        markerEdge.isFromGoodSplit = true;
        decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].isFromGoodSplit = true;
        // Go to the next marker edge
        markerEdgeIndex++;
      }
    }
  }
}

// Perform BFS to find the nodes that are contracted into the marker node in the decomposition tree
// We then assign the same cluster ID to all nodes that are contracted into the marker node and return this cluster ID
// NB: The marker edge is directed, i.e. we find the nodes of one side of the split that are contracted into the marker node
NodeID CanonicalDecomposition::find_nodes_contracted_to_marker_node_in_decomposition_tree(mt_kahypar::parallel::scalable_vector<ClusterID> &clusterID,
                                                                                          const ClusterIndex sourceDecompositionIndex,
                                                                                          const MarkerEdge markerEdge)
{
  // Store whether a decomposition element has been visited during the BFS
  std::vector<bool> hasVisitedDecompositionElement(numDecompositionElements, false);
  hasVisitedDecompositionElement[sourceDecompositionIndex] = true;
  hasVisitedDecompositionElement[markerEdge.targetDecompositionIndex] = true;
  // Initialize the queue for the BFS
  std::queue<int> queue;
  queue.push(markerEdge.targetDecompositionIndex);
  // Store the first node ID that is contracted into the marker node
  // NB: We use the node ID of the first seen node as the representative node ID for the marker node
  NodeID firstSeenNodeID = std::numeric_limits<NodeID>::max();
  // Perform the BFS
  while (!queue.empty())
  {
    ClusterIndex currentDecompositionIndex = queue.front();
    queue.pop();

    // Make sure that the current decomposition element is not deleted
    assert(!decompositionTree[currentDecompositionIndex].isDeleted);

    // Get the included nodes of the current decomposition element
    std::vector<NodeID> &currentIncludedNodes = decompositionTree[currentDecompositionIndex].includedNodes;

    if (!currentIncludedNodes.empty())
    {
      // Set the first seen node ID if it is not set yet
      if (firstSeenNodeID == std::numeric_limits<NodeID>::max())
        firstSeenNodeID = currentIncludedNodes[0];
      // Loop over all included nodes of the current decomposition element and contract them into the marker node
      for (NodeID nodeID : currentIncludedNodes)
        clusterID[nodeID] = firstSeenNodeID;
    }

    // Loop over all marker edges of the current decomposition element
    for (const MarkerEdge &currentMarkerEdge : decompositionTree[currentDecompositionIndex].markerEdges)
      // Add the target decomposition element to the queue if it has not been visited yet
      if (!hasVisitedDecompositionElement[currentMarkerEdge.targetDecompositionIndex])
      {
        queue.push(currentMarkerEdge.targetDecompositionIndex);
        hasVisitedDecompositionElement[currentMarkerEdge.targetDecompositionIndex] = true;
      }
  }
  return firstSeenNodeID;
}

// Check if the hypergraph is a solid polygon
// A solid polygon is a hypergraph that consists of a GRAPH cycle (including all nodes) where each edge has the same weight A and a hyperedge containing all nodes with weight B.
// NB: If A = 0 (i.e. there is no graph cycle) OR if B = 0 (i.e. there is no hyperedge containing all nodes), then the hypergraph is still a solid polygon.
bool CanonicalDecomposition::is_solid_polygon(StaticHypergraph &hypergraph) const
{
  NodeIndex numNodes = hypergraph.initialNumNodes();
  EdgeIndex numEdges = hypergraph.initialNumEdges();

  // First check if the hypergraph only has a single hyperedge
  if (numEdges == 1)
    // If so, the hypergraph is a solid polygon if the hyperedge contains all nodes
    return (hypergraph.edgeSize(0) == numNodes);

  // Otherwise, the hypergraph is only a solid polygon if it contains a graph cycle (with numNodes edges) and maybe a hyperedge containing all nodes
  if (numEdges != numNodes && numEdges != numNodes + 1)
    return false;

  // Store whether the hypergraph should have a hyperedge containing all nodes
  bool shouldFindHyperedgeContainingAllNodes = (numEdges == numNodes + 1);

  // Make sure that the node degrees fit the requirements of a solid polygon
  for (NodeID nodeID : hypergraph.nodes())
    if (hypergraph.nodeDegree(nodeID) != 2 + (shouldFindHyperedgeContainingAllNodes ? 1 : 0))
      return false;

  // Make sure that the edges of the graph cycle have the same weight and that the hyperedge containing all nodes exists (if necessary)
  EdgeWeight edgeWeightInGraphCycle = std::numeric_limits<EdgeWeight>::max();
  for (EdgeID edgeID : hypergraph.edges())
  {
    // Check if the hyperedge is a graph edge
    if (hypergraph.edgeSize(edgeID) == 2)
    {
      // If so, check if the weight is the same as the weight of the first graph edge
      if (edgeWeightInGraphCycle == std::numeric_limits<EdgeWeight>::max())
        edgeWeightInGraphCycle = hypergraph.edgeWeight(edgeID);
      else if (hypergraph.edgeWeight(edgeID) != edgeWeightInGraphCycle)
        return false;
    }
    // Check if the hyperedge contains all nodes
    else if (hypergraph.edgeSize(edgeID) == numNodes)
    {
      // If so, we do not have a solid polygon if we should not find a hyperedge containing all nodes
      if (!shouldFindHyperedgeContainingAllNodes)
        return false;
      // We no longer need to find a hyperedge containing all nodes (i.e. if there is another, we would conclude that the hypergraph is not a solid polygon)
      shouldFindHyperedgeContainingAllNodes = false;
    }
    else
    {
      // If the hyperedge is neither a graph edge nor a hyperedge containing all nodes (if necessary), then the hypergraph is not a solid polygon
      return false;
    }
  }
  // Make sure that the graph edges form a single cycle by starting at a node and reaching all other nodes exactly once using BFS on the graph edges
  NodeID startNodeID = 0;
  std::vector<bool> isVisited(numNodes, false);
  NodeIndex numVisitedNodes = 0;
  std::queue<int> queue;
  queue.push(startNodeID);
  isVisited[startNodeID] = true;
  numVisitedNodes++;
  while (!queue.empty())
  {
    NodeID nodeID = queue.front();
    queue.pop();

    // Loop over all incident edges of the node
    for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
    {
      // Ignore hyperedges that are not of size two
      if (hypergraph.edgeSize(edgeID) != 2)
        continue;

      // Get the iterator range of the pins of the edge
      auto pinsIterator = hypergraph.pins(edgeID).begin();
      // Get the target of the hyperedge of size two
      NodeID targetID = 0;
      do
      {
        targetID = *(pinsIterator++);
      } while ((targetID == nodeID || !hypergraph.nodeIsEnabled(targetID)) && pinsIterator != hypergraph.pins(edgeID).end());

      // If the target node has not been visited yet, mark it as visited and add it to the queue
      if (!isVisited[targetID])
      {
        isVisited[targetID] = true;
        numVisitedNodes++;
        queue.push(targetID);
      }
    }
  }
  // If we have visited all nodes, then the hypergraph is a solid polygon
  return (numVisitedNodes == numNodes);
}

// Create the hypercactus from the canonical decomposition
// NB: The weight of a node in the hypercactus indicates how many nodes from the original hypergraph it represents
mt_kahypar_hypergraph_t CanonicalDecomposition::create_hypercactus(StaticHypergraph &hypergraph, std::vector<ClusterID> &hypercactusMapping, bool allowOutputFloatEdgeWeights)
{
  // Make sure that we have at least one decomposition element (i.e. there exists a prime decomposition)
  assert(numDecompositionElements > 0 && numDecompositionElements == decompositionTree.size());
  // Make sure that the hypercactus mapping has the same size as the number of nodes in the hypergraph
  assert(hypercactusMapping.size() == hypergraph.initialNumNodes());

  // Get the initial number of nodes in the hypergraph
  const NodeID initialNumNodes = hypergraph.initialNumNodes();

  // Initialize the cluster ID vector (only necessary if a decomposition element is a solid polygon with a graph cycle and we must rebuild it)
  mt_kahypar::parallel::scalable_vector<ClusterID> clusterID(initialNumNodes);

  // Initialize the adjacency list of the hypercactus and the list of the weight of the edges and of the nodes
  // NB: We use upper bounds for the number of hypercactus nodes, edges and pins, since we will only know the exact numbers later
  //     For the number of hypercactus nodes, we use numDecompositionElements to account for the auxiliary star nodes that are created for each prime decomposition element
  //     For the number of hypercactus edges, we know that each orignal node accounts for at most one hypercactus edge and each marker node for at most two
  //     For the number of hypercactus pins, we know that we have the most pins if solid polygons are converted to graph cycles, i.e. we only have hyperedges of size 2 and thus twice as many pins as hyperedges
  //     The space that was allocated too much will be simply ignored by MT-KaHyPar.
  NodeIndex upperBoundNumHypercactusNodes = initialNumNodes + numMarkers + numDecompositionElements;
  EdgeIndex upperBoundNumHypercactusEdges = initialNumNodes + 2 * numMarkers;
  NodeIndex upperBoundNumHypercactusPins = 2 * upperBoundNumHypercactusEdges;

  // If we do not allow floating hyperedge weights, then we might keep the hyperedge containing all nodes in the decomposition element (if it is a solid polygon)
  // In the worst case, each decomposition element is a solid polygon with such a hyperedge containing all nodes
  if (!allowOutputFloatEdgeWeights)
  {
    upperBoundNumHypercactusEdges += numDecompositionElements;
    upperBoundNumHypercactusPins += initialNumNodes;
  }

  std::unique_ptr<size_t[]> edgesAdjIndices = std::make_unique<size_t[]>(upperBoundNumHypercactusEdges + 1);
  std::unique_ptr<mt_kahypar_hyperedge_id_t[]> edgesAdjList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(upperBoundNumHypercactusPins);
  std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edgesWeights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(upperBoundNumHypercactusEdges);
  std::unique_ptr<mt_kahypar_hypernode_weight_t[]> nodesWeights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(upperBoundNumHypercactusNodes);

  // Store the number of hypercactus nodes, edges and pins
  // NB: We know that the hypercactus has at least one node per marker
  NodeIndex numHypercactusNodes = 0;
  EdgeIndex numHypercactusEdges = 0;
  NodeIndex numHypercactusPins = 0;

  // Store whether we already have seen the marker node
  std::vector<bool> hasSeenMarkerNode(initialNumMarkers, false);

  // Loop over all decomposition elements
  for (ClusterIndex i = 0; i < numDecompositionElements; ++i)
  {
    // Skip the decomposition element if it is deleted
    if (decompositionTree[i].isDeleted)
      continue;

    // Check if the decomposition element is a solid polygon
    if (decompositionTree[i].isSolidPolygon)
    {
      // Check if the decomposition element has a graph cycle
      if (decompositionTree[i].hasGraphCycle)
      {
        // If so, we must rebuild the hypergraph of the solid polygon to determine the order of the nodes in the cycle
        // Reset the cluster ID vector to contain 0, 1, 2, ..., initialNumNodes-1
        std::iota(clusterID.begin(), clusterID.end(), 0);

        // Loop over the marker edges of the decomposition element to find which original nodes are contracted into the marker node via BFS
        // NB: We also store the ID of the marker node (i.e. the first node ID that is contracted into the marker node) in the markerNodeIDs vector
        std::vector<NodeID> markerNodeIDs(decompositionTree[i].markerEdges.size());
        for (EdgeIndex j = 0; j < decompositionTree[i].markerEdges.size(); ++j)
          markerNodeIDs[j] = find_nodes_contracted_to_marker_node_in_decomposition_tree(clusterID, i, decompositionTree[i].markerEdges[j]);

        // Contract the hypergraph with the cluster ID vector to obtain the hypergraph of the solid polygon
        StaticHypergraph solidPolygon = hypergraph.contract(clusterID);
        NodeIndex numNodesInSolidPolygon = solidPolygon.initialNumNodes();
        EdgeIndex numEdgesInSolidPolygon = solidPolygon.initialNumEdges();

        // Store for each node in the graph cycle of the solid polygon either the index of the original/included node or the index of the marker node it represents
        // The second value of the pair indicates whether it represents a marker node (true) or an original/included node (false)
        // This is necessary to transfer the order of the nodes in the cycle to the original/included nodes and the marker nodes
        std::vector<std::pair<NodeIndex, bool>> includedNodeOrMarkerNodeForCycleNode(numNodesInSolidPolygon);
        for (NodeIndex j = 0; j < decompositionTree[i].includedNodes.size(); ++j)
          includedNodeOrMarkerNodeForCycleNode[clusterID[decompositionTree[i].includedNodes[j]]] = {j, false};
        for (EdgeIndex j = 0; j < decompositionTree[i].markerEdges.size(); ++j)
          includedNodeOrMarkerNodeForCycleNode[clusterID[markerNodeIDs[j]]] = {j, true};

        // Start with the first node in the cycle
        NodeID currentNodeID = 0;
        // Store the new node ID of the start of the cycle (so that we can close the cycle later)
        NodeID startNewNodeID = std::numeric_limits<NodeID>::max();
        // Store the previous node in the graph cycle, both the node ID in the solid polygon and the new node ID in the hypercactus
        NodeID previousNodeID = std::numeric_limits<NodeID>::max();
        NodeID previousNewNodeID = std::numeric_limits<NodeID>::max();
        // Stop as soon as we reach the start node for the second time
        while (previousNodeID == std::numeric_limits<NodeID>::max() || currentNodeID != 0)
        {
          NodeID newNodeID;
          // Check if the node represents a marker node
          if (includedNodeOrMarkerNodeForCycleNode[currentNodeID].second)
          {
            newNodeID = get_new_node_id_of_marker_in_hypercactus(decompositionTree[i].markerEdges[includedNodeOrMarkerNodeForCycleNode[currentNodeID].first], hasSeenMarkerNode, numHypercactusNodes);
            nodesWeights[newNodeID] = 0;
          }
          else
          {
            newNodeID = numHypercactusNodes++;
            hypercactusMapping[decompositionTree[i].includedNodes[includedNodeOrMarkerNodeForCycleNode[currentNodeID].first]] = newNodeID;
            nodesWeights[newNodeID] = 1;
          }

          // Check if we have a previous node in the graph cycle
          if (previousNewNodeID != std::numeric_limits<NodeID>::max())
          {
            // Add the edge between the previous node and the current node
            edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
            edgesWeights[numHypercactusEdges] = minEdgeCut / 2;
            numHypercactusEdges++;
            edgesAdjList[numHypercactusPins++] = previousNewNodeID;
            edgesAdjList[numHypercactusPins++] = newNodeID;
          }
          // Otherwise store the new node ID of the start of the cycle
          else
            startNewNodeID = newNodeID;

          // Loop over all incident edges of the node
          for (EdgeID edgeID : solidPolygon.incidentEdges(currentNodeID))
          {
            // Ignore hyperedges that are not of size two
            if (solidPolygon.edgeSize(edgeID) != 2)
              continue;

            // Get the iterator range of the pins of the edge
            auto pinsIterator = solidPolygon.pins(edgeID).begin();
            // Get the target of the hyperedge of size two
            NodeID targetID = 0;
            do
            {
              targetID = *(pinsIterator++);
            } while ((targetID == currentNodeID || !solidPolygon.nodeIsEnabled(targetID)) && pinsIterator != solidPolygon.pins(edgeID).end());

            // If the target node is not the previous node in the graph cycle, we can continue
            if (targetID != previousNodeID)
            {
              previousNodeID = currentNodeID;
              previousNewNodeID = newNodeID;
              currentNodeID = targetID;
              break;
            }
          }
        }

        // Close the cycle by adding the edge between the last node and the start node
        edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
        edgesWeights[numHypercactusEdges] = minEdgeCut / 2;
        numHypercactusEdges++;
        edgesAdjList[numHypercactusPins++] = previousNewNodeID;
        edgesAdjList[numHypercactusPins++] = startNewNodeID;

        // Check if the mincut is odd and if there exists a hyperedge containing all nodes
        // If so, the hyperedge containing all nodes has an odd weight and must be kept to avoid floating point hyperedge weights
        // We assign a weight of 1 to the hyperedge containing all nodes and and a weight of minEdgeCut / 2 (integer division) to the cycle edges
        if (!allowOutputFloatEdgeWeights && minEdgeCut % 2 == 1 && numEdgesInSolidPolygon == numNodesInSolidPolygon + 1)
        {
          edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
          edgesWeights[numHypercactusEdges] = 1;
          numHypercactusEdges++;
          NodeIndex oldNumHypercactusPins = numHypercactusPins;
          for (NodeIndex i = 0; i < numNodesInSolidPolygon; ++i)
          {
            // Retrieve the new node ID by looking at the target of the just added graph edges
            NodeID newNodeID = edgesAdjList[oldNumHypercactusPins - 1 - 2 * i];
            edgesAdjList[numHypercactusPins++] = newNodeID;
          }
        }
      }
      else
      {
        // Otherwise, we have a single hyperedge that contains all nodes and marker nodes
        edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
        edgesWeights[numHypercactusEdges] = minEdgeCut;
        numHypercactusEdges++;
        // Assign the original nodes as pins of the hyperedge
        for (NodeID nodeID : decompositionTree[i].includedNodes)
        {
          NodeID newNodeID = numHypercactusNodes++;
          hypercactusMapping[nodeID] = newNodeID;
          nodesWeights[newNodeID] = 1;
          edgesAdjList[numHypercactusPins++] = newNodeID;
        }
        // Assign the marker nodes as pins of the hyperedge
        // NB: Each marker node persists in the hypercactus
        for (MarkerEdge &markerEdge : decompositionTree[i].markerEdges)
        {
          NodeID newNodeID = get_new_node_id_of_marker_in_hypercactus(markerEdge, hasSeenMarkerNode, numHypercactusNodes);
          nodesWeights[newNodeID] = 0;
          edgesAdjList[numHypercactusPins++] = newNodeID;
        }
      }
    }
    // Otherwise, the decomposition element is prime, i.e. we must transform it into a star structure
    else
    {
      // Observation: Since the decomposition element is prime, it must have at least one node with a trivial mincut.
      //              If numDecompositionElements > 1, then each decomposition element must have at least one marker node.
      //              If numDecompositionElements == 1, then there is no split, which means that the mincut must be trivial.

      // Create the auxiliary star node for the decomposition element
      NodeID auxiliaryStarNodeID = numHypercactusNodes++;
      nodesWeights[auxiliaryStarNodeID] = 0;
      // Assign the original nodes as pins of the hyperedge
      for (NodeID nodeID : decompositionTree[i].includedNodes)
      {
        // Compute the cut when isolating the node
        CutValue isolationCut = 0;
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
          if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
            isolationCut += hypergraph.edgeWeight(edgeID);
        // If the isolating cut is equal to the mincut, then the original node persists in the hypercactus
        if (isolationCut == minEdgeCut)
        {
          NodeID newNodeID = numHypercactusNodes++;
          hypercactusMapping[nodeID] = newNodeID;
          nodesWeights[newNodeID] = 1;
          // Add the edge between the original node and the auxiliary star node
          edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
          edgesWeights[numHypercactusEdges] = minEdgeCut;
          numHypercactusEdges++;
          edgesAdjList[numHypercactusPins++] = auxiliaryStarNodeID;
          edgesAdjList[numHypercactusPins++] = newNodeID;
        }
        // Otherwise, the original node is merged into the auxiliary star node
        else
        {
          hypercactusMapping[nodeID] = auxiliaryStarNodeID;
          nodesWeights[auxiliaryStarNodeID] += 1;
        }
      }

      // Assign the marker nodes as pins of the hyperedge
      // NB: Each marker node persists in the hypercactus
      for (MarkerEdge &markerEdge : decompositionTree[i].markerEdges)
      {
        NodeID newNodeID = get_new_node_id_of_marker_in_hypercactus(markerEdge, hasSeenMarkerNode, numHypercactusNodes);
        nodesWeights[newNodeID] = 0;
        // Add the edge between the marker node and the auxiliary star node
        edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;
        edgesWeights[numHypercactusEdges] = minEdgeCut;
        numHypercactusEdges++;
        edgesAdjList[numHypercactusPins++] = auxiliaryStarNodeID;
        edgesAdjList[numHypercactusPins++] = newNodeID;
      }
    }
  }

  // Add the sentinel that points to the index behind the last pin
  edgesAdjIndices[numHypercactusEdges] = numHypercactusPins;

  // Create the hypercactus
  return mt_kahypar_create_hypergraph(DETERMINISTIC,
                                      numHypercactusNodes,
                                      numHypercactusEdges,
                                      edgesAdjIndices.get(),
                                      edgesAdjList.get(),
                                      edgesWeights.get(),
                                      nodesWeights.get());
}

// Return for a marker edge the new node ID of the marker node in the hypercactus
NodeID CanonicalDecomposition::get_new_node_id_of_marker_in_hypercactus(MarkerEdge &markerEdge, std::vector<bool> &hasSeenMarkerNode, NodeIndex &numHypercactusNodes)
{
  if (hasSeenMarkerNode[markerEdge.markerID])
    return markerEdge.representativeNodeID;

  NodeID newNodeID = numHypercactusNodes++;
  hasSeenMarkerNode[markerEdge.markerID] = true;
  markerEdge.representativeNodeID = newNodeID;
  decompositionTree[markerEdge.targetDecompositionIndex].markerEdges[markerEdge.markerEdgeIndexAtTarget].representativeNodeID = newNodeID;
  return newNodeID;
}

// Print the decomposition tree
void CanonicalDecomposition::print_decomposition_tree(const bool isCanonical) const
{
  std::cout << "===================================================" << std::endl;
  if (isCanonical)
    std::cout << "############# CANONICAL DECOMPOSITION #############" << std::endl;
  else
    std::cout << "############### PRIME DECOMPOSITION ###############" << std::endl;
  std::cout << "===================================================" << std::endl;
  for (ClusterIndex i = 0; i < numDecompositionElements; ++i)
  {
    if (decompositionTree[i].isDeleted)
      continue;
    if (i > 0)
      std::cout << std::endl;
    std::cout << "decomposition_element \t\t\t" << i + 1 << std::endl;
    if (verbose)
    {
      std::cout << "included_nodes \t\t\t\t";
      for (NodeID nodeID : decompositionTree[i].includedNodes)
        std::cout << nodeID + 1 << " ";
      std::cout << std::endl;
    }
    else
      std::cout << "num_included_nodes \t\t\t" << decompositionTree[i].includedNodes.size() << std::endl;
    std::cout << "marker_edges \t\t\t\t";
    if (decompositionTree[i].markerEdges.empty())
      std::cout << "/" << std::endl;
    else
    {
      for (const MarkerEdge &markerEdge : decompositionTree[i].markerEdges)
        if (verbose)
          std::cout << "(" << markerEdge.markerID + 1 << ", " << markerEdge.representativeNodeID + 1 << ", "
                    << markerEdge.targetDecompositionIndex + 1 << ", " << markerEdge.markerEdgeIndexAtTarget + 1 << ", "
                    << markerEdge.isFromGoodSplit << ") ";
        else
          std::cout << markerEdge.markerID + 1 << " ";
      std::cout << std::endl;
    }
  }
}