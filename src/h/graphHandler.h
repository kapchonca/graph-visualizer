#ifndef LAB4_H_GRAPHHANDLER_H_
#define LAB4_H_GRAPHHANDLER_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bmpHandler.h"
#include "graphMath.h"
#include "graphVisualizator.h"

using vert_p = std::shared_ptr<Vertex>;

// Represents a graph
class Graph {
 public:
  // Constructor reads data from a file and initializes the graph.
  Graph(const std::string& file_path);

  void DrawGraph(const std::string& filename);

 private:
  // Purpose: Perform Breadth-First Search (BFS) on the graph starting from a specified vertex.
  // Input: the starting vertex for BFS.
  // Output: Updates the distances field of the starting vertex with BFS traversal information.
  std::unordered_map<Vertex*, int> BFS(Vertex* start, int max_depth) const;

  // Function to find a k-center approximation of the graph using a greedy approach.
  // Input: k - The number of centers to select
  // Output: Set of vertices representing the selected k-centers
  std::unordered_set<std::shared_ptr<Vertex>> kCenter(long unsigned k) const;

  // Assigns random x- and y-coordinates to verticeso on a plane
  void RandomLayout() const;

  // Finds all neighboors up to a certain depth
  void FormNeighbourhood(Vertex* chip_chrome_and_the_mono_tones, int radius);

  // Сalculates the best location of vertices in small groups
  void LocalLayout(Vertex* p, int radius);

  // Сalculates the best location of vertices in the entire graph
  void GlobalLayout();

  std::vector<std::shared_ptr<Vertex>> vertices_;
  GraphMath graph_math;
  BMPWriter write;
  int vertex_num_;
  int edge_num_;
  const int kEdgeLen = 40;
  const int kLocalRadius = 7;
  const int kIterations = 4;
  const int kRatio = 3;
  const int kMinSize = 2;
};

#endif  // LAB4_H_GRAPHHANDLER_H_