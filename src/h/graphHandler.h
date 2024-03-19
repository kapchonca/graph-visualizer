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

// Structure representing a vertex in the graph that stores adjacent vertices and distances to them
struct Vertex {
  Vertex(const int ind) : kVertInd(ind){};
  const int kVertInd;
  std::vector<std::weak_ptr<Vertex>> adjacent_vertices;
  std::unordered_map<Vertex*, int> distances;
  double x, y;
  std::unordered_set<Vertex*> neighboorhood;
};

// Represents a graph
class Graph {
 public:
  // Constructor reads data from a file and initializes the graph.
  Graph(const std::string file_path);

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

  // Calculates derivative with respect to given coordinate of energy function
  double OneVariableDerivative(Vertex* parameter, char with_respect);

  // Calculates an auxiliary value that helps to identify the vertex that is most deviated from its optimal placement
  double CalculateDelta(Vertex* parameter);

  // Calculates a coefficient for linear equations with respect to given coordinates
  double TwoVariablesDerivative(Vertex* parameter, char with_respect1,
                                char with_respect2);

  // Distance between two vertices on a plane
  double EuclideanDistance(Vertex* v, Vertex* u, float power);

  // Linear equation to adjust placement of a vertex
  void SolveLinearEquations(Vertex* p);

  // Сalculates the best location of vertices in small groups
  void LocalLayout(Vertex* p, int radius);

  // Сalculates the best location of vertices in the entire graph
  void GlobalLayout();

  // Moves vertices to optimal coordinates for image size optimization
  std::pair<int, int> MoveCoordinates();

  void RoundCoords(Vertex* v);

  std::vector<std::vector<int8_t>> GetPixels(std::pair<int, int> dimensions);

  std::vector<std::vector<int8_t>> DrawNumbers(
      std::vector<std::vector<int8_t>> image_data);

  std::vector<std::shared_ptr<Vertex>>
      vertices_;  // public for now for easier debugging
  int vertex_num_;

 private:
  int edge_num_;
  const int kEdgeLen = 40;
  const int kLocalRadius = 7;
  const int kIterations = 4;
  const int kRatio = 3;
  const int kMinSize = 2;
};

#endif  // LAB4_H_GRAPHHANDLER_H_