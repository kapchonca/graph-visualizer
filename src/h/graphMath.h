#ifndef LAB4_H_GRAPHMATH_H_
#define LAB4_H_GRAPHMATH_H_

#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
// #include "graphHandler.h"

// Structure representing a vertex in the graph that stores adjacent vertices and distances to them
struct Vertex {
  Vertex(const int ind) : kVertInd(ind){};
  const int kVertInd;
  std::vector<std::weak_ptr<Vertex>> adjacent_vertices;
  std::unordered_map<Vertex*, int> distances;
  double x, y;
  std::unordered_set<Vertex*> neighboorhood;
};

class GraphMath {
 public:
  //   GraphMath(std::vector<std::shared_ptr<Vertex>> vertices)
  //       : vertices_(vertices){};
  double EuclideanDistance(Vertex* v, Vertex* u, float power);
  double OneVariableDerivative(Vertex* parameter, char with_respect);
  double CalculateDelta(Vertex* parameter);
  double TwoVariablesDerivative(Vertex* parameter, char with_respect1,
                                char with_respect2);
  void SolveLinearEquations(Vertex* p);

 private:
  const int kEdgeLen = 40;
  //   std::vector<std::shared_ptr<Vertex>> vertices_;
};

#endif  // LAB4_H_GRAPHMATH_H_