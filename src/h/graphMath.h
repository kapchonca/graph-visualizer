#ifndef LAB4_H_GRAPHMATH_H_
#define LAB4_H_GRAPHMATH_H_

#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Structure representing a vertex in the graph that stores adjacent vertices and distances to them
struct Vertex {
  Vertex(const int ind) : kVertInd(ind){};
  const int kVertInd;
  std::vector<Vertex*> adjacent_vertices;
  std::unordered_map<Vertex*, int> distances;
  double x, y;
  std::vector<Vertex*> neighboorhood;
};

class GraphMath {
 public:
  // Distance between two vertices on a plane
  double EuclideanDistance(Vertex* v, Vertex* u, float power);

  // Calculates derivative with respect to given coordinate of energy function
  double OneVariableDerivative(Vertex* parameter, char with_respect);

  // Calculates an auxiliary value that helps to identify the vertex that is most deviated from its optimal placement
  double CalculateDelta(Vertex* parameter);

  // Calculates a coefficient for linear equations with respect to given coordinates
  double TwoVariablesDerivative(Vertex* parameter, char with_respect1,
                                char with_respect2);

  // Linear equation to adjust placement of a vertex
  void SolveLinearEquations(Vertex* p);

 private:
  const int kEdgeLen = 40;
};

#endif  // LAB4_H_GRAPHMATH_H_