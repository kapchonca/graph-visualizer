#ifndef LAB4_H_GRAPHVISUALIZATOR_H_
#define LAB4_H_GRAPHVISUALIZATOR_H_

#include <limits>
#include <memory>
#include <vector>

#include "graphHandler.h"

class Visualizator {
 public:
  Visualizator(std::vector<std::shared_ptr<Vertex>> vertices)
      : vertices_(vertices){};
  std::vector<std::vector<int8_t>> GetPixels(std::pair<int, int> dimensions);
  std::vector<std::vector<int8_t>> DrawNumbers(
      std::vector<std::vector<int8_t>> image_data);
  void RoundCoords(Vertex* v);
  std::pair<int, int> MoveCoordinates();

 private:
  const int kEdgeLen = 40;

  std::vector<std::shared_ptr<Vertex>> vertices_;
};
#endif  //LAB4_H_GRAPHVISUALIZATOR_H_
