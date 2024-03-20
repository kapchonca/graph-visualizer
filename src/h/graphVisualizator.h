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
  // Return image data
  std::vector<std::vector<int8_t>> GetPixels(
      const std::pair<int, int>& dimensions);
  // Moves vertices to optimal coordinates for image size optimization and returns their coords
  std::pair<int, int> GetCoordinates();

 private:
  // Round coordinates to the closest intereger values
  void RoundCoords(Vertex* v);

  // Add numbers to the image
  void DrawNumbers();

  // Draw a circle
  void DrawVertex(Vertex* v);

  // Draw all edges connected to the vertex
  void DrawEdges(Vertex* v);

  const int kEdgeLen = 40;
  std::vector<std::vector<int8_t>> image_data_;

  std::vector<std::shared_ptr<Vertex>> vertices_;
};
#endif  //LAB4_H_GRAPHVISUALIZATOR_H_
