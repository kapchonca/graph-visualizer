#include "h/bmpHandler.h"
#include "h/graphHandler.h"
#include "h/graphVisualizator.h"

int main() {
  const std::string filename = "binary_tree";
  Graph graph("examples/" + filename + ".txt");
  graph.GlobalLayout();
  Visualizator vis(graph.vertices_);
  std::pair<int, int> dimensions = vis.MoveCoordinates();
  BMPWriter write;
  write.SetDimensions(dimensions.first, dimensions.second);
  write.imageData_ = vis.GetPixels(dimensions);
  write.save("output/" + filename + ".bmp");

  // BMPWriter bmpWriter;
  // bmpWriter.SetDimensions(8, 16);
  // std::vector<std::vector<int8_t>> pixels = bmpWriter.Read("./sprites/2.bmp");
  // bmpWriter.imageData_ = pixels;
  // bmpWriter.save("huy.bmp");

  // std::cout << dimensions.first << dimensions.second;
  // bmpWriter.SetDimensions(dimensions.first, dimensions.second);
  // for (int x = 0; x < dimensions.first / 2; x += 1) {
  //   for (int y = 0; y < dimensions.second / 2; y += 1) {
  //     bmpWriter.setPixel(x, y, 0, 120, 0);
  //   }
  // }
  // bmpWriter.save("output.bmp");

  // BMPWriter reader;
  // std::vector<std::vector<int8_t>> pixels =
  //     reader.Read("../src/sprites/circle.bmp");

  // const std::string file = "graph";
  // Graph graph(file);
  // graph.GlobalLayout();
  // auto c = graph.MoveCoordinates();
  // std::cout << "width: " << c.first << " height: " << c.second << std::endl;
  // for (auto v : graph.vertices_) {
  //   std::cout << "vertex index: " << v->kVertInd + 1 << " x: " << v->x
  //             << " y: " << v->y << std::endl;
  // }

  // graph.RandomLayout();
  // for (auto v : graph.vertices_) {
  //   v->distances = graph.BFS(v, graph.vertex_num_);
  // std::cout << v->x_coord << ' ' << v->y_coord << std::endl;
  // std::cout << "main vertex: " << v->kVertInd + 1 << std::endl;
  // for (const auto& n : v->distances) {
  //   std::cout << "vertex: " << n.first->kVertInd + 1
  //             << "distance: " << n.second;
  // }
  // std::cout << std::endl;
  // auto c = graph.kCenter(3);

  // for (auto x : c) {
  //   std::cout << x->kVertInd + 1;
  // }
  // }
  // auto c = graph.kCenter(3);
  // graph.FormNeighbourhood(c);

  // for (auto n : c) {
  //   std::cout << "center ind: " << n->kVertInd + 1 << "neighbours: ";
  //   for (auto y : n->neighboorhood) {
  //     std::cout << y->kVertInd + 1 << ' ';
  //   }
  //   std::cout << "delta:" << graph.CalculateDelta(n) << std::endl;
  // }
  // std::cout << graph.CalculateEnergy();
  return 0;
}
