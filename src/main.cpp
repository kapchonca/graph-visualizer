#include "h/bmpHandler.h"
#include "h/graphHandler.h"

int main() {
  const int width = 800;
  const int height = 600;

  BMPWriter bmpWriter("output.bmp", width, height);

  for (int x = 0; x < width / 2; x += 1) {
    for (int y = 0; y < height / 2; y += 1) {
      bmpWriter.setPixel(x, y, 0, 120, 0);
    }
  }

  bmpWriter.save();

  // const std::string file = "graph";
  // Graph graph(file);
  // for (auto v : graph.vertices_) {
  //   graph.BFS(v);
  //   std::cout << "main vertex: " << v->kVertInd + 1 << std::endl;
  //   for (const auto& n : v->distances) {
  //     std::cout << "depth: " << n.first << " neighbours: {";
  //     for (auto k : n.second) {
  //       std::cout << k.lock()->kVertInd + 1 << ',';
  //     }
  //     std::cout << "}\t";
  //   }
  //   std::cout << std::endl;
  // }
  return 0;
}
