#include "h/bmpHandler.h"
#include "h/graphHandler.h"

int main() {
  // const int width = 800;
  // const int height = 600;

  // BMPWriter bmpWriter("output.bmp", width, height);

  // for (int x = 0; x < width / 2; x += 1) {
  //   for (int y = 0; y < height / 2; y += 1) {
  //     bmpWriter.setPixel(x, y, 0, 120, 0);
  //   }
  // }

  // bmpWriter.save();

  const std::string file = "graph";
  Graph graph(file);
  graph.RandomLayout();
  for (auto v : graph.vertices_) {
    v->distances = graph.BFS(v, graph.vertex_num_);
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
  }
  auto c = graph.kCenter(3);
  graph.FormNeighbourhood(c);

  for (auto n : c) {
    std::cout << "center ind: " << n->kVertInd + 1 << "neighbours: ";
    for (auto y : n->neighboorhood) {
      std::cout << y->kVertInd + 1 << ' ';
    }
    std::cout << "delta:" << graph.CalculateDelta(n) << std::endl;
  }
  // std::cout << graph.CalculateEnergy();
  return 0;
}
