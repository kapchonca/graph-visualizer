#include <chrono>
#include "h/bmpHandler.h"
#include "h/graphHandler.h"
#include "h/graphVisualizator.h"

int main() {
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  const std::string filename = "100v";
  Graph graph("examples/" + filename + ".txt");
  graph.DrawGraph(filename);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time difference = "
            << std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                     begin)
                   .count()
            << "[µs]" << std::endl;
  return 0;
}
