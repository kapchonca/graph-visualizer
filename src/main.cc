#include "h/bmpHandler.h"
#include "h/graphHandler.h"
#include "h/graphVisualizator.h"

int main() {
  const std::string filename = "sun";
  Graph graph("examples/" + filename + ".txt");
  graph.DrawGraph(filename);
  return 0;
}
