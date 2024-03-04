#include "h/graphHandler.h"

Graph::Graph(const std::string file_path) {
  std::ifstream input_file(file_path);

  if (!input_file.is_open()) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }

  input_file >> vertex_num_ >> edge_num_;

  // Resize the vertices vector and initialize vertices.
  vertices_.resize(vertex_num_, nullptr);
  for (int i = 0; i < vertex_num_; ++i) {
    std::shared_ptr<Vertex> v = std::make_shared<Vertex>(i);
    vertices_[i] = v;
  }

  // Read edges and establish connections between vertices.
  int num1, num2;
  while (input_file >> num1 >> num2) {
    vertices_[num1 - 1]->adjacent_vertices.push_back(vertices_[num2 - 1]);
    vertices_[num2 - 1]->adjacent_vertices.push_back(vertices_[num1 - 1]);
  }

  input_file.close();
}

void Graph::PrintData() const {
  for (std::shared_ptr<Vertex> v : vertices_) {
    std::cout << v->kVertInd << std::endl;
    for (std::weak_ptr<Vertex> adj : v->adjacent_vertices) {
      std::cout << adj.lock()->kVertInd << ' ';
    }
    std::cout << std::endl;
  }
}

void Graph::BFS(std::shared_ptr<Vertex> start) const {
  std::vector<bool> visited(vertex_num_, false);
  std::queue<std::pair<std::shared_ptr<Vertex>, int>> q;
  std::unordered_map<int, std::vector<std::weak_ptr<Vertex>>> traversal;

  visited[start->kVertInd] = true;
  q.push({start, 0});

  while (!q.empty()) {
    std::shared_ptr<Vertex> currentVertex = q.front().first;
    int currentDepth = q.front().second;
    q.pop();
    traversal[currentDepth].push_back(currentVertex);

    // Explore neighbors and update distances.
    for (std::weak_ptr<Vertex> neighbor : currentVertex->adjacent_vertices) {
      if (!visited[neighbor.lock()->kVertInd]) {
        visited[neighbor.lock()->kVertInd] = true;
        q.push({neighbor.lock(), currentDepth + 1});
      }
    }
  }

  start->distances = traversal;
}
