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
  std::queue<std::pair<Vertex*, int>> q;
  std::unordered_map<Vertex*, int> traversal;

  visited[start->kVertInd] = true;
  q.push({start.get(), 0});

  while (!q.empty()) {
    Vertex* currentVertex = q.front().first;
    int currentDepth = q.front().second;
    q.pop();
    traversal[currentVertex] = currentDepth;

    // Explore neighbors and update distances.
    for (std::weak_ptr<Vertex> neighbor : currentVertex->adjacent_vertices) {
      if (!visited[neighbor.lock()->kVertInd]) {
        visited[neighbor.lock()->kVertInd] = true;
        q.push({neighbor.lock().get(), currentDepth + 1});
      }
    }
  }
  traversal.erase(start.get());
  start->distances = traversal;
}

// Greedy 2-Approximation Algorithm for k-Center
std::unordered_set<std::shared_ptr<Vertex>> Graph::kCenter(
    long unsigned k) const {
  std::unordered_set<std::shared_ptr<Vertex>> centers;
  std::unordered_set<std::shared_ptr<Vertex>> remainingVertices;

  // for (int i = 0; i < vertices; ++i) {
  //     remainingVertices.insert(i);
  // }
  remainingVertices.insert(vertices_.begin(), vertices_.end());
  // Select the first center arbitrarily
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> distribution(
      0, vertex_num_ - 1);

  auto firstVertex = vertices_[distribution(rng)];
  centers.insert(firstVertex);
  std::cout << "inserted vertex " << firstVertex->kVertInd + 1 << std::endl;
  remainingVertices.erase(firstVertex);

  while (centers.size() < k) {
    std::shared_ptr<Vertex> farthestVertex = nullptr;
    int maxDistance = -1;

    // Find the vertex farthest from the current set of centers
    for (auto v : remainingVertices) {
      int minDistance = std::numeric_limits<int>::max();

      for (auto center : centers) {
        minDistance = std::min(minDistance, v->distances[center.get()]);
      }

      if (minDistance > maxDistance) {
        maxDistance = minDistance;
        farthestVertex = v;
      }
    }

    // Add the farthest vertex as a new center
    centers.insert(farthestVertex);
    std::cout << "inserted vertex " << farthestVertex->kVertInd + 1
              << std::endl;
    remainingVertices.erase(farthestVertex);
  }

  return centers;
}