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

std::unordered_map<Vertex*, int> Graph::BFS(Vertex* start,
                                            int max_depth) const {
  std::vector<bool> visited(vertex_num_, false);
  std::queue<std::pair<Vertex*, int>> q;
  std::unordered_map<Vertex*, int> traversal;

  visited[start->kVertInd] = true;
  q.push({start, 0});

  while (!q.empty() && q.front().second < max_depth) {
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
  traversal.erase(start);
  return traversal;
}

// Greedy 2-Approximation Algorithm for k-Center
std::unordered_set<std::shared_ptr<Vertex>> Graph::kCenter(
    long unsigned k) const {
  std::unordered_set<std::shared_ptr<Vertex>> centers;
  std::unordered_set<std::shared_ptr<Vertex>> remainingVertices;

  remainingVertices.insert(vertices_.begin(), vertices_.end());
  // Select the first center arbitrarily
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> distribution(
      0, vertex_num_ - 1);

  auto firstVertex = vertices_[distribution(rng)];
  centers.insert(firstVertex);
  // std::cout << "inserted vertex " << firstVertex->kVertInd + 1 << std::endl;
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
    // std::cout << "inserted vertex " << farthestVertex->kVertInd + 1
    //           << std::endl;
    remainingVertices.erase(farthestVertex);
  }

  return centers;
}

void Graph::RandomLayout() const {
  std::random_device dvc;
  std::mt19937 rng(dvc());
  std::uniform_real_distribution<> distribution(0, vertex_num_ * kEdgeLen);

  for (auto v : vertices_) {
    v->x = distribution(rng);
    v->y = distribution(rng);
    // std::cout << "vertex index: " << v->kVertInd + 1
    //           << " initial placement: " << v->x << ' ' << v->y << std::endl;
  }
}

void Graph::FormNeighbourhood(Vertex* center, int radius) {
  center->neighboorhood.clear();
  std::unordered_map<Vertex*, int> traversal = BFS(center, radius);
  for (auto pair : traversal) {
    center->neighboorhood.insert(pair.first);
  }
}

void Graph::LocalLayout(Vertex* p, int radius) {
  double max_delta;
  Vertex* v_to_adjust;
  for (unsigned long i = 0; i < kIterations * vertices_.size(); ++i) {
    max_delta = 0;
    for (auto v : p->neighboorhood) {
      FormNeighbourhood(v, radius);
      double delta = graph_math.CalculateDelta(v);
      if (delta > max_delta) {
        max_delta = delta;
        v_to_adjust = v;
      }
    }
    graph_math.SolveLinearEquations(v_to_adjust);
  }
}

void Graph::GlobalLayout() {

  for (auto v : vertices_) {
    v->distances = BFS(v.get(), vertex_num_);
  }

  RandomLayout();

  int radius;
  int k = kMinSize;

  // GraphMath graph_math(vertices_);

  while (k <= vertex_num_) {
    auto c = kCenter(k);

    int maxDistance = -1;

    // Find the vertex farthest from the current set of centers
    // graph diameter
    for (auto v : c) {
      int minDistance = std::numeric_limits<int>::max();

      for (auto n : c) {
        if (n != v) {
          minDistance = std::min(minDistance, v->distances[n.get()]);
        }
      }

      if (minDistance > maxDistance) {
        maxDistance = minDistance;
      }
    }
    radius = kLocalRadius * maxDistance;
    for (auto center : c) {
      FormNeighbourhood(center.get(), radius);
      LocalLayout(center.get(), radius);
    }
    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<> distribution(0, 1);
    for (auto v : vertices_) {
      v->x += distribution(gen);
      v->y += distribution(gen);
    }
    k *= kRatio;
  }
}
