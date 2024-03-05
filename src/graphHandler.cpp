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

std::unordered_map<Vertex*, int> Graph::BFS(std::shared_ptr<Vertex> start,
                                            int max_depth) const {
  std::vector<bool> visited(vertex_num_, false);
  std::queue<std::pair<Vertex*, int>> q;
  std::unordered_map<Vertex*, int> traversal;

  visited[start->kVertInd] = true;
  q.push({start.get(), 0});

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
  traversal.erase(start.get());
  return traversal;
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

void Graph::RandomLayout() const {
  // int estimated_img_size = 30 * vertex_num_;
  std::random_device dvc;
  std::mt19937 rng(dvc());
  std::uniform_real_distribution<> distribution(0, vertex_num_ * 30);

  for (auto v : vertices_) {
    v->x_coord = distribution(rng);
    v->y_coord = distribution(rng);
  }
}

// double Graph::CalculateEnergy() {
//   double energy = 0;
//   for (auto v : vertices_) {  // clang-format off
//     for (auto neighbor : v->neighboorhood) {
//       energy += pow(std::sqrt(pow(v->x_coord - neighbor->x_coord, 2) +
//                               pow(v->y_coord - neighbor->y_coord, 2)) -
//                         30 * v->distances[neighbor], 2) /
//                         v->distances[neighbor];
//     }  // clang-format on
//   }
//   return energy;
// }

void Graph::FormNeighbourhood(
    std::unordered_set<std::shared_ptr<Vertex>> centers) {
  for (auto center : centers) {
    std::unordered_map<Vertex*, int> traversal = BFS(center, k);
    for (auto pair : traversal) {
      center->neighboorhood.insert(pair.first);
    }
  }
}

double Graph::CalculateXDerivative(std::shared_ptr<Vertex> parameter) {
  double first_derivative;
  for (auto neighbor : parameter->neighboorhood) {
    first_derivative +=
        ((parameter->x_coord - neighbor->x_coord) -
         ((30 * parameter->distances[neighbor]) *
              (parameter->x_coord - neighbor->x_coord) /
              std::sqrt(std::pow(parameter->x_coord - neighbor->x_coord, 2)) +
          pow(parameter->y_coord - neighbor->y_coord, 2))) /
        parameter->distances[neighbor];
  }
  return first_derivative;
}

double Graph::CalculateYDerivative(std::shared_ptr<Vertex> parameter) {
  double second_derivative;
  for (auto neighbor : parameter->neighboorhood) {
    second_derivative +=
        ((parameter->y_coord - neighbor->y_coord) -
         ((30 * parameter->distances[neighbor]) *
              (parameter->y_coord - neighbor->y_coord) /
              std::sqrt(std::pow(parameter->x_coord - neighbor->x_coord, 2)) +
          pow(parameter->y_coord - neighbor->y_coord, 2))) /
        parameter->distances[neighbor];
  }
  return second_derivative;
}

double Graph::CalculateDelta(std::shared_ptr<Vertex> parameter) {
  double first_derivative = CalculateXDerivative(parameter);
  double second_derivative = CalculateYDerivative(parameter);
  double delta =
      std::sqrt(std::pow(first_derivative, 2) + std::pow(second_derivative, 2));
  return delta;
}
