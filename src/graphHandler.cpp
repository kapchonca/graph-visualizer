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
    std::cout << "vertex index: " << v->kVertInd + 1
              << " initial placement: " << v->x << ' ' << v->y << std::endl;
  }
}

void Graph::FormNeighbourhood(Vertex* center, int radius) {
  center->neighboorhood.clear();
  std::unordered_map<Vertex*, int> traversal = BFS(center, radius);
  for (auto pair : traversal) {
    center->neighboorhood.insert(pair.first);
  }
}

double Graph::EuclideanDistance(Vertex* v, Vertex* u, float power) {
  return std::pow((std::pow(v->x - u->x, 2)) + pow(v->y - u->y, 2), power);
}

double Graph::CalculateXDerivative(Vertex* parameter) {
  double first_derivative = 0.0;
  for (auto neighbor : parameter->neighboorhood) {
    first_derivative += ((parameter->x - neighbor->x) -
                         ((kEdgeLen * parameter->distances[neighbor]) *
                          (parameter->x - neighbor->x) /
                          EuclideanDistance(parameter, neighbor, 0.5))) /
                        parameter->distances[neighbor];
  }
  return first_derivative;
}

double Graph::CalculateYDerivative(Vertex* parameter) {
  double second_derivative = 0.0;
  for (auto neighbor : parameter->neighboorhood) {
    second_derivative += ((parameter->y - neighbor->y) -
                          (kEdgeLen * parameter->distances[neighbor]) *
                              (parameter->y - neighbor->y) /
                              EuclideanDistance(parameter, neighbor, 0.5)) /
                         parameter->distances[neighbor];
  }
  return second_derivative;
}

double Graph::CalculateDelta(Vertex* parameter) {
  double first_derivative = CalculateXDerivative(parameter);
  double second_derivative = CalculateYDerivative(parameter);
  double delta =
      std::sqrt(std::pow(first_derivative, 2) + std::pow(second_derivative, 2));
  return delta;
}

double Graph::CalculateX_XDerivative(Vertex* parameter) {
  double x_x_derivative = 0.0;
  for (auto n : parameter->neighboorhood) {
    x_x_derivative += (1 - (kEdgeLen * parameter->distances[n] *
                            std::pow(parameter->y - n->y, 2)) /
                               EuclideanDistance(parameter, n, 1.5)) /
                      parameter->distances[n];
  }
  return x_x_derivative;
}

double Graph::CalculateX_YDerivative(Vertex* parameter) {
  double x_y_derivative = 0.0;
  for (auto n : parameter->neighboorhood) {
    x_y_derivative += (kEdgeLen * parameter->distances[n] *
                       (parameter->y - n->y) * (parameter->x - n->x)) /
                      EuclideanDistance(parameter, n, 1.5) /
                      parameter->distances[n];
  }
  return x_y_derivative;
}

double Graph::CalculateY_YDerivative(Vertex* parameter) {
  double y_y_derivative = 0.0;
  for (auto n : parameter->neighboorhood) {
    y_y_derivative += (1 - (kEdgeLen * parameter->distances[n] *
                            std::pow(parameter->x - n->x, 2)) /
                               EuclideanDistance(parameter, n, 1.5)) /
                      parameter->distances[n];
  }
  return y_y_derivative;
}

void Graph::SolveLinearEquations(Vertex* p) {

  // Coefficients of the linear equations
  double a_1 = CalculateX_XDerivative(p);
  double b_1 = CalculateX_YDerivative(p);
  double c_1 = -CalculateXDerivative(p);

  double a_2 = CalculateX_YDerivative(p);
  double b_2 = CalculateY_YDerivative(p);
  double c_2 = -CalculateYDerivative(p);

  double delta_x = 0;
  double delta_y = 0;

  if (a_1 != 0) {  // Normalize first coeff
    b_1 /= a_1;
    c_1 /= a_1;
    a_1 = 1.0;

    if (a_2 != 0) {  // Normalize first coeff
      b_2 /= a_2;
      c_2 /= a_2;
      a_2 = 1.0;

      // Subtract the 2nd equation from the 1st equation
      a_1 = 0.0;
      b_1 -= b_2;
      c_1 -= c_2;

      if (b_1 != 0) {
        delta_y = c_1 / b_1;
        delta_x = c_2 - delta_y * b_2;
      }

    } else if (b_2 != 0) {
      delta_y = c_2 / b_2;
      delta_x = (c_1 - delta_y * b_1) / a_1;
    }

  } else if (b_1 != 0) {
    delta_y = c_1 / b_1;

    if (a_2 != 0) {
      delta_x = (c_2 - delta_y * b_2) / a_2;
    }
  }

  p->x += delta_x;
  p->y += delta_y;
}

void Graph::LocalLayout(Vertex* p, int radius) {
  double max_delta;
  Vertex* v_to_adjust;
  for (unsigned long i = 0; i < kIterations * vertices_.size(); ++i) {
    max_delta = 0;
    for (auto v : p->neighboorhood) {
      FormNeighbourhood(v, radius);
      double delta = CalculateDelta(v);
      if (delta > max_delta) {
        max_delta = delta;
        v_to_adjust = v;
      }
    }
    SolveLinearEquations(v_to_adjust);
  }
}

void Graph::GlobalLayout() {

  for (auto v : vertices_) {
    v->distances = BFS(v.get(), vertex_num_);
  }

  RandomLayout();

  int radius;
  int k = kMinSize;
  while (k <= vertex_num_) {
    auto c = kCenter(k);

    int maxDistance = -1;

    // Find the vertex farthest from the current set of centers
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

void Graph::RoundCoords(Vertex* v) {
  v->x = std::round(v->x);
  v->y = std::round(v->y);
}

std::pair<int, int> Graph::MoveCoordinates() {

  double min_x = std::numeric_limits<int>::max();
  double min_y = std::numeric_limits<int>::max();

  double max_x = std::numeric_limits<int>::min();
  double max_y = std::numeric_limits<int>::min();

  for (auto v : vertices_) {
    RoundCoords(v.get());
    min_x = std::min(min_x, v->x);
    min_y = std::min(min_y, v->y);
    max_x = std::max(max_x, v->x);
    max_y = std::max(max_y, v->y);
  }

  for (auto v : vertices_) {
    v->x += kEdgeLen - min_x;
    v->y += kEdgeLen - min_y;
  }

  return std::pair<int, int>(
      {max_x + 2 * kEdgeLen - min_x, max_y + 2 * kEdgeLen - min_y});
}