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

double Graph::EuclideanDistance(Vertex* v, Vertex* u, float power) {
  return std::pow((std::pow(v->x - u->x, 2)) + pow(v->y - u->y, 2), power);
}

double Graph::OneVariableDerivative(Vertex* parameter, char with_respect) {
  double coordinate_p = (with_respect == 'x') ? parameter->x : parameter->y;
  double derivative = 0.0;
  for (auto neighbor : parameter->neighboorhood) {
    double coordinate_n = (with_respect == 'x') ? neighbor->x : neighbor->y;
    derivative += ((coordinate_p - coordinate_n) -
                   (kEdgeLen * parameter->distances[neighbor]) *
                       (coordinate_p - coordinate_n) /
                       EuclideanDistance(parameter, neighbor, 0.5)) /
                  parameter->distances[neighbor];
  }
  return derivative;
}

double Graph::CalculateDelta(Vertex* parameter) {
  double first_derivative = OneVariableDerivative(parameter, 'x');
  double second_derivative = OneVariableDerivative(parameter, 'y');
  double delta =
      std::sqrt(std::pow(first_derivative, 2) + std::pow(second_derivative, 2));
  return delta;
}

double Graph::TwoVariablesDerivative(Vertex* parameter, char with_respect1,
                                     char with_respect2) {
  double derivative = 0.0;
  double coordinate_p1 = (with_respect1 == 'x') ? parameter->x : parameter->y;
  double coordinate_p2 = (with_respect2 == 'x') ? parameter->x : parameter->y;
  for (auto n : parameter->neighboorhood) {
    double coordinate_n1 = (with_respect1 == 'x') ? n->x : n->y;
    double coordinate_n2 = (with_respect2 == 'x') ? n->x : n->y;
    double intermed_value =
        (kEdgeLen * parameter->distances[n] * (coordinate_p1 - coordinate_n1) *
         (coordinate_p2 - coordinate_n2)) /
        EuclideanDistance(parameter, n, 1.5);
    if (!(with_respect1 == 'x' && with_respect2 == 'y')) {
      intermed_value = 1 - intermed_value;
    }
    derivative += intermed_value / parameter->distances[n];
  }
  return derivative;
}

void Graph::SolveLinearEquations(Vertex* p) {

  // Coefficients of the linear equations
  double a_1 = TwoVariablesDerivative(p, 'y', 'y');
  double b_1 = TwoVariablesDerivative(p, 'x', 'y');
  double c_1 = -OneVariableDerivative(p, 'x');

  double a_2 = TwoVariablesDerivative(p, 'x', 'y');
  double b_2 = TwoVariablesDerivative(p, 'x', 'x');
  double c_2 = -OneVariableDerivative(p, 'y');

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

std::vector<std::vector<int8_t>> Graph::GetPixels(
    std::pair<int, int> dimensions) {
  std::vector<std::vector<int8_t>> image_data;

  const int width = dimensions.first;
  const int height = dimensions.second;

  image_data.resize(height, std::vector<int8_t>(3 * width, 255));
  image_data = DrawNumbers(image_data);
  for (auto v : vertices_) {
    std::vector<std::pair<int, int>> circlePoints;

    for (int i = v->x - 7; i <= v->x + 7; i++) {
      for (int j = v->y - 7; j <= v->y + 7; j++) {
        if (std::pow(i - v->x, 2) + std::pow(j - v->y, 2) <= 7 * 7) {
          image_data[j][3 * i] = 0;
          image_data[j][3 * i + 1] = 0;
          image_data[j][3 * i + 2] = 0;
        }
      }
    }

    for (auto u : v->adjacent_vertices) {

      int x1 = v->x;
      int y1 = v->y;
      int x2 = u.lock()->x;
      int y2 = u.lock()->y;
      // Function to draw a line using Bresenham's algorithm
      int dx = std::abs(x2 - x1);
      int dy = std::abs(y2 - y1);

      int sx = (x1 < x2) ? 1 : -1;
      int sy = (y1 < y2) ? 1 : -1;

      int err = dx - dy;

      while (true) {
        // Set the pixel at (x1, y1) to 0
        image_data[y1][3 * x1] = 0;
        image_data[y1][3 * x1 + 1] = 0;
        image_data[y1][3 * x1 + 2] = 0;

        // Check if we've reached the end point
        if (x1 == x2 && y1 == y2) {
          break;
        }

        int e2 = 2 * err;

        // Move along the x-axis
        if (e2 > -dy) {
          err -= dy;
          x1 += sx;
        }

        // Move along the y-axis
        if (e2 < dx) {
          err += dx;
          y1 += sy;
        }
      }
    }
  }

  return image_data;
}

std::vector<std::vector<int8_t>> Graph::DrawNumbers(
    std::vector<std::vector<int8_t>> image_data) {

  BMPWriter number_writer;
  int num_height = 16;
  int num_width = 8;
  number_writer.SetDimensions(num_width, num_height);
  int padding_x = 24;
  int padding_y = 0;
  int dig_count;
  for (auto v : vertices_) {
    dig_count = 1;
    std::string num = std::to_string(v->kVertInd);
    for (auto& digit : num) {
      std::string filename(1, digit);
      filename = "./sprites/" + filename + ".bmp";
      std::vector<std::vector<int8_t>> digit_pixels =
          number_writer.Read(filename);
      for (int y = 0; y < num_height; ++y) {
        for (int x = 0; x < num_width; ++x) {
          image_data[padding_y + y + v->y][padding_x * dig_count +
                                           3 * (x + v->x)] = digit_pixels[y][x];
          image_data[padding_y + y + v->y]
                    [padding_x * dig_count + 3 * (x + v->x) + 1] =
                        digit_pixels[y][x];
          image_data[padding_y + y + v->y]
                    [padding_x * dig_count + 3 * (x + v->x) + 2] =
                        digit_pixels[y][x];
        }
      }
      ++dig_count;
    }
  }
  return image_data;
}