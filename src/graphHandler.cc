#include "h/graphHandler.h"

Graph::Graph(const std::string& file_path) {
  std::ifstream input_file(file_path);

  if (!input_file.is_open()) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }

  input_file >> vertex_num_ >> edge_num_;

  // Resize the vertices vector and initialize vertices.
  vertices_.resize(vertex_num_, nullptr);
  for (int i = 0; i < vertex_num_; ++i) {
    vert_p v = std::make_shared<Vertex>(i);
    vertices_[i] = v;
  }

  // Read edges and establish connections between vertices.
  int num1, num2;
  int counter = 0;
  while (input_file >> num1 >> num2) {
    if (num1 < 1 || num1 > vertex_num_ || num2 < 1 || num2 > vertex_num_) {
      throw std::runtime_error("Invalid vertex index in edge definition");
    }

    if (counter == edge_num_) {
      throw std::runtime_error("Invalid number of rows in input file");
    }

    if (num1 != num2) {
      vertices_[num1 - 1]->adjacent_vertices.emplace_back(
          vertices_[num2 - 1].get());
      vertices_[num2 - 1]->adjacent_vertices.emplace_back(
          vertices_[num1 - 1].get());
    }

    ++counter;
  }

  if (counter < edge_num_) {
    throw std::runtime_error("Invalid number of rows in input file");
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
    Vertex* current_vertex = q.front().first;
    int current_depth = q.front().second;
    q.pop();
    traversal[current_vertex] = current_depth;

    // Explore neighbors and update distances.
    for (Vertex* neighbor : current_vertex->adjacent_vertices) {
      if (!visited[neighbor->kVertInd]) {
        visited[neighbor->kVertInd] = true;
        q.push({neighbor, current_depth + 1});
      }
    }
  }
  traversal.erase(start);
  return traversal;
}

// Greedy 2-Approximation Algorithm for k-Center
std::unordered_set<vert_p> Graph::kCenter(long unsigned k) const {
  std::unordered_set<vert_p> centers;
  std::unordered_set<vert_p> remaining_vertices;

  remaining_vertices.insert(vertices_.begin(), vertices_.end());
  // Select the first center arbitrarily
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> distribution(
      0, vertex_num_ - 1);

  vert_p first_vertex = vertices_[distribution(rng)];
  centers.insert(first_vertex);
  remaining_vertices.erase(first_vertex);

  while (centers.size() < k) {
    vert_p farthest_vertex = nullptr;
    int max_distance = -1;

    // Find the vertex farthest from the current set of centers
    for (vert_p v : remaining_vertices) {
      int min_distance = std::numeric_limits<int>::max();

      for (vert_p center : centers) {
        min_distance = std::min(min_distance, v->distances[center.get()]);
      }

      if (min_distance > max_distance) {
        max_distance = min_distance;
        farthest_vertex = v;
      }
    }

    // Add the farthest vertex as a new center
    centers.insert(farthest_vertex);
    remaining_vertices.erase(farthest_vertex);
  }

  return centers;
}

void Graph::RandomLayout() const {
  std::random_device dvc;
  std::mt19937 rng(dvc());
  std::uniform_real_distribution<> distribution(0, vertex_num_ * kEdgeLen);
  int solo_x = kEdgeLen;
  int solo_y = kEdgeLen;

  for (vert_p v : vertices_) {
    if (v->adjacent_vertices.size() <= 1) {
      v->x = solo_x;
      v->y = solo_y;
      if (solo_x >= vertex_num_ * kEdgeLen / 10 - kEdgeLen) {
        solo_y += kEdgeLen;
        solo_x = kEdgeLen;
      } else {
        solo_x += kEdgeLen;
      }
      continue;
    }
    v->x = distribution(rng);
    v->y = distribution(rng);
  }
}

void Graph::FormNeighbourhood(Vertex* center, int radius) {
  center->neighboorhood.clear();
  for (int i = 1; i < radius; ++i) {
    for (Vertex* v : shortest_paths_[center][i]) {
      center->neighboorhood.emplace_back(v);
    }
  }
}

void Graph::LocalLayout(Vertex* p, int radius) {
  if (p->adjacent_vertices.size() == 0) {
    return;
  }
  double max_delta;
  Vertex* v_to_adjust = nullptr;
  for (unsigned long i = 0; i < kIterations * vertices_.size(); ++i) {
    max_delta = 0;
    for (Vertex* v : p->neighboorhood) {
      FormNeighbourhood(v, radius);
      double delta = graph_math.CalculateDelta(v);
      if (delta > max_delta) {
        max_delta = delta;
        v_to_adjust = v;
      }
    }
    if (v_to_adjust != nullptr) {
      graph_math.SolveLinearEquations(v_to_adjust);
    }
  }
}

void Graph::ShortestPaths(Vertex* start) {
  std::vector<bool> visited(vertex_num_, false);
  std::queue<std::pair<Vertex*, int>> q;

  visited[start->kVertInd] = true;
  q.push({start, 0});

  while (!q.empty()) {
    Vertex* current_vertex = q.front().first;
    int current_depth = q.front().second;
    q.pop();
    shortest_paths_[start][current_depth].insert(current_vertex);

    // Explore neighbors and update distances.
    for (Vertex* neighbor : current_vertex->adjacent_vertices) {
      if (!visited[neighbor->kVertInd]) {
        visited[neighbor->kVertInd] = true;
        q.push({neighbor, current_depth + 1});
      }
    }
  }
}

void Graph::GlobalLayout() {

  for (vert_p v : vertices_) {
    v->distances = BFS(v.get(), vertex_num_);
    ShortestPaths(v.get());
  }

  RandomLayout();

  int radius;
  int k = kMinSize;

  while (k <= vertex_num_) {
    std::unordered_set<vert_p> c = kCenter(k);

    int max_distance = -1;

    // Find the vertex farthest from the current set of centers
    // graph diameter
    for (vert_p v : c) {
      int min_distance = std::numeric_limits<int>::max();

      for (vert_p n : c) {
        if (n != v) {
          min_distance = std::min(min_distance, v->distances[n.get()]);
        }
      }

      if (min_distance > max_distance) {
        max_distance = min_distance;
      }
    }
    radius = kLocalRadius * max_distance;
    for (vert_p center : c) {
      FormNeighbourhood(center.get(), radius);
      LocalLayout(center.get(), radius);
    }
    std::cout
        << k
        << " centers are processesed successfully\nproceed to next iteration\n"
        << std::flush;
    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<> distribution(0, 1);
    for (vert_p v : vertices_) {
      v->x += distribution(gen);
      v->y += distribution(gen);
    }
    k *= kRatio;
  }
}

void Graph::DrawGraph(const std::string& filename) {

  GlobalLayout();

  Visualizator visual(vertices_);
  std::pair<int, int> dimensions = visual.GetCoordinates();

  write.SetDimensions(dimensions.first, dimensions.second);
  write.SetImageData(visual.GetPixels(dimensions));
  write.save("output/" + filename + ".bmp");
}