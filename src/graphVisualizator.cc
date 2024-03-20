#include "h/graphVisualizator.h"

void Visualizator::RoundCoords(Vertex* v) {
  v->x = std::round(v->x);
  v->y = std::round(v->y);
}

std::pair<int, int> Visualizator::GetCoordinates() {

  double min_x = std::numeric_limits<int>::max();
  double min_y = std::numeric_limits<int>::max();

  double max_x = std::numeric_limits<int>::min();
  double max_y = std::numeric_limits<int>::min();

  for (vert_p v : vertices_) {
    RoundCoords(v.get());
    min_x = std::min(min_x, v->x);
    min_y = std::min(min_y, v->y);
    max_x = std::max(max_x, v->x);
    max_y = std::max(max_y, v->y);
  }

  for (vert_p v : vertices_) {
    v->x += kEdgeLen - min_x;
    v->y += kEdgeLen - min_y;
  }

  return std::pair<int, int>(
      {max_x + 2 * kEdgeLen - min_x, max_y + 2 * kEdgeLen - min_y});
}

std::vector<std::vector<int8_t>> Visualizator::GetPixels(
    const std::pair<int, int>& dimensions) {

  const int width = dimensions.first;
  const int height = dimensions.second;

  image_data_.resize(height, std::vector<int8_t>(3 * width, 255));
  DrawNumbers();
  for (vert_p v : vertices_) {

    DrawVertex(v.get());
    DrawEdges(v.get());
  }

  return image_data_;
}

void Visualizator::DrawVertex(Vertex* v) {
  for (int i = v->x - 7; i <= v->x + 7; i++) {
    for (int j = v->y - 7; j <= v->y + 7; j++) {
      if (std::pow(i - v->x, 2) + std::pow(j - v->y, 2) <= 7 * 7) {
        image_data_[j][3 * i] = 0;
        image_data_[j][3 * i + 1] = 0;
        image_data_[j][3 * i + 2] = 0;
      }
    }
  }
}

void Visualizator::DrawEdges(Vertex* v) {
  for (Vertex* u : v->adjacent_vertices) {

    int x1 = v->x;
    int y1 = v->y;
    int x2 = u->x;
    int y2 = u->y;
    // Function to draw a line using Bresenham's algorithm
    int dx = std::abs(x2 - x1);
    int dy = std::abs(y2 - y1);

    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;

    int err = dx - dy;

    while (true) {
      // Set the pixel at (x1, y1) to 0
      image_data_[y1][3 * x1] = 0;
      image_data_[y1][3 * x1 + 1] = 0;
      image_data_[y1][3 * x1 + 2] = 0;

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

void Visualizator::DrawNumbers() {

  BMPWriter number_writer;
  int num_height = 16;
  int num_width = 8;
  number_writer.SetDimensions(num_width, num_height);

  int padding_x = 24;
  std::vector<std::vector<std::vector<int8_t>>> digits(10);

  for (int digit = 0; digit < 10; ++digit) {
    std::string filename = std::to_string(digit);
    filename = "./sprites/" + filename + ".bmp";
    digits[digit] = number_writer.Read(filename);
  }

  int dig_count;
  for (vert_p v : vertices_) {
    dig_count = 1;
    std::string num = std::to_string(v->kVertInd + 1);
    for (char digit : num) {
      int digit_s = digit - '0';
      for (int y = 0; y < num_height; ++y) {
        for (int x = 0; x < num_width; ++x) {
          image_data_[y + v->y][padding_x * dig_count + 3 * (x + v->x)] =
              digits[digit_s][y][x];
          image_data_[y + v->y][padding_x * dig_count + 3 * (x + v->x) + 1] =
              digits[digit_s][y][x];
          image_data_[y + v->y][padding_x * dig_count + 3 * (x + v->x) + 2] =
              digits[digit_s][y][x];
        }
      }
      ++dig_count;
    }
  }
}