#include "h/bmpHandler.h"

BMPWriter::BMPWriter(const std::string& filename, int width, int height) {
  file_.open(filename, std::ios::out | std::ios::binary);

  if (!file_.is_open()) {
    throw std::runtime_error("Error: Could not open file for writing.");
    // return;
  }

  header_.signature = 0x4D42;
  header_.fileSize = 0;
  header_.reserved1 = 0;
  header_.reserved2 = 0;
  header_.dataOffset = sizeof(BMPHeader);
  header_.headerSize = sizeof(BMPHeader) - 14;
  header_.width = width;
  header_.height = height;
  header_.planes = 1;
  header_.bitsPerPixel = 24;
  header_.compression = 0;
  header_.imageSize = 0;
  header_.xPixelsPerMeter = 2835;
  header_.yPixelsPerMeter = 2835;
  header_.colorsUsed = 0;
  header_.colorsImportant = 0;

  imageData_.resize(height, std::vector<int8_t>(3 * width, 0));
}

BMPWriter::~BMPWriter() {
  file_.close();
}

void BMPWriter::setPixel(int x, int y, int8_t red, int8_t green, int8_t blue) {
  if (x >= 0 && x < header_.width && y >= 0 && y < header_.height) {
    int index = 3 * x;
    imageData_[y][index] = blue;
    imageData_[y][index + 1] = green;
    imageData_[y][index + 2] = red;
  }
}

bool BMPWriter::save() {
  if (!file_.is_open()) {
    std::cerr << "Error: File not open for writing." << std::endl;
    return false;
  }

  header_.imageSize = header_.height * ((header_.width * 3 + 3) / 4 * 4);
  header_.fileSize = header_.dataOffset + header_.imageSize;
  file_.write(reinterpret_cast<char*>(&header_), sizeof(BMPHeader));

  for (int i = header_.height - 1; i >= 0; --i) {
    file_.write(reinterpret_cast<char*>(imageData_[i].data()),
                3 * header_.width);
    for (int j = 0; j < (4 - (3 * header_.width) % 4) % 4; ++j) {
      file_.put(0);
    }
  }

  std::cout << "BMP image saved successfully." << std::endl;
  return true;
}
