#include "h/bmpHandler.h"

bool BMPWriter::save(const std::string& filename) {
  std::ofstream file_(filename);
  if (!file_.is_open()) {
    std::cerr << "Error: File not open for writing." << std::endl;
    return false;
  }
  header_.bitsPerPixel = 24;
  header_.imageSize = header_.height * ((header_.width * 3 + 3) / 4 * 4);
  header_.fileSize = header_.dataOffset + header_.imageSize;
  file_.write(reinterpret_cast<char*>(&header_), sizeof(BMPHeader));

  for (int i = 0; i < header_.height; ++i) {
    file_.write(reinterpret_cast<char*>(imageData_[i].data()),
                3 * header_.width);
    for (int j = 0; j < (4 - (3 * header_.width) % 4) % 4; ++j) {
      file_.put(0);
    }
  }

  file_.write(reinterpret_cast<char*>(imageData_.data()), header_.height);
  std::cout << "BMP image saved successfully." << std::endl;
  file_.close();
  return true;
}

std::vector<std::vector<int8_t>> BMPWriter::Read(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Error: Could not open file for reading.");
  }
  file.read(reinterpret_cast<char*>(&header_), sizeof(BMPHeader));

  int bytesToSkip = header_.dataOffset - sizeof(BMPHeader);

  file.seekg(bytesToSkip, std::ios::cur);

  int rowSize = calculateRowSize();
  int paddedSize = rowSize * header_.height;

  std::vector<int8_t> temp_img_data(paddedSize, 0);

  file.read(reinterpret_cast<char*>(temp_img_data.data()), paddedSize);

  file.close();

  for (int h = 0; h < header_.height; ++h) {
    for (int w = 0; w < header_.width; ++w) {
      imageData_[h][w] = temp_img_data[(h * header_.width + w)];
    }
  }

  return imageData_;
}

int BMPWriter::calculateRowSize() {
  int bytesPerPixel = header_.bitsPerPixel / 8;
  int rowSizeBytes = ((header_.width * bytesPerPixel + 3) / 4) * 4;
  return rowSizeBytes;
}

void BMPWriter::SetDimensions(int width, int height) {
  header_.width = width;
  header_.height = height;
  imageData_.resize(height,
                    std::vector<int8_t>(header_.bitsPerPixel / 8 * width, 0));
}

void BMPWriter::SetImageData(std::vector<std::vector<int8_t>> imageData) {
  imageData_ = imageData;
}