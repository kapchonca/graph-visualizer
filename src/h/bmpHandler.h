#ifndef LAB4_H_BMPHANDLER_H_
#define LAB4_H_BMPHANDLER_H_

#include <fstream>
#include <iostream>
#include <vector>

#pragma pack(push, 1)

struct BMPHeader {
  int16_t signature = 0x4D42;
  int32_t fileSize = 0;
  int16_t reserved1 = 0;
  int16_t reserved2 = 0;
  int32_t dataOffset = sizeof(BMPHeader);
  int32_t headerSize = sizeof(BMPHeader) - 14;
  int32_t width = 0;
  int32_t height = 0;
  int16_t planes = 1;
  int16_t bitsPerPixel = 24;
  int32_t compression = 0;
  int32_t imageSize = 0;
  int32_t xPixelsPerMeter = 2835;
  int32_t yPixelsPerMeter = 2835;
  int32_t colorsUsed = 0;
  int32_t colorsImportant = 0;
};

#pragma pack(pop)

class BMPWriter {
 public:
  // Write data to a bmp image
  bool save(const std::string& filename);

  // Read image data from the given file
  std::vector<std::vector<int8_t>> Read(const std::string& filename);

  // Calculates row size with padding
  int calculateRowSize();

  // Sets dimensions of the output image
  void SetDimensions(int width, int height);

  // Sets image data for output image
  void SetImageData(std::vector<std::vector<int8_t>> image_data);

 private:
  std::vector<std::vector<int8_t>> image_data_;
  BMPHeader header_;
};

#endif  // LAB4_H_BMPHANDLER_H_