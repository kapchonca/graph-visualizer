#ifndef LAB4_H_BMPHANDLER_H_
#define LAB4_H_BMPHANDLER_H_

#include <fstream>
#include <iostream>
#include <vector>

#pragma pack(push, 1)

struct BMPHeader {
  int16_t signature;
  int32_t fileSize;
  int16_t reserved1;
  int16_t reserved2;
  int32_t dataOffset;
  int32_t headerSize;
  int32_t width;
  int32_t height;
  int16_t planes;
  int16_t bitsPerPixel;
  int32_t compression;
  int32_t imageSize;
  int32_t xPixelsPerMeter;
  int32_t yPixelsPerMeter;
  int32_t colorsUsed;
  int32_t colorsImportant;
};

#pragma pack(pop)

class BMPWriter {
 public:
  BMPWriter();

  bool save(const std::string filename);
  std::vector<std::vector<int8_t>> Read(const std::string& filename);
  int calculateRowSize(BMPHeader header);
  void SetDimensions(int width, int height);
  std::vector<std::vector<int8_t>> imageData_;

 private:
  BMPHeader header_;
};

#endif  // LAB4_H_BMPHANDLER_H_