#pragma once

#include <stdint.h>

class pSort {

public:

  typedef struct {
    uint32_t x;
    char a, b;
    int8_t payload[58];
  } dataType;

  void init();
  void close();
  void sort(dataType *data, int32_t n);
};
