#pragma once

#include <stdexcept>

#define ASSURE(condition, msg)                                                 \
  do {                                                                         \
    if (!(condition))                                                          \
      throw std::logic_error(msg);                                             \
  } while (false)
