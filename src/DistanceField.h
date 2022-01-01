#pragma once
//------------------------------------------------------------------------------
#include "VoxelVolume.h"
#include <optional>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------

struct DistanceField : public VoxelVolume<float> {
  DistanceField() : distanceFieldCreated(false) {}

  template <typename T>
  void
  create_distance_field(Vector3l const &s, const char *filename,
                        std::optional<float> const &isoValue = std::nullopt);

  template <typename T>
  void
  create_distance_field(VoxelVolume<T> const &voxelVolume,
                        std::optional<float> const &isoValue = std::nullopt);

  void calculate_porosity() const;

  bool distanceFieldCreated;
}; //------------------------------------------------------------------------------
} // namespace fred
