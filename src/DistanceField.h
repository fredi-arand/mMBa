#pragma once
//------------------------------------------------------------------------------
#include "VoxelVolume.h"
#include <optional>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
struct DistanceField : public VoxelVolume<float> {

  template <typename T>
  static DistanceField
  create(Vector3l const &s, const char *filename,
         std::optional<float> const &isoValue = std::nullopt);

  template <typename T>
  static DistanceField
  create(VoxelVolume<T> const &voxelVolume,
         std::optional<float> const &isoValue = std::nullopt);

  void calculate_porosity() const;

private:
  DistanceField() = default;
};
//------------------------------------------------------------------------------
} // namespace fred
