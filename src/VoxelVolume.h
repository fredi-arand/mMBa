#pragma once
//------------------------------------------------------------------------------
#include <Eigen/Dense>
#include <vector>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using Vector3l = Eigen::Vector3<long>;
//------------------------------------------------------------------------------
template <typename T> struct VoxelVolume {

  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;

  VoxelVolume() : s(0, 0, 0), spacing(0, 0, 0) {}

  void resize(Vector3l const &s, T const &value = T());

  // standard functions to translate between ID and Coord
  Vector3l vxID_to_vx(size_t vxID) const {
    return Vector3l(vxID % s(0), (vxID / s(0)) % s(1), vxID / (s(0) * s(1)));
  }

  // use this to always get vxID even for illegal values of vx
  size_t vx_to_vxID(Vector3l const &vx) const;

  // import volume from raw file
  void import_raw_volume(Vector3l const &s, const char *filename);

  // export volume as pgm stacks
  void export_pgm_stacks(const char *foldername) const;
  void export_raw(const char *filename) const;

  // export for gnuplot
  void export_stack_for_gp(long stackID, const char *filename) const;

  /// Convenience functions to save typing
  /// @{
  const_reference operator[](size_t vxID) const { return data[vxID]; }
  reference operator[](size_t vxID) { return data[vxID]; }

  const_reference operator[](Vector3l const &x) const {
    return data[vx_to_vxID(x)];
  }

  const_reference operator()(long x0, long x1, long x2) const {
    return data[vx_to_vxID({x0, x1, x2})];
  }

  std::vector<T> const &operator()() const { return data; }
  std::vector<T> &operator()() { return data; }
  /// @}

  std::vector<T> data;
  Vector3l s, spacing;
};
//------------------------------------------------------------------------------
struct MorphologyValue {
  enum State : uint32_t { INIT, ENCLOSED, THROAT, BACKGROUND };
  uint32_t state : 2;
  uint32_t parentId : 30;
};
//------------------------------------------------------------------------------
static_assert(sizeof(MorphologyValue) == 4);
//------------------------------------------------------------------------------
} // namespace fred
