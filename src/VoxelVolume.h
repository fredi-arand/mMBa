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

  VoxelVolume()
      : s(Vector3l(0, 0, 0)), spacing(Vector3l(1, 1, 1)), hasValues(false) {}

  // standard functions to translate between ID and Coord
  Vector3l vxID_to_vx(size_t vxID) const {
    return Vector3l(vxID % s(0), (vxID / s(0)) % s(1), vxID / (s(0) * s(1)));
  }

  // use this to always get vxID even for illegal values of vx
  size_t vx_to_vxID(Vector3l const &vx) const;

  void set_spacing_and_voxelValues_from_s() {
    // quickly set spacing and voxelValues for given s
    spacing << 1, s(0), s(0) * s(1);
    voxelValues.resize(s.cast<size_t>().prod());
  }

  // switch x and y
  void switch_xy();

  // import volume from raw file
  void import_raw_volume(Vector3l const &s, const char *filename);

  // binarize volume
  T binarize_volume(VoxelVolume<bool> &isFict, VoxelVolume<bool> &isPhys) const;
  void binarize_volume(VoxelVolume<bool> &isFict, VoxelVolume<bool> &isPhys,
                       T isoValue) const;

  // create from multiplication of two other volumes
  void create_from_two_volumes(VoxelVolume<T> const &volA,
                               VoxelVolume<T> const &volB);

  // export volume as pgm stacks
  void export_pgm_stacks(const char *foldername) const;
  void export_raw(const char *filename) const;

  // export for gnuplot
  void export_stack_for_gp(long stackID, const char *filename) const;

  T operator[](size_t vxID) const { return voxelValues[vxID]; }

  T operator()(Vector3l x) const { return this->operator[](vx_to_vxID(x)); }
  T operator()(long x0, long x1, long x2) const {
    return this->operator()(Vector3l(x0, x1, x2));
  }

  size_t size() { return voxelValues.size(); }

  // voxelValues: value of voxel (i,j,k) is stored at voxelValues(spacing(2)*k +
  // spacing(1)*j + i) s:           size of voxel grid spacing:     spacing in
  // 0/1/2-direction
  std::vector<T> voxelValues;
  Vector3l s, spacing;
  bool hasValues;
};
//------------------------------------------------------------------------------
} // namespace fred
