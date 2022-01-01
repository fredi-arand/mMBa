#pragma once

#include "DistanceField.h"
#include "VoxelVolume.h"
#include <map>
#include <vector>

#ifdef ENABLE_GNU_PARALLEL
#include <parallel/algorithm>
#endif

#define PI 3.14159265358979323846
namespace fred {

struct PoreMorphology {

  uint32_t initValue{0}, enclosedValue{uint32_t(1) << 30},
      throatValue{uint32_t(2) << 30}, backgroundValue{uint32_t(3) << 30};
  uint32_t parentCounter{0};
  float smallPadding{0.5};
  float epsilon{0.2f};

  bool parallelFlag{true};

  // Silin, Patzek (2006): Pore space morphology analysis using maximal
  // inscribed spheres

  PoreMorphology(DistanceField const &distanceField)
      : poreMorphologyCreated(false), throatsReduced(false),
        distanceFieldP(&distanceField) {}

  void create_pore_morphology() {
    //      auto rMax = *max_element(distanceFieldP->voxelValues.begin(),
    //      distanceFieldP->voxelValues.end());
    //      create_pore_morphology(max(rMax/20.0,1.0), 0.0);
    create_pore_morphology(1.0, 0.0);
  }
  //    void
  //    create_pore_morphology_single(){create_pore_morphology_single(1.0,0.0);}

  void create_pore_morphology(float rMinMaster, float rMinBall);
  //    void create_pore_morphology_single(float rMinMaster, float rMinBall);

  void update_neighbors_flood(size_t const &voxelIndex_i);
  void update_neighbors_box(const size_t &voxelIndex_i);
  bool quick_neighbor_check(size_t i);

  void export_ppm_stacks(const char *foldername);

  void reduce_throat_volume();

  void export_statistics(const char *foldername);

  void merge_pores(float throatRatio);

  //    void export_effective_radius_statistics();

  template <typename T>
  void export_histogram(std::vector<T> const &values,
                        const char *filename) const;

  template <typename T>
  void export_histogram(std::vector<T> const &values, const char *filename,
                        double minValue, double maxValue,
                        size_t nrOfBins) const;

  void create_legacy_volumes(VoxelVolume<uint32_t> &morphologyVolume,
                             VoxelVolume<uint8_t> &stateVolume);

  //    void alter_morphology_highest_connected_pore();

  uint32_t get_flag(uint32_t const &morphologyValue) const;
  uint32_t get_parent(uint32_t const &morphologyValue) const;

  VoxelVolume<uint32_t> morphologyVolume;

  std::map<uint32_t, size_t> parentToVoxelIndex;

  bool poreMorphologyCreated, throatsReduced;
  DistanceField const *distanceFieldP;
};
//------------------------------------------------------------------------------
} // namespace fred
