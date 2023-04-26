#pragma once

#include "definitions.h"
#include "voxelvolume.h"
#include <map>

namespace morphologyStatistics{


  std::map<PoreID,Pore> extract_pores(const VoxelVolume<PoreID> &morphologyVolume,
                                      VoxelVolume<uint8_t> const & stateVolume, bool includeBoundaryPores);

  std::map<ThroatID,Throat> extract_throats(VoxelVolume<PoreID> const & morphologyVolume,
                                            VoxelVolume<uint8_t> const & stateVolume);

  std::map<PoreID,float> calculate_pore_volumes(std::map<PoreID,Pore> const & pores);

  std::map<PoreID,float> calculate_effective_radii(const std::map<PoreID, float> &poreVolumes);

  std::map<ThroatID,float> calculate_throat_volumes(const std::map<ThroatID, Throat> &throats);

  float volume_of_(Pore const & poreOrThroat);

  void cout_void_fractions(const std::map<PoreID, Pore> &pores,
                           const std::map<ThroatID, Throat> &throats,
                           const VoxelVolume<float> &distanceField, float r_small_large);

  Eigen::Vector3f center_of_mass(Pore const & pore,
                                 VoxelVolume<float> const & distanceField);

  Eigen::Matrix3f covariance_matrix(Pore const & pore,
                                    Eigen::Vector3f const & centerOfMass,
                                    VoxelVolume<float> const & distanceField);

  void logarithmic_histogram_r_eff(std::string filename,
                                         const std::map<PoreID, float> &poreRadii, long bins);

  std::map<PoreID,std::vector<PoreID> > pore_network(std::map<PoreID,Pore> const & pores,
                                                     std::map<ThroatID,Throat> const & throats);

  void average_connectivity_vs_effective_radius(std::string filename,
                                                std::map<PoreID,float> const & poreRadii,
                                                std::map<PoreID,std::vector<PoreID> > const & poreNetwork,
                                                long bins);

  void pore_projections(std::string directory,
                        std::map<PoreID,Pore> const & pores,
                        VoxelVolume<float> const & distanceField);

  void seperate_pores_by_size(float r_eff,
                              std::map<PoreID,Pore> & smallPores,
                              std::map<PoreID,Pore> & largePores,
                              std::map<PoreID,Pore> const & pores,
                              VoxelVolume<float> const & distanceField );

  Eigen::Matrix3d cout_average_principal_component_analysis(std::map<PoreID,Pore> const & pores,
                                                 VoxelVolume<float> const & distanceField);

  void compare_anisotropy(std::map<PoreID,Pore> const & pores,
                          VoxelVolume<float> const & distanceField,
                          std::string filename);

  void histogramm_r_eff(std::string filename,
                        std::map<PoreID,float> const & poreRadii,
                        long bins);

  void mle_normal(std::map<PoreID,float> const & poreRadii, float samplingFactor);

  float effective_radius(float volume);

  void material_voxels_in_direction(int direction, VoxelVolume<float> const & distanceField,
                                    std::string filename);

  VoxelVolume<uint8_t> create_voronoi_seeds(std::map<PoreID,Pore> const & pores, const VoxelVolume<float> &distanceField);

  std::vector<VoxelIndex> find_largest_inscribed_sphere_centers(auto const & pores,
                                                                VoxelVolume<float> const & distanceField);

  void nodes_vs_ligaments_distribution(std::map<PoreID,Pore> const & pores,
                                       std::map<PoreID,std::vector<PoreID> > const & poreNetwork,
                                       VoxelVolume<float> const & distanceField);

}
