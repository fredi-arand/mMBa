#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <Eigen/Dense>
#include <map>
#include <set>
#include <vector>
#include <utility>


using PoreID = uint32_t;
using ThroatID = std::set<PoreID>;
using VoxelIndex = size_t;
using VoxelPosition = Eigen::Vector3i;
using FloatPosition = Eigen::Vector3f;

using Pore = std::vector<VoxelIndex>;
using Throat = Pore;

uint32_t const INIT_VALUE = 0;
uint32_t const ENCLOSED = 1;
uint32_t const THROAT = 2;
uint32_t const BACKGROUND = 3;

#endif // DEFINITIONS_H
