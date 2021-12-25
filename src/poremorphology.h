#pragma once

#include "distancefield.h"
#include "voxelvolume.h"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stdlib.h>
#include <time.h>
#include <tuple>
#include <utility>
#ifdef ENABLE_GNU_PARALLEL
#include <parallel/algorithm>
#endif

#define PI 3.14159265358979323846

using namespace std;
using namespace std::chrono;
using namespace Eigen;

typedef Matrix<uint8_t, 3, 1> Vector3ui8;

struct PoreMorphology {

  uint32_t initValue{0}, enclosedValue{uint32_t(1) << 30},
      throatValue{uint32_t(2) << 30}, backgroundValue{uint32_t(3) << 30};
  uint32_t parentCounter{0};
  float smallPadding{0.5};
  float epsilon{0.2f};

  bool parallelFlag{true};

  // Silin, Patzek (2006): Pore space morphology analysis using maximal insribed
  // spheres

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

  void export_ppm_stacks(string foldername);

  void reduce_throat_volume();

  void export_statistics(string foldername);

  void merge_pores(float throatRatio);

  //    void export_effective_radius_statistics();

  template <typename T>
  void export_histogram(vector<T> const &values, string filename) const;

  template <typename T>
  void export_histogram(vector<T> const &values, string filename,
                        double minValue, double maxValue,
                        size_t nrOfBins) const;

  void create_legacy_volumes(VoxelVolume<uint32_t> &morphologyVolume,
                             VoxelVolume<uint8_t> &stateVolume);

  //    void alter_morphology_highest_connected_pore();

  uint32_t get_flag(uint32_t const &morphologyValue) const;
  uint32_t get_parent(uint32_t const &morphologyValue) const;

  VoxelVolume<uint32_t> morphologyVolume;

  map<uint32_t, size_t> parentToVoxelIndex;

  bool poreMorphologyCreated, throatsReduced;
  DistanceField const *distanceFieldP;
};

//*************************************************************************************************

bool PoreMorphology::quick_neighbor_check(size_t i) {

  bool ignoreBall = false;

  Vector3i x_i = morphologyVolume.vxID_to_vx(i);

  bool hasNeighbor = false;
  bool hasTwoNeighbors = false;
  uint32_t neighborLabel = 0;
  for (int K = -1; K <= 1; ++K)
    for (int J = -1; J <= 1; ++J)
      for (int I = -1; I <= 1; ++I) {
        Vector3i x_j = x_i + Vector3i(I, J, K);
        uint32_t m_j = morphologyVolume(x_j);

        uint32_t m_j_state = get_flag(m_j);
        uint32_t m_j_label = get_parent(m_j);

        if (m_j_state == backgroundValue || m_j_state == throatValue ||
            m_j_label == 0)
          continue;

        if (!hasNeighbor) {
          hasNeighbor = true;
          neighborLabel = m_j_label;
          continue;
        }

        if (m_j_label != neighborLabel)
          hasTwoNeighbors = true;
      }

  if (hasTwoNeighbors) {
    morphologyVolume[i] = throatValue;
    ignoreBall = true;
    return ignoreBall;
  }

  if (hasNeighbor) {
    morphologyVolume[i] = neighborLabel;
  }

  return ignoreBall;
}

//*************************************************************************************************

uint32_t PoreMorphology::get_parent(uint32_t const &morphologyValue) const {
  static uint32_t morphologyReader =
      (~uint32_t(0)) - (uint32_t(1) << 31) - (uint32_t(1) << 30);
  return morphologyReader & morphologyValue;
}

//*************************************************************************************************

uint32_t PoreMorphology::get_flag(uint32_t const &morphologyValue) const {
  static uint32_t flagReader = (uint32_t(1) << 31) + (uint32_t(1) << 30);
  return flagReader & morphologyValue;
}

//*************************************************************************************************

void PoreMorphology::create_legacy_volumes(
    VoxelVolume<uint32_t> &_morphologyVolume,
    VoxelVolume<uint8_t> &stateVolume) {
  auto const &s = morphologyVolume.s;
  auto const &spacing = morphologyVolume.spacing;
  _morphologyVolume.s = s;
  _morphologyVolume.spacing = spacing;
  _morphologyVolume.voxelValues.clear();
  _morphologyVolume.voxelValues.resize(s.cast<size_t>().prod(), 0);
  stateVolume.s = s;
  stateVolume.spacing = spacing;
  stateVolume.voxelValues.clear();
  stateVolume.voxelValues.resize(s.cast<size_t>().prod(), 0);

  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n) {
    uint32_t flag = get_flag(morphologyVolume[n]);
    uint32_t parent = get_parent(morphologyVolume[n]);

    stateVolume[n] = uint8_t(flag >> 30);
    _morphologyVolume[n] = uint32_t(parent);
  }
}

//*************************************************************************************************

void PoreMorphology::merge_pores(float throatRatio) {

  cout << "\nMerging Pores with maximum throat ratio > " << throatRatio
       << " with larger Pores...\n";

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  set<set<uint32_t>> throats;
  map<set<uint32_t>, float> maxThroatRadii;
  map<uint32_t, uint32_t> changeSet;
  // and: map<uint32_t,size_t> parentToVoxelIndex;

  auto const &s = morphologyVolume.s;

  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    if (get_flag(morphologyVolume[voxelIndex]) != throatValue)
      continue;

    set<uint32_t> throat;
    Vector3i voxelCoordinate = morphologyVolume.vxID_to_vx(voxelIndex);

    for (int K = -1; K <= 1; ++K)
      for (int J = -1; J <= 1; ++J)
        for (int I = -1; I <= 1; ++I) {
          Vector3i neighborCoordinate = voxelCoordinate + Vector3i(I, J, K);
          uint32_t neighborMorphology = morphologyVolume(neighborCoordinate);
          if (get_flag(neighborMorphology) != enclosedValue)
            continue;

          throat.insert(get_parent(neighborMorphology));
        }

    if (throat.size() < 2) {
      cout << endl;
      for (int K = -1; K <= 1; ++K) {
        for (int J = -1; J <= 1; ++J) {
          for (int I = -1; I <= 1; ++I) {
            Vector3i neighborCoordinate = voxelCoordinate + Vector3i(I, J, K);
            uint32_t neighborMorphology = morphologyVolume(neighborCoordinate);
            uint32_t neighborFlag = get_flag(neighborMorphology);
            cout << (neighborFlag >> 30) << " "
                 << get_parent(neighborMorphology) << "    ";
          }
          cout << endl;
        }
        cout << endl;
      }
      cout << endl;
      return;
    }

    throats.insert(throat);

    if (maxThroatRadii.count(throat) == 0) {
      maxThroatRadii[throat] = (*distanceFieldP)[voxelIndex];
      continue;
    }

    if (maxThroatRadii[throat] < (*distanceFieldP)[voxelIndex])
      maxThroatRadii[throat] = (*distanceFieldP)[voxelIndex];
  }

  // sort ascending
  // http://stackoverflow.com/questions/8736997/using-lambdas-in-maps
  map<size_t, set<set<uint32_t>>,
      function<bool(size_t const &, size_t const &)>>
  poreCentersToThroats([&](size_t const &i, size_t const &j) -> bool {
    return (*distanceFieldP)[i] == (*distanceFieldP)[j]
               ? i < j
               : (*distanceFieldP)[i] < (*distanceFieldP)[j];
  });

  for (auto const &throat : throats)
    for (auto const &parent : throat) {
      size_t voxelIndex = parentToVoxelIndex.at(parent);
      poreCentersToThroats[voxelIndex].insert(throat);
    }

  for (auto poreCenterToThroats = poreCentersToThroats.begin();
       poreCenterToThroats != poreCentersToThroats.end();
       ++poreCenterToThroats) {

    // current pore
    size_t poreVoxelIndex = poreCenterToThroats->first;
    uint32_t parent = get_parent(morphologyVolume[poreVoxelIndex]);
    set<set<uint32_t>> connectedThroats = poreCenterToThroats->second;

    // current radius
    float rPore = (*distanceFieldP)[poreVoxelIndex];

    // pore without throats --> continue
    if (connectedThroats.size() == 0)
      continue;

    // find throat with maximal radius
    float maxThroatRadius = 0.0;
    set<uint32_t> maxThroat;
    for (auto const &throat : connectedThroats)
      if (maxThroatRadii.at(throat) > maxThroatRadius) {
        maxThroatRadius = maxThroatRadii.at(throat);
        maxThroat = throat;
      }

    //    cout << endl << maxThroatRadius/rPore << endl;

    if (maxThroatRadius / rPore <= throatRatio)
      continue;

    // if largest throat is larger than the throat ratio, change pore to
    // the largest one it is connected to via the throat

    if (maxThroat.size() < 2) {
      cout << endl << *maxThroat.begin() << endl;
      return;
    }

    uint32_t mergeToParent = 0;
    float parentRadius = 0.0;
    for (auto const &neighborParent : maxThroat) {
      if (parentToVoxelIndex.count(neighborParent) == 0)
        cout << endl << neighborParent << endl;

      if (neighborParent != parent &&
          (*distanceFieldP)[parentToVoxelIndex.at(neighborParent)] >
              parentRadius) {
        parentRadius = (*distanceFieldP)[parentToVoxelIndex.at(neighborParent)];
        mergeToParent = neighborParent;
      }
    }

    // 1: add to changeset
    changeSet[parent] = mergeToParent;
    if (parent == 0 || mergeToParent == 0 || parent == mergeToParent) {
      cout << endl;
      for (auto const &neighborParent : maxThroat) {
        cout << neighborParent << " ";
      }
      cout << endl << parent << " " << mergeToParent;
      return;
    }

    // 2: update other throats
    for (set<uint32_t> const &throat : connectedThroats) {
      set<uint32_t> newThroat = throat;
      newThroat.erase(parent);
      newThroat.insert(mergeToParent);

      auto throatAndRadius = maxThroatRadii.find(throat);
      float throatRadius = throatAndRadius->second;
      maxThroatRadii.erase(throatAndRadius);

      if (newThroat.size() > 1)
        maxThroatRadii[newThroat] = throatRadius;

      // 3: update poreCentersToThroats
      for (uint32_t parentAtThroat : newThroat) {
        size_t poreCenter = parentToVoxelIndex[parentAtThroat];
        poreCentersToThroats[poreCenter].erase(throat);

        if (newThroat.size() > 1)
          poreCentersToThroats[poreCenter].insert(newThroat);
      }
    }

    // 4: remove from parentToVoxelIndex
    parentToVoxelIndex.erase(parent);
    //    cout << "\nmerged: " << parent  << " --> " << mergeToParent;
  }

  //  size_t counter = 0;‚
  for (auto firstChange = changeSet.begin(); firstChange != changeSet.end();
       ++firstChange) {
    vector<uint32_t> from_parents(1, firstChange->first);
    uint32_t to_parent = firstChange->second;

    while (changeSet.count(to_parent) != 0) {

      if (to_parent == changeSet.at(to_parent)) {
        cout << endl << to_parent << "-->" << changeSet[to_parent] << endl;
        return;
      }

      //      cout << to_parent << "-->" <<  changeSet.at(to_parent) << endl;

      from_parents.push_back(to_parent);
      to_parent = changeSet.at(to_parent);
    }

    if (from_parents.size() == 1)
      continue;

    for (auto const &from_parent : from_parents)
      changeSet[from_parent] = to_parent;
  }

  // use changeSet on morphologyVolume
  for (uint32_t &morphologyValue : morphologyVolume.voxelValues) {
    uint32_t flag = get_flag(morphologyValue);
    if (flag != enclosedValue)
      continue;

    uint32_t parent = get_parent(morphologyValue);

    if (changeSet.count(parent) == 0)
      continue;

    morphologyValue = changeSet[parent] + enclosedValue;
  }

  // look for "ghost" throats
  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    if (get_flag(morphologyVolume[voxelIndex]) != throatValue)
      continue;

    set<uint32_t> throat;
    Vector3i voxelCoordinate = morphologyVolume.vxID_to_vx(voxelIndex);

    for (int K = -1; K <= 1; ++K)
      for (int J = -1; J <= 1; ++J)
        for (int I = -1; I <= 1; ++I) {
          Vector3i neighborCoordinate = voxelCoordinate + Vector3i(I, J, K);
          uint32_t neighborMorphology = morphologyVolume(neighborCoordinate);
          if (get_flag(neighborMorphology) == enclosedValue)
            throat.insert(get_parent(neighborMorphology));
        }

    if (throat.size() > 1)
      continue;

    morphologyVolume[voxelIndex] = *(throat.begin()) + enclosedValue;
  }

  cout << "\nRemoved Pores: " << changeSet.size() << endl;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}

//*************************************************************************************************

void PoreMorphology::export_ppm_stacks(string foldername) {
  //  srand(time(NULL));
  srand(0);

  //  if(!poreMorphologyCreated)
  //  {cout << "\nWARNING: Can't export ppm stacks, no Pore Morphology
  //  created!\n"; return;}

  cout << "\nExporting as ppm stacks to " << foldername << " ...\n";

  VoxelVolume<Vector3ui8> colorVolume;
  colorVolume.s = morphologyVolume.s;
  colorVolume.spacing = morphologyVolume.spacing;
  colorVolume.voxelValues.resize(colorVolume.s.cast<size_t>().prod(),
                                 Vector3ui8(0, 0, 0));

  vector<size_t> colorShuffle;
  map<size_t, size_t> masterPoreToColor;

  size_t dummyCounter = 0;
  for (auto const &parentAndVoxelIndex : parentToVoxelIndex) {
    uint32_t parent = parentAndVoxelIndex.first;
    masterPoreToColor[parent] = dummyCounter;
    colorShuffle.push_back(dummyCounter);
    ++dummyCounter;
  }

  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  //  unsigned seed = 0;
  shuffle(colorShuffle.begin(), colorShuffle.end(),
          std::default_random_engine(seed));

#pragma omp parallel for
  for (size_t n = 0; n < morphologyVolume.voxelValues.size(); ++n)
    if (get_flag(morphologyVolume[n]) != backgroundValue) {
      if (get_flag(morphologyVolume[n]) == throatValue) {
        colorVolume[n] = Vector3ui8(127, 127, 127);
      } else {
        size_t colorID = get_parent(morphologyVolume[n]);

        colorID = colorShuffle[masterPoreToColor[colorID]];

        size_t r = (min(colorID, colorShuffle.size() - colorID) * 512) /
                   colorShuffle.size();
        while (r > 255)
          --r;

        size_t dummyMasterID =
            (colorID + colorShuffle.size() / 3) % colorShuffle.size();
        size_t g =
            (min(dummyMasterID, colorShuffle.size() - dummyMasterID) * 512) /
            colorShuffle.size();
        while (g > 255)
          --g;

        dummyMasterID =
            (colorID + (colorShuffle.size() * 2) / 3) % colorShuffle.size();
        size_t b =
            (min(dummyMasterID, colorShuffle.size() - dummyMasterID) * 512) /
            colorShuffle.size();
        while (b > 255)
          --b;

        Vector3ui8 someColor(r, g, b);

        if (colorShuffle.size() == 1)
          someColor << 0, 0, 255;

        colorVolume[n] = someColor;
      }
    }

  for (int k = 0; k < colorVolume.s(2); ++k) {
    vector<uint8_t> currImage(colorVolume.s(1) * colorVolume.s(0) * 3, 0);

#pragma omp parallel for
    for (int j = 0; j < colorVolume.s(1); ++j) {
      auto pxIt = currImage.begin() +
                  3 * (colorVolume.spacing(1) * (colorVolume.s(1) - 1 - j));

      for (auto vxIt = colorVolume.voxelValues.begin() +
                       colorVolume.spacing.dot(Vector3i(0, j, k));
           vxIt != colorVolume.voxelValues.begin() +
                       colorVolume.spacing.dot(Vector3i(0, j + 1, k));
           ++vxIt, pxIt += 3) {
        *pxIt = (*vxIt)(0);
        *(pxIt + 1) = (*vxIt)(1);
        *(pxIt + 2) = (*vxIt)(2);
      }
    }

    char numberBuffer[64];
    sprintf(numberBuffer, "%06i", k);
    ofstream myfile(foldername + "stack" + numberBuffer + ".ppm");
    myfile << "P6\n"
           << colorVolume.s(0) << " " << colorVolume.s(1) << endl
           << 255 << endl;
    myfile.write((const char *)currImage.data(), currImage.size());
  }
}

//*************************************************************************************************

void PoreMorphology::reduce_throat_volume() {

  if (!poreMorphologyCreated) {
    cout << "\nCreate Pore Morphology first!\n";
    return;
  }

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  cout << "\nReducing Throat Volume ...\n";

  Vector3i const &s = morphologyVolume.s;

  VoxelVolume<uint8_t> throatVoxelVolume;
  throatVoxelVolume.s = s;
  throatVoxelVolume.spacing = morphologyVolume.spacing;
  throatVoxelVolume.voxelValues.resize(s.cast<size_t>().prod(), 0);

  vector<size_t> throatVoxelsToSeperate;
  throatVoxelsToSeperate.reserve(
      static_cast<size_t>(sqrt(s.cast<float>().prod())));
  for (size_t voxelIndex = 0; voxelIndex < morphologyVolume.voxelValues.size();
       ++voxelIndex)
    if (get_flag(morphologyVolume[voxelIndex]) == throatValue) {
      throatVoxelsToSeperate.push_back(voxelIndex);
      throatVoxelVolume[voxelIndex] = 1;
    }

  vector<vector<size_t>> throatsAndConnectedVoxels;
  throatsAndConnectedVoxels.reserve(sqrt(throatVoxelsToSeperate.size()));

  while (throatVoxelsToSeperate.size() != 0) {
    //    cout << endl << throatVoxels.size();
    if (throatVoxelVolume[throatVoxelsToSeperate.back()] == 0) {
      throatVoxelsToSeperate.pop_back();
      continue;
    }

    size_t floodFillSeed = throatVoxelsToSeperate.back();
    vector<size_t> floodFillRegion(1, floodFillSeed);
    floodFillRegion.reserve(throatVoxelsToSeperate.size());

    throatVoxelVolume[floodFillSeed] = 0;

    vector<size_t> floodFillStack = floodFillRegion;
    floodFillStack.reserve(throatVoxelsToSeperate.size());

    while (floodFillStack.size() != 0) {
      Vector3i coordinate = morphologyVolume.vxID_to_vx(floodFillStack.back());
      floodFillStack.pop_back();

      for (int K = -1; K <= 1; ++K)
        for (int J = -1; J <= 1; ++J)
          for (int I = -1; I <= 1; ++I) {
            Vector3i neighborCoordinate = coordinate + Vector3i(I, J, K);
            size_t neighborIndex =
                morphologyVolume.vx_to_vxID(neighborCoordinate);
            if (throatVoxelVolume[neighborIndex] == 0)
              continue;

            floodFillRegion.push_back(neighborIndex);
            floodFillStack.push_back(neighborIndex);
            throatVoxelVolume[neighborIndex] = 0;
          }
    }

    throatsAndConnectedVoxels.push_back(floodFillRegion);
  }

  if (!parallelFlag)
    omp_set_num_threads(1);

#pragma omp parallel for
  for (size_t throatID = 0; throatID < throatsAndConnectedVoxels.size();
       ++throatID) {

    //    cout << endl << throatID << endl;

    set<size_t, function<bool(size_t const &, size_t const &)>> throatVoxels(
        [&](size_t const &a, size_t const &b) -> bool {
          return (*distanceFieldP)[a] == (*distanceFieldP)[b]
                     ? (a > b)
                     : (*distanceFieldP)[a] > (*distanceFieldP)[b];
        });

    //    cout << "\nset defined\n";

    vector<size_t> &throatVoxelVector = throatsAndConnectedVoxels[throatID];

    //    cout << "\nvector size: " << throatVoxelVector.size();
    //    cout << endl;

    sort(throatVoxelVector.begin(), throatVoxelVector.end(),
         [&](size_t const &a, size_t const &b) {
           return (*distanceFieldP)[a] == (*distanceFieldP)[b]
                      ? (a > b)
                      : (*distanceFieldP)[a] > (*distanceFieldP)[b];
         });

    //    cout << "\nsorted\n";

    throatVoxels.insert(throatVoxelVector.begin(), throatVoxelVector.end());

    //    cout << "\nThroat Voxels inserted\n";

    // shrink connected region according to watershed logic
    auto indexIterator = throatVoxels.begin();
    while (indexIterator != throatVoxels.end()) {

      size_t vxID = *indexIterator;

      //      cout << endl << vxID << endl;

      Vector3i coordinate = morphologyVolume.vxID_to_vx(vxID);
      if (get_flag(morphologyVolume(coordinate)) != throatValue) {
        cout << "\n!\n";
      }

      //    cout << endl << (*distanceFieldP)[vxID] << endl;

      // check for neighbours
      uint32_t neighbourValue = 0;
      bool neighborFound = false;

      // assume that throat will be changed
      bool changeThroat = true;
      bool hasThroatNeighbor = false;

      // check all neighbours
      for (int K = (coordinate(2) == 0 ? 0 : -1);
           changeThroat &&
           K <= (coordinate(2) == morphologyVolume.s(2) - 1 ? 0 : 1);
           ++K)
        for (int J = (coordinate(1) == 0 ? 0 : -1);
             changeThroat &&
             J <= (coordinate(1) == morphologyVolume.s(1) - 1 ? 0 : 1);
             ++J)
          for (int I = (coordinate(0) == 0 ? 0 : -1);
               changeThroat &&
               I <= (coordinate(0) == morphologyVolume.s(0) - 1 ? 0 : 1);
               ++I) {
            if (K == 0 && J == 0 && I == 0)
              continue;

            Vector3i checkVx = coordinate + Vector3i(I, J, K);
            size_t checkVxID = checkVx.cast<size_t>().dot(
                morphologyVolume.spacing.cast<size_t>());

            // ignore if checkVx has throatValue or belongs to background
            if (get_flag(morphologyVolume[checkVxID]) == backgroundValue)
              continue;

            if (get_flag(morphologyVolume[checkVxID]) == throatValue) {
              hasThroatNeighbor = true;
              continue;
            }

            // try to change current throat to value of first checkVx which
            // belongs to a pore
            uint32_t otherPoreID = get_parent(morphologyVolume[checkVxID]);

            if (!neighborFound) {
              neighbourValue = otherPoreID;
              neighborFound = true;
              continue;
            }
            // if throat value is going to be changed, and same pore is hit, do
            // nothing
            else if (otherPoreID == neighbourValue)
              continue;
            // don't change throat if two different neighbouring pores
            else
              changeThroat = false;
          }

      //      if(!hasThroatNeighbor && !neighborFound)
      //      {
      //#pragma omp critical
      //        cout << "\nbad! removing throat voxel enclosed by material.\n";
      //        throatVoxels.erase(indexIterator);
      //        indexIterator = throatVoxels.begin();
      //        morphologyVolume[vxID] = backgroundValue;
      //        continue;
      //      }

      // only throat voxels and material in vicinity. try next
      if (!neighborFound) {
        ++indexIterator;
        continue;
      }

      // current throat voxel will be processed. remove from connected voxels
      throatVoxels.erase(indexIterator);
      indexIterator = throatVoxels.begin();

      // two different pores as neighbours. remove from connected voxels
      if (!changeThroat) { /*cout << "\nactual throat voxel\n";*/
        continue;
      }

      // change current throat voxel to neighbouring value
      morphologyVolume(coordinate) = neighbourValue + enclosedValue;
      //        cout << "\nwatersheded throat voxel\n";
    }
  }

  throatsReduced = true;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}

//*************************************************************************************************

void PoreMorphology::create_pore_morphology(float rMinMaster, float rMinBall) {

  if (!distanceFieldP->distanceFieldCreated) {
    cout << "\nWARNING: Can't create Pore Morphology, no Distance Field "
            "created!\n";
    return;
  }

  DistanceField const &distanceField = *distanceFieldP;
  Vector3i const &s = distanceField.s;

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  cout << "\nCreating Pore Morphology:\n";

  morphologyVolume.s = distanceFieldP->s;
  morphologyVolume.spacing = distanceFieldP->spacing;

  morphologyVolume.voxelValues.clear();
  morphologyVolume.voxelValues.resize(morphologyVolume.s.cast<size_t>().prod(),
                                      backgroundValue);

  // each voxel in the void space is its own master
  size_t voidVoxels = 0;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    if (distanceField[n] > rMinBall) {
      ++voidVoxels;
      morphologyVolume[n] = initValue;
    }

  if (voidVoxels == 0) {
    cout << "\nnothing to do\n";
    return;
  }

  vector<size_t> processingOrder;
  processingOrder.clear();
  processingOrder.reserve(distanceField.voxelValues.size() / 16);

  float r_max;
#ifdef ENABLE_GNU_PARALLEL
  r_max = *(__gnu_parallel::max_element(distanceField.voxelValues.begin(),
                                        distanceField.voxelValues.end()));
#else
  r_max = *(max_element(distanceField.voxelValues.begin(),
                        distanceField.voxelValues.end()));
#endif

  float r_infimum = r_max;

  while (r_infimum != 0.0) {
    if (r_max <= 2.0 || r_max <= 2.0 * rMinBall)
      r_infimum = rMinBall;
    else
      r_infimum = r_max / 2.0;

    //    r_infimum = 0.0;

    processingOrder.clear();
    for (size_t index = 0; index < s.prod(); ++index)
      if (distanceField[index] > r_infimum && distanceField[index] <= r_max &&
          (get_flag(morphologyVolume[index]) == initValue))
        processingOrder.push_back(index);

    //    processingOrder.shrink_to_fit();

    cout << "Pores: " << parentToVoxelIndex.size() << endl;
    cout << scientific << r_infimum << " < r <= " << r_max << endl;

#ifdef ENABLE_GNU_PARALLEL
    __gnu_parallel::sort(processingOrder.begin(), processingOrder.end(),
                         [&](size_t const &i, size_t const &j) {
                           return distanceField[i] == distanceField[j]
                                      ? i < j // for reproducibility
                                      : distanceField[i] > distanceField[j];
                         });
#else
    sort(processingOrder.begin(), processingOrder.end(),
         [&](size_t const &i, size_t const &j) {
           return distanceField[i] == distanceField[j]
                      ? i < j // for reproducibility
                      : distanceField[i] > distanceField[j];
         });
#endif

    if (processingOrder.size() == 0) {
      r_max = r_infimum;
      continue;
    }

    //    size_t progressCounter = 0;
    //    size_t forLoopCounter=0;
    for (auto const &voxelIndex_i : processingOrder) {

      float const &r_i = distanceField[voxelIndex_i];

      //      while(progressCounter <=
      //      (forLoopCounter*100)/processingOrder.size())
      //      {
      //        cout << scientific << progressCounter << " %,\tr=" << r_i
      //             << ",\tpores: " << parentToVoxelIndex.size() << endl;
      //        ++progressCounter;
      //      }
      //      ++forLoopCounter;

      //    if(roundedR_i<omp_get_num_threads())
      //      omp_set_num_threads(1);

      uint32_t &morphologyValue_i = morphologyVolume[voxelIndex_i];
      uint32_t flag_i = get_flag(morphologyValue_i);
      uint32_t parent_i = get_parent(morphologyValue_i);

      // cases: throat, enclosed
      if (flag_i != initValue)
        continue;

      if (1) {
        // ignored voxels due to unlucky inclusion with epsilon
        if (parent_i == 0)
          if (quick_neighbor_check(voxelIndex_i))
            continue;

        // morphologyValue_i may be changed by quick neighbor check
        flag_i = get_flag(morphologyValue_i);
        parent_i = get_parent(morphologyValue_i);
      }

      // case: not allowed to be parent
      if (r_i < rMinMaster && flag_i == initValue && parent_i == 0)
        continue;

      // case: parent.
      if (parent_i == 0) {
        ++parentCounter;
        parentToVoxelIndex[parentCounter] = voxelIndex_i;
        morphologyValue_i = parentCounter;
        parent_i = get_parent(morphologyValue_i);
      }

      // ball always encloses itself. Morphology is fixed at this point.
      morphologyValue_i += enclosedValue;
      flag_i = get_flag(morphologyValue_i);

      // check and update neighborhood
      update_neighbors_box(voxelIndex_i);
      //    update_neighbors_flood(voxelIndex_i);
    }

    r_max = r_infimum;
  }

  //  omp_set_num_threads(0);

  // count changed voxels
  size_t ignoredVoxels = 0;
  for (auto &morphologyValue : morphologyVolume.voxelValues)
    if (get_flag(morphologyValue) == initValue) {
      morphologyValue = backgroundValue;
      ++ignoredVoxels;
    }

  poreMorphologyCreated = true;

  cout << "\nIgnored Void Voxel Fraction: " << float(ignoredVoxels) / voidVoxels
       << endl;
  cout << "Pores: " << parentToVoxelIndex.size() << endl;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}

//*************************************************************************************************

void PoreMorphology::update_neighbors_flood(size_t const &voxelIndex_i) {
  auto const &s = morphologyVolume.s;
  auto const &distanceField = *distanceFieldP;

  Vector3i const voxelCoordinate_i = morphologyVolume.vxID_to_vx(voxelIndex_i);

  uint32_t const &morphologyValue_i = morphologyVolume[voxelIndex_i];
  uint32_t parent_i = get_parent(morphologyValue_i);

  float const &r_i = distanceField[voxelIndex_i];
  long const roundedR_i = floor(r_i);
  float const &r_i_squared = r_i * r_i;

  vector<Vector3i> floodStack;
  floodStack.reserve(pow(2 * roundedR_i + 1, 3));
  floodStack.push_back(voxelCoordinate_i);

  VoxelVolume<uint8_t> processedVoxels;
  processedVoxels.s =
      Vector3i(2 * roundedR_i + 1, 2 * roundedR_i + 1, 2 * roundedR_i + 1);
  processedVoxels.spacing = Vector3i(
      1, processedVoxels.s(0), processedVoxels.s(0) * processedVoxels.s(1));
  processedVoxels.voxelValues.resize(processedVoxels.s.prod(), 0);
  processedVoxels(roundedR_i, roundedR_i, roundedR_i) = 1;

  for (size_t n = 0; n != floodStack.size(); ++n) {

    Vector3i const voxelCoordinate_j_old = floodStack[n];

    for (long neighborIndex = 0; neighborIndex < 6; ++neighborIndex) {
      Vector3i voxelCoordinate_j = voxelCoordinate_j_old;
      voxelCoordinate_j(neighborIndex / 2) += neighborIndex % 2 ? 1 : -1;

      Vector3i d_voxel_ij = voxelCoordinate_j - voxelCoordinate_i;

      float r_ij_squared = d_voxel_ij.cast<float>().squaredNorm();

      if (r_ij_squared > r_i_squared || (voxelCoordinate_j.array() < 0).any() ||
          (voxelCoordinate_j.array() >= s.array()).any())
        continue;

      if ((d_voxel_ij.array().abs() > roundedR_i).any() ||
          processedVoxels(d_voxel_ij +
                          Vector3i(roundedR_i, roundedR_i, roundedR_i)))
        continue;

      processedVoxels(d_voxel_ij +
                      Vector3i(roundedR_i, roundedR_i, roundedR_i)) = true;

      size_t const voxelIndex_j =
          morphologyVolume.vx_to_vxID(voxelCoordinate_j);

      uint32_t &morphologyValue_j = morphologyVolume[voxelIndex_j];
      uint32_t flag_j = get_flag(morphologyValue_j);
      uint32_t parent_j = get_parent(morphologyValue_j);

      if (flag_j != initValue && flag_j != throatValue)
        continue;

      float const &r_j = distanceField[voxelIndex_j];

      if (r_j > r_i)
        continue;

      floodStack.push_back(voxelCoordinate_j);

      float r_ij = sqrt(r_ij_squared);

      // change to child if possible
      if (parent_j == 0) {
        morphologyValue_j = parent_i;
        parent_j = parent_i;
        //        cout << endl << "child";
      }

      if (parent_j == parent_i) {

        // try to enclose
        if (r_ij + r_j <= r_i + 0.2 * r_j) {
          morphologyValue_j += enclosedValue;
          flag_j = enclosedValue;
          //          cout << endl << "enclosed";
        }

        continue;
      }

      // some value other than the current master has been written
      // --> mark as throat
      morphologyValue_j += throatValue;
      flag_j = throatValue;
      //      cout << endl << "throat";
    }
  }
}

//*************************************************************************************************

void PoreMorphology::update_neighbors_box(size_t const &voxelIndex_i) {

  auto const &s = morphologyVolume.s;
  auto const &distanceField = *distanceFieldP;

  Vector3i const voxelCoordinate_i = morphologyVolume.vxID_to_vx(voxelIndex_i);

  uint32_t const &morphologyValue_i = morphologyVolume[voxelIndex_i];
  uint32_t parent_i = get_parent(morphologyValue_i);

  float const &r_i = distanceField[voxelIndex_i];
  float const r_i_padded = r_i + 0.5;
  long const roundedR_i_padded = floor(r_i_padded);
  float const &r_i_padded_squared = r_i_padded * r_i_padded;

  for (long K = -roundedR_i_padded; K <= roundedR_i_padded; ++K)
    for (long J = -roundedR_i_padded; J <= roundedR_i_padded; ++J)
      for (long I = -roundedR_i_padded; I <= roundedR_i_padded; ++I) {
        if (K == 0 && J == 0 && I == 0)
          continue;

        Vector3i const voxelCoordinate_j =
            voxelCoordinate_i + Vector3i(I, J, K);
        if ((voxelCoordinate_j.array() < 0).any() ||
            (voxelCoordinate_j.array() >= s.array()).any())
          continue;

        float r_ij_squared =
            (voxelCoordinate_i - voxelCoordinate_j).cast<float>().squaredNorm();
        if (r_ij_squared > r_i_padded_squared)
          continue;

        size_t const voxelIndex_j =
            morphologyVolume.vx_to_vxID(voxelCoordinate_j);

        uint32_t &morphologyValue_j = morphologyVolume[voxelIndex_j];
        uint32_t flag_j = get_flag(morphologyValue_j);
        uint32_t parent_j = get_parent(morphologyValue_j);

        if (flag_j != initValue)
          continue;

        float const &r_j = distanceField[voxelIndex_j];
        float const r_j_padded = r_j + 0.5;
        if (r_j_padded > r_i_padded) { /*cout << "\nblub\n";*/
          continue;
        }

        float r_ij = sqrt(r_ij_squared);

        // change to slave if possible
        if (parent_j == 0) {
          morphologyValue_j = parent_i;
          parent_j = parent_i;
        }

        if (parent_j == parent_i) {

          // try to enclose
          if (r_ij + r_j <= r_i + epsilon * r_j) {
            morphologyValue_j += enclosedValue;
            flag_j = enclosedValue;
          }

          continue;
        }

        // some value other than the current master has been written
        // --> mark as throat
        morphologyValue_j += throatValue;
        flag_j = enclosedValue;
      }
}
