#include "PoreMorphology.h"
#include "Diagnose.h"
#include "Parallel.h"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <set>
#include <tuple>
#include <utility>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using namespace std;
using namespace std::chrono;
using namespace Eigen;
//------------------------------------------------------------------------------
void PoreMorphology::set_from_voxel_neighborhood(size_t i) {
  MorphologyValue &m_i = morphologyVolume[i];
  ASSURE(m_i.parentId == 0, "");
  Vector3l x_i = morphologyVolume.vxID_to_vx(i);
  for (int K = -1; K <= 1; ++K)
    for (int J = -1; J <= 1; ++J)
      for (int I = -1; I <= 1; ++I) {
        // Neighbor j.
        Vector3l x_j = x_i + Vector3l(I, J, K);
        MorphologyValue m_j = morphologyVolume[x_j];

        if (m_j.state == MorphologyValue::BACKGROUND ||
            m_j.state == MorphologyValue::THROAT || m_j.parentId == 0)
          continue; // Neighbor j does not influence i.

        if (m_i.parentId == 0) {
          m_i.parentId = m_j.parentId;
        } else if (m_i.parentId != m_j.parentId) {
          m_i.state = MorphologyValue::THROAT;
          return;
        }
      }
}
//------------------------------------------------------------------------------
void PoreMorphology::create_legacy_volumes(
    VoxelVolume<uint32_t> &_morphologyVolume,
    VoxelVolume<uint8_t> &stateVolume) {
  Vector3l const &s = morphologyVolume.s;
  _morphologyVolume.resize(s);
  stateVolume.resize(s);

  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n) {
    uint32_t flag = morphologyVolume[n].state;
    uint32_t parent = morphologyVolume[n].parentId;

    stateVolume[n] = static_cast<uint8_t>(flag);
    _morphologyVolume[n] = parent;
  }
}
//------------------------------------------------------------------------------
namespace {
//------------------------------------------------------------------------------
enum class SortOrder { Ascending, Descending };
//------------------------------------------------------------------------------
template <SortOrder order> class DistanceFieldCompare {
public:
  explicit DistanceFieldCompare(const DistanceField &d) : d_(d) {}

  bool operator()(size_t i, size_t j) const {
    switch (order) {
    case SortOrder::Ascending:
      return d_[i] == d_[j] ? i < j : d_[i] < d_[j];
    case SortOrder::Descending:
      return d_[i] == d_[j] ? i > j : d_[i] > d_[j];
    }
  }

private:
  const DistanceField &d_;
};
//------------------------------------------------------------------------------
} // namespace
//------------------------------------------------------------------------------
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
    if (morphologyVolume[voxelIndex].state != MorphologyValue::THROAT)
      continue;

    set<uint32_t> throat;
    Vector3l voxelCoordinate = morphologyVolume.vxID_to_vx(voxelIndex);

    for (int K = -1; K <= 1; ++K)
      for (int J = -1; J <= 1; ++J)
        for (int I = -1; I <= 1; ++I) {
          Vector3l neighborCoordinate = voxelCoordinate + Vector3l(I, J, K);
          MorphologyValue neighborMorphology =
              morphologyVolume[neighborCoordinate];
          if (neighborMorphology.state != MorphologyValue::ENCLOSED)
            continue;

          throat.insert(neighborMorphology.parentId);
        }

    if (throat.size() < 2) {
      cout << endl;
      for (int K = -1; K <= 1; ++K) {
        for (int J = -1; J <= 1; ++J) {
          for (int I = -1; I <= 1; ++I) {
            Vector3l neighborCoordinate = voxelCoordinate + Vector3l(I, J, K);
            MorphologyValue neighborMorphology =
                morphologyVolume[neighborCoordinate];
            uint32_t neighborFlag = neighborMorphology.state;
            cout << neighborFlag << " " << neighborMorphology.parentId
                 << "    ";
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
      maxThroatRadii[throat] = distanceField[voxelIndex];
      continue;
    }

    if (maxThroatRadii[throat] < distanceField[voxelIndex])
      maxThroatRadii[throat] = distanceField[voxelIndex];
  }

  // sort ascending
  DistanceFieldCompare<SortOrder::Ascending> cmp(distanceField);
  map<size_t, set<set<uint32_t>>, decltype(cmp)> poreCentersToThroats(cmp);

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
    uint32_t parent = morphologyVolume[poreVoxelIndex].parentId;
    set<set<uint32_t>> connectedThroats = poreCenterToThroats->second;

    // current radius
    float rPore = distanceField[poreVoxelIndex];

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
          distanceField[parentToVoxelIndex.at(neighborParent)] > parentRadius) {
        parentRadius = distanceField[parentToVoxelIndex.at(neighborParent)];
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
  for (MorphologyValue &morphologyValue : morphologyVolume()) {
    uint32_t flag = morphologyValue.state;
    if (flag != MorphologyValue::ENCLOSED)
      continue;

    uint32_t parent = morphologyValue.parentId;

    if (changeSet.count(parent) == 0)
      continue;

    morphologyValue.state = MorphologyValue::ENCLOSED;
    morphologyValue.parentId = changeSet[parent];
  }

  // look for "ghost" throats
  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    if (morphologyVolume[voxelIndex].state != MorphologyValue::THROAT)
      continue;

    set<uint32_t> throat;
    Vector3l voxelCoordinate = morphologyVolume.vxID_to_vx(voxelIndex);

    for (int K = -1; K <= 1; ++K)
      for (int J = -1; J <= 1; ++J)
        for (int I = -1; I <= 1; ++I) {
          Vector3l neighborCoordinate = voxelCoordinate + Vector3l(I, J, K);
          MorphologyValue neighborMorphology =
              morphologyVolume[neighborCoordinate];
          if (neighborMorphology.state == MorphologyValue::ENCLOSED)
            throat.insert(neighborMorphology.parentId);
        }

    if (throat.size() > 1)
      continue;

    morphologyVolume[voxelIndex] = {MorphologyValue::ENCLOSED,
                                    *(throat.begin())};
  }

  cout << "\nRemoved Pores: " << changeSet.size() << endl;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}
//------------------------------------------------------------------------------
void PoreMorphology::export_ppm_stacks(const char *foldername) {
  //  srand(time(NULL));
  srand(0);

  //  if(!poreMorphologyCreated)
  //  {cout << "\nWARNING: Can't export ppm stacks, no Pore Morphology
  //  created!\n"; return;}

  cout << "\nExporting as ppm stacks to " << foldername << " ...\n";

  using Vector3ui8 = Vector3<uint8_t>;

  VoxelVolume<Vector3ui8> colorVolume;
  colorVolume.resize(morphologyVolume.s, Vector3ui8(0, 0, 0));

  vector<size_t> colorShuffle;
  map<size_t, size_t> parentPoreToColor;

  size_t dummyCounter = 0;
  for (auto const &parentAndVoxelIndex : parentToVoxelIndex) {
    uint32_t parent = parentAndVoxelIndex.first;
    parentPoreToColor[parent] = dummyCounter;
    colorShuffle.push_back(dummyCounter);
    ++dummyCounter;
  }

  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  //  unsigned seed = 0;
  shuffle(colorShuffle.begin(), colorShuffle.end(),
          std::default_random_engine(seed));

  parallelFor(morphologyVolume().size(), [&](size_t n) {
    if (morphologyVolume[n].state != MorphologyValue::BACKGROUND) {
      if (morphologyVolume[n].state == MorphologyValue::THROAT) {
        colorVolume[n] = Vector3ui8(127, 127, 127);
      } else {
        size_t colorID = morphologyVolume[n].parentId;

        colorID = colorShuffle[parentPoreToColor[colorID]];

        size_t r = (min(colorID, colorShuffle.size() - colorID) * 512) /
                   colorShuffle.size();
        while (r > 255)
          --r;

        size_t dummyParentID =
            (colorID + colorShuffle.size() / 3) % colorShuffle.size();
        size_t g =
            (min(dummyParentID, colorShuffle.size() - dummyParentID) * 512) /
            colorShuffle.size();
        while (g > 255)
          --g;

        dummyParentID =
            (colorID + (colorShuffle.size() * 2) / 3) % colorShuffle.size();
        size_t b =
            (min(dummyParentID, colorShuffle.size() - dummyParentID) * 512) /
            colorShuffle.size();
        while (b > 255)
          --b;

        Vector3ui8 someColor(r, g, b);

        if (colorShuffle.size() == 1)
          someColor << 0, 0, 255;

        colorVolume[n] = someColor;
      }
    }
  });

  for (int k = 0; k < colorVolume.s(2); ++k) {
    vector<uint8_t> currImage(colorVolume.s(1) * colorVolume.s(0) * 3, 0);

    parallelFor(colorVolume.s(1), [&](long j) {
      auto pxIt = currImage.begin() +
                  3 * (colorVolume.spacing(1) * (colorVolume.s(1) - 1 - j));

      for (auto vxIt = colorVolume().begin() +
                       colorVolume.spacing.dot(Vector3l(0, j, k));
           vxIt != colorVolume().begin() +
                       colorVolume.spacing.dot(Vector3l(0, j + 1, k));
           ++vxIt, pxIt += 3) {
        *pxIt = (*vxIt)(0);
        *(pxIt + 1) = (*vxIt)(1);
        *(pxIt + 2) = (*vxIt)(2);
      }
    });

    char numberBuffer[64];
    snprintf(numberBuffer, 64, "%06i", k);
    ofstream myFile(string(foldername) + "stack" + numberBuffer + ".ppm");
    myFile << "P6\n"
           << colorVolume.s(0) << " " << colorVolume.s(1) << endl
           << 255 << endl;
    myFile.write((const char *)currImage.data(), currImage.size());
  }
}
//------------------------------------------------------------------------------
void PoreMorphology::reduce_throat_volume() {

  if (!poreMorphologyCreated) {
    cout << "\nCreate Pore Morphology first!\n";
    return;
  }

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  cout << "\nReducing Throat Volume ...\n";

  Vector3l const &s = morphologyVolume.s;

  VoxelVolume<uint8_t> throatVoxelVolume;
  throatVoxelVolume.resize(s, 0);

  vector<size_t> throatVoxelsToSeparate;
  throatVoxelsToSeparate.reserve(
      static_cast<size_t>(sqrt(s.cast<float>().prod())));
  for (size_t voxelIndex = 0; voxelIndex < morphologyVolume().size();
       ++voxelIndex)
    if (morphologyVolume[voxelIndex].state == MorphologyValue::THROAT) {
      throatVoxelsToSeparate.push_back(voxelIndex);
      throatVoxelVolume[voxelIndex] = 1;
    }

  vector<vector<size_t>> throatsAndConnectedVoxels;
  throatsAndConnectedVoxels.reserve(sqrt(throatVoxelsToSeparate.size()));

  while (throatVoxelsToSeparate.size() != 0) {
    //    cout << endl << throatVoxels.size();
    if (throatVoxelVolume[throatVoxelsToSeparate.back()] == 0) {
      throatVoxelsToSeparate.pop_back();
      continue;
    }

    size_t floodFillSeed = throatVoxelsToSeparate.back();
    vector<size_t> floodFillRegion(1, floodFillSeed);
    floodFillRegion.reserve(throatVoxelsToSeparate.size());

    throatVoxelVolume[floodFillSeed] = 0;

    vector<size_t> floodFillStack = floodFillRegion;
    floodFillStack.reserve(throatVoxelsToSeparate.size());

    while (floodFillStack.size() != 0) {
      Vector3l coordinate = morphologyVolume.vxID_to_vx(floodFillStack.back());
      floodFillStack.pop_back();

      for (int K = -1; K <= 1; ++K)
        for (int J = -1; J <= 1; ++J)
          for (int I = -1; I <= 1; ++I) {
            Vector3l neighborCoordinate = coordinate + Vector3l(I, J, K);
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

  DistanceFieldCompare<SortOrder::Descending> cmp(distanceField);

  parallelFor(throatsAndConnectedVoxels.size(), [&](size_t throatID) {
    //    cout << endl << throatID << endl;

    set<size_t, decltype(cmp)> throatVoxels(cmp);

    //    cout << "\nset defined\n";

    vector<size_t> &throatVoxelVector = throatsAndConnectedVoxels[throatID];

    //    cout << "\nvector size: " << throatVoxelVector.size();
    //    cout << endl;

    sort(throatVoxelVector.begin(), throatVoxelVector.end(), cmp);

    //    cout << "\nsorted\n";

    throatVoxels.insert(throatVoxelVector.begin(), throatVoxelVector.end());

    //    cout << "\nThroat Voxels inserted\n";

    // shrink connected region according to watershed logic
    auto indexIterator = throatVoxels.begin();
    while (indexIterator != throatVoxels.end()) {

      size_t vxID = *indexIterator;

      //      cout << endl << vxID << endl;

      Vector3l coordinate = morphologyVolume.vxID_to_vx(vxID);
      if (morphologyVolume[coordinate].state != MorphologyValue::THROAT) {
        cout << "\n!\n";
      }

      //    cout << endl << distanceField[vxID] << endl;

      // check for neighbors
      uint32_t neighbourValue = 0;
      bool neighborFound = false;

      // assume that throat will be changed
      bool changeThroat = true;
      // bool hasThroatNeighbor = false;

      // check all neighbors
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

            Vector3l checkVx = coordinate + Vector3l(I, J, K);
            size_t checkVxID = checkVx.cast<size_t>().dot(
                morphologyVolume.spacing.cast<size_t>());

            // ignore if checkVx has MorphologyValue::THROAT or belongs to
            // background
            if (morphologyVolume[checkVxID].state ==
                MorphologyValue::BACKGROUND)
              continue;

            if (morphologyVolume[checkVxID].state == MorphologyValue::THROAT) {
              // hasThroatNeighbor = true;
              continue;
            }

            // try to change current throat to value of first checkVx which
            // belongs to a pore
            uint32_t otherPoreID = morphologyVolume[checkVxID].parentId;

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
      // #pragma omp critical
      //        cout << "\nbad! removing throat voxel enclosed by material.\n";
      //        throatVoxels.erase(indexIterator);
      //        indexIterator = throatVoxels.begin();
      //        morphologyVolume[vxID] = MorphologyValue::BACKGROUND;
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

      // two different pores as neighbors. remove from connected voxels
      if (!changeThroat) { /*cout << "\nactual throat voxel\n";*/
        continue;
      }

      // change current throat voxel to neighbouring value
      morphologyVolume[morphologyVolume.vx_to_vxID(coordinate)] = {
          MorphologyValue::ENCLOSED, neighbourValue};
    }
  });

  throatsReduced = true;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}
//------------------------------------------------------------------------------
void PoreMorphology::create_pore_morphology(float rMinParent, float rMinBall) {

  Vector3l const &s = distanceField.s;

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  cout << "\nCreating Pore Morphology:\n";

  morphologyVolume.resize(s, {MorphologyValue::BACKGROUND, 0});

  // each voxel in the void space is its own parent
  size_t voidVoxels = 0;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    if (distanceField[n] > rMinBall) {
      ++voidVoxels;
      morphologyVolume[n] = {MorphologyValue::INIT, 0};
    }

  cout << "\nVoid space fraction: "
       << double(voidVoxels) / distanceField.s.prod() << "\n";

  if (voidVoxels == 0) {
    cout << "\nnothing to do\n";
    return;
  }

  vector<size_t> processingOrder;
  processingOrder.clear();
  processingOrder.reserve(distanceField().size() / 16);

  float r_max;
  r_max = *(max_element(distanceField().begin(), distanceField().end()));

  VoxelVolume<float> skeletonVolume;
  float skeletonValue;
  bool visualizeSkeleton = !exportSkeletonPath.empty();
  if (visualizeSkeleton) {
    skeletonValue = 2 * log2(1.f + r_max);
    skeletonVolume = distanceField;
    for (auto &v : skeletonVolume.data)
      v = log2(1.f + v);
  }

  float r_infimum = r_max;

  while (r_infimum != 0.0) {
    if (r_max <= 2.0 || r_max <= 2.0 * rMinBall)
      r_infimum = rMinBall;
    else
      r_infimum = r_max / 2.0;

    //    r_infimum = 0.0;

    processingOrder.clear();
    for (size_t index = 0; index < s.cast<size_t>().prod(); ++index)
      if (distanceField[index] > r_infimum && distanceField[index] <= r_max &&
          (morphologyVolume[index].state == MorphologyValue::INIT))
        processingOrder.push_back(index);

    //    processingOrder.shrink_to_fit();

    cout << "Pores: " << parentToVoxelIndex.size() << endl;
    cout << scientific << r_infimum << " < r <= " << r_max << endl;

    parallelSort(processingOrder.begin(), processingOrder.end(),
                 DistanceFieldCompare<SortOrder::Descending>(distanceField));

    if (processingOrder.size() == 0) {
      r_max = r_infimum;
      continue;
    }

    //    size_t progressCounter = 0;
    //    size_t forLoopCounter=0;
    for (size_t i : processingOrder) {

      float const &r_i = distanceField[i];

      //      while(progressCounter <=
      //      (forLoopCounter*100)/processingOrder.size())
      //      {
      //        cout << scientific << progressCounter << " %,\tr=" << r_i
      //             << ",\tpores: " << parentToVoxelIndex.size() << endl;
      //        ++progressCounter;
      //      }
      //      ++forLoopCounter;

      MorphologyValue &m_i = morphologyVolume[i];

      // cases: throat, enclosed
      if (m_i.state != MorphologyValue::INIT)
        continue;

      bool i_from_direct_neighborhood = false;
      if (1) {
        // Due to the epsilon inclusion criterion, it might happen that i was
        // not found by a larger maximal ball j, since j could have been
        // enclosed by an even larger ball k
        if (m_i.parentId == 0) {
          set_from_voxel_neighborhood(i);
          i_from_direct_neighborhood =
              (m_i.state != MorphologyValue::THROAT) && (m_i.parentId != 0);
        }
      }

      // case: not allowed to be parent
      if (r_i < rMinParent && m_i.state == MorphologyValue::INIT &&
          m_i.parentId == 0)
        continue;

      // Basically, one can assume that all voxels which contribute to the
      // morphology are skeleton voxels.
      if (visualizeSkeleton && !i_from_direct_neighborhood)
        skeletonVolume[i] = skeletonValue;

      // case: parent.
      if (m_i.parentId == 0) {
        ++parentCounter;
        parentToVoxelIndex[parentCounter] = i;
        m_i.parentId = parentCounter;
      }

      // ball always encloses itself. Morphology is fixed at this point.
      m_i.state = MorphologyValue::ENCLOSED;

      // check and update neighborhood
      update_neighbors(i);
      //    update_neighbors_flood(i);
    }

    r_max = r_infimum;
  }

  // count changed voxels
  size_t ignoredVoxels = 0;
  for (auto &morphologyValue : morphologyVolume())
    if (morphologyValue.state == MorphologyValue::INIT) {
      morphologyValue.state = MorphologyValue::BACKGROUND;
      ++ignoredVoxels;
    }

  poreMorphologyCreated = true;

  cout << "\nIgnored Void Voxel Fraction: "
       << double(ignoredVoxels) / voidVoxels << endl;
  cout << "Pores: " << parentToVoxelIndex.size() << endl;

  if (visualizeSkeleton)
    skeletonVolume.export_pgm_stacks(exportSkeletonPath.c_str());

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;
}
//------------------------------------------------------------------------------
void PoreMorphology::update_neighbors(size_t i) {

  Vector3l const &s = morphologyVolume.s;
  Vector3l const x_i = morphologyVolume.vxID_to_vx(i);

  MorphologyValue const &m_i = morphologyVolume[i];

  float const &r_i = distanceField[i];
  float const r_i_padded = r_i + 0.5;
  long const roundedR_i_padded = floor(r_i_padded);
  float const r_i_padded_squared = r_i_padded * r_i_padded;

  for (long K = -roundedR_i_padded; K <= roundedR_i_padded; ++K)
    for (long J = -roundedR_i_padded; J <= roundedR_i_padded; ++J)
      for (long I = -roundedR_i_padded; I <= roundedR_i_padded; ++I) {
        if (K == 0 && J == 0 && I == 0)
          continue;

        Vector3l const x_j = x_i + Vector3l(I, J, K);
        if ((x_j.array() < 0).any() || (x_j.array() >= s.array()).any())
          continue;

        float d_ij_squared = (x_i - x_j).cast<float>().squaredNorm();
        if (d_ij_squared > r_i_padded_squared)
          continue;

        size_t const j = morphologyVolume.vx_to_vxID(x_j);

        MorphologyValue &m_j = morphologyVolume[j];

        if (m_j.state != MorphologyValue::INIT)
          continue;

        float const &r_j = distanceField[j];
        float const r_j_padded = r_j + 0.5;
        if (r_j_padded > r_i_padded) { /*cout << "\nblub\n";*/
          continue;
        }

        float d_ij = sqrt(d_ij_squared);

        // update parent if applicable
        if (m_j.parentId == 0) {
          m_j.parentId = m_i.parentId;
        }

        if (m_j.parentId == m_i.parentId) {

          // try to enclose
          if (d_ij + r_j <= r_i + .9f) {
            m_j.state = MorphologyValue::ENCLOSED;
          }

          continue;
        }

        // some value other than the current parent has been written
        // --> mark as throat
        m_j.state = MorphologyValue::THROAT;
      }
}
//------------------------------------------------------------------------------
} // namespace fred
