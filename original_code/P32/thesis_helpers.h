#pragma once

#include <iostream>
#include <fstream>
#include "distancefield.h"
#include <set>
#include "poremorphology.h"

namespace thesis_helpers {

//*************************************************************************************************

  void gnuplot_distance_field_and_maximal_balls( string const & folderName,
                                                 DistanceField const & distanceField )
  {

    ofstream file_gp(folderName + "distance_field.txt");
    set<pair <float, pair <unsigned,unsigned> > > maximalBalls;

    Vector3i const & s = distanceField.s;

    srand(0);
    map<size_t,float> color;
    for(size_t n=0; n<s.prod(); ++n)
      color[n] = float(rand())/RAND_MAX;

    for(unsigned i=0; i<s(0); ++i)
      for(unsigned j=0; j<s(1); ++j){
        file_gp << i << " " << j << " " << distanceField(i,j,0)
                << " " << color.at(distanceField.vx_to_vxID(Vector3i(i,j,0))) << endl;
        if(distanceField(i,j,0) == 0.0)
          continue;
        maximalBalls.insert( make_pair( distanceField(i,j,0), make_pair(s(1)-j,s(0)-i) ) );
      }

    ofstream file2_gp(folderName + "maximal_balls.txt");
    for(auto it=maximalBalls.end(); it--!=maximalBalls.begin();){
      float d = it->first;
      unsigned j = (it->second).first;
      j = s(1)-j;
      unsigned i = (it->second).second;
      i = s(0)-i;
      file2_gp << i << " " << j << " " << d
               << " " << color.at(distanceField.vx_to_vxID(Vector3i(i,j,0))) << endl;
    }

  }

//*************************************************************************************************

  template<typename T>
  void gnuplot_volume_to_images( string const & folderName,
                                 VoxelVolume<T> const & voxelVolume )
  {

    Vector3i const & s = voxelVolume.s;

    for(unsigned k=0; k<s(2); ++k)
    {
      ofstream myFile(folderName + "stack_" + to_string(k) + ".txt");
      for(unsigned j=0; j<s(1); ++j)
        for(unsigned i=0; i<s(0); ++i)
        {
          myFile << i << " " << j << " " << float(voxelVolume(i,j,k)) << endl;
        }
    }
  }

//*************************************************************************************************

  uint32_t get_parent(uint32_t const & morphologyValue)
  {
    static uint32_t morphologyReader = (~0)-(1<<31)-(1<<30);
    return morphologyReader & morphologyValue;
  }

//*************************************************************************************************

  uint32_t get_flag(uint32_t const & morphologyValue)
  {
    static uint32_t flagReader = (1<<31)+(1<<30);
    return flagReader & morphologyValue;
  }

//*************************************************************************************************

  uint32_t initValue{0},
           enclosedValue{ uint32_t(1) << 30 },
           throatValue{ uint32_t(1) << 31 },
           backgroundValue{ (uint32_t(1) << 31) + (uint32_t(1) << 30) };

  unsigned fileCounter=0;

//*************************************************************************************************

  void gnuplot_palette_file(string fileName, PoreMorphology const & poreMorphology)
  {

    auto const & morphologyVolume = poreMorphology.morphologyVolume;
    Vector3i const & s = morphologyVolume.s;
    auto const & parentToVoxelIndex = poreMorphology.parentToVoxelIndex;

    srand(0);
    map<size_t,float> color;
    for(size_t n=0; n<s.prod(); ++n)
      color[n] = float(rand())/RAND_MAX;


    ofstream image_palette_file(fileName);
    for(size_t voxelIndex=0; voxelIndex<s.prod(); ++voxelIndex)
    {
      uint32_t morphologyValue = morphologyVolume[voxelIndex];
      uint32_t flag = get_flag(morphologyValue);
      uint32_t parent = get_parent(morphologyValue);
      Vector3i position = morphologyVolume.vxID_to_vx(voxelIndex);

      if(flag == backgroundValue)
      {
        image_palette_file << position.transpose() << " -1.0\n";
        continue;
      }

      if(flag == throatValue)
      {
        image_palette_file << position.transpose() << " -2.0 \n";
        continue;
      }

      if( parentToVoxelIndex.count(parent) == 0 )
      {
        image_palette_file << position.transpose() << " -3.0 \n";
        continue;
      }

      size_t parentIndex = parentToVoxelIndex.at(parent);
      float x = color[ parentIndex ];

      image_palette_file << position.transpose() << " " << x << endl;
    }
  }

//*************************************************************************************************

  void update_neighbors_box(DistanceField const & distanceField,
                            VoxelVolume<uint32_t> & morphologyVolume,
                            size_t const & voxelIndex_i,
                            map<uint32_t,size_t> const & parentToVoxelIndex)
  {

    auto const & s = morphologyVolume.s;

    double padding=0.5;

    Vector3i const voxelCoordinate_i = morphologyVolume.vxID_to_vx(voxelIndex_i);

    uint32_t const & morphologyValue_i = morphologyVolume[voxelIndex_i];
    uint32_t parent_i = get_parent(morphologyValue_i);

    float const & r_i = distanceField[voxelIndex_i]+padding;
    long const roundedR_i = floor(r_i);
    float const & r_i_squared = r_i*r_i;

    for(long K=-roundedR_i; K<=roundedR_i; ++K)
      for(long J=-roundedR_i; J<=roundedR_i; ++J)
        for(long I=-roundedR_i; I<=roundedR_i; ++I)
        {
          if(K==0 && J==0 && I==0)
            continue;

          Vector3i const voxelCoordinate_j = voxelCoordinate_i + Vector3i(I,J,K);
          if((voxelCoordinate_j.array()<0).any() ||
             (voxelCoordinate_j.array()>=s.array()).any() )
            continue;

          float r_ij_squared = (voxelCoordinate_i-voxelCoordinate_j).cast<float>().squaredNorm();
          if(r_ij_squared>r_i_squared)
            continue;

          size_t const voxelIndex_j = morphologyVolume.vx_to_vxID(voxelCoordinate_j);

          uint32_t & morphologyValue_j = morphologyVolume[voxelIndex_j];
          uint32_t flag_j = get_flag(morphologyValue_j);
          uint32_t parent_j = get_parent(morphologyValue_j);

          if(flag_j != initValue)
            continue;

          float const & r_j = distanceField[voxelIndex_j]+padding;
          if(r_j > r_i)
          { /*cout << "\nblub\n";*/ continue;}

          float r_ij = sqrt(r_ij_squared);

          // change to slave if possible
          if( parent_j == 0 )
          {
            morphologyValue_j = parent_i;
            parent_j = parent_i;
          }

          if(parent_j == parent_i)
          {

            // try to enclose
            if(r_ij+r_j <= r_i+.2*r_j )
            {
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



    // paint

    ofstream activeBall("thesis/steps/"+to_string(fileCounter)+"ball.txt");
    ofstream allBalls("thesis/steps/"+to_string(fileCounter)+"all.txt");

    srand(0);
    map<size_t,float> color;
    for(size_t n=0; n<s.prod(); ++n)
      color[n] = float(rand())/RAND_MAX;

    vector<size_t> processingOrder;

    for(size_t index=0; index<s.prod(); ++index)
      if(distanceField[index]>0)
        processingOrder.push_back(index);

    __gnu_parallel::sort(processingOrder.begin(), processingOrder.end(),
                         [&](size_t const & i, size_t const & j){
      return distanceField[i]==distanceField[j]
          ? i<j // for reproducibility
          : distanceField[i]>distanceField[j];
    });

    for(size_t n=0; n<processingOrder.size(); ++n)
    {
      size_t voxelIndex = processingOrder[n];
      uint32_t morphologyValue = morphologyVolume[voxelIndex];
      uint32_t flag = get_flag(morphologyValue);
      uint32_t parent = get_parent(morphologyValue);
      Vector3i position = morphologyVolume.vxID_to_vx(voxelIndex);

      if(flag == initValue )
      {
        if(parentToVoxelIndex.count(parent) == 0)
          allBalls << position.transpose()
                   << " " << distanceField[voxelIndex]
                      << " " << -2.0 << endl;
        else
        {
          auto parentIndex = parentToVoxelIndex.at(parent);
          allBalls << position.transpose()
                   << " " << distanceField[voxelIndex]
                      << " " << -1.0+color.at(parentIndex)  << endl;
        }

        continue;
      }

      if(flag == throatValue)
      {
        allBalls << position.transpose()
                 << " " << distanceField[voxelIndex]
                 << " " << 2.0 << endl;
        continue;
      }

      auto parentIndex = parentToVoxelIndex.at(parent);
      allBalls << position.transpose()
               << " " << distanceField[voxelIndex]
               << " " << color.at(parentIndex) << endl;
    }

    activeBall << voxelCoordinate_i.transpose()
               << " " << distanceField[voxelIndex_i]
               << " " << color.at(parentToVoxelIndex.at(parent_i)) << endl;

    vector<Vector3i> colors2;
    for(size_t n=0; n<s.prod(); ++n)
      colors2.push_back(Vector3i(rand()%255,rand()%255,rand()%255));


    ofstream image_file("thesis/steps/" + to_string(fileCounter) + "img.txt");

    for(size_t voxelIndex=0; voxelIndex<s.prod(); ++voxelIndex)
    {
      uint32_t morphologyValue = morphologyVolume[voxelIndex];
      uint32_t flag = get_flag(morphologyValue);
      uint32_t parent = get_parent(morphologyValue);
      Vector3i position = morphologyVolume.vxID_to_vx(voxelIndex);

      if(flag == backgroundValue)
      {
        image_file << position.transpose() << " 255 255 255 0\n";
        continue;
      }

      if(flag == throatValue)
      {
        image_file << position.transpose() << " 127 127 127 255\n";
        continue;
      }

      if( parentToVoxelIndex.count(parent) == 0 )
      {
        image_file << position.transpose() << " 0 0 0 255\n";
        continue;
      }

      image_file << position.transpose() << " "
                 << colors2[ parentToVoxelIndex.at(parent) ].transpose() << " 255\n";

    }


// in fig10.gp, "set palette rgb 33,13,10" is used.
// gnuplot: "show palette rgbformulae" explains how values from x=0 to 1 are mapped:
// 33 -> r = |2*x-0.5|
// 13 -> g = sin(180*x)
// 10 -> b = cos(90*x)

    ofstream image_palette_file("thesis/steps/" + to_string(fileCounter) + "img_palette.txt");
    for(size_t voxelIndex=0; voxelIndex<s.prod(); ++voxelIndex)
    {
      uint32_t morphologyValue = morphologyVolume[voxelIndex];
      uint32_t flag = get_flag(morphologyValue);
      uint32_t parent = get_parent(morphologyValue);
      Vector3i position = morphologyVolume.vxID_to_vx(voxelIndex);

      if(flag == backgroundValue)
      {
        image_palette_file << position.transpose() << " -1.0\n";
        continue;
      }

      if(flag == throatValue)
      {
        image_palette_file << position.transpose() << " -2.0 \n";
        continue;
      }

      if( parentToVoxelIndex.count(parent) == 0 )
      {
        image_palette_file << position.transpose() << " -3.0 \n";
        continue;
      }

      size_t parentIndex = parentToVoxelIndex.at(parent);
      float x = color[ parentIndex ];

      image_palette_file << position.transpose() << " " << x << endl;
    }


    ++fileCounter;
  }

//*************************************************************************************************

  void mb_step_by_step(DistanceField const & distanceField,
                       float rMinMaster,
                       float rMinBall){

    uint32_t parentCounter{0};
    map<uint32_t,size_t> parentToVoxelIndex;


    Vector3i const & s = distanceField.s;

    cout << "\nCreating Pore Morphology:\n";

    VoxelVolume<uint32_t> morphologyVolume;

    morphologyVolume.s = distanceField.s;
    morphologyVolume.spacing = distanceField.spacing;

    morphologyVolume.voxelValues.clear();
    morphologyVolume.voxelValues.resize(morphologyVolume.s.prod(),backgroundValue);

    // each voxel in the void space is its own master
    size_t voidVoxels=0;
    for(size_t n=0; n<s.prod(); ++n)
      if(distanceField[n]>rMinBall)
      {
        ++voidVoxels;
        morphologyVolume[n] = initValue;
      }

    if(voidVoxels == 0)
    {
      cout << "\nnothing to do\n";
      return;
    }


    vector<size_t> processingOrder;
    processingOrder.clear(); processingOrder.reserve(distanceField.voxelValues.size()/16);

    float r_max = *(__gnu_parallel::max_element(distanceField.voxelValues.begin(),
                                                distanceField.voxelValues.end()));

    float r_infimum = r_max;

    while(r_infimum != 0.0)
    {
      if(r_max <= 2.0 ||
         r_max <= 2.0*rMinBall)
        r_infimum = rMinBall;
      else
        r_infimum = r_max/2.0;

  //    r_infimum = 0.0;

      processingOrder.clear();
      for(size_t index=0; index<s.prod(); ++index)
        if(distanceField[index] > r_infimum &&
           distanceField[index] <=r_max &&
           get_flag(morphologyVolume[index]) == initValue)
          processingOrder.push_back(index);

  //    processingOrder.shrink_to_fit();

      cout << "Pores: " << parentToVoxelIndex.size() << endl;
      cout << scientific << r_infimum << " < r <= " << r_max << endl;

      __gnu_parallel::sort(processingOrder.begin(), processingOrder.end(),
                           [&](size_t const & i, size_t const & j){
        return distanceField[i]==distanceField[j]
            ? i<j // for reproducibility
            : distanceField[i]>distanceField[j];
      });

      if(processingOrder.size()==0)
      {r_max = r_infimum; continue;}

      size_t progressCounter = 0;
      size_t forLoopCounter=0;
      for(auto const & voxelIndex_i : processingOrder)
      {

        float const & r_i = distanceField[voxelIndex_i];

  //      while(progressCounter <= (forLoopCounter*100)/processingOrder.size())
  //      {
  //        cout << scientific << progressCounter << " %,\tr=" << r_i
  //             << ",\tpores: " << parentToVoxelIndex.size() << endl;
  //        ++progressCounter;
  //      }
  //      ++forLoopCounter;

        //    if(roundedR_i<omp_get_num_threads())
        //      omp_set_num_threads(1);

        uint32_t & morphologyValue_i = morphologyVolume[voxelIndex_i];
        uint32_t flag_i = get_flag(morphologyValue_i);
        uint32_t parent_i = get_parent(morphologyValue_i);

        // cases: throat, enclosed
        if(flag_i != initValue)
          continue;

        // case: not allowed to be parent
        if(r_i<rMinMaster && flag_i==initValue && parent_i==0)
          continue;

        // case: parent.
        if(parent_i == 0)
        {
          ++parentCounter;
          parentToVoxelIndex[parentCounter] = voxelIndex_i;
          morphologyValue_i = parentCounter;
          parent_i = get_parent(morphologyValue_i);
        }

        // ball always encloses itself. Morphology is fixed at this point.
        morphologyValue_i += enclosedValue;
        flag_i = get_flag(morphologyValue_i);

        // check and update neighborhood
        update_neighbors_box(distanceField, morphologyVolume, voxelIndex_i,
                             parentToVoxelIndex);
        //    update_neighbors_flood(voxelIndex_i);
      }

      r_max = r_infimum;
    }


  //  omp_set_num_threads(0);

    // count changed voxels
    size_t ignoredVoxels=0;
    for(auto & morphologyValue : morphologyVolume.voxelValues)
      if( get_flag(morphologyValue)==initValue )
      {morphologyValue = backgroundValue; ++ignoredVoxels;}

    cout << "\nIgnored Void Voxel Fraction: " << float(ignoredVoxels)/voidVoxels << endl;
    cout << "Pores: " << parentToVoxelIndex.size() << endl;

  }

}
