#include "morphologystatistics.h"
#include <algorithm>
#include "logarithmic_bins.h"
#include <parallel/algorithm>

using namespace std;
using namespace Eigen;

namespace morphologyStatistics{

  map<PoreID,Pore> extract_pores(VoxelVolume<PoreID> const & morphologyVolume,
                                 VoxelVolume<uint8_t> const & stateVolume,
                                 bool includeBoundaryPores)
  {

    cout << "\nExtracting pores ...\n";

    map<PoreID,Pore> pores;

    for(VoxelIndex n=0; n<stateVolume.s.cast<size_t>().prod(); ++n)
    {

      if(stateVolume[n] == BACKGROUND ||
         stateVolume[n] == THROAT )
      {continue;}

      PoreID poreID = morphologyVolume[n];

      (pores[poreID]).push_back(n);
    }

    if(includeBoundaryPores)
      return pores;

    // remove boundary pores
    auto const & s = stateVolume.s;
    set<PoreID> removedPores;
    for(long sideIndex=0; sideIndex<6; ++sideIndex)
    {
      long axisIndex = sideIndex/2;
      long isTopSide = sideIndex%2;

      VoxelPosition vx(0,0,0);

      vx(axisIndex) = (s(axisIndex)-1)*isTopSide;

      Vector2i iterateDimensions((axisIndex+1)%3,(axisIndex+2)%3);
      for(vx(iterateDimensions(0))=0; vx(iterateDimensions(0))<s(iterateDimensions(0)); ++vx(iterateDimensions(0)))
        for(vx(iterateDimensions(1))=0; vx(iterateDimensions(1))<s(iterateDimensions(1)); ++vx(iterateDimensions(1)))
        {
          PoreID poreID = morphologyVolume(vx);
          pores.erase(poreID);
          removedPores.insert(poreID);
        }
    }

    cout << "\nRemoved boundary pores:\n" << removedPores.size() << endl;

    return pores;

  }

//*************************************************************************************************

  map<ThroatID,Throat> extract_throats(const VoxelVolume<PoreID> &morphologyVolume,
                                       VoxelVolume<uint8_t> const & stateVolume)
  {

    cout << "\nExtracting throats ...\n";

    map<ThroatID,Throat> throats;

    for(VoxelIndex n=0; n<stateVolume.s.prod(); ++n)
    {
      if( stateVolume[n] != THROAT )
        continue;

      VoxelPosition throatCoordinate = stateVolume.to_vx(n);
      ThroatID throatID;

      for(long K=-1; K<=1; ++K)
        for(long J=-1; J<=1; ++J)
          for(long I=-1; I<=1; ++I)
          {
            VoxelPosition neighborCoordinate = throatCoordinate + VoxelPosition(I,J,K);
            VoxelIndex neighborIndex = stateVolume.to_vxID(neighborCoordinate);

            if( stateVolume[neighborIndex] == THROAT ||
                stateVolume[neighborIndex] == BACKGROUND )
              continue;

            throatID.insert( morphologyVolume[neighborIndex] );
          }

      if(throatID.size()<=1)
        cout << "\nWARNING: Throat Voxel does not connect at least two Pores\n";

      throats[throatID].push_back(n);

    }

    return throats;

  }

//*************************************************************************************************

  float effective_radius(float volume)
  {
    return pow(0.75*volume/M_PI, 1.0/3.0);
  }

//*************************************************************************************************

  std::map<PoreID,float> calculate_pore_volumes(std::map<PoreID,Pore> const & pores)
  {

    cout << "\nCalculating pore volumes ...\n";
    map<PoreID,float> poreVolumes;

    for(auto const & IDAndPore : pores)
      poreVolumes[IDAndPore.first] = volume_of_(IDAndPore.second);

    return poreVolumes;

  }

//*************************************************************************************************

  std::map<PoreID,float> calculate_effective_radii(std::map<PoreID,float> const & poreVolumes)
  {
    cout << "\nCalculating effective radii ...\n";
    map<PoreID,float> poreRadii = poreVolumes;

    for(auto & IDAndRadius : poreRadii )
      IDAndRadius.second = effective_radius(IDAndRadius.second);

    return poreRadii;
  }

//*************************************************************************************************

  std::map<ThroatID,float> calculate_throat_volumes(const std::map<ThroatID, Throat> &throats)
  {
    cout << "\nCalculating Throat Volumes ...\n";
    map<ThroatID,float> throatVolumes;

    for(auto const & IDAndThroat : throats)
      throatVolumes[IDAndThroat.first] = volume_of_(IDAndThroat.second);

    return throatVolumes;
  }

//*************************************************************************************************

  float volume_of_(Pore const & poreOrThroat)
  {

    double volume=0.0;
    for( auto const & index : poreOrThroat )
      volume += 1.0;

    return volume;

  }


//*************************************************************************************************

  void cout_void_fractions(std::map<PoreID,Pore> const & pores,
                           const std::map<ThroatID, Throat> &throats,
                           VoxelVolume<float> const & distanceField,
                           float r_small_large)
  {

    Vector3i const & s = distanceField.s;

    double smallPoreVolume=0.0;
    double largePoreVolume=0.0;
    double voidSpaceVolume=0.0;

    for(auto const & poreIDAndPore : pores)
    {
      double poreVolume = volume_of_(poreIDAndPore.second);
      double poreRadius = effective_radius(poreVolume);
      if(poreRadius<r_small_large)
        smallPoreVolume += poreVolume;
      else
        largePoreVolume += poreVolume;
    }

    for(auto const & throatIDAndThroat : throats)
    {
      double throatVolume=volume_of_(throatIDAndThroat.second);
      for(auto const & pore : throatIDAndThroat.first)
      {
        double poreVolume = volume_of_(pores.at(pore));
        double poreRadius = effective_radius(poreVolume);
        if( poreRadius<r_small_large )
          smallPoreVolume += throatVolume/throatIDAndThroat.first.size();
        else
          largePoreVolume += throatVolume/throatIDAndThroat.first.size();
      }
    }

    for(auto const & value : distanceField.voxelValues)
      if(value>0.0)
        voidSpaceVolume += 1.0;

    cout << scientific << "\nsmall pore porosity " << smallPoreVolume/s.prod() << endl;
    cout << scientific << "\nlarge pore porosity " << largePoreVolume/s.prod() << endl;
    cout << scientific << "\nactual porosity " << voidSpaceVolume/s.prod() << endl << endl;


//    size_t ignoredCounter=0;
//    for(size_t n=0; n<stateVolume.voxelValues.size();++n)
//      if(stateVolume[n]==BACKGROUND && distanceField[n]>0)
//        ++ignoredCounter;




  }

//*************************************************************************************************


  Vector3f center_of_mass(Pore const & pore,
                          VoxelVolume<float> const & distanceField)
  {

    Vector3d center_of_mass(Vector3d::Zero());

    for( auto const & index : pore )
      center_of_mass += distanceField.to_vx(index).cast<double>();

    center_of_mass /= double(pore.size());

    return center_of_mass.cast<float>();

  }

//*************************************************************************************************

  Matrix3f covariance_matrix(Pore const & pore,
                             Vector3f const & centerOfMass,
                             VoxelVolume<float> const & distanceField)
  {

    Matrix3d covarianceMatrix(Matrix3d::Zero());

    auto const centerOfMassDouble = centerOfMass.cast<double>();

    for( auto const & index : pore )
    {
      Vector3d displacement = distanceField.to_vx(index).cast<double>() - centerOfMassDouble;
      covarianceMatrix += displacement * displacement.transpose();
    }

    covarianceMatrix /= double(pore.size());

    return covarianceMatrix.cast<float>();

//    // solve for eigenvalues
//    SelfAdjointEigenSolver<Matrix3f> eigensolver(principalComponents);
//    if( eigensolver.info() != Success )
//      cout << "\nNo Eigenvalues found. Shouldn't happen.\n";
//    Vector3f eigenValues = eigensolver.eigenvalues();
//    Matrix3f eigenVectors = eigensolver.eigenvectors();


  }

//*************************************************************************************************

  void logarithmic_histogram_r_eff(std::string filename,
                                         std::map<PoreID,float> const & poreRadii,
                                         long bins)
  {

    auto poreMinMax = minmax_element(poreRadii.begin(), poreRadii.end(),
                                     [](pair<PoreID,float> const & eleA,
                                        pair<PoreID,float> const & eleB){
                                        return eleA.second < eleB.second; } );
    float minRadius = (*poreMinMax.first).second;
    float maxRadius = (*poreMinMax.second).second;

    if(minRadius == maxRadius)
    {cout << "\nNo statistics: min radius is max radius!\n"; return;}

    cout << "\nExporting logarithmic histogram: Counts vs effective voxel radius\n";

    ofstream myfile(filename);

    myfile << "# bin_minimum \t bin_supremum \t bin_height\n";

    map<float,long> binsAndCounts;

    for(int binID=0; binID<bins; ++binID)
      binsAndCounts[ minRadius*pow(maxRadius/minRadius,float(binID)/float(bins)) ] = 0;

    for(auto const & IDAndRadius : poreRadii)
      ++((*(--binsAndCounts.upper_bound(IDAndRadius.second))).second);

    for(auto itBinAndCount = binsAndCounts.begin(); itBinAndCount != binsAndCounts.end();
        ++ itBinAndCount)
    {
      auto itNextBin = itBinAndCount;
      ++itNextBin;

      myfile << ((*itBinAndCount).first) << "\t"
             << ( itNextBin == binsAndCounts.end()
                               ? (maxRadius)
                               : ((*itNextBin).first) ) << "\t"
             << (*itBinAndCount).second << endl;
    }

  }

//*************************************************************************************************

  void histogramm_r_eff(std::string filename,
                        std::map<PoreID,float> const & poreRadii,
                        long bins)
  {
    auto poreMinMax = minmax_element(poreRadii.begin(), poreRadii.end(),
                                     [](pair<PoreID,float> const & eleA,
                                        pair<PoreID,float> const & eleB){
                                        return eleA.second < eleB.second; } );
    float minRadius = ((*poreMinMax.first).second);
    float maxRadius = ((*poreMinMax.second).second);

    if(minRadius == maxRadius)
    {cout << "\nNo statistics: min volume is max volume!\n"; return;}

    cout << "\nExporting histogram: Counts vs effective voxel radius\n";

    ofstream myfile(filename);

    myfile << "# bin_minimum \t bin_supremum \t bin_height\n";

    map<float,long> binsAndCounts;

    for(int binID=0; binID<bins; ++binID)
      binsAndCounts[ minRadius+(maxRadius-minRadius)*float(binID)/float(bins) ] = 0;

    for(auto const & IDAndVolume : poreRadii)
      ++((*(--binsAndCounts.upper_bound((IDAndVolume.second)))).second);

    for(auto itBinAndCount = binsAndCounts.begin(); itBinAndCount != binsAndCounts.end();
        ++ itBinAndCount)
    {
      auto itNextBin = itBinAndCount;
      ++itNextBin;

      myfile << (*itBinAndCount).first << "\t"
             << ( itNextBin == binsAndCounts.end()
                               ? maxRadius
                               : (*itNextBin).first )  << "\t"
             << (*itBinAndCount).second << endl;
    }
  }

//*************************************************************************************************

  std::map<PoreID,std::vector<PoreID> > pore_network(std::map<PoreID,Pore> const & pores,
                                                     std::map<ThroatID,Throat> const & throats)
  {
    cout << "\nCreating pore to pore network...\n";

    map<PoreID,vector<PoreID> > poreNetwork;

    for( auto const & poreIDAndVoxels : pores )
      poreNetwork[poreIDAndVoxels.first] = vector<PoreID>();

    for( auto const & throatIDAndVoxels : throats )
    {
      auto connectedPores = throatIDAndVoxels.first;

      if(connectedPores.size()<=1){
        cout << "\nNon connecting Throat, this shouldn't happen\n";
        return map<PoreID,vector<PoreID> >();
      }

      for(auto itPoreA = connectedPores.begin(); itPoreA!=connectedPores.end(); ++itPoreA)
      {
        auto itPoreB = itPoreA;
        ++itPoreB;
        if(itPoreB == connectedPores.end())
          break;

        for(; itPoreB != connectedPores.end(); ++itPoreB)
        {
          if( poreNetwork.find(*itPoreA) != poreNetwork.end() )
            (poreNetwork.at(*itPoreA)).push_back(*itPoreB);
          if( poreNetwork.find(*itPoreB) != poreNetwork.end() )
            (poreNetwork.at(*itPoreB)).push_back(*itPoreA);
        }
      }

    }

    return poreNetwork;

  }

//*************************************************************************************************

  void average_connectivity_vs_effective_radius(std::string filename,
                                                std::map<PoreID,float> const & poreRadii,
                                                const std::map<PoreID, std::vector<PoreID> > &poreNetwork,
                                                long bins)
  {

    auto poreMinMax = minmax_element(poreRadii.begin(), poreRadii.end(),
                                     [](pair<PoreID,float> const & eleA,
                                        pair<PoreID,float> const & eleB){
                                        return eleA.second < eleB.second; } );

    float minRadius = (*poreMinMax.first).second;
    float maxRadius = (*poreMinMax.second).second;

    if(minRadius == maxRadius)
    {cout << "\nNo statistics: min radius is max radius!\n"; return;}

    cout << "\nExporting: Average Connectivity vs effective voxel radius ...\n";

    ofstream myfile(filename);

    myfile << "# bin_minimum \t bin_supremum \t average_connectivity\n";

    map<float,vector<unsigned> > binsAndConnectivities;

    for(int binID=0; binID<bins; ++binID)
      binsAndConnectivities[ minRadius*pow(maxRadius/minRadius,float(binID)/float(bins)) ]
          = vector<unsigned>();

    for(auto const & IDAndVolume : poreRadii)
    {
      PoreID const & poreID = IDAndVolume.first;
      float const & poreVolume = IDAndVolume.second;

      auto binIterator = --binsAndConnectivities.upper_bound(poreVolume);

      unsigned poreConnectivity = (poreNetwork.at(poreID)).size();

      (*binIterator).second.push_back(poreConnectivity);
    }

    for(auto itBinAndConnectivities = binsAndConnectivities.begin(); itBinAndConnectivities != binsAndConnectivities.end();
        ++ itBinAndConnectivities)
    {
      auto itNextBin = itBinAndConnectivities;
      ++itNextBin;

      myfile << ((*itBinAndConnectivities).first) << "\t"
             << ( itNextBin == binsAndConnectivities.end()
                               ? (maxRadius)
                               : ((*itNextBin).first) ) << "\t";

      float averageConnectivity = 0.0;
      for(auto const & connectivity : (*itBinAndConnectivities).second)
        averageConnectivity += connectivity;

      if((*itBinAndConnectivities).second.size() != 0)
        averageConnectivity /= float( (*itBinAndConnectivities).second.size() );

      myfile << averageConnectivity << endl;
    }


    return;


    // Extra part: connectivities among pores with r_eff > 6.0
    ofstream myfile2(filename+".special");

    myfile2 << "# bin_minimum \t bin_supremum \t average_connectivity_r_eff_larger_6.0\n";

    binsAndConnectivities.clear();

    for(int binID=0; binID<bins; ++binID)
      binsAndConnectivities[ minRadius*pow(maxRadius/minRadius,float(binID)/float(bins)) ]
          = vector<unsigned>();

    for(auto const & IDAndVolume : poreRadii)
    {
      PoreID const & poreID = IDAndVolume.first;
      float const & poreVolume = IDAndVolume.second;

      if((poreVolume) < 6.0)
        continue;

      auto binIterator = --binsAndConnectivities.upper_bound(poreVolume);

      unsigned poreConnectivity = 0;

      for( auto const & connectedPoreID : poreNetwork.at(poreID) )
        if( (poreRadii.at(connectedPoreID)) >= 6.0 )
          ++poreConnectivity;

      (*binIterator).second.push_back(poreConnectivity);
    }

    for(auto itBinAndConnectivities = binsAndConnectivities.begin(); itBinAndConnectivities != binsAndConnectivities.end();
        ++ itBinAndConnectivities)
    {
      auto itNextBin = itBinAndConnectivities;
      ++itNextBin;

      myfile2 << ((*itBinAndConnectivities).first) << "\t"
             << ( itNextBin == binsAndConnectivities.end()
                               ? (maxRadius)
                               : ((*itNextBin).first) ) << "\t";

      float averageConnectivity = 0.0;
      for(auto const & connectivity : (*itBinAndConnectivities).second)
        averageConnectivity += connectivity;

      if((*itBinAndConnectivities).second.size() != 0)
        averageConnectivity /= float( (*itBinAndConnectivities).second.size() );

      myfile2 << averageConnectivity << endl;
    }



  }

//*************************************************************************************************

  void pore_projections(std::string directory,
                        std::map<PoreID,Pore> const & pores,
                        VoxelVolume<float> const & distanceField)
  {

    Vector3i const & s = distanceField.s;

    cout << "\nCreating projections of pores ...\n";

    // centers of mass

    vector<Vector3f> centersOfMass;

    for(auto const & IDAndPore : pores)
      centersOfMass.push_back( center_of_mass(IDAndPore.second,distanceField) );

    // upper and lower bound of pores

    vector<pair<Vector3i,Vector3i> > corners000and111;
    for(auto const & IDAndPore : pores)
    {
      auto const & pore = IDAndPore.second;
      Vector3i corner000 = s;
      Vector3i corner111 = Vector3i::Zero();
      for( auto const & index : pore )
      {
        Vector3i coordinate = distanceField.to_vx(index);
        for(int dim=0; dim<3; ++dim)
        {
          if(coordinate(dim) < corner000(dim))
            corner000(dim) = coordinate(dim);
          if(coordinate(dim) > corner111(dim))
            corner111(dim) = coordinate(dim);
        }
      }
      corners000and111.push_back(make_pair(corner000,corner111));
    }

    // max size of images
    Vector3i imageXYZ(0,0,0);
    for(long n=0; n<centersOfMass.size(); ++n)
    {
      Vector3i s;
      for(int dim=0; dim<3; ++dim)
      {
        s(dim) = 2*max( round(centersOfMass[n](dim))-corners000and111[n].first(dim),
                        corners000and111[n].second(dim)-round(centersOfMass[n](dim)) ) + 1;

        if( s(dim) > imageXYZ(dim) )
          imageXYZ(dim) = s(dim);
      }
    }

    auto maxCoeffXyz = imageXYZ.maxCoeff();
    imageXYZ = Vector3i(maxCoeffXyz,maxCoeffXyz,maxCoeffXyz);
    Vector3i imageCenter = imageXYZ.array()/2;

    // project pores to planes
    vector<VoxelVolume<float> > planeProjections(3);

    for(int dim=0; dim<3; ++dim)
    {
      planeProjections[dim].s = imageXYZ;
      planeProjections[dim].s(dim) = 1;
      planeProjections[dim].voxelValues.resize(planeProjections[dim].s.prod(),1.0);
    }

    long n=0;
    for( auto const & IDAndPore : pores )
    {

      auto currentPoreProjections = planeProjections;
      for(int dim=0; dim<3; ++dim)
      {
        currentPoreProjections[dim].voxelValues.clear();
        currentPoreProjections[dim].voxelValues.resize(currentPoreProjections[dim].s.prod(),0.0);
      }

      Vector3i centerOfMass(round(centersOfMass[n](0)),
                            round(centersOfMass[n](1)),
                            round(centersOfMass[n](2))  );

      ++n;

      for(auto const & voxelIndex : IDAndPore.second)
      {
        Vector3i coordinate = distanceField.to_vx(voxelIndex);
        Vector3i imageCoordinate = coordinate-centerOfMass+imageCenter;

        for(int dim=0; dim<3; ++dim)
        {
          Vector3i planeCoordinate = imageCoordinate;
          planeCoordinate(dim) = 1;
          auto planeIndex = currentPoreProjections[dim].to_vxID(planeCoordinate);
          currentPoreProjections[dim][planeIndex] = 1.0;
        }
      }

      for(int dim=0; dim<3; ++dim)
        for(long n=0; n<currentPoreProjections[dim].s.prod(); ++n)
        {
          planeProjections[dim][n] += currentPoreProjections[dim][n];
        }
    }

    for(int dim=0; dim<3; ++dim)
      for(long n=0; n<planeProjections[dim].s.prod(); ++n )
        planeProjections[dim][n] = log(planeProjections[dim][n]);

    for(int dim=0; dim<3; ++dim)
    {
      ofstream myfile(directory + "projection" + to_string(dim) + ".mat" );
      long n=0;
      for(long dim0=0; dim0<maxCoeffXyz; ++dim0)
      {
        for(long dim1=0; dim1<maxCoeffXyz; ++dim1)
          myfile << planeProjections[dim][n++] << "\t";
        myfile << endl;
      }
    }

  }

//*************************************************************************************************

  void seperate_pores_by_size(float r_eff,
                              std::map<PoreID,Pore> & smallPores,
                              std::map<PoreID,Pore> & largePores,
                              std::map<PoreID,Pore> const & pores,
                              VoxelVolume<float> const & distanceField)
  {

    cout << "\nSeperating small from large pores\n";

    smallPores.clear();
    largePores.clear();

    for(auto const & IDAndPore : pores)
    {
      auto const & ID = IDAndPore.first;
      auto const & pore = IDAndPore.second;

      float volume = volume_of_(pore);
      float r_eff_pore = effective_radius(volume);

      if(r_eff_pore < r_eff)
        smallPores[ID] = pore;
      else
        largePores[ID] = pore;
    }

  }

//*************************************************************************************************

  Matrix3d cout_average_principal_component_analysis(std::map<PoreID,Pore> const & pores,
                                                 VoxelVolume<float> const & distanceField)
  {

    cout << "\nPrincipal componenent anaylsis ...\n";

    // add covariance matrices
    Matrix3d SigmaAverage = Matrix3d::Zero();
    for( auto const & IDAndPore : pores )
    {
      auto const & pore = IDAndPore.second;
      Matrix3f SigmaPore = covariance_matrix(pore,
                                              center_of_mass(pore,distanceField),
                                              distanceField);
      SigmaAverage += SigmaPore.cast<double>();
    }

    SigmaAverage /= double(pores.size());

    // solve for eigenvalues
    SelfAdjointEigenSolver<Matrix3d> eigensolver(SigmaAverage);
    if( eigensolver.info() != Success )
      cout << "\nNo Eigenvalues found. Shouldn't happen.\n";
    Vector3d eigenValues = eigensolver.eigenvalues();
    Matrix3d eigenVectors = eigensolver.eigenvectors();

    Vector3d sigmas = eigenValues.array().sqrt();

    cout << "\nSigma on principle axis:\n" << sigmas.transpose() << endl;
    cout << "\nPrinciple axis directions:\n" << eigenVectors << endl;

    cout << "\nRatio of sigmas: "
         << sigmas(0)/sigmas(1) << " : "
         << 1 << " : "
         << sigmas(2)/sigmas(1) << endl;

    return SigmaAverage;

  }

//*************************************************************************************************

  void compare_anisotropy(std::map<PoreID,Pore> const & pores,
                          VoxelVolume<float> const & distanceField,
                          std::string filename)
  {


    // see function cout_average_principal_component_analysis
    // do same for each bin




    auto poreRadii = calculate_effective_radii(calculate_pore_volumes(pores));
    auto poreMinMax = minmax_element(poreRadii.begin(), poreRadii.end(),
                                     [](pair<PoreID,float> const & eleA,
                                        pair<PoreID,float> const & eleB){
                                        return eleA.second < eleB.second; } );
    float minRadius = (*poreMinMax.first).second;
    float maxRadius = (*poreMinMax.second).second;


    size_t bins = sqrt(pores.size());
    LogarithmicBins<double> logBins(minRadius,maxRadius,bins);

    vector<Matrix3d> cumulativeSigmas(bins,Matrix3d::Zero());
    vector<size_t> binCounter(bins,0);

    auto radiusIterator = poreRadii.begin();
    for(auto const & ele : pores)
    {
      size_t bin = logBins.to_bin((*radiusIterator).second);
      cumulativeSigmas[bin] += covariance_matrix(ele.second,
                                                 center_of_mass(ele.second,distanceField),
                                                 distanceField).cast<double>();

      ++radiusIterator;
    }

    ofstream myfile(filename);
    myfile << "#bin_start\tbin_end\tsig_3/sig_2\tv3_x\tv3_y\tv3_z\tsig_1/sig_2\tv2_x\tv2_y\tv2_z\tv1_x\tv1_y\tv1_z\n";
    for(size_t n=bins; --n>0; )
    {
      auto Sigma = cumulativeSigmas[n];

      if((Sigma.array() == 0).all())
        continue;

      SelfAdjointEigenSolver<Matrix3d> eigensolver(Sigma);
      if( eigensolver.info() != Success )
        cout << "\nNo Eigenvalues found. Shouldn't happen.\n";
      Vector3d eigenValues = eigensolver.eigenvalues();
      Matrix3d eigenVectors = eigensolver.eigenvectors();

      double sigmaRatio = sqrt(eigenValues(2)/eigenValues(1));

      myfile << logBins.to_value(n) << "\t"
             << logBins.to_value(n+1) << "\t"
             << sigmaRatio << "\t"
             << eigenVectors.col(2).transpose() << "\t"
             << sqrt(eigenValues(0)/eigenValues(1)) << "\t"
             << eigenVectors.col(1).transpose() << "\t"
             << eigenVectors.col(0).transpose() << endl;
    }


  }

//*************************************************************************************************

  void mle_normal(std::map<PoreID,float> const & poreRadii, float samplingFactor)
  {

    cout << "\nEstimating parameters of normal distribution ...\n";

    // https://en.wikipedia.org/wiki/Log-normal_distribution#Maximum_likelihood_estimation_of_parameters
    double muHat = 0.0;
    for(auto const & poreAndRadius : poreRadii)
      muHat += poreAndRadius.second*samplingFactor;
    muHat /= poreRadii.size();

    double sigmaHatSquared = 0.0;
    for(auto const & poreAndRadius : poreRadii)
      sigmaHatSquared += pow(poreAndRadius.second*samplingFactor-muHat,2);
    sigmaHatSquared /= poreRadii.size()-1;

    cout << "\nEstimated mu:\n" << muHat << endl;
    cout << "\nEstimated sigma:\n" << sqrt(sigmaHatSquared) << endl;
  }

//*************************************************************************************************

  vector<VoxelIndex> find_largest_inscribed_sphere_centers(auto const & pores,
                                                           VoxelVolume<float> const & distanceField)
  {

    vector<VoxelIndex> indicesAtSphereCenters;
    indicesAtSphereCenters.reserve(pores.size());

    for( auto const & poreID_and_voxelIndices : pores )
    {
      vector<VoxelIndex> const & voxelIndices = poreID_and_voxelIndices.second;

      VoxelIndex d_max_index = *(max_element(
                                 voxelIndices.begin(),
                                 voxelIndices.end(),
                                 [&](VoxelIndex const & i, VoxelIndex const & j){
                                   return distanceField[i] == distanceField[j]
                                          ? i<j
                                          : distanceField[i] < distanceField[j];} ));

      indicesAtSphereCenters.push_back(d_max_index);
    }

    return indicesAtSphereCenters;
  }

//*************************************************************************************************

  void nodes_vs_ligaments_distribution(std::map<PoreID,Pore> const & pores,
                                       std::map<PoreID, std::vector<PoreID> > const & poreNetwork,
                                       VoxelVolume<float> const & distanceField)
  {
    auto poreCenters = find_largest_inscribed_sphere_centers(pores, distanceField);

    int64_t ligamentCounter=0;
    int64_t nodeCounter=0;

    auto pore_iterator = pores.begin();
    auto network_iterator = poreNetwork.begin();
    auto center_iterator = poreCenters.begin();
    for( ; pore_iterator != pores.end() ; ++pore_iterator,
                                          ++network_iterator,
                                          ++center_iterator )
    {
      int connectedPores = (network_iterator->second).size();
      if(connectedPores < 3)
      {
        ligamentCounter += (pore_iterator->second).size();
        continue;
      }

//      nodeCounter += (pore_iterator->second).size();
//      continue;

      float r = distanceField[*center_iterator];
      float r_squared = r*r;
      VoxelPosition x_i = distanceField.to_vx(*center_iterator);

      for(VoxelIndex const & j : pore_iterator->second )
      {
        VoxelPosition x_j = distanceField.to_vx(j);
        if((x_i-x_j).squaredNorm() <= r_squared)
          ++nodeCounter;
        else
          ++ligamentCounter;
      }
    }

    cout << scientific << "\nligament fraction: "
         << double(ligamentCounter)/double(ligamentCounter+nodeCounter)
         << "\nnode fraction: "
         << double(nodeCounter)/double(ligamentCounter+nodeCounter) << endl;
  }

//*************************************************************************************************

  VoxelVolume<uint8_t> create_voronoi_seeds( std::map<PoreID,Pore> const & pores,
                                             VoxelVolume<float> const & distanceField)
  {

    VoxelPosition const & s = distanceField.s;

    VoxelVolume<uint8_t> seedVolume; seedVolume.s = s;
    seedVolume.voxelValues.resize(s.prod(),0);

    auto seedCenters = find_largest_inscribed_sphere_centers(pores, distanceField);

    for( VoxelIndex const & seedCenter : seedCenters )
      seedVolume[seedCenter] = 255;

    return seedVolume;

  }

//*************************************************************************************************

  void material_voxels_in_direction(int direction, VoxelVolume<float> const & distanceField,
                                    std::string filename)
  {
    Vector3i const & s = distanceField.s;

    vector<size_t> voxelsInPlane(s(direction),0);

#pragma omp parallel for
    for(int dir0 = 0; dir0<s(direction); ++dir0)
    {
      int direction1 = (direction+1)%3;
      int direction2 = (direction+2)%3;
      for(int dir1 = 0; dir1<s(direction1); ++dir1)
        for(int dir2 = 0; dir2<s(direction2); ++dir2)
          if( distanceField(dir0,dir1,dir2) <= 0 )
            ++voxelsInPlane[dir0];
    }

    ofstream myFile(filename);
    for(size_t counter : voxelsInPlane)
      myFile << counter << endl;

  }

}
