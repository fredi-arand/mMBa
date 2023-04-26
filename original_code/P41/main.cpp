#include <iostream>
#include "voxelio.h"
#include <string>
#include "definitions.h"
#include "morphologystatistics.h"
#include <map>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include "histogram.h"

using namespace std;
using namespace voxelIO;
using namespace morphologyStatistics;
using namespace Eigen;

void experiment1(); // pore radius distribution of R
void experiment2(); // pore radius distribution of R_large
void experiment3(); // void volume fractions
void experiment4(); // projections
void experiment5(); // principal components
void experiment6(); // reciprocal space
void experiment7(); // data for anisotropy fits
void experiment8(); // data for size fits
void experiment9(); // artificial foam: large pores
void experiment10(); // artificial foam: all pores
void experiment11(); // count material
void experiment12(); // create voronoi seeds
void experiment13(); // measure material at nodes and ligaments for sample
void experiment14(); // measure material at nodes and ligaments for model
void experiment15(); // close small pores of sample
void experiment16(); // statistics of modified model
void experiment17(); // thesis: mMBa validation data set
void experiment18(); // thesis: hcp verification
void experiment18a(); // thsis: addendum hcp verification
void experiment19(); // thesis: last verification
void experiment20(); // thesis: last verification 2
void experiment21(); // thesis: carbon foam: 8 experiments
void experiment22(); // thesis: carbon foam: large region
void experiment23(); // thesis: anisotropy
void experiment23a(); // thesis: addendum anisotropy of model
void experiment24(); // thesis: relative ratios porosity
void experiment25(); // thesis: model: 8 experiments
void experiment26(); // thesis: model: large region

//*************************************************************************************************

int main(int argc, char *argv[])
{

  cout << "\nPlease edit main.cpp to change program behaviour.\n";

  experiment23a();

  return 0;

}

//*************************************************************************************************

void experiment26()
{
  Vector3i s(1000,1000,1000);
  int bins=4*50;

  VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                             "input/8_subvolumes_model/large_morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                             "input/8_subvolumes_model/large_stateVolume.raw",s);

  bool includeBoundaryPores=false;
  auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
  auto radii = calculate_effective_radii(calculate_pore_volumes(pores));

  Histogram<double> histogram( log10(0.1), log10(1000.0), bins );
  for( auto const & id_radius : radii )
    ++histogram( log10(2.0*id_radius.second) );

  ofstream myfile("model_subvolumes_results/large_region.txt");
  for(int bin=0; bin<histogram.bins; ++bin)
    myfile << histogram.to_value(bin) << " "
           << histogram.to_value(bin+1) << " "
           << histogram.counts[bin] << endl;
}

//*************************************************************************************************

void experiment25()
{
  int bins = 4*50;

  Histogram<double> cumulative_histogram( log10(0.1), log10(1000.0), bins );
  for(int i=0; i<8; ++i)
  {
    Vector3i s(1000,1000,1000);
//    const VoxelVolume<float> distanceField = to_voxelVolume<float>(
//                                               "input/8_subvolumes/"+to_string(i)+"distanceField.raw",s);
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                               "input/8_subvolumes_model/"+to_string(i)+"morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                               "input/8_subvolumes_model/"+to_string(i)+"stateVolume.raw",s);

    bool includeBoundaryPores=false;
    auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
    auto radii = calculate_effective_radii(calculate_pore_volumes(pores));

    Histogram<double> histogram( log10(0.1), log10(1000.0), bins );
    for( auto const & id_radius : radii )
    {
      ++histogram( log10(id_radius.second) );
      ++cumulative_histogram( log10(id_radius.second) );
    }

    ofstream myfile("model_subvolumes_results/"+to_string(i)+".txt");
    for(int bin=0; bin<histogram.bins; ++bin)
      myfile << histogram.to_value(bin) << " "
             << histogram.to_value(bin+1) << " "
             << histogram.counts[bin] << endl;
  }

  ofstream myfile("model_subvolumes_results/cumulative.txt");
  for(int bin=0; bin<cumulative_histogram.bins; ++bin)
    myfile << cumulative_histogram.to_value(bin) << " "
           << cumulative_histogram.to_value(bin+1) << " "
           << cumulative_histogram.counts[bin] << endl;


}

//*************************************************************************************************

void experiment24()
{
  Vector3i const s(1000,1000,1000);
  int64_t smallPoreVoxels = 0;
  int64_t largePoreVoxels = 0;
  int64_t voidVoxels = 0;

  VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
        "input/8_subvolumes/0morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume = to_voxelVolume<uint8_t>(
                                       "input/8_subvolumes/0stateVolume.raw",s);

  bool includeBoundaryPores=true;
  auto pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  map<PoreID,Pore> smallPores, largePores;
  seperate_pores_by_size(8.0, smallPores, largePores, pores, VoxelVolume<float>());

  for(auto const & ID_and_voxels : smallPores)
      smallPoreVoxels += ID_and_voxels.second.size();

  for(auto const & ID_and_voxels : largePores)
      largePoreVoxels += ID_and_voxels.second.size();

  for (auto const & value : stateVolume.voxelValues)
    if(value != BACKGROUND)
      ++ voidVoxels;

  cout << endl << "small pore voxels: " << smallPoreVoxels;
  cout << endl << "large pore voxels: " << largePoreVoxels;
  cout << endl << "void voxels: " << voidVoxels << endl;


}

//*************************************************************************************************

void experiment23()
{
  Vector3i const s(1000,1000,1000);
  int bins=4*50;
  double hist_min = log10(0.1);
  double hist_max = log10(1000.0);
  bool includeBoundaryPores=false;
  ofstream myFile("../../source/thesis/chapter_2/fig35.txt");

  Histogram<double> h_combined(hist_min, hist_max, bins);
  vector<Matrix3f> average_covariance_matrices( bins, Matrix3f::Zero() );

  // handle R_large volume first
  {
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                             "input/8_subvolumes/large_morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                               "input/8_subvolumes/large_stateVolume.raw",s);

    map<PoreID,Pore> largePores;
    {
      auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
      map<PoreID,Pore> smallPores;
      seperate_pores_by_size( 5.0, smallPores, largePores, pores, VoxelVolume<float>() );
    }

    VoxelVolume<float> distanceFieldDummy;
    distanceFieldDummy.s = s;
    for( auto const & id_and_indices : largePores )
    {
      auto const & indices = id_and_indices.second;

      double r_eff = 2.0*effective_radius( indices.size() );
      ++h_combined(log10(r_eff));

      average_covariance_matrices[ h_combined.to_bin(log10(r_eff)) ] +=
          4.0*covariance_matrix(indices, center_of_mass(indices,distanceFieldDummy),
                                distanceFieldDummy);
    }

  }

  // Proceed with R 1..8,
  for(int i=0; i<8; ++i)
  {
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                               "input/8_subvolumes/"+to_string(i)+"morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                               "input/8_subvolumes/"+to_string(i)+"stateVolume.raw",s);

    map<PoreID,Pore> smallPores;
    {
      auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
      map<PoreID,Pore> largePores;
      seperate_pores_by_size( 10.0, smallPores, largePores, pores, VoxelVolume<float>() );
    }

    VoxelVolume<float> distanceFieldDummy;
    distanceFieldDummy.s = s;
    for( auto const & id_and_indices : smallPores )
    {
      auto const & indices = id_and_indices.second;

      double r_eff = effective_radius( indices.size() );
      ++h_combined(log10(r_eff));

      average_covariance_matrices[ h_combined.to_bin(log10(r_eff)) ] +=
              covariance_matrix(indices, center_of_mass(indices,distanceFieldDummy),
                                distanceFieldDummy);
    }

  }

  myFile << "# log10(r_i)" << "\t"
         << "log10(r_{i+1})" << "\t"
         << "Counts" << "\t"
         << "lambda_1" << "\t"
         << "lambda_2" << "\t"
         << "lambda_3" << "\t"
         << "v_11" << "\t"
         << "v_21" << "\t"
         << "v_31" << "\t"
         << "v_12" << "\t"
         << "v_22" << "\t"
         << "v_32" << "\t"
         << "v_13" << "\t"
         << "v_23" << "\t"
         << "v_33" << "\n";

  for(long n=0; n<average_covariance_matrices.size(); ++n)
  {
    myFile << h_combined.to_value(n) << "\t"
           << h_combined.to_value(n+1) << "\t"
           << h_combined.counts[n] << "\t";

    if( h_combined.counts[n] != 0 )
      average_covariance_matrices[n] /= h_combined.counts[n];

    auto const & Sigma = average_covariance_matrices[n];

    SelfAdjointEigenSolver<Matrix3f> eigensolver(Sigma);
    if( eigensolver.info() != Success )
      cout << "\nNo Eigenvalues found. Shouldn't happen.\n";
    Vector3f eigenValues = eigensolver.eigenvalues();
    Matrix3f eigenVectors = eigensolver.eigenvectors();

    myFile << (eigenValues.array().abs()).transpose() << "\t";
    myFile << eigenVectors.col(0).transpose() << "\t"
           << eigenVectors.col(1).transpose() << "\t"
           << eigenVectors.col(2).transpose() << "\n";
  }


}

//*************************************************************************************************

void experiment23a()
{
  Vector3i const s(1000,1000,1000);
  size_t bins=4*50;
  double hist_min = log10(0.1);
  double hist_max = log10(1000.0);
  bool includeBoundaryPores=false;
  ofstream myFile("/Users/arand/uni/source/thesis/chapter_3/fig19.txt");

  Histogram<double> h_combined(hist_min, hist_max, bins);
  vector<Matrix3f> average_covariance_matrices( bins, Matrix3f::Zero() );

  // handle R_large volume first
  {
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                             "/Users/arand/from_vgpc32/build/P41/input/8_subvolumes_model/large_morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                             "/Users/arand/from_vgpc32/build/P41/input/8_subvolumes_model/large_stateVolume.raw",s);

    map<PoreID,Pore> largePores;
    {
      auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
      cout << "\nTotal pores: " << pores.size() << endl;
      map<PoreID,Pore> smallPores;
      seperate_pores_by_size( 5.f/.69f, smallPores, largePores, pores, VoxelVolume<float>() );
    }

    VoxelVolume<float> distanceFieldDummy;
    distanceFieldDummy.s = s;
    for( auto const & id_and_indices : largePores )
    {
      auto const & indices = id_and_indices.second;

      double r_eff = 2.0*effective_radius( indices.size() );
      ++h_combined(log10(r_eff));

      average_covariance_matrices[ h_combined.to_bin(log10(r_eff)) ] +=
          4.0*covariance_matrix(indices, center_of_mass(indices,distanceFieldDummy),
                                distanceFieldDummy);
    }

  }

  // Proceed with R 1..8,
  for(int i=0; i<8; ++i)
  {
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                               "/Users/arand/from_vgpc32/build/P41/input/8_subvolumes_model/"+to_string(i)+"morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                               "/Users/arand/from_vgpc32/build/P41/input/8_subvolumes_model/"+to_string(i)+"stateVolume.raw",s);

    map<PoreID,Pore> smallPores;
    {
      auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
      cout << "\nTotal pores: " << pores.size() << endl;
      map<PoreID,Pore> largePores;
      seperate_pores_by_size( 10.f/.69f, smallPores, largePores, pores, VoxelVolume<float>() );
    }

    VoxelVolume<float> distanceFieldDummy;
    distanceFieldDummy.s = s;
    for( auto const & id_and_indices : smallPores )
    {
      auto const & indices = id_and_indices.second;

      double r_eff = effective_radius( indices.size() );
      ++h_combined(log10(r_eff));

      average_covariance_matrices[ h_combined.to_bin(log10(r_eff)) ] +=
              covariance_matrix(indices, center_of_mass(indices,distanceFieldDummy),
                                distanceFieldDummy);
    }

  }

  myFile << "# log10(r_i)" << "\t"
         << "log10(r_{i+1})" << "\t"
         << "Counts" << "\t"
         << "lambda_1" << "\t"
         << "lambda_2" << "\t"
         << "lambda_3" << "\t"
         << "v_11" << "\t"
         << "v_21" << "\t"
         << "v_31" << "\t"
         << "v_12" << "\t"
         << "v_22" << "\t"
         << "v_32" << "\t"
         << "v_13" << "\t"
         << "v_23" << "\t"
         << "v_33" << "\n";

  for(long n=0; n<average_covariance_matrices.size(); ++n)
  {
    myFile << h_combined.to_value(n) << "\t"
           << h_combined.to_value(n+1) << "\t"
           << h_combined.counts[n] << "\t";

    if( h_combined.counts[n] != 0 )
      average_covariance_matrices[n] /= h_combined.counts[n];

    auto const & Sigma = average_covariance_matrices[n];

    SelfAdjointEigenSolver<Matrix3f> eigensolver(Sigma);
    if( eigensolver.info() != Success )
      cout << "\nNo Eigenvalues found. Shouldn't happen.\n";
    Vector3f eigenValues = eigensolver.eigenvalues();
    Matrix3f eigenVectors = eigensolver.eigenvectors();

    myFile << (eigenValues.array().abs()).transpose() << "\t";
    myFile << eigenVectors.col(0).transpose() << "\t"
           << eigenVectors.col(1).transpose() << "\t"
           << eigenVectors.col(2).transpose() << "\n";
  }


}

//*************************************************************************************************

void experiment22()
{
  Vector3i s(1000,1000,1000);
  int bins=4*50;

  VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                             "input/8_subvolumes/large_morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                             "input/8_subvolumes/large_stateVolume.raw",s);

  bool includeBoundaryPores=false;
  auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
  auto radii = calculate_effective_radii(calculate_pore_volumes(pores));

  Histogram<double> histogram( log10(0.1), log10(1000.0), bins );
  for( auto const & id_radius : radii )
    ++histogram( log10(2.0*id_radius.second) );

  ofstream myfile("carbonfoam_subvolumes_results/large_region.txt");
  for(int bin=0; bin<histogram.bins; ++bin)
    myfile << histogram.to_value(bin) << " "
           << histogram.to_value(bin+1) << " "
           << histogram.counts[bin] << endl;
}

//*************************************************************************************************

void experiment21()
{
  int bins = 4*50;

  Histogram<double> cumulative_histogram( log10(0.1), log10(1000.0), bins );
  for(int i=0; i<8; ++i)
  {
    Vector3i s(1000,1000,1000);
//    const VoxelVolume<float> distanceField = to_voxelVolume<float>(
//                                               "input/8_subvolumes/"+to_string(i)+"distanceField.raw",s);
    VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>(
                                               "input/8_subvolumes/"+to_string(i)+"morphologyVolume.raw",s);
    VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>(
                                               "input/8_subvolumes/"+to_string(i)+"stateVolume.raw",s);

    bool includeBoundaryPores=false;
    auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);
    auto radii = calculate_effective_radii(calculate_pore_volumes(pores));

    Histogram<double> histogram( log10(0.1), log10(1000.0), bins );
    for( auto const & id_radius : radii )
    {
      ++histogram( log10(id_radius.second) );
      ++cumulative_histogram( log10(id_radius.second) );
    }

    ofstream myfile("carbonfoam_subvolumes_results/"+to_string(i)+".txt");
    for(int bin=0; bin<histogram.bins; ++bin)
      myfile << histogram.to_value(bin) << " "
             << histogram.to_value(bin+1) << " "
             << histogram.counts[bin] << endl;
  }

  ofstream myfile("carbonfoam_subvolumes_results/cumulative.txt");
  for(int bin=0; bin<cumulative_histogram.bins; ++bin)
    myfile << cumulative_histogram.to_value(bin) << " "
           << cumulative_histogram.to_value(bin+1) << " "
           << cumulative_histogram.counts[bin] << endl;


}

//*************************************************************************************************

void experiment20()
{

  const Vector3i s(128,128,128);
  const VoxelVolume<float> distanceField = to_voxelVolume<float>(   "input/verification2/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>("input/verification2/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>( "input/verification2/stateVolume.raw",s);

  bool includeBoundaryPores = true;
  auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);

  set<PoreID> boundaryPores;
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
        boundaryPores.insert(poreID);
      }
  }

  map<float,vector<FloatPosition> > radii_centers;
  ifstream spheres_file("../P57/centers_radii.txt");
  {
    float x0,x1,x2,r;
    while( spheres_file >> x0 >> x1 >> x2 >> r )
    {
      FloatPosition x(x0/8.0,x1/8.0,x2/8.0);
      r /= 8.0;
      radii_centers[r].push_back(x);
    }
  }

  map<float,int> is_material;
  map<float,int> is_boundary;
  map<PoreID,vector<float> > pore_to_r;
  map<float,int> is_merged;
  map<float,int> is_detected;

  for( auto const & r_centers : radii_centers )
  {
    float r = r_centers.first;
    for( auto const & x : r_centers.second )
    {

      VoxelPosition x_int = x.array().floor().cast<int>();

      if(distanceField(x_int)==0.0)
      {
        ++is_material[r];
        continue;
      }

      if( stateVolume(x_int)==THROAT )
      {
        ++is_merged[r];
        continue;
      }

      PoreID m = morphologyVolume(x_int);
//      if(boundaryPores.count(m)==1)
//      {
//        ++is_boundary[r];
//        continue;
//      }

      if(m==0)
        cout << "\nblub\n";

      pore_to_r[m].push_back(r);
    }
  }

  for( auto const & m_rs : pore_to_r )
  {
    vector<float> const & rs = m_rs.second;
    PoreID m = m_rs.first;
    if( boundaryPores.count(m)==1)
      ++is_boundary[rs.back()];
    else
      ++is_detected[rs.back()];

    for(auto it = --(m_rs.second.end()); it-- != m_rs.second.begin();  )
      ++is_merged[*it];
  }

  cout << "\nMaterial:\n";
  for( auto const & r_count : is_material )
    cout << r_count.first << " " << r_count.second << endl;

  cout << "\nBoundary:\n";
  for( auto const & r_count : is_boundary )
    cout << r_count.first << " " << r_count.second << endl;

  cout << "\nMerged:\n";
  for( auto const & r_count : is_merged )
    cout << r_count.first << " " << r_count.second << endl;

  cout << "\nDetected:\n";
  for( auto const & r_count : is_detected )
    cout << r_count.first << " " << r_count.second << endl;

  map<int,map<float,int> > gp_table;
  for( auto const & r_count : is_material )
    gp_table[0][r_count.first] = r_count.second;
  for( auto const & r_count : is_boundary )
    gp_table[1][r_count.first] = r_count.second;
  for( auto const & r_count : is_merged )
    gp_table[2][r_count.first] = r_count.second;
  for( auto const & r_count : is_detected )
    gp_table[3][r_count.first] = r_count.second;

  for(int line=0; line<4; ++line)
  {
    cout << "\nsome text";
    for(float r=0.5; r!=32.0; r*=2.0)
      cout << " & " << gp_table[line][r];
    cout << " \\\\ \\hline";
  }
  cout << endl;

  cout << "some text";
  for(float r=0.5; r!=32.0; r*=2.0)
  {
    int sum=0;
    for(int line=0; line<4; ++line)
      sum += gp_table[line][r];
    cout << " & " << sum;
  }
  cout << endl;

}

//*************************************************************************************************

void experiment19()
{
  const Vector3i s(128,128,128);
  const VoxelVolume<float> distanceField = to_voxelVolume<float>(   "input/verification2/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume = to_voxelVolume<PoreID>("input/verification2/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =       to_voxelVolume<uint8_t>( "input/verification2/stateVolume.raw",s);

  bool includeBoundaryPores = false;
  auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);

  auto radii = calculate_effective_radii(
    calculate_pore_volumes( pores )
  );

  float r_min = 10000, r_max=0.0;
  for(auto const & id_radius : radii)
  {
    float r = id_radius.second;
    if(r>r_max)
      r_max = r;
    if(r<r_min)
      r_min = r;
  }

  Histogram<float> histogram(log2(.5/sqrt(2.0)),log2(16.0*sqrt(2.0)),
                             6);
  for( auto const & id_radius : radii )
    ++histogram( log2(id_radius.second) );

  ofstream myfile("histogram_validation.txt");
  for(int n=0; n<histogram.bins; ++n)
    myfile << histogram.to_value(n) << " "
           << histogram.to_value(n+1) << " "
           << histogram.counts[n] << endl;

  Histogram<float> histogram2(log2(.5/sqrt(2.0)),log2(16.0*sqrt(2.0)),
                             5*6);
  for( auto const & id_radius : radii )
    ++histogram2( log2(id_radius.second) );

  ofstream myfile2("histogram_validation2.txt");
  for(int n=0; n<histogram2.bins; ++n)
    myfile2 << histogram2.to_value(n) << " "
            << histogram2.to_value(n+1) << " "
            << histogram2.counts[n] << endl;
}

//*************************************************************************************************

void experiment18a()
{
  const Vector3i s(512,512,512);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/hcp/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/hcp/stateVolume.raw",s);

  bool includeBoundaryPores = false;
  map<PoreID,Pore> pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);


  auto poreRadii = calculate_effective_radii(
                     calculate_pore_volumes(pores)
                   );

  Histogram<float> histogram(13.44-1.2, 13.44+1.2, 48);
  for(auto const & indexAndRadius : poreRadii)
    ++histogram(indexAndRadius.second);

  ofstream myfile("hcp_pore_radii_r_eff.txt");
  for(int n=0; n<histogram.bins; ++n)
  {
      myfile << histogram.to_value(n) << " "
             << histogram.to_value(n+1) << " "
             << histogram.counts[n] << endl;
  }

}

//*************************************************************************************************

void experiment18()
{
  const Vector3i s(512,512,512);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/hcp/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/hcp/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/hcp/stateVolume.raw",s);

  bool includeBoundaryPores = false;
  map<PoreID,Pore> pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);

//  {
//    auto pore_radii = calculate_effective_radii(
//                        calculate_pore_volumes(pores)
//                        );
//    for(auto const ele : pore_radii)
//      cout << ele.second << endl;

//    return;
//  }

  vector<VoxelIndex> sphereCenters = find_largest_inscribed_sphere_centers(pores,distanceField);
  vector<float> poreRadii;
  for( auto const & index : sphereCenters )
    poreRadii.push_back(distanceField[index]);

  sort( poreRadii.begin(), poreRadii.end() );
  Histogram<float> histogram(13.44-1.2, 13.44+1.2, 48);
  for(auto const & r : poreRadii)
    ++histogram(r);

  ofstream myfile("hcp_pore_radii.txt");
  for(int n=0; n<histogram.bins; ++n)
    myfile << histogram.to_value(n) << " "
           << histogram.to_value(n+1) << " "
           << histogram.counts[n] << endl;

  ofstream myfile3("hcp_pore_centers.txt");
  for(auto const & index : sphereCenters )
    myfile3 << distanceField.to_vx(index).transpose() << endl;


  map<ThroatID,Throat> throats = extract_throats(morphologyVolume,stateVolume);

  // remove boundary throats!
  vector<ThroatID> remove_these_throats;
  for(auto const & throatID_throatVoxels : throats)
    for(auto const & index : throatID_throatVoxels.second)
    {
      VoxelPosition position = distanceField.to_vx(index);
      if((position.array() == 0).any()
         || (position.array() == s.array()-1).any() )
      {
        remove_these_throats.push_back(throatID_throatVoxels.first);
        break;
      }
    }

  for(auto const & throatID : remove_these_throats)
    throats.erase(throatID);

//  vector<VoxelIndex> throatCenters = find_largest_inscribed_sphere_centers(throats,distanceField);
  vector<VoxelIndex> throatCenters;
  throatCenters.reserve(throats.size());

  for( auto const & poreID_and_voxelIndices : throats )
  {
    vector<VoxelIndex> const & voxelIndices = poreID_and_voxelIndices.second;

    VoxelIndex d_max_index = *(max_element(
                               voxelIndices.begin(),
                               voxelIndices.end(),
                               [&](VoxelIndex const & i, VoxelIndex const & j){
                                 return distanceField[i] == distanceField[j]
                                        ? i<j
                                        : distanceField[i] < distanceField[j];} ));

    throatCenters.push_back(d_max_index);
  }
  vector<float> throatRadii;
  for( auto const & index : throatCenters )
    throatRadii.push_back(distanceField[index]);

  sort( throatRadii.begin(), throatRadii.end() );

  Histogram<float> histogram_throat(4.098-1.2, 4.098+1.2, 48);
  for(auto const & r : throatRadii)
    ++histogram_throat(r);

  ofstream myfile2("hcp_throat_radii.txt");
  for(int n=0; n<histogram_throat.bins; ++n)
    myfile2 << histogram_throat.to_value(n) << " "
           << histogram_throat.to_value(n+1) << " "
           << histogram_throat.counts[n] << endl;

  ofstream myfile4("hcp_throat_centers.txt");
  for(auto const & index : throatCenters )
    myfile4 << distanceField.to_vx(index).transpose() << endl;




}

//*************************************************************************************************

void experiment17()
{
  const Vector3i s(500,500,500);

  for(float epsilon=0.0; epsilon < 0.21; epsilon += 0.05)
  {
    const VoxelVolume<float> distanceField =
        to_voxelVolume<float>("input/thesis_verification/distanceField_"
                              + to_string(int(20*epsilon)) + ".raw",s);
    VoxelVolume<PoreID> morphologyVolume =
        to_voxelVolume<PoreID>("input/thesis_verification/morphologyVolume_"
                                 + to_string(int(20*epsilon))+".raw",s);
    VoxelVolume<uint8_t> stateVolume =
        to_voxelVolume<uint8_t>("input/thesis_verification/stateVolume_"
                                + to_string(int(20*epsilon))+".raw",s);

    bool includeBoundaryPores = true;
    auto pores = extract_pores(morphologyVolume,stateVolume, includeBoundaryPores);

    vector<VoxelIndex> sphereCenters = find_largest_inscribed_sphere_centers(pores, distanceField);

    ofstream myFile( "verification_thesis/centers_" + to_string(int(20*epsilon)) + ".txt" );
    for(VoxelIndex const & i : sphereCenters)
      myFile << stateVolume.to_vx(i).transpose() << endl;
  }
}

//*************************************************************************************************

void experiment16()
{
  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/digitalFoam_modified/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/digitalFoam_modified/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/digitalFoam_modified/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  auto pores = extract_pores(morphologyVolume,stateVolume,includeBoundaryPores);

  logarithmic_histogram_r_eff("histogram_r_eff_eroded_dilated.txt",
                              calculate_effective_radii(
                                calculate_pore_volumes(pores)),
                              sqrt(pores.size()));
}

//*************************************************************************************************

void experiment15()
{
  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  VoxelVolume<uint8_t> voxelVolume =
      to_voxelVolume<uint8_t>("../P32/carbonFoam_denoised_connected_750_750_750.raw",s);

  bool includeBoundaryPores = true;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  auto throats = extract_throats(morphologyVolume,stateVolume);
  auto poreRadii = calculate_effective_radii(calculate_pore_volumes(pores));
  float smallPoreRadius = 7.0;

  for(auto const & poreID_and_voxelIndices : pores )
  {
    PoreID const & poreID = poreID_and_voxelIndices.first;
    if(poreRadii.at(poreID)>=smallPoreRadius)
      continue;

    auto const & voxelIndices = poreID_and_voxelIndices.second;
    for( auto const & voxelIndex : voxelIndices )
      voxelVolume[voxelIndex] = 255;
  }

  for(auto const & throatID_and_voxelIndices : throats)
  {
    ThroatID const & throatID = throatID_and_voxelIndices.first;
    bool connectedToSmallPore = false;
    for(auto const & poreID : throatID)
    {
      if(poreRadii.at(poreID)>=smallPoreRadius)
        continue;

      connectedToSmallPore = true;
    }

    if(!connectedToSmallPore)
      continue;

    auto const & voxelIndices = throatID_and_voxelIndices.second;
    for( auto const & voxelIndex : voxelIndices )
      voxelVolume[voxelIndex] = 255;
  }

  for(VoxelIndex n=0; n<s.prod(); ++n)
    if(distanceField[n] > 0.0
       && (stateVolume[n] == BACKGROUND
           || stateVolume[n] == INIT_VALUE ) )
      voxelVolume[n] = 255;

  to_file(voxelVolume, "../../volumedata/ITV_carbon_foam/sample_without_small_pores/sample_without_small_pores_750_750_750_uint8.raw");

}

//*************************************************************************************************

void experiment14()
{
  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/material_model/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/material_model/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/material_model/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  map<ThroatID,Throat> throats = extract_throats(morphologyVolume,stateVolume);

  auto poreNetwork = pore_network(pores,throats);

  nodes_vs_ligaments_distribution(pores,poreNetwork,distanceField);
}

//*************************************************************************************************

void experiment13()
{
  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/material_sample/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/material_sample/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/material_sample/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  map<ThroatID,Throat> throats = extract_throats(morphologyVolume,stateVolume);

  auto poreNetwork = pore_network(pores,throats);

  nodes_vs_ligaments_distribution(pores,poreNetwork,distanceField);
}

//*************************************************************************************************

void experiment12()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/digitalFoam/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/digitalFoam/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/digitalFoam/stateVolume.raw",s);

  bool includeBoundaryPores = true;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);

  VoxelVolume<uint8_t> voronoiSeeds = create_voronoi_seeds(pores, distanceField);

  voxelIO::to_file(voronoiSeeds,"voronoiSeeds_750_750_750_uint8.raw");

}

//*************************************************************************************************

void experiment11()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/digitalFoam/distanceField.raw",s);


  material_voxels_in_direction(0, distanceField, "digitalFoam_material_voxels_x");

  distanceField = to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  material_voxels_in_direction(0, distanceField, "actualFoam_material_voxels_x");
}

//*************************************************************************************************

void experiment10()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/digitalFoam/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/digitalFoam/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/digitalFoam/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  logarithmic_histogram_r_eff("artificial_pores_logarithmic.txt",
                              calculate_effective_radii(calculate_pore_volumes(pores)),
                              sqrt(pores.size()));

  map<PoreID,Pore> smallPores, largePores;
  seperate_pores_by_size(7.0,smallPores,largePores,pores,distanceField);

  cout << endl << "Small Pores: " << smallPores.size() << endl;
  cout << endl << "Large Pores: " << largePores.size() << endl;

}

//*************************************************************************************************

void experiment9()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/digitalFoam_a/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/digitalFoam_a/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/digitalFoam_a/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  logarithmic_histogram_r_eff("artificial_large_pores_logarithmic.txt",
                              calculate_effective_radii(calculate_pore_volumes(pores)),
                              sqrt(pores.size()));


}

//*************************************************************************************************

void experiment8()
{
  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  map<PoreID,Pore> smallPores,largePores;

  seperate_pores_by_size(7.0, smallPores,largePores,pores,distanceField);
  histogramm_r_eff("pores_smaller_7um.txt",
                   calculate_effective_radii(calculate_pore_volumes(smallPores)),
                   sqrt(smallPores.size()));

  distanceField =
      to_voxelVolume<float>("input/carbonFoamDownsampled/merged0.8/distanceField.raw",s);
  morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoamDownsampled/merged0.8/morphologyVolume.raw",s);
  stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoamDownsampled/merged0.8/stateVolume.raw",s);


  pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  seperate_pores_by_size(3.5, smallPores,largePores,pores,distanceField);

  logarithmic_histogram_r_eff("pores_larger_7um_downsampled_log.txt",
                              calculate_effective_radii(calculate_pore_volumes(largePores)),
                              sqrt(largePores.size()));

}

//*************************************************************************************************

void experiment7()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoamDownsampled/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoamDownsampled/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoamDownsampled/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  map<PoreID,Pore> smallPores,largePores;

  seperate_pores_by_size(3.5, smallPores,largePores,pores,distanceField);

  compare_anisotropy(largePores,distanceField,"anisotropy_comparison_7.0-large_pores_downsampled.txt");


}

//*************************************************************************************************

void experiment6()
{
  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  map<PoreID,Pore> smallPores,mediumPores,largePores;
  seperate_pores_by_size(5.0, smallPores, largePores, pores, distanceField);
  seperate_pores_by_size(2.0, largePores, mediumPores, smallPores, distanceField);

  Matrix3d Sigma_a =  cout_average_principal_component_analysis(mediumPores, distanceField);
  SelfAdjointEigenSolver<Matrix3d> eigensolver_a(Sigma_a);
  Matrix3d V_a = eigensolver_a.eigenvectors();

  seperate_pores_by_size(12.0, smallPores, largePores, pores, distanceField);
  seperate_pores_by_size(9.0, largePores, mediumPores, smallPores, distanceField);

  Matrix3d Sigma_b =  cout_average_principal_component_analysis(mediumPores, distanceField);
  SelfAdjointEigenSolver<Matrix3d> eigensolver_b(Sigma_b);
  Matrix3d V_b = eigensolver_b.eigenvectors();

  cout << endl << acos((V_a.col(2).transpose())*V_b.col(2))/M_PI*180.0 << endl;

  cout << endl << acos((V_a.col(0).transpose())*V_b.col(2))/M_PI*180.0 << endl;
}

//*************************************************************************************************

void experiment5()
{
  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  map<PoreID,Pore> smallPores,largePores;
  seperate_pores_by_size(2.0, smallPores, largePores, pores, distanceField);

  compare_anisotropy(largePores, distanceField, "anisotropy_comparison_2.0-large_pores.txt");
}

//*************************************************************************************************

void experiment4()
{

  const Vector3i s(750,750,750);
  VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoamDownsampled/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoamDownsampled/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoamDownsampled/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  map<PoreID,Pore> smallPores,largePores;
  seperate_pores_by_size(3.5, smallPores, largePores, pores, distanceField);

  pore_projections("./",largePores,distanceField);


  distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);

  seperate_pores_by_size(7.0, smallPores,largePores,pores,distanceField);

  pore_projections("./small_pore_projections/",smallPores,distanceField);


}

//*************************************************************************************************

void experiment3()
{

  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = true;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;
  map<ThroatID,Throat> throats = extract_throats(morphologyVolume,stateVolume);

  float r_small_large = 7.0;
  cout_void_fractions(pores,throats,distanceField,r_small_large);

}

//*************************************************************************************************

void experiment2()
{

  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoamDownsampled/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoamDownsampled/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoamDownsampled/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  logarithmic_histogram_r_eff("poreRadiusDistributionDownsampledLog.txt",
                              calculate_effective_radii(
                                calculate_pore_volumes(pores)),
                              sqrt(pores.size()));

}

//*************************************************************************************************

void experiment1()
{

  const Vector3i s(750,750,750);
  const VoxelVolume<float> distanceField =
      to_voxelVolume<float>("input/carbonFoam/merged0.8/distanceField.raw",s);
  VoxelVolume<PoreID> morphologyVolume =
      to_voxelVolume<PoreID>("input/carbonFoam/merged0.8/morphologyVolume.raw",s);
  VoxelVolume<uint8_t> stateVolume =
      to_voxelVolume<uint8_t>("input/carbonFoam/merged0.8/stateVolume.raw",s);

  bool includeBoundaryPores = false;

  map<PoreID,Pore> pores = extract_pores(morphologyVolume, stateVolume, includeBoundaryPores);
  cout << endl << "Pores: " << pores.size() << endl;

  logarithmic_histogram_r_eff("poreRadiusDistributionLog.txt",
                              calculate_effective_radii(
                                calculate_pore_volumes(pores)),
                              sqrt(pores.size()));

  map<PoreID,Pore> smallPores, largePores;
  seperate_pores_by_size(7.0,smallPores,largePores,pores,distanceField);

  cout << endl << "Small Pores: " << smallPores.size() << endl;
  cout << endl << "Large Pores: " << largePores.size() << endl;
}
