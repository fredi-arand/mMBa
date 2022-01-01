#include "VoxelVolume.h"
#include "distancefield.h"
#include "examplespaper.h"
#include "importhomberg.h"
#include "poremorphology.h"
#include "thesis_helpers.h"
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace fred;

void carbonFoamPaper1();
void carbonFoamPaper2();
void carbonFoamPaper3();
void carbonFoamPaper4();
void carbonFoamPaper5();
void carbonFoamPaper5a();
void thesis1();
void thesis2();
void krakowska1();
void carbonFoamPaper6(); // analyze material of sample without small pores
void carbonFoamPaper7(); // analyze material of model withtout small pores
void carbonFoamPaper8(); // modified model
void thesis3();          // mMBa verification
void thesis4();          // mMBa verification (see cageo article)
void thesis5();          // mMBa verification (I promise, last one)
void thesis6();          // 8 segmentations of subregions
void thesis7();          // downsampled large region
void thesis8();          // 8 segmentations (model)
void thesis9();          // downsampled large region (model)
void berea();

int main(int, char **) {

  berea();
  return 0;
}
//------------------------------------------------------------------------------
void berea() {
  Vector3l s(400, 400, 400);
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(s, "../../Berea.raw");
  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);
  poreMorphology.export_ppm_stacks("../../output/stacks/");
}
//------------------------------------------------------------------------------

void thesis9() {
  Vector3l s(1000, 1000, 1000);

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      s,
      "/Users/arand/from_vgpc32/build/P62/input_morphology/"
      "full_foam_downsampled.raw",
      4.0);
  PoreMorphology poreMorphology(distanceField);

  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);
  //  poreMorphology.export_ppm_stacks("output/stacks/");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
  morphologyVolume.export_raw("/Users/arand/from_vgpc32/build/P41/input/"
                              "8_subvolumes_model/large_morphologyVolume.raw");
  stateVolume.export_raw("/Users/arand/from_vgpc32/build/P41/input/"
                         "8_subvolumes_model/large_stateVolume.raw");
}

//------------------------------------------------------------------------------

void thesis8() {
  Vector3l s(1000, 1000, 1000);

  for (int i = 0; i < 8; ++i) {
    DistanceField distanceField;
    distanceField.create_distance_field<uint8_t>(
        s,
        (string(
             "/Users/arand/from_vgpc32/build/P62/input_morphology/full_foam_") +
         to_string(i) + ".raw")
            .c_str(),
        4.0);
    PoreMorphology poreMorphology(distanceField);

    poreMorphology.create_pore_morphology(0.0, 0.0);
    poreMorphology.reduce_throat_volume();
    poreMorphology.merge_pores(0.8);

    if (i == 0)
      poreMorphology.export_ppm_stacks("output/stacks/");

    VoxelVolume<uint32_t> morphologyVolume;
    VoxelVolume<uint8_t> stateVolume;
    poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
    string pathBase =
        string("/Users/arand/from_vgpc32/build/P41/input/8_subvolumes_model/") +
        to_string(i);
    morphologyVolume.export_raw((pathBase + "morphologyVolume.raw").c_str());
    stateVolume.export_raw((pathBase + "stateVolume.raw").c_str());
  }
}

//------------------------------------------------------------------------------

void thesis7() {
  Vector3l s(1000, 1000, 1000);

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      s,
      "../../volumedata/ITV_carbon_foam/"
      "carbonFoam2_1000_1000_1000_uint8_downsampled.raw",
      144.5);
  PoreMorphology poreMorphology(distanceField);

  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);
  //  poreMorphology.export_ppm_stacks("output/stacks/");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
  morphologyVolume.export_raw(
      "../P41/input/8_subvolumes/large_morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/8_subvolumes/large_stateVolume.raw");
}

//------------------------------------------------------------------------------

void thesis6() {
  Vector3l s(1000, 1000, 1000);

  for (int i = 0; i < 8; ++i) {
    DistanceField distanceField;
    distanceField.create_distance_field<uint8_t>(
        s,
        (string("../../volumedata/ITV_carbon_foam/subvolumes_2000/vol") +
         to_string(i) + "_1000_1000_1000_uint8_denoised_connected.raw")
            .c_str(),
        144.5);
    PoreMorphology poreMorphology(distanceField);

    poreMorphology.create_pore_morphology(0.0, 0.0);
    poreMorphology.reduce_throat_volume();
    poreMorphology.merge_pores(0.8);

    if (i == 0)
      poreMorphology.export_ppm_stacks("output/stacks/");

    VoxelVolume<uint32_t> morphologyVolume;
    VoxelVolume<uint8_t> stateVolume;
    poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

    string pathBase = string("../P41/input/8_subvolumes/") + to_string(i);
    morphologyVolume.export_raw((pathBase + "morphologyVolume.raw").c_str());
    stateVolume.export_raw((pathBase + "stateVolume.raw").c_str());
  }
}

//------------------------------------------------------------------------------

void thesis5() {
  Vector3l s(128, 128, 128);
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      s, "../P57/artificial_128_128_128_uint8.raw", 127.5);
  PoreMorphology poreMorphology(distanceField);

  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);
  //  poreMorphology.export_ppm_stacks("output/stacks/");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
  distanceField.export_raw("../P41/input/verification2/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/verification2/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/verification2/stateVolume.raw");
}

//------------------------------------------------------------------------------

void thesis4() {
  Vector3l s(512, 512, 512);
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      s, "../P55/hcp_512_512_512_uint8.raw", 127.5);
  PoreMorphology poreMorphology(distanceField);

  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);
  poreMorphology.export_ppm_stacks("output/stacks/");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
  distanceField.export_raw("../P41/input/hcp/distanceField.raw");
  morphologyVolume.export_raw("../P41/input/hcp/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/hcp/stateVolume.raw");
}

//------------------------------------------------------------------------------

void thesis3() {

  Vector3l s(500, 500, 500);
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      s, "../P54/random_balls_downsampled_500_500_500_uint8_t.raw", 127.5);

  for (float epsilon = 0.0; epsilon <= 0.21; epsilon += 0.05) {
    PoreMorphology poreMorphology(distanceField);
    poreMorphology.epsilon = epsilon;
    poreMorphology.create_pore_morphology(1.0, 0.0);
    poreMorphology.reduce_throat_volume();
    poreMorphology.merge_pores(0.8);

    //    poreMorphology.export_ppm_stacks("output/stacks/");

    VoxelVolume<uint32_t> morphologyVolume;
    VoxelVolume<uint8_t> stateVolume;
    poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

    string pathBegin = "../P41/input/thesis_verification/";
    string pathEnd = to_string(int(20 * epsilon)) + ".raw";
    distanceField.export_raw((pathBegin + "distanceField_" + pathEnd).c_str());
    morphologyVolume.export_raw(
        (pathBegin + "morphologyVolume_" + pathEnd).c_str());
    stateVolume.export_raw((pathBegin + "stateVolume_" + pathEnd).c_str());
  }
}

//------------------------------------------------------------------------------

void carbonFoamPaper8() {
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(750, 750, 750),
      "../../volumedata/ITV_carbon_foam/artificial_large_and_small_pores/"
      "eroded_dilated_gauss_750_750_750_8bit.raw",
      127.5);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  distanceField.export_raw(
      "../P41/input/digitalFoam_modified/distanceField.raw");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  morphologyVolume.export_raw(
      "../P41/input/digitalFoam_modified/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/digitalFoam_modified/stateVolume.raw");

  //  poreMorphology.export_ppm_stacks("output/stacks/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper7() {
  Vector3l const s(750, 750, 750);

  VoxelVolume<uint8_t> voxelVolume;
  voxelVolume.import_raw_volume(
      s, "../../volumedata/ITV_carbon_foam/artificial_large_pores/"
         "750cubed_gauss_8bit.raw");

  float threshold = 127.5;
  for (auto &value : voxelVolume.voxelValues)
    value = 255 - value;

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(voxelVolume, 255.0 - threshold);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  //  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw("../P41/input/material_model/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/material_model/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/material_model/stateVolume.raw");

  //  poreMorphology.export_ppm_stacks("output/stacks/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper6() {
  Vector3l const s(750, 750, 750);

  VoxelVolume<uint8_t> voxelVolume;
  //  voxelVolume.import_raw_volume(s,
  //  "carbonFoam_denoised_connected_750_750_750.raw");
  voxelVolume.import_raw_volume(
      s, "../../volumedata/ITV_carbon_foam/sample_without_small_pores/"
         "sample_without_small_pores_750_750_750_uint8.raw");

  float threshold = 119.5;
  for (auto &value : voxelVolume.voxelValues)
    value = 255 - value;

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(voxelVolume, 255.0 - threshold);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  //  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw("../P41/input/material_sample/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/material_sample/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/material_sample/stateVolume.raw");

  poreMorphology.export_ppm_stacks("output/stacks/");
}

//------------------------------------------------------------------------------

void krakowska1() {
  DistanceField distanceField;
  distanceField.create_distance_field<float>(
      Vector3l(1500, 1500, 400), "../../volumedata/paulina/12888_1.vol", 23.0);
  //  distanceField.create_distance_field<uint8_t>(
  //        Vector3l(750,750,750),"carbonFoam_denoised_connected_750_750_750.raw",120);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  poreMorphology.export_ppm_stacks("output/stacks/");
}

//------------------------------------------------------------------------------

void thesis2() {

  VoxelVolume<uint8_t> coverageRepresentation0;
  coverageRepresentation0.import_raw_volume(
      Vector3l(128, 128, 1),
      "../../volumedata/miniVolume3/blub_128_128_1_uint8.raw");

  gnuplot_volume_to_images("thesis/other_", coverageRepresentation0);

  VoxelVolume<uint8_t> coverageRepresentation;
  coverageRepresentation.import_raw_volume(
      Vector3l(64, 64, 1),
      "../../volumedata/miniVolume3/blub_64_64_1_uint8.raw");

  gnuplot_volume_to_images("thesis/", coverageRepresentation);

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(64, 64, 1),
      "../../volumedata/miniVolume3/blub_64_64_1_uint8.raw", 127.5);

  gnuplot_distance_field_and_maximal_balls("thesis/", distanceField);

  // maximal balls, step by step... see thesis_helpers.h
  mb_step_by_step(distanceField, 1.0, 0.0);

  // use algorithm as is

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  gnuplot_palette_file("thesis/reduced_throats.txt", poreMorphology);
  poreMorphology.merge_pores(0.8);
  gnuplot_palette_file("thesis/merged.txt", poreMorphology);
}

//------------------------------------------------------------------------------

void thesis1() {

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(8, 8, 1), "../../source/thesis/chapter_2/image_8_8_1_uint8.raw",
      8.0);

  gnuplot_distance_field_and_maximal_balls("../../source/thesis/chapter_2/",
                                           distanceField);
}

//------------------------------------------------------------------------------

void carbonFoamPaper5() {
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(750, 750, 750),
      "../../volumedata/ITV_carbon_foam/artificial_large_and_small_pores/"
      "750cubed_8bit.raw",
      127.5);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  distanceField.export_raw("../P41/input/digitalFoam/distanceField.raw");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  morphologyVolume.export_raw("../P41/input/digitalFoam/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/digitalFoam/stateVolume.raw");

  poreMorphology.export_ppm_stacks(
      "output/morphology_carbonFoam_denoised_connected/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper5a() {
  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(750, 750, 750),
      "../../volumedata/ITV_carbon_foam/artificial_large_and_small_pores/"
      "artificial_foam_downsampled.raw",
      127.5);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  distanceField.export_raw("../P41/input/digitalFoam_a/distanceField.raw");

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  morphologyVolume.export_raw(
      "../P41/input/digitalFoam_a/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/digitalFoam_a/stateVolume.raw");

  //  poreMorphology.export_ppm_stacks("output/morphology_carbonFoam_denoised_connected/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper4() {

  DistanceField distanceField;
  distanceField.create_distance_field<float>(Vector3l(750, 750, 750),
                                             "../P42/upsampled2.raw", 120.0);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw(
      "../P41/input/carbonFoamUpsampled2/merged0.8/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/carbonFoamUpsampled2/merged0.8/morphologyVolume.raw");
  stateVolume.export_raw(
      "../P41/input/carbonFoamUpsampled2/merged0.8/stateVolume.raw");

  poreMorphology.export_ppm_stacks(
      "output/morphology_carbonFoam_denoised_connected/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper3() {

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(750, 750, 750),
      "../P39/carbonFoam1500_downsampled_denoised_connected", 120.0);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw(
      "../P41/input/carbonFoamDownsampled/merged0.8/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/carbonFoamDownsampled/merged0.8/morphologyVolume.raw");
  stateVolume.export_raw(
      "../P41/input/carbonFoamDownsampled/merged0.8/stateVolume.raw");

  //  poreMorphology.export_ppm_stacks("output/morphology_carbonFoam_denoised_connected/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper2() {

  DistanceField distanceField;
  distanceField.create_distance_field<float>(Vector3l(750, 750, 750),
                                             "../P42/upsampled.raw", 120.0);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw(
      "../P41/input/carbonFoamUpsampled/merged0.8/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/carbonFoamUpsampled/merged0.8/morphologyVolume.raw");
  stateVolume.export_raw(
      "../P41/input/carbonFoamUpsampled/merged0.8/stateVolume.raw");

  poreMorphology.export_ppm_stacks(
      "output/morphology_carbonFoam_denoised_connected/");
}

//------------------------------------------------------------------------------

void carbonFoamPaper1() {

  DistanceField distanceField;
  distanceField.create_distance_field<uint8_t>(
      Vector3l(750, 750, 750), "carbonFoam_denoised_connected_750_750_750.raw",
      120);

  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(1.0, 0);
  poreMorphology.reduce_throat_volume();

  //  distanceField.export_raw("../P41/input/carbonFoam/notMerged/distanceField.raw");
  //  poreMorphology.morphologyVolume.export_raw("../P41/input/carbonFoam/notMerged/morphologyVolume.raw");
  //  poreMorphology.stateVolume.export_raw("../P41/input/carbonFoam/notMerged/stateVolume.raw");

  poreMorphology.merge_pores(0.8);

  VoxelVolume<uint32_t> morphologyVolume;
  VoxelVolume<uint8_t> stateVolume;
  poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);

  distanceField.export_raw(
      "../P41/input/carbonFoam/merged0.8/distanceField.raw");
  morphologyVolume.export_raw(
      "../P41/input/carbonFoam/merged0.8/morphologyVolume.raw");
  stateVolume.export_raw("../P41/input/carbonFoam/merged0.8/stateVolume.raw");

  poreMorphology.export_ppm_stacks(
      "output/morphology_carbonFoam_denoised_connected/");
}
