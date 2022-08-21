#include "RunFromYaml.h"
#include "DistanceField.h"
#include "PoreMorphology.h"
#include "VoxelVolume.h"
#include <optional>
#include <string>
#include <yaml-cpp/yaml.h>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
static DistanceField createDistanceFile(const std::string &voxelFormat,
                                        const Vector3l &s, const char *path,
                                        const std::optional<float> &isoValue) {
  if (voxelFormat.compare("u8") == 0)
    return DistanceField::create<uint8_t>(s, path, isoValue);
  else if (voxelFormat.compare("f32") == 0)
    return DistanceField::create<float>(s, path);
  return DistanceField::create<float>(s, nullptr, isoValue);
}
//------------------------------------------------------------------------------
void runFromYaml(const char *yamlPath) {

  YAML::Node config = YAML::LoadFile(yamlPath);
  YAML::Node cIv = config["input volume"];
  auto path = cIv["path"].as<std::string>();
  YAML::Node cIvR = cIv["resolution"];
  auto toLong = [&cIvR](const char *str) { return cIvR[str].as<long>(); };
  Vector3l s(toLong("width"), toLong("height"), toLong("depth"));
  std::optional<float> isoValue;
  YAML::Node cIvIv = cIv["iso value"];
  // Not happy with using this, but not performance relevant
  try {
    isoValue = cIvIv.as<float>();
  } catch (const YAML::BadConversion &) {
  }
  auto voxelFormat = cIv["voxel format"].as<std::string>();

  DistanceField distanceField =
      createDistanceFile(voxelFormat, s, path.c_str(), isoValue);
  PoreMorphology poreMorphology(distanceField);
  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  if (YAML::Node outputVolumesPath = config["output volumes path"]) {
    VoxelVolume<uint32_t> morphologyVolume;
    VoxelVolume<uint8_t> stateVolume;
    poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
    auto thePath = outputVolumesPath.as<std::string>();
    morphologyVolume.export_raw((thePath + "/morphologyVolume.raw").c_str());
    stateVolume.export_raw((thePath + "/stateVolume.raw").c_str());
  }

  if (YAML::Node visualizationPath = config["visualization path"]) {
    auto thePath = visualizationPath.as<std::string>() + '/';
    poreMorphology.export_ppm_stacks(thePath.c_str());
  }
}
//------------------------------------------------------------------------------
} // namespace fred
