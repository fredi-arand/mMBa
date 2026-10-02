#include "CommandLine.h"
#include "DistanceField.h"
#include "PoreMorphology.h"
#include "VoxelVolume.h"
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
static const char *usage =
    R"(Usage: mMBa [options] <input.raw>

Analyzes the void space of a raw voxel volume.

Options:
  --size <WxHxD>            number of voxels in x, y, z, e.g. 400x400x400
  --format <u8|f32>         voxel format: u8 (default) or f32
  --iso <value>             voxels below this value are void space
                            (default: average of minimum and maximum value)
  --output-volumes <dir>    write morphologyVolume.raw and stateVolume.raw
  --visualization <dir>     write the morphology as ppm stacks
  --skeleton <dir>          write the skeleton as pgm stacks
  -h, --help                show this help
)";
//------------------------------------------------------------------------------
struct Options {
  std::string inputPath;
  Vector3l size = Vector3l::Zero();
  std::string format = "u8";
  std::optional<float> isoValue;
  std::string outputVolumesPath;
  std::string visualizationPath;
  std::string skeletonPath;
};
//------------------------------------------------------------------------------
static Vector3l parseSize(const std::string &str) {
  Vector3l s;
  size_t begin = 0;
  for (int i = 0; i < 3; ++i) {
    size_t end = i < 2 ? str.find('x', begin) : str.size();
    if (end == std::string::npos)
      throw std::invalid_argument("--size must look like 400x400x400");
    s(i) = std::stol(str.substr(begin, end - begin));
    begin = end + 1;
  }
  if ((s.array() <= 0).any())
    throw std::invalid_argument("--size must be positive");
  return s;
}
//------------------------------------------------------------------------------
// Returns nullopt if the program should exit (e.g. after printing help)
static std::optional<Options> parse(int argc, char **argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto value = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::invalid_argument(arg + " requires a value");
      return argv[++i];
    };

    if (arg == "-h" || arg == "--help") {
      std::cout << usage;
      return std::nullopt;
    } else if (arg == "--size")
      o.size = parseSize(value());
    else if (arg == "--format")
      o.format = value();
    else if (arg == "--iso")
      o.isoValue = std::stof(value());
    else if (arg == "--output-volumes")
      o.outputVolumesPath = value();
    else if (arg == "--visualization")
      o.visualizationPath = value();
    else if (arg == "--skeleton")
      o.skeletonPath = value();
    else if (arg.rfind("-", 0) == 0)
      throw std::invalid_argument("unknown option " + arg);
    else if (o.inputPath.empty())
      o.inputPath = arg;
    else
      throw std::invalid_argument("more than one input file given");
  }

  if (o.inputPath.empty())
    throw std::invalid_argument("no input file given");
  if (o.size.isZero())
    throw std::invalid_argument("--size is required");
  if (o.format != "u8" && o.format != "f32")
    throw std::invalid_argument("--format must be u8 or f32");
  return o;
}
//------------------------------------------------------------------------------
static void run(const Options &o) {
  const char *path = o.inputPath.c_str();
  DistanceField distanceField =
      o.format == "u8"
          ? DistanceField::create<uint8_t>(o.size, path, o.isoValue)
          : DistanceField::create<float>(o.size, path, o.isoValue);
  PoreMorphology poreMorphology(distanceField);
  poreMorphology.exportSkeletonPath = o.skeletonPath;

  poreMorphology.create_pore_morphology(0.0, 0.0);
  poreMorphology.reduce_throat_volume();
  poreMorphology.merge_pores(0.8);

  if (!o.outputVolumesPath.empty()) {
    VoxelVolume<uint32_t> morphologyVolume;
    VoxelVolume<uint8_t> stateVolume;
    poreMorphology.create_legacy_volumes(morphologyVolume, stateVolume);
    morphologyVolume.export_raw(
        (o.outputVolumesPath + "/morphologyVolume.raw").c_str());
    stateVolume.export_raw((o.outputVolumesPath + "/stateVolume.raw").c_str());
  }

  if (!o.visualizationPath.empty())
    poreMorphology.export_ppm_stacks((o.visualizationPath + '/').c_str());
}
//------------------------------------------------------------------------------
int runFromCommandLine(int argc, char **argv) {
  std::optional<Options> options;
  try {
    options = parse(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "mMBa: " << e.what() << "\n\n" << usage;
    return 1;
  }
  if (options)
    run(*options);
  return 0;
}
//------------------------------------------------------------------------------
} // namespace fred
