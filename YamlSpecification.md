# Specification for yaml file

A yaml file is necessary to set up your program.

## Version 0.1

```yaml
version: 0.1
input volume:
  path: PATH
  voxel format: VOXEL_FORMAT
  resolution:
    width: WIDTH
    height: HEIGHT
    depth: DEPTH
  iso value: ISO_VALUE
output volumes path: OUTPUT_VOLUMES_PATH
visualization path: VISUALIZATION_PATH
```

`PATH` specifies the path to the volume you wish to analyze.
The volume must store voxels directly as bytes.  
Example: `/Users/BallsExpert/mmba/Berea.raw`

`VOXEL_FORMAT`: might take one of the following values:  
`u8`: unsigned integer (8 bit)  
`f32`: floating point (32 bit)  

`WIDTH`, `HEIGHT`, `DEPTH`: Number of voxels in x, y, z direction.  
Example: `256`

`ISO_VALUE`: A value that distinguishes void space from material.
Values smaller than `ISO_VALUE` belong to void space and will be analyzed.
Note that you might also write `automatic`.
Then, `ISO_VALUE` is set to the average of minimum and maximum voxel value.  
Example: `0.5`

(optional)
`OUTPUT_VOLUMES_PATH` specifies a path to store output volumes which can be
analyzed with another program (ask me for it!).
It must be an existing directory.
After the algorithm finishes, it will contain two files:
`MorphologyVolume.raw` and `StateVolume.raw`.  
Example: `/Users/BallsExpert/tmp/BereaMorphology/`

(optional)
`VISUALIZATION_PATH` specifies a path to store ppm files that visualize the
morphology.
It must be an existing directory.
After the algorithm finishes, it will contain `ppm` files which represent stacks
in z-direction, starting with `stack000000.ppm`
