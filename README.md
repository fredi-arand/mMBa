# Modified Maximal Ball Algorithm (mMBa)

Algorithm which was used for my
[PhD thesis](https://archiv.ub.uni-heidelberg.de/volltextserver/26476/).
It is a slightly improved version of the algorithm presented in
[Computers & Geosciences](https://www.sciencedirect.com/science/article/pii/S0098300416305180).

Being a physicist, I had to learn programming, so excuse the mess 😉

## First steps

Typically, you might want to "see something". In order to do so, download the
[Berea sandstone](https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/berea-sandstone/)
sample.

Next, have a look at `demo.sh` (and maybe `src/main.cpp`). If you think that it
will probably work out, run

```bash
./demo.sh
```

Next, you'll need to create a yaml file, for example `berea.yaml`:

```yaml
version: 0.1
input volume:
  path: PATH/TO/Berea.raw
  voxel format: u8
  resolution:
    width: 400
    height: 400
    depth: 400
  iso value: automatic
output volumes path: PATH/TO/volumes
visualization path: PATH/TO/visualization
```

See also [YamlSpecification.md](./YamlSpecification.md).

Finally, run

```bash
build/bin/mMBA PATH/TO/berea.yaml
```

### Troubleshooting

[Contact me 😉](mailto:fredi.arand@gmail.com)
