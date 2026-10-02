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

Next, install [CMake](https://cmake.org) and [Eigen](https://eigen.tuxfamily.org),
e.g.

```bash
brew install cmake eigen                # macOS
sudo apt install cmake libeigen3-dev    # Debian/Ubuntu
```

Then build:

```bash
./build.sh
```

Finally, run

```bash
build/mMBa PATH/TO/Berea.raw --size 400x400x400 --visualization PATH/TO/visualization
```

See `build/mMBa --help` for all options. To keep track of your runs, put the
commands in a shell script.

### Tab completion (zsh)

zsh can complete the options from `--help`; add this to your `~/.zshrc`:

```zsh
compdef _gnu_generic mMBa
```

### Original code

The more or less undocumented original code can be found here:
[original_code](./original_code/)

### Troubleshooting

[Contact me 😉](mailto:fredi.arand@gmail.com)
