# j-eri-regent

Standalone Regent code that calculates two-electron repulsion integrals (ERIs) to form the J matrix, a building block in many quantum chemistry methods. Regent enables straightforward execution on multiple architectures: CPUs, GPUs, and multiple nodes.

Authors: K. Grace Johnson, Ellis Hoag, Seema Mirchandaney, Alex Aiken, Todd J. Martinez

## Setup

Running Regent code requires installing the Legion programming system. The dependencies of Legion are:
- Linux, macOS, or another Unix
- A C++ 98 (or newer) compiler (GCC, Clang, Intel, or PGI) and GNU Make
- Python (for Regent and build, profiling, and debugging tools)
- *Optional*: CUDA 5.0 or newer (for NVIDIA GPUs)
- *Optional*: GASNet (for networking, see installation instructions here: https://legion.stanford.edu/gasnet/)

### Legion Installation
Optional arguments for compiling with GPU support (CUDA) and multi-node support (GASNet) are commented. 

```bash
git clone https://gitlab.com/StanfordLegion/legion.git -b hijack_registration_hack
export LEGION_SRC=$PWD/legion
export CC=<path to icc or gcc, etc.>
export CXX=<path to icc or g++, etc.>
#export GASNET=<path to gasnet install directory>             # Optional, multi-node
#export CONDUIT=ibv                                           # Optional, multi-node (also set as `mpi`)
#export GPU_ARCH=maxwell                                      # Optional, GPU (set to desired arch.)
#export USE_CUDA=1                                            # Optional, GPU
#export CUDA_BIN_PATH=$CUDA_HOME                              # Optional, GPU
$LEGION_SRC/language/scripts/setup_env.py --cmake  \
    --terra-url https://github.com/StanfordLegion/terra.git \
    --terra-branch luajit2.1 \
    --extra="-DCMAKE_INSTALL_PREFIX=$LEGION_INSTALL_PATH" \
    --extra="-DCMAKE_BUILD_TYPE=Release" \
    --extra="-DLegion_HIJACK_CUDART=OFF" \
#   --extra="--with-gasnet ${GASNET}" \                       # Optional, multi-node
#   --extra="--cuda" \                                        # Optional, GPU
export REGENT=$LEGION_SRC/language/regent.py
alias regent=$REGENT
```

Further information on Legion installation can be found here: https://legion.stanford.edu/starting/

## Running

### Running a pre-compiled version

### Notes on angular momentum and compilation time

Be sure to select the appropriate angular momentum using the `-L [S|P|D|F|G]` option. This will tell Lua to produce the correct number of Regent tasks. Higher angular momentums need more and larger kernels which can take a longer time to compile to CUDA code. The number of J kernels needed is <code>(2L-1)<sup>2</sup></code>.

| Angular Momentum | Number of J Kernels | Memory     | Compilation wall-time |
|:----------------:|:-------------------:|:----------:|:---------------------:|
| S = 1            | 1                   | Negligible | < 1 Minute            |
| P = 2            | 9                   | 2 GB       | 2 Minutes             |
| D = 3            | 25                  | > 4 GB     | > 5 Minutes           |
| F = 4            | 49                  | > 7 GB     | > 7 Minutes           |
| G = 5            | 81                  | > 31 GB    | > 1 Hour              |


### Submodule

Normally this project (eri-regent) is a submodule of the larger TeraChem project and is located at production/regintbox/src/eri-regent.  
When this is the case you should build Legion using the scripts in production/scripts/xstream, such as build_legion_k80.bash.

### Standalone

If you are building eri-regent in a standalone configuration follow these instructions:

When more than 1 GB of memory is needed, you must build Legion with `luajit2.1`.
Instructions for building on Ubuntu Linux:

```bash
cd $LEGION_DIR/language
./install.py --cmake --terra-url https://github.com/StanfordLegion/terra.git --terra-branch luajit2.1
make -C build install
export REGENT=/usr/local/bin/regent.py
alias regent=$REGENT
```

Instructions for building on xstream:
```bash
cd $LEGION_DIR/language
module load GCC/4.9.2-binutils-2.25  
module load OpenMPI/1.8.5
module load Python/3.6.0
module load CMake/3.5.2
module load CUDA/8.0.61
module load LLVM/3.7.0
export CONDUIT=ibv
export USE_CUDA=1
export CUDA=$CUDA_HOME
export CC=gcc
export CXX=g++
export LEGION_SRC=$HOME/work/legion
export LEGION_INSTALL_PATH=$HOME/work/legion_install
./scripts/setup_env.py --cmake  --terra-url https://github.com/StanfordLegion/terra.git --terra-branch luajit2.1 --extra=-DCMAKE_INSTALL_PREFIX=$LEGION_INSTALL_PATH
export REGENT=$LEGION_DIR/language/regent.py
alias regent=$REGENT
```

## Building

Use the Makefile to compile and run inside C++. This will generate a header file and a library for the eri tasks so they can be called within C++. The Makefile assumes the `RG_MAX_MOMENTUM` environment variable has been set. If you want to compile for a new `RG_MAX_MOMENTUM` then you need to run `make rg.clean` before the environment variable will affect the build.

```bash
cd eri-regent
export RG_MAX_MOMENTUM=P
make
```

## Running and Testing
Run with Regent using `top_jfock.rg` or `top_kfock.rg` inside `src/` for testing. Note that running eri-regent with this method does not require you to run `make`.

```bash
cd eri-regent/src
# To run JFock algorithm
regent top_jfock.rg -L P -i tests/integ/h2o -v tests/integ/h2o/output.dat
# To run KFock algorithm
regent top_kfock.rg -L S -i tests/integ/h2 -v tests/integ/h2/kfock_output.dat
# Use option `-fflow 0` to compile eri-regent faster
```

To test eri-regent with C++, compile the test program inside `src/tests/cpp` after building eri-regent.
```bash
cd eri-regent/src/tests/cpp
make
```

This will produce a binary inside `eri-regent/build`.
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LEGION_DIR/language/build/lib/
cd eri-regent
# To run JFock algorithm
build/eri_regent_test -i src/tests/integ/h2o -a jfock
# To run KFock algorithm
build/eri_regent_test -i src/tests/integ/h2o -a kfock
```

### Higher Angular Momentum Systems

Be sure to select the appropriate angular momentum using the `-L [S|P|D|F|G]` option. This will tell Lua to produce the correct number of Regent tasks. Higher angular momentums need more and larger kernels which can take a long time to compile to Cuda code. The number of J kernels needed is <code>(2L-1)<sup>2</sup></code> and the number of K kernels needed is <code>L<sup>2</sup> * (L<sup>2</sup> + 1) / 2</code>.

| Angular Momentum | Number of J Kernels | Number of K Kernels | Memory     | Wall-time   |
|:----------------:|:-------------------:|:-------------------:|:----------:|:-----------:|
| S = 1            | 1                   | 1                   | Negligible | < 1 Minute  |
| P = 2            | 9                   | 10                  | 2 GB       | 2 Minutes   |
| D = 3            | 25                  | 45                  | > 4 GB     | > 5 Minutes |
| F = 4            | 49                  | 136                 | > 7 GB     | > 7 Minutes |
| G = 5            | 81                  | 325                 | > 31 GB    | > 1 Hour    |

## Testing with Python3
First compile the test program in `eri-regent/src/tests/cpp`, then you can use `python3` to run the binary on all test inputs.
```bash
python scripts/test.py
python scripts/test_boys.py
```
