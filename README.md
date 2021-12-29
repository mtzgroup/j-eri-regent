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

### Pull j-eri-regent code
```bash
git clone https://github.com/GraceJohnson/j-eri-regent.git
#TODO: change this to mtzgroup when finalized
```

## Running

```bash
export LEGION_SRC=<path to dir with Legion build>/legion
export REGENT=$LEGION_SRC/language/regent.py
alias regent=$REGENT
```
Run with Regent using `top_jfock.rg` which contains the top level task, the Regent equivalent of main in C/C++.
```bash
cd j-eri-regent/src
# To run J matrix algorithm:
regent top_jfock.rg -L P -i ../tests/h2o -v ../tests/h2o/output.dat
# To partition tasks and run in parallel on 2 GPUs:
regent top_jfock.rg -L P -i ../tests/h2o -v ../tests/h2o/output.dat -p 2 -ll:gpu 2
```
Note that this command will both compile and execute the Regent code.

#### Options
- `-L [S|P|D|F|G]` specifies the max angular momentum. Compiler will generate all kernels up to and including those containing `L`.
- `-i` specifies path to directory containing input files (see below)
- `-v` verify output with reference data in this file (see below)
- `-p` specifies the number of partitions of each integral task. Default is 1.
- `-ll:gpu` directive passed to Legion specifying the number of GPUs to parallelize across.
- `-ll:cpu` directive passed to Legion specifying the number of CPUs to parallelize across. See all Legion command-line flags here: https://legion.stanford.edu/starting/
- `-fflow 0` compiles kernels faster
- `-h` print usage (including these and more options) and exit

#### Tests
The `tests` directory contains 5 tests on different systems with different max angular momentum:
- `h2`: one hydrogen molecule, 6-311g basis. L = S
- `h2o`: one water molecule, 6-311g basis. L = P
- `co2`: one carbon dioxide molecule, 6-311g basis. L = P
- `dna-pair`: one DNA base pair, 6-31g basis. L = P
- `fe`: one iron atom, 6-31g basis. L = D

Each test has sample data generated from TeraChem on the first SCF iteration of an RHF calculation in the same manner as described in the paper.
- `bras.dat` coordinates (x,y,x) and Gaussian basis information (eta, C) of the bra in the Hermite basis
- `kets.dat` coordinates (x,y,x) and Gaussian basis information (eta, C) of the ket in the Hermite basis  and corresponding density value(s)
- `parameters.dat` parameters from TeraChem input, only `thredp` is currrently used (for Schwartz screening)
- `output.dat` output data (in Hermite basis) generated from TeraChem used to verify the Regent calculation

New input tests can be generated for any systems/angular momenta by conforming to the file formats in the `.dat` files (i.e. separated by angular momentum group, written in hex, labeled, and in the Hermite basis).

### Code structure
- `src`: contains the top level task in`top_jfock.rg`, the driver for kernel generation and execution in `jfock.rg`, region specifications in `fields.rg`, and helper functions in `helper.rg`.
- `src/utils`: code to read and parse inputs
- `src/md`: code to compute intergrals using the McMurchie-Davidson algorithm

### Notes on angular momentum and compilation time

Be sure to select the appropriate angular momentum using the `-L [S|P|D|F|G]` option. This will tell Lua to produce the correct number of Regent tasks. Higher angular momenta require more and larger kernels which can take a longer time to compile to CUDA code. The number of J kernels needed is <code>(2L-1)<sup>2</sup></code>.

| Max Angular Momentum | Number of J Kernels | Memory     | Compilation wall-time |
|:--------------------:|:-------------------:|:----------:|:---------------------:|
| S = 0                | 1                   | Negligible | < 1 Minute            |
| P = 1                | 9                   | 2 GB       | 2 Minutes             |
| D = 2                | 25                  | > 4 GB     | > 5 Minutes           |
| F = 3                | 49                  | > 7 GB     | > 7 Minutes           |
| G = 4                | 81                  | > 31 GB    | > 1 Hour              |
