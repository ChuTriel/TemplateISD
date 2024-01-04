# TemplateISD

Implementation of a modified Prange that exploits a known error weight distribution.

## Requirements and Dependencies

For now the png package is probably needed as [M4RI](https://bitbucket.org/malb/m4ri) uses it internally. However, because nowhere on their site it is listed as a dependency, it is also possible that the package is only used when installed on the system.
Furthermore OpenMP is needed.
```bash
sudo apt install autoconf automake make cmake libpng-dev libomp-dev
```

To install OpenMP on MacOs, execute the following command:
```bash
brew install libomp
```

## Setup and Build

To setup and build the project, start with:

```bash
git clone git@github.com:ChuTriel/TemplateISD.git
./setup.sh
```

setup.sh is needed to setup [M4RI](https://bitbucket.org/malb/m4ri) for our use. <br>
To build the project after it is sucessfully prepared, type the following commands and replace \<target\> with what you need:
```bash
mkdir build && cd build
cmake ..
make <target>
```

Valid targets are:

| target | Description |
| --- | --- |
| main | Single-threaded prange (OpenMP is not required to be installed to build this target) |
| prange | Multi-threaded prange (Requires OpenMP to be installed) |
| dumer | Multi-threaded dumer (Requires OpenMP to be installed) |
| bench | Benchmark for prange and dumer |

For example, if you want to build and run the multi-threaded prange variant, the whole command chain looks like this:
```bash
mkdir build && cd build
cmake ..
make prange
./prange
```

## Additional Information

### Main

Nothing special to mention here.

### Prange

```bash
./prange <nr_threads>
```

You can optionally specify the number of threads by providing it as a cli argument. If no argument is provided then the maximum possible number of threads is chosen.

### Dumer

```bash
./dumer <nr_threads>
```

You can optionally specify the number of threads by providing it as a cli argument. If no argument is provided then the maximum possible number of threads is chosen. The parallelization is performed permutation-wise for now. Perhaps nested threading where the birthday decoding is also parallelized will be enabled in the future (using a second cli argument).

