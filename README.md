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

For example, if you want to build and run the multi-threaded prange variant, the whole command chain looks like this:
```bash
mkdir build && cd build
cmake ..
make prange
./prange
```

