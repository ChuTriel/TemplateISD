# TemplateISD

Implementation of a modified Prange that exploits a known error weight distribution.

## Requirements and Dependencies

For now the png package is probably needed as [M4RI](https://bitbucket.org/malb/m4ri) uses it internally. However, because nowhere on their site it is listed as a dependency, it is also possible that the package is only used when installed on the system.
```bash
sudo apt install libpng-dev
```
## Setup and Build

To setup and build the project, start with:

```bash
git clone git@github.com:ChuTriel/TemplateISD.git
./setup.sh
```

setup.sh is needed to setup [M4RI](https://bitbucket.org/malb/m4ri) for our use. <br>
To build the project after it is sucessfully prepared, type the following commands:

```bash
mkdir build && cd build
cmake ..
make main
```

You can run the program with
```bash
./main
```