# MicMac

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Requirements: Clang (or g++)
```bash
sudo apt install g++-5 clang-3.6
```

## Get MicMac up and running on Linux

### How to download and build

To get MicMac from a machine connected to the internet, type in a terminal:
```bash
git clone https://github.com/ThibaultLatrille/MicMac.git
```

This should create a folder called `MicMac` (the MicMac root folder). You must go there before compiling MicMac:

```bash
cd MicMac
```

Then, to build MicMac simply run:

```bash
make
```

### How to run MicMac

Basic usage for MicMac is (from the MicMac root folder):

```bash
bin/neutral --help
bin/stabilizing --help
```
