# phyloGQuanti

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Requirements: Clang (or g++)
```bash
sudo apt install g++-5 clang-3.6
```

## Get phyloGQuanti up and running on Linux

### How to download and build

To get phyloGQuanti from a machine connected to the internet, type in a terminal:
```bash
git clone https://github.com/ThibaultLatrille/phyloGQuanti.git
```

This should create a folder called `phyloGQuanti` (the phyloGQuanti root folder). You must go there before compiling phyloGQuanti:

```bash
cd phyloGQuanti
```

Then, to build phyloGQuanti simply run:

```bash
make
```

### How to run phyloGQuanti

Basic usage for phyloGQuanti is (from the phyloGQuanti root folder):

```bash
bin/neutral --help
bin/stabilizing --help
```
