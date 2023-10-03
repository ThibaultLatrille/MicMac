**A test of diversifying selection for a trait from within and between species variations**\
Thibault Latrille, Mélodie Bastian, Théo Gaboriau, Nicolas Salamin\
_bioRxiv_\
[doi.org/10.1101/2023.10.02.559886](https://doi.org/10.1101/2023.10.02.559886)

**Compiled binaries and instructions for BayesCode are available
at [github.com/ThibaultLatrille/bayescode](https://github.com/ThibaultLatrille/bayescode)**

# MicMac

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).

The experiments are meant to run on Linux/Unix/MacOS operating systems.

If problems and/or questions are encountered, feel free
to [open issues](https://github.com/ThibaultLatrille/MicMac/issues).

## 0. Local copy

Clone the repository and `cd` to the dir.

```
git clone https://github.com/ThibaultLatrille/MicMac
cd MicMac
```

## 1. Installation

### General dependencies

Install python3 packages

```
sudo apt install -qq -y python3-dev python3-pip
pip3 install snakemake scipy numpy matplotlib pandas ete3 bio statsmodels --user
```

Install [BayesCode](https://github.com/ThibaultLatrille/bayescode) to obtain the executable `nodetraits`
and `readnodetraits`:

```
conda install -c bioconda -c conda-forge bayescode
```

Alternatively, you can compile BayesCode from source in the `utils` folder:

```bash
mkdir utils
cd utils
git clone https://github.com/ThibaultLatrille/BayesCode
cd BayesCode
make tiny
``` 

## 2. Run empirical analysis

In folder `data_empirical` run `snakemake`:

```
cd data_empirical
snakemake -j 8 -k
```

## 3. Run simulated analysis

Compile simulation program (requires `g++`, `make` and `cmake`):

```
cd simulator
make release
cd ..
```

In folder `data_simulated` run one the experiment with the configuration file:

```
cd data_simulated
python3 simulated_experiment.py -c constant_pop_size.yaml -j 32
```

This will create a folder `data_simulated/constant_pop_size` with the configuration file `constant_pop_size.yaml` and
run the experiment with 32 threads.

## Add features or debug in the python scripts

You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille8/MicMac/pulls)

## Licence

The MIT License (MIT)

Copyright (c) 2023 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
