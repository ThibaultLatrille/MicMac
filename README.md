**Detecting diversifying selection for a trait from within and between-species genotypes and phenotypes**,\
Thibault Latrille, Melodie Bastian, Theo Gaboriau, Nicolas Salamin,\
_Journal of Evolutionary Biology_, voae084,\
2024,\
[doi.org/10.1093/jeb/voae084](https://doi.org/10.1093/jeb/voae084)

**Compiled binaries and instructions for BayesCode are available at [github.com/ThibaultLatrille/bayescode](https://github.com/ThibaultLatrille/bayescode)**

# MicMac

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.\
The experiments can either run on a local computer or in a cluster configuration (slurm).\
The experiments are meant to run on Linux/Unix/MacOS operating systems.

## Repository structure

The repository is split into six main directories, some of which have subdirectories.\
Within each directory is a `README.md` file which summarizes the purpose of that directory.\
The main directories are as follows:

### **`data_empirical`**
This folder contains the snakemake pipeline (`Snakefile`) and data to run the empirical analysis shown in table 1.\
The details are described in the `README.md` file in the `data_empirical` folder.

### **`data_simulated`**
This folder contains the snakemake pipeline (`Snakefile`) and data to run simulated analysis shown in figure 3 and S2.\
The details are described in the `README.md` file in the `data_simulated` folder.

### **`manuscripts`**
This folder contains the different versions of the manuscript (.tex), supplementary material (.tex) and figures (.pdf).

### **`scripts`**
This folder contains the python scripts used by the snakemake pipelines to run the analysis (see `data_empirical` and `data_simulated`).

### **`simulator`**
This folder contains the simulation program to generate the data for the simulated analysis, required to run the experiments in the `data_simulated` folder.\
The details to compile the simulation program and its usage are described in the `README.md` file in the `simulator` folder.

### **`utils`**
This folder contains the [BayesCode](https://github.com/ThibaultLatrille/bayescode) software, necessary to run the analysis in the `data_simulated` and `data_empirical` folders.\
The details to install BayesCode and its usage are described in the `README.md` file in the `utils` folder.

## Installation and dependencies

### Python and snakemake dependencies

Install python3, snakemake and the necessary python packages.

```
sudo apt install -qq -y python3-dev python3-pip
pip3 install snakemake scipy numpy matplotlib pandas ete3 bio statsmodels --user
```

### BayesCode and simulator dependencies
The [BayesCode](https://github.com/ThibaultLatrille/bayescode) software is necessary to run the analysis in the `data_simulated` and `data_empirical` folders.\
See the `README.md` file in the `utils` folder to install BayesCode.

The simulator program is necessary to run the experiments in the `data_simulated` folder.\
See the `README.md` file in the `simulator` folder to compile the simulation program.

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
