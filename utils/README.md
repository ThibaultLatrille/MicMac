### **`utils`**
This folder contains the [BayesCode](https://github.com/ThibaultLatrille/bayescode) software, necessary to run the analysis in the `data_simulated` and `data_empirical` folders.

Install BayesCode to obtain the executable `nodetraits` and `readnodetraits` with conda:
```
conda install -c bioconda -c conda-forge bayescode
```

Alternatively, you can compile BayesCode from source code (copied from https://github.com/ThibaultLatrille/bayescode) in the `utils/BayesCode` folder.\
Compiling BayesCode requires the buildsystem `make` and `cmake` and a c++ compiler (`g++` or `clang`).\
Compiling is then done by running the following commands in the `utils/BayesCode` folder:
```bash
cd BayesCode
make tiny
```
This will create the executable `nodetraits` and `readnodetraits` in the `BayesCode/bin` folder.

The file `bayescode_nodetraits.pdf` contains for documentation for `nodetraits` and `readnodetraits` used in this study.\
The full documentation of BayesCode can be found in the supplementary materials at https://github.com/ThibaultLatrille/bayescode/wiki.
