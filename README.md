# MolOPTGA

[//]: # (Badges)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is a repository of MolOPTGA(Multi-molecular properties Optimization with Genetic Algorithm), MolOPTGA can use 
the Genetic Algorithm and Pareto Front to optimize the molecular properties, including vina score, plogp and sascore.

### Installation

Run:

`conda env create -f environment.yml`

This will take care of installing all required dependencies.

The only required dependency is the latest Conda package manager, which you can download with the Anaconda Python distribution [here](https://www.anaconda.com/distribution/).

### Running

You can use MolOPTGA to optimize the plogp, vina score and sascore of the molecule.

`python main.py -i example.txt -r pro.pdbqt -c folder1 -p folder2 -o folder3 -s ./src/components/QuickVina2-GPU-2-1-CUDA -c_x 32.068 -c_y 16.86 -c_z 131.12  -ps 40 -cs 400 -t 10`
