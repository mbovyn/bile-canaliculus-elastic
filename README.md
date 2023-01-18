# Elastic model of bile canaliculus

This repository hosts the code and results for the paper "Hepatocyte apical bulkheads provide a mechanical means to oppose bile pressure" by Bebelman and Bovyn et al. JCB 2023.

We use a elastic model of the bile canaliculus, implemented in the [Surface Evolver](http://facstaff.susqu.edu/brakke/evolver/evolver.html) programming language, to show that apical bulkheads increase the ability of bile canaliculi to hold pressure.

The repo contains:
- A comparison of Surface Evolver's implementation of elasticity to theory (test of elastic sphere)
- Code generated to create and manipulate the the model of the bile canaliculus (elastic_BC_code)
- A simple example which runs the code and Surface Evolver (simple_example)
- Code which generates the results of the paper and the code used to analyze the results (MB & MB et al. JCB 2023)

## Basics

[Surface Evolver](http://facstaff.susqu.edu/brakke/evolver/evolver.html) is used to find shapes of the bile canalicus. You'll need it to run the code.

Julia is used to process output from Surface Evolver and plot results.

Mathematica is used for generating an initial geometry and performing parameter inference.
