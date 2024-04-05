# Instances-and-Formulations-Parallel-Machine-Scheduling

## Parallel Machine Scheduling Instances
This repository contains 12 instances based on real-world data.

## Parallel Machine Scheduling Formulations
We also provide the formulations for the parallel machine scheduling problem tested in the paper:

> *Carrilho, L. M., Oliveira, F., & Hamacher, S. (2024).*
> A novel exact formulation for parallel machine scheduling problems.
> Computers & Chemical Engineering.
> **DOI:** [10.1016/j.compchemeng.2024.108649](https://doi.org/10.1016/j.compchemeng.2024.108649)

## Shared files
- **instance**: Instance struct, containing parameters of the instance problem.
- **param**: Parameter struct, defining how to execute the model.
- **variable_condition**: Variable condition struct, defining how the general condition of variables in the model was established.
- **create_variable_condition**: Functions to create the variable condition, dependent on models.
- **deterministic_models**: Functions to construct deterministic models.
- **auxiliary_functions**: Auxiliary functions.
- **main**: Script to run model optimization.

Feel free to explore the data and formulations provided here!
