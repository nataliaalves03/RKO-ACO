# RKO-ACO for QMC-VSBPP

[![Paper](https://img.shields.io/badge/Paper-Computers_%26_Operations_Research-blue)](https://www.sciencedirect.com/science/article/pii/S0305054826001000?via%3Dihub)

This repository contains the official implementation and data for the paper **"Random-key optimizer and linearization for the quadratic multiple constraints variable-sized bin packing problem"**, accepted in *Computers & Operations Research* (2026).

This work introduces **RKO-ACO**, a continuous-domain Ant Colony Optimization algorithm integrated into the Random-Key Optimizer framework, and a linearized mathematical model to solve the challenging QMC-VSBPP.

## 📌 Citation

If you use this code or the mathematical models in your research, please cite our paper:

```bibtex
@article{SANTOS2026107482,
  title = {Random-key optimizer and linearization for the quadratic multiple constraints variable-sized bin packing problem},
  journal = {Computers & Operations Research},
  volume = {192},
  pages = {107482},
  year = {2026},
  issn = {0305-0548},
  doi = {https://doi.org/10.1016/j.cor.2026.107482},
  url = {https://www.sciencedirect.com/science/article/pii/S0305054826001000},
  author = {Natalia Alves Santos and Marlon Jeske and Antônio Augusto Chaves},
  keywords = {Bin packing problem, Random-key optimizer, Ant colony optimization, Metaheuristics, Mathematical programming}
}
```

## 🏗️ Framework Lineage & Acknowledgements

This repository is an increment built upon the **Random-Key Optimizer (RKO) Framework v2.0**. 
* **Original Framework Repository:** [RKO-solver/RKO_Cpp_v2.0](https://github.com/RKO-solver/RKO_Cpp_v2.0)
* **Base Framework Reference:** Chaves, A.A., et al. *A Random-Key Optimizer for Combinatorial Optimization*. ([paper](https://link.springer.com/article/10.1007/s10732-025-09568-z))

While the base RKO framework provides a robust architecture for combinatorial optimization, this repository specifically provides the custom **Continuous Ant Colony Optimization (ACO)** implementation, the adaptive Q-learning parameter control tailored for this context, and the problem-specific decoders for the **QMC-VSBPP**.

## 🚀 Scope and Features

This code extends the RKO to solve the **Quadratic Multiple Constraints Variable-Sized Bin Packing Problem (QMC-VSBPP)**, featuring:
* Multiple capacity dimensions and heterogeneous bin types.
* A linearized mathematical model enabling exact solvers (like Gurobi) to compute strong lower bounds.
* **RKO-ACO**: Adaptive continuous-domain ACO enhanced with Q-learning parameter control and efficient local search.

## 💻 Running the Algorithm

The algorithm is written in C++ (C++20) and requires the OpenMP paradigm to be enabled for parallel execution.

1. **Enter the Program directory:** `cd Program`
2. **Compile the code:** `make rebuild`
   *Alternatively, compile via terminal:*
   `g++ -std=c++20 -o runTest main.cpp -O3 -fopenmp`
3. **Run the RKO-ACO:** `./runTest ../Instances/QMC-VSBPP/instance_name.txt T` 
   *(Where `T` is the maximum running time in seconds).*
   *In Windows:* `runTest.exe ../Instances/QMC-VSBPP/instance_name.txt T`

*Note: You must ensure the folders `Instances/QMC-VSBPP` (containing the benchmark instances) and `Results` (where output files are written) exist in the root directory.*

## 📂 Code Structure

* **SPECIFIC_CODE:**
    * **Problem_QMC.h**: Contains the data structure for the QMC-VSBPP, the data reading functions, and the specific decoder.
* **GENERAL_CODE:**
    * **/MH**: Contains the metaheuristic mechanisms, including the newly introduced ACO.
    * **/Main**: Contains the main function and shared variables.
    * **Data.h / Output.h**: Data structures and output handlers for statistical analysis.

## ⚙️ Configuration (`config_tests.conf`)

The `config_tests.conf` file controls the execution parameters. Each metaheuristic runs in a separate thread. Available methods include SA, ILS, VNS, BRKGA, BRKGA-CS, PSO, GA, LNS, GRASP, IPR, and **ACO**.

**Key configuration flags:**
* `MAXRUNS`: Maximum number of runs.
* `debug`: Execution mode (0 for test, 1 for debug).
* `control`: Parameter configuration mode (0 for offline tuning via `ParametersOffline.txt`, 1 for online Q-Learning configuration via `ParametersOnline.txt`).
* `strategy`: Local search strategy (1 for first improvement, 2 for best improvement).
* `restart`: Percentage of max time at which restart triggers.
* `sizePool`: Size of the elite solution pool.
