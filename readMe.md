
# Deep learning linking mechanistic models to single-cell transcriptomics data reveals transcriptional bursting in response to DNA damage
![image](https://github.com/cellfateTX/DeepTX/blob/master/trainInferSSA/logo/deepTXlogo.jpg)
Code and data for paper  "Deep learning linking mechanistic models to single-cell transcriptomics data reveals transcriptional bursting in response to DNA damage".

## System Requirements

The Julia code used for training and inference of the model comes with versions.

    Julia Version 1.7.1
    Commit ac5cc99908 (2021-12-22 19:35 UTC)
    Platform Info:
    OS: Windows (x86_64-w64-mingw32)
    CPU: AMD Ryzen 5 4600H with Radeon Graphics
    WORD_SIZE: 64
    LIBM: libopenlibm
    LLVM: libLLVM-12.0.1 (ORCJIT, znver2)

The R code of data analysis and visualization for investigating the burst kinetic in the plotsAndAnalysis folder are all R scripts with versions.

    R version 4.2.1 (2022-06-23 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)

The matlab version information of getStatsData.m script in TraingInferSSA folder.

    MATLAB Version: 9.10.0.1649659 (R2021a) Update 1
    Operating System: Microsoft Windows 10 Professional Version 10.0 (Build 19044)
    Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

## Directory
### TraingInferSSA
Employ the Stochastic Simulation Algorithm (SSA) algorithm to generate data and train a neural network capable of solving the mechanism model. Evaluate the performance of both solving and inference of the neural network model on the test set.
#### Workflow
* Generate data for the TX model, comprising model parameters and their corresponding stationary distribution. (simulation.jl)
* Compute statistics, including mean, variance, and bimodal coefficient, relevant to the model parameters. (getStatsData.m)
* Partition the dataset into training, validation, and test sets for neural network training. (generateNNset.jl)
* Train the neural network model. (trainSolver.jl)
* Validate statistics and distribution errors between SSA results and neural network predictions on the test set. (validateStatsDist.jl)
* Employ the neural network model to infer dynamic parameters from the test set data. (inferTestSet.jl)
* Compute the posterior distribution for the dynamic parameters of the test set. (posteriorDistribution.jl)
* Conduct stability analysis of model inference. (inferenceRobust.jl)

### InferObservedData
The trained neural network is utilized to infer potential kinetic model parameters from single-cell sequencing data. The analysis pipeline in the three folders are analogous. InferIdUDMSO serves as an illustrative example to elucidate the process of code execution involves the systematic processing and inference of observed data.
#### Workflow
* Perform preprocessing of single-cell sequencing (scRNA-seq) data. (data_preprocess.R)
* Utilize the trained neural network to infer scRNA-seq data. (TX_inferrer.jl)
* Utilize inferred parameters to perform SSA simulation and obtain the true distribution of inferred parameters. (generate_distribution.jl)
* Compute the KL divergence between the inferred distribution and the distribution of scRNA-seq to assess the accuracy of model inference. (filter_gene_by_kl.R)
### PlotsAndAnalysis
The folder contains analyses including differential analysis and gene enrichment analysis of the burst kinetics inferred from scRNA-seq data, along with visualization. The names of the code files align with the serial numbers of the figure in the manuscript.
