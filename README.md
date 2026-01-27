This repository is created for the WCCI2026 competition  of dynamic multi-objective optimisation, There are two tracks of competition, briefly described as follows:  
Track 1: Dynamic Unconstrained Multi-Objective Optimisation  
Track 2: Dynamic Constrained Multi-Objective Optimisation  

All the benchmark functions have been implemented in MATLAB code^[1]^, your competition results can be submitted as a brief technical report.   
Please send your results directly to Dr Xiaozhong Yu (xzyu@ smail.xtu.edu.cn)  

More details can be found in "https://zoujuan1.github.io/#page-top"



# Search space

| Problems | $x_i=1$ | $x_i=2,..n$ |
| -------- | ------- | ----------- |
| RDP1     | [0,1]   | [-1,1]      |
| RDP2     | [0,1]   | [-1,1]      |
| RDP3     | [1,4]   | [0,1]       |
| RDP4     | [0,1]   | [-1,1]      |
| RDP5     | [0,1]   | [-1,1]      |
| RDP6     | [0,1]   | [-1,1]      |
| RDP7     | [0,1]   | [-1,1]      |
| RDP8     | [0,1]   | [-1,1]      |
| RDP9     | [0,1]   | [-1,1]      |
| RDP10    | [0,1]   | [-1,1]      |
| RDC1     | [0,1]   | [0,1]       |
| RDC2     | [0,1]   | [-1,1]      |
| RDC3     | [0,1]   | [0,1]       |
| RDC4     | [0,1]   | [0,1]       |
| RDC5     | [0,1]   | [-1,1]      |
| RDC6     | [0,1]   | [-1,1]      |
| RDC7     | [0,1]   | [-1,1]      |
| RDC8     | [0,1]   | [-1,1]      |
| RDC9     | [0,1]   | [0,1]       |



# Parameter Setting

- Population Size:100
- Number of variables: 10
- Frequency of change $\tau_t$: 10 (fast changing environments), 20 (slow changing environments).
- Severity of change $n_t$: 5 (severe changing environments), 10 (moderate changing environments).
- Number of changes: 30.
- Number of independent runs: 20
- Stopping criterion: a maximum number of 100(30$\tau_t$+50) fitness evaluations, where 500 fitness evaluations are given before the first environmental change occurs.
- Metrics: MIGD、MHV^[2]^

# Result Submission

It is expected that competition results can be submitted in tables in a format exemplified in Table 1. However, other ways of result presentation are also acceptable. Please do make sure your result is of high readability for submission, and multiple types of results shown in Table are clearly recorded, including the mean and standard deviation of the MIGD/MHV values of each test instance.

| Problem | $(\tau_t,n_t)$                        | MIGD(mean(std.))   | MHV(mean(std.))    |
| ------- | ------------------------------------- | ------------------ | ------------------ |
| RDP1    | 10,5 <br> 10,10 <br> 20,5  <br> 20,10 | 1.234E-2(1.234E-3) | 1.234E-2(1.234E-3) |
| RDP2    |                                       |                    |                    |
| ……      |                                       |                    |                    |
| RDP10   | 10,5 <br> 10,10 <br> 20,5 <br> 20,10  |                    |                    |

# Reference

[1] Zou J, Hou Z, Jiang S, et al. Knowledge Transfer With Mixture Model in Dynamic Multi-Objective Optimization[J]. IEEE Transactions on Evolutionary Computation, 2025.

[2] Yu X, Zheng J, Hu Y, et al. An adaptive response algorithm based on dual-space detection for dynamic multiobjective optimization[J]. Swarm and Evolutionary Computation, 2025, 98: 102092.
