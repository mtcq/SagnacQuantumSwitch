## Code to accompany: [Demonstration of a quantum SWITCH in a Sagnac configuration](https://arxiv.org/abs/xxxx.xxxx)

#### Teodor Strömberg, Peter Schiansky, Robert W. Peterson, Marco Túlio Quintino, and Philip Walther


This is a repository for the article [Demonstration of a quantum SWITCH in a Sagnac configuration](https://arxiv.org/abs/xxxx.xxxx)

 This code requires:
- [cvx](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.

This repository consists of:

- [MaxWorstCase.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/MaxWorstCase.m):
Script finds the maximal worst case success probability for Sequential.

- [MaxSuccessProb_com_anticom_task.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/MaxSuccessProb_com_anticom_task.m):
Script finds the maximal success probabilities for uniform average channel discrimination for Parallel, Sequential, and general strategies with standard SDP methods.

- [MaxSuccessProb_com_anticom_task_symbolic.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/MaxSuccessProb_com_anticom_task_symbolic.m):
Script finds and verifies the upper bounds for sequential strategy for uniform average channel discrimination in a rigorous computer assisted proof way.

- [anglecalc.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/anglecalc.py):
Script that finds the waveplate angles to implement the SU(2) unitary U using the reciprocal polarization gadget 

- [DefineSets_com_anticom.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/DefineSets_com_anticom.m):
Script defines the sets considered in the discrimination task.

- [DefineSets_com_anticom_symbolic.m](https://github.com/mtcq/SagnacQuantumSwitch/blob/main/DefineSets_com_anticom_symbolic.m):
Script defines the sets considered in the discrimination task as symbolic variables

- [auxiliary_functions](https://github.com/mtcq/SagnacQuantumSwitch/tree/main/auxiliary_functions):
Folder with auxiliary functions.
