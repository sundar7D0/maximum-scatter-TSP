# maximum-scatter-TSP
Six heuristics and their smoothed analysis for NP-hard Maximum Scatter Traveling Salesman Problem (MSTSP)

In this work, we describe six algorithms for MSTSP with improved formulations of prior work that enhance their real-world efficacy. Further, we perform experimental studies motivated by smoothed analysis to comprehensively evaluate these algorithms in terms of run-time, quality and stability. This work, done under the guidance of [Prof. Raghavendra Rao](http://www.cse.iitm.ac.in/profile.php?arg=MTU=), was published as paper (lead author) at [SIAM-ALENEX22](https://www.siam.org/conferences/cm/conference/alenex22).

The six algorithms that we describe in this work are:

1. Naive Greedy Algorithm
2. Naive Weave Algorithm
3. Hoffmann Weave Algorithm
4. Dirac Algorithm
5. Pure 2-opt Algorithm
6. Randomized 2-opt Algorithm

Our benchmarking experiments can be broadly split into three categories:

* Closeness of algorithm predictions to the scatter bound
* Deviation of maximum scatter predictions under perturbation
* Variation in the runtime of algorithms

## Key Contributions:

1. We observe that the Naive Greedy algorithm is very fast and easy-to-implement baseline for MSTSP.
2. We present the Naive Weave and Hoffmann Weave algorithms which introduce an improved formulation of the work by Arkin and Hoffmann to extend their usability to non-regular grids.
3. We introduce Pure 2-opt and Randomised 2-opt as very close approximation algorithms for the MSTSP.
4. We use a real-world dataset augmented using five graph perturbations and evaluated with three edge cost metrics to perform a comprehensive perturbation analysis of the algorithms and compare results on three critical performance measures, namely, the quality, runtime and stability of the algorithms.

## Using the code
A single directory contains code of 6 different algorithms and another data. It can be easily ported to run with public-cloud like GCP, AWS, etc.

## Requirements
* C++11

## Additional resources
1. [Paper](https://epubs.siam.org/doi/abs/10.1137/1.9781611977042.13)

