### Code
There are 5 code files in `code_files` folder:
- `main.cpp`: Main file to run
- `naive_greedy.cpp`: Contains code for naive greedy algorithm
- `weave.cpp`: Contains code for naive and Hoffmann weave algorithms
- `dirac.cpp`: Contains code for Dirac algorithm
- `_2opt.cpp`: Contains code for pure 2-opt and randomized 2-opt algorithms

**To run the code**
- Place the `TSP_data_files` folder inside the `code_files` folder
- cd to the `code_files` folder
- Run: `g++ main.cpp dirac.cpp _2opt.cpp weave.cpp naive_greedy.cpp -o main.out`. Use '-std=c++11' compiler option if required. 

### Data
The `TSP_data_files` folder provides graph input files used for experiments performed in the paper and are derived from the TSPLIB library.
