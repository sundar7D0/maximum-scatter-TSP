#ifndef COMMON_H
#define COMMON_H

#include <math.h>
#include <random>
#include <dirent.h>	
#include <iostream>
#include <iterator>
#include <sys/types.h>
#include <bits/stdc++.h>
#include <unordered_map>
#include <unordered_set>
	
using namespace std;

void print_path(vector<int>);
vector<int> rand_ham_path(int);
int get_index(vector<int>, int);
vector<int> dirac(vector<vector<float>>, float*, float*);
float naive_greedy(vector<vector<float>> &a);
vector<int> _2opt_step(vector<int>, int, int, int, int);
float scatter(vector<vector<float>>, vector<int>, int, int*);
float naive_weave(vector<vector<float>>, float);
float hoffmann_weave(vector<vector<float>>, float);
vector<int> _2opt(vector<int>, vector<vector<float>>, int*, float*, int);

#endif