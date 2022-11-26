#include "common.h"

vector<int> _2opt(vector<int> ham_path, vector<vector<float>> adj_mat, int* MAX_ITER, float *scatter_value, int randomize=1)
{// returns an approximate max scatter based on the 2-opt algorithm
	int n=adj_mat.size(), min_index, found=1, p, q, r, s;  //'p'=V_i, 'q'=V_{i+1}, 'r'=V_j, 's'=V_{j+1} => 'pq', 'rs' are edges of current ham path 
	if (randomize)
		ham_path=rand_ham_path(n);  //randomly initialize an hamiltonian path
	float min_dist;
	*MAX_ITER=0;  //can be used to limit the below infinite while
	while(true)  //loop until no further '_2opt_step' can be done on current ham path
	{
		*MAX_ITER+=1;
		found=0;
		min_dist=scatter(adj_mat,ham_path,n,&min_index);  //min_index = (small) vertex id of current scatter edge
		for(int i=0;i<n;i++)
		{
			p=ham_path[i];
			q=ham_path[(i+1)%n];
			r=ham_path[min_index];
			s=ham_path[(min_index+1)%n];
			if (adj_mat[p][r]>min_dist && adj_mat[q][s]>min_dist)  //a 2-opt step can be performed to maximize current scatter value
			{
				ham_path=_2opt_step(ham_path,i,(i+1)%n,min_index,(min_index+1)%n);  //perform a 2opt step => 'pr', 'qs' become new edges
				found=1;
				break;
			}
		}
		if (found==1)
			continue;
		else  //no further '_2opt_step' on current ham path
		{
			*scatter_value=scatter(adj_mat,ham_path,n,&min_index);
			return (ham_path);
		}
	}
	*scatter_value=scatter(adj_mat,ham_path,n,&min_index);
	return (ham_path);
}