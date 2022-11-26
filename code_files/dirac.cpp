#include "common.h"

vector<int> dirac(vector<vector<float> > adj_mat, float *scatter_value, float *scatter_bound)  //(adj_mat: square adjacency matrix with non-negative entries)
{  //returns an approximate max scatter based on Arkin's maxmin_1 heuristic and Dirac's theorem for hamiltonian paths 
    int n=adj_mat.size(), S_i, p, q, min_node;
    float min=INFINITY, median, temp_dist;
    vector<int> ham_path=rand_ham_path(n);  //random hamiltonian path initialized; indexed from 0->(n-1)
    vector<float> to_sort(n);
    vector<vector<int> > E_(n);    
    if(n<2)  //not a valid graph
    {
        *scatter_value=*scatter_bound=0;
        return (ham_path);
    }
    for(int i=0;i<n;i++)  //to validate 'adj_mat' & find 'scatter_bound', 'min' (minimum of medians / threshold)
    {
        if(adj_mat[i].size()==n)
        {
            to_sort=adj_mat[i];
            sort(to_sort.begin(),to_sort.end());
            median=to_sort[ceil((n+1)/2)];
            if(min>median)
                min=median;
            if(*scatter_bound>to_sort[n-2])
                *scatter_bound=to_sort[n-2];
        }
        else  //in-valid adjacency matrix 
        {
            cout<<"unexpected non-square adjacency-matrix, exiting!";
            exit(0);
        }
    }
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            if(adj_mat[i][j]>=min)
                E_[i].push_back(j);  //storing all vertices to which a vertex has an edge length not less than threshold 
    for(int i=0;i<n;i++)
    {
        p=ham_path[i],q=ham_path[(i+1)%n];
        temp_dist=adj_mat[p][q];  //an edge of current path
        if(temp_dist<min)  //current edge being less than threshold
        {
            vector<int> S_p(E_[p].size()), S_q(E_[q].size()), S_pq(n);
            for(int j=0;j<E_[q].size();j++)  //find all vertices {'s'} such that edge 'qs' is not less than threshold
            {
                S_q[j]=get_index(ham_path,E_[q][j])-1;
                if(S_q[j]==-1)
                    S_q[j]+=n;
            }
            for(int j=0;j<E_[p].size();j++)  //find all vertices {'r'} such that edge 'pr' is not less than threshold
                S_p[j]=get_index(ham_path,E_[p][j]);
            sort(S_p.begin(),S_p.end());sort(S_q.begin(),S_q.end());  //sort sets to reduce the complexity of finding intersection
            set_intersection(S_p.begin(),S_p.end(),S_q.begin(),S_q.end(),S_pq.begin());  //find one set of adjacent vertices: ('r','s')
            S_i=S_pq[0];  //certainly will yield a non-empty intersection value as a result of Dirac's theorem
            ham_path=_2opt_step(ham_path,i,(i+1)%n,S_i,(S_i+1)%n);  //perform a 2opt step => 'pr', 'qs' become new edges
        }
    }
    *scatter_value=scatter(adj_mat,ham_path,n,&min_node);
    return (ham_path);
}