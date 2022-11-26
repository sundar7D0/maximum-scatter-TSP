#include "common.h"

float naive_greedy(vector<vector<float> > &adj_mat) // Computes max scatter using naive greedy approach
// Logic: Start from a random node, visit the farthest unvisited node at each iteration until all nodes are visited.
{
    int dim, tmp_dim;
    dim = tmp_dim = adj_mat.size();

    int visited_nodes[dim]; // vector of all nodes (0-not visited, 1-visited)
    while(tmp_dim--)
        visited_nodes[tmp_dim]=0;  // initially all nodes are unvisited
    
    int curr_node = rand()%dim; // choose a random node to start
    int first_node = curr_node;
    visited_nodes[curr_node] = 1;

    float max_dist = 0;
    int tmp_node = 0;
    float max_scatter = INFINITY; // initialize max_scatter to a large value

    for(int i=0; i<dim-1; i++)
    {
        max_dist = 0;
        int j=0;
        for(j=0; j<dim; j++)
        {
            if(visited_nodes[j]==0 && adj_mat[curr_node][j]>=max_dist) 
            // From the current node, look for the farthest unvisited node and make an edge to visit it
            {
                tmp_node = j;
                max_dist = adj_mat[curr_node][j];
            }
        }
        curr_node = tmp_node;
        visited_nodes[curr_node] = 1; // Set the next node to visited
        if(max_dist<max_scatter)
            max_scatter = max_dist; // Update max scatter if a smaller edge is found
    }

    max_scatter = min(max_scatter, adj_mat[curr_node][first_node]); // Update max scatter if the last edge (returning to 1st node) is smallest
    return max_scatter;
}
