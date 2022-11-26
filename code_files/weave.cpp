#include<iostream>
using namespace std;
#include <unordered_map>
#include<bits/stdc++.h>

float euclidean_dist(float, float, float, float);

float euclidean_pair_dist(pair<float, float> p1, pair<float, float> p2) // Euclidean distance of a pair
{
//    return sqrt(pow(p2.first-p1.first, 2)+pow(p2.second-p1.second, 2));  //euclidean dist
    return abs(p2.first-p1.first)+abs(p2.second-p1.second);  //manhattan dist
//  return ceil(sqrt(pow(p2.first-p1.first, 2)+pow(p2.second-p1.second, 2)));  //ceil-2d dist

}

void transpose(vector<vector<int> > &b) // Transpose a matrix b
{
    if (b.size() == 0)
        return;

    vector<vector<int> > trans_vec(b[0].size(), vector<int>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    b = trans_vec; 
}

// We first enforce a fixed number of points per column and make a grid. Later, we transpose the grid,
// so that rows become columns and vice-versa

vector<vector<int> > make_grid(vector<vector<float> > nodes, float fact)
{
    int dim = nodes.size();

    // Convert the dimesion to a perfect square
    int num_rows, num_cols, max_pts_per_col;
    max_pts_per_col = floor(sqrt(dim))*fact;

    // ****** Remove the last few points to ensure each column has equal no. of points ******
    dim = dim - dim % max_pts_per_col;

    vector<pair<float, int> > col_idx;

    for(int i=0; i<dim; i++)
    {
        col_idx.push_back(make_pair(nodes[i][1], i));
    }

    sort(col_idx.begin(), col_idx.end());

    vector<pair<float, pair<int, int> > > row_idx;

    for (int i=0; i<dim; i++)
    {
        row_idx.push_back(make_pair(nodes[col_idx[i].second][0], make_pair(col_idx[i].second, floor(i/max_pts_per_col))));
    }

    sort(row_idx.begin(), row_idx.end());
    num_cols = 1 + floor((dim-1)/max_pts_per_col);

    vector<vector<int> > grid_positions;
    vector<int> curr_row(num_cols);
    fill(curr_row.begin(), curr_row.end(), -1);

    for(int i=0; i<dim; i++)
    {
        int curr_idx = row_idx[i].second.second;
        if(curr_row[curr_idx]==-1)
        {
            curr_row[curr_idx] = row_idx[i].second.first;
        }
        else
        {
            grid_positions.push_back(curr_row); //grid_positions stores a copy of curr_row
            fill(curr_row.begin(), curr_row.end(), -1);
            curr_row[curr_idx] = row_idx[i].second.first;
        }
    }

    grid_positions.push_back(curr_row);
    num_rows = grid_positions.size();
    return grid_positions;
}


vector<int> arkin_1D_seq(int num_cols) // To generate the 1D permuations that gives max scatter by Arkin's algorithm
{
    int k = floor(num_cols/2);
    vector<int> idx_sequence(num_cols+1, 0);

    if(num_cols%2==1) // If no. of points in 1D sequence is odd
    {
        for(int i=0; i<=2*k-2; i+=2)
        {
            idx_sequence[i] = int(i/2);
            idx_sequence[i+1] = k+int(i/2)+1;
        }
        idx_sequence[2*k] = k;
    }

    else // If no. of points in 1D sequence is odd
    {
        if(k%2==0) 
        {
            idx_sequence[0] = 0;
            int i=1;
            while(i<=k-3)
            {
                idx_sequence[i] = idx_sequence[i-1]+(k+1);
                i++;
                idx_sequence[i] = idx_sequence[i-1]-(k-1);
                i++;
            }
            idx_sequence[i++] = num_cols-1;
            idx_sequence[i++] = k-1;

            while(i<num_cols)
            {
                idx_sequence[i] = idx_sequence[i-1]+(k-1);
                i++;
                idx_sequence[i] = idx_sequence[i-1]-(k+1);
                i++;
            }
        }

        if(k%2==1)
        {
            idx_sequence[0] = 0;
            int i=1;
            while(i<=k-2)
            {
                idx_sequence[i] = idx_sequence[i-1]+(k+1);
                i++;
                idx_sequence[i] = idx_sequence[i-1]-(k-1);
                i++;
            }

            idx_sequence[i++] = num_cols-1;
            while(i<num_cols)
            {
                idx_sequence[i] = idx_sequence[i-1]-(k+1);
                i++;
                idx_sequence[i] = idx_sequence[i-1]+(k-1);
                i++;
            }
        }
    }
    idx_sequence[idx_sequence.size()-1] = 0;
    return idx_sequence;
}


float find_path_dist(vector<pair<float, float> > grid_row) // Find max scatter value for a given 1D sequence using Arkin's
{
    int num_cols = grid_row.size();

    // if(num_cols<=3)
    //     cout<<"Encountered a too small row (even)!! "<<num_cols<<endl;

    if(num_cols<3 && num_cols%2==1) // When a row has just 1 point, Arkin will return 0 which must be ignored
    {
        // cout<<"Encountered a too small row (odd)!! "<<num_cols<<endl;
        return pow(2,16)-1;
    }    
    
    vector<int> idx_sequence = arkin_1D_seq(num_cols);

    float max_scatter = pow(2, 16)-1;
    float edge_dist;

    for(int i=0; i<num_cols-1; i++) 
    // Do not measure the length of the edge connecting the last point to the first point, hence num_cols-1
    {
        edge_dist = euclidean_pair_dist(grid_row[idx_sequence[i]], grid_row[idx_sequence[i+1]]);
        max_scatter = min(edge_dist, max_scatter);
    }
    return max_scatter;
}

float naive_weave(vector<vector<float> > nodes, float fact)
// Logic: Use Arkin's 1D algorithm to visit all points within a row. Then, move to the next row.
// Within a row, only travel between grids that have a point in them.
{
    vector<vector<int> > grid_positions = make_grid(nodes, fact);
    transpose(grid_positions);
    int num_rows = grid_positions.size();
    int num_cols = grid_positions[0].size();
    vector<float> row_max_scatter;
    vector<pair<float, float> > pts_connecting_rows;

    for(int i=0; i<num_rows; i++)
    {
        vector<pair<float, float> > grid_row;
        for(int j=0; j<num_cols; j++)
        {
            if(grid_positions[i][j]!=-1)
                grid_row.push_back(make_pair(nodes[grid_positions[i][j]][0], nodes[grid_positions[i][j]][1]));
        }

        row_max_scatter.push_back(find_path_dist(grid_row)); // All points within a row - each point as a coord pair
        pts_connecting_rows.push_back(grid_row[0]); // 1st visited point (by Arkin's order) within a row
        pts_connecting_rows.push_back(grid_row[floor(grid_row.size()/2)]); // Last visited point (by Arkin's order) within a row
    }

    float max_scatter = *min_element(row_max_scatter.begin(), row_max_scatter.end());
    float tmp_max_scatter = max_scatter;

    // Distances in travelling between the rows
    float edge_dist;
    for(int i=1; i<2*num_rows-1; i+=2)
    {
        edge_dist = euclidean_pair_dist(pts_connecting_rows[i], pts_connecting_rows[i+1]);
        max_scatter = min(max_scatter, edge_dist);
    }
    // To join the 1st and last visited points of the whole 2D tour
    edge_dist = euclidean_pair_dist(pts_connecting_rows[0], pts_connecting_rows[2*num_rows-1]);
    
    max_scatter = min(max_scatter, edge_dist);

    return max_scatter;
}

float hoffmann_pair_dist(vector<vector<pair<float, float> > > grid_valid_pts, int idx1, int idx2, int idx3, vector<int> idx_sequence)
// Utility function for Hoffmann weave. To find min edge cost within a path along a Hoffmann pair/triplet
// Set idx3 to -1 for Hoffmann pair and to the 3rd row index for triplet
{
    float row_min_dist = pow(2, 16)-1;
    int a, b;
    float edge_dist;
    for(int i=0; i<grid_valid_pts[0].size()-1; i++)
    {        
        // For the case of i%2==0 or i%3==0:
        a = idx1;
        b = idx2;

        if(idx3==-1 && i%2==1)
        {
            a = idx2;
            b = idx1;
        }

        else if(idx3!=-1)
        {
            if(i%3==1)
            {
                a = idx2;
                b = idx3;
            }
            else if(i%3==2)
            {
                a = idx3;
                b = idx1;
            }
        }

        edge_dist = euclidean_pair_dist(grid_valid_pts[a][idx_sequence[i]], grid_valid_pts[b][idx_sequence[i+1]]);
        row_min_dist = min(row_min_dist, edge_dist);
    }
    return row_min_dist;
}

float hoffmann_weave(vector<vector<float> > nodes, float fact)
// Logic: Perform 
{
    vector<vector<int> > grid_positions = make_grid(nodes, fact);
    // Transpose rows and columns since Hoffmann weave gives better results
    // when Arkin 1D is run on a fixed (large) set of points
    transpose(grid_positions); 
    int num_rows = grid_positions.size();
    int num_cols = grid_positions[0].size();

    int num_invalid_pts = 0; // Find grid locations where no point is present
    for(int i=0;i<num_cols; i++)
    {
        if(grid_positions[0][i]==-1)
            num_invalid_pts+=1;
    }

    int num_pts_per_row = num_cols - num_invalid_pts; // Number of points per row (fixed due to transposition, but in different subset of columns)
    vector<int> idx_sequence = arkin_1D_seq(num_pts_per_row); // Get desired ordering of points for Arkin 1D
    rotate(idx_sequence.begin(), idx_sequence.begin()+1, idx_sequence.end()); // Column ordering for Hoffmann requires left rotation by 1 step

    vector<float> row_max_scatter;
    vector<pair<float, float> > pts_connecting_rows;
    vector<vector<pair<float, float> > > grid_valid_pts;
    vector<vector<pair<float, float> > > grid_aux_valid_pts;

    float max_scatter = pow(2, 16)-1;
    int orig_num_rows_odd = 0;
    int k = floor(num_pts_per_row/2);
    int t = floor(num_rows/2);

    for(int i=0; i<num_rows; i++)
    {
        vector<pair<float, float> > grid_row;
        
        for(int j=0; j<num_cols; j++)
        {
            if(grid_positions[i][j]!=-1) // Add coordinates of valid points in each row to a vector
            {
                grid_row.push_back(make_pair(nodes[grid_positions[i][j]][0], nodes[grid_positions[i][j]][1]));
            }
        }
        grid_valid_pts.push_back(grid_row);
        // Collect all points for the case of odd number of rows except those that form the triplet
        if(num_rows%2==1 && i!=0 && i!=t && i!=num_rows-1)
            grid_aux_valid_pts.push_back(grid_row); 
    }

    for(int i=0; i<grid_valid_pts.size(); i++) // Sanity check to ensure all rows have equal no. of points
    {
        if(grid_valid_pts[i].size()!=grid_valid_pts[0].size())
            cout<<i<<" "<<grid_valid_pts[i].size()<<" "<<grid_valid_pts[0].size()<<endl;
    }

    if(num_rows%2==1)
    {
        orig_num_rows_odd = 1;
        
        int ordering[6] = {0, t, num_rows-1, 0, t, num_rows-1};
        
        // Traversing triplets 0, floor(num_rows/2), num_rows-1

        // In each traversal, 1/3rd of the points from each row of the triplet are covered
        max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, 0, t, num_rows-1, idx_sequence));
        // Connect the last point from the previous traversal to the 1st point in the next traversal
        max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[ordering[2+num_pts_per_row%3]][0], grid_valid_pts[num_rows-1][k+1]));
        max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, num_rows-1, 0, t, idx_sequence));
        max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[ordering[1+num_pts_per_row%3]][0], grid_valid_pts[t][k+1]));
        max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, t, num_rows-1, 0, idx_sequence));
        
        if(num_rows>3) 
        {
            // Edge to the 1st point visited by Hoffmann weave among the remaining N-3 rows
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[ordering[0+num_pts_per_row%3]][0], grid_valid_pts[1][k+1]));
            // Last edge of the entire tour including all points
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[num_rows-2][0], grid_valid_pts[0][k+1]));
        }
        else // If there are only 3 rows in total, make the return path and return the max scatter
        {
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[ordering[0+num_pts_per_row%3]][0], grid_valid_pts[0][k+1]));
            return max_scatter;
        }

        grid_valid_pts = grid_aux_valid_pts; // All points except those in the triplet
        
        // Reset parameters
        num_rows = num_rows-3;
        t = floor(num_rows/2);
        if(num_rows!=grid_valid_pts.size())
            cout<<"Issue here 1 "<<num_rows<<" "<<grid_valid_pts.size()<<endl;
    }

    // *** The number of rows is now even in all cases (since triplet was removed for odd case) ***

    if(num_pts_per_row%2==1)
    // When number of pts per row is odd, the row idx of starting and ending pts of 1 call to hoffmann_pair_dist are the same.
    {
        for(int i=0; i<t; i++)
        {
            
            // (0, k+1) -> (t, 2) -> .... -> (0,0)
            max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, i, i+t, -1, idx_sequence));
            // Connect the last point to the 1st point of the next path
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[i][0], grid_valid_pts[i+t][k+1]));
            // (t, k+1) -> (0, 2) -> .... -> (t,0)
            max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, i+t, i, -1, idx_sequence));
        }
        if(!orig_num_rows_odd) // For even number of total rows, find cost of last returning edge of entire tour
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[num_rows-1][0], grid_valid_pts[0][k+1]));
    }

    else
    // When number of pts per row is even, the row idx of starting and ending pts of 1 call to hoffmann_pair_dist are opposite.
    {
        // Visit rows in pairs of (0,t,0,t,..), (1,t+1,1,t+1,...), ...
        for(int i=0; i<t; i++)
        {
            max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, i, i+t, -1, idx_sequence));
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[i+t][0], grid_valid_pts[i+1][k+1]));
        }
        max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[t-1][0], grid_valid_pts[t][k+1]));
        
        // Visit unvisited points in pairs of (t,0,t,0,...), (t+1,1,t+1,1,...), ...
        for(int i=t; i<num_rows; i++)
        {
            max_scatter = min(max_scatter, hoffmann_pair_dist(grid_valid_pts, i, i-t, -1, idx_sequence));
            if(i<num_rows-1)
                max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[i-t][0], grid_valid_pts[i+1][k+1]));
        }
        if(!orig_num_rows_odd) // For even number of total rows, find cost of last returning edge of entire tour
            max_scatter = min(max_scatter, euclidean_pair_dist(grid_valid_pts[num_rows-1][0], grid_valid_pts[0][k+1]));
    }

    return max_scatter;
}
