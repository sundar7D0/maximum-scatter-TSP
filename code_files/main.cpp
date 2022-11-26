#include "common.h"

string strip(string &in)  //strips leading and trailing whitespaces of a string
{
    string::const_iterator b = in.begin(), e = in.end();
    while (isspace(*b))
        ++b;
    while (isspace(*(e-1)) && b!=e)
        --e;
    return string(b,e);
}

float euclidean_dist(float x1, float y1, float x2, float y2)  //calculates euclidean distance
{
    return sqrt(pow(x2-x1, 2)+pow(y2-y1, 2));
}

float manhattan_dist(float x1, float y1, float x2, float y2)  //calculates Manhattan distance
{
    return abs(x2-x1)+abs(y2-y1);
}

float max_2D(float x1, float y1, float x2, float y2)  //calculates max 2D distance
{
    return max(abs(x2-x1), abs(y2-y1));
}

float ceil_2D(float x1, float y1, float x2, float y2)  //calculates ceil 2D distance
{
    return ceil(euclidean_dist(x1, y1, x2, y2));
}

vector<vector<float>> gaussian_noise(vector<vector<float>> nodes, default_random_engine* generator, float mean, float std_dev)  //N(mean,std_dev)
{
    cout<<"gaussian noise with: "<<"(mean, std_dev)=("<<mean<<','<<std_dev<<')'<<endl;
    normal_distribution<double> dist(mean, std_dev);  //std object for gaussian noise generation
    int n=nodes.size();
    for (auto& point : nodes)  //add Gaussian noise to (x,y) coordinates
    {
        point[0] = point[0] + dist(*generator);
        point[1] = point[1] + dist(*generator);
    }
    return nodes;
}

vector<vector<float>> uniform_noise(vector<vector<float>> nodes, default_random_engine* generator, float a, float b)  //U(a,b)
{
    cout<<"uniform noise with: "<<"(a,b)=("<<a<<','<<b<<')'<<endl;
    std::uniform_real_distribution<> dist(a, b);  //std object for uniform noise generation
    int n=nodes.size();
    for (auto& point : nodes)  //add uniform noise to (x,y) coordinates
    {
        point[0] = point[0] + dist(*generator);
        point[1] = point[1] + dist(*generator);
    }
    return nodes;
}

bool validate_file(string filename)
{
    string curr_line; 
    ifstream MyReadFile(filename);  //reads the file
    unordered_map<string, string> file_specs;  //gets the tsplib file specs from the first few lines into a dictionary
    while(true)
    {
        getline (MyReadFile, curr_line);
        int idx = curr_line.find(":");
        if (idx<curr_line.length())
        {   
            string s1 = curr_line.substr(0,idx);
            string s2 = curr_line.substr(idx+1);
            file_specs[strip(s1)] = strip(s2);
        }
        else    
            break;
    }
    if (file_specs["EDGE_WEIGHT_TYPE"]=="EUC_2D" || file_specs["EDGE_WEIGHT_TYPE"]=="ATT" || file_specs["EDGE_WEIGHT_TYPE"]=="GEO")
        return true;
    else
        return false;
}

vector<string> find_valid_files()  //finds tsp data files in the desired data format 
{
    string fname;
    vector<string>  valid_fnames;
    unordered_map<string, string> file_specs;
    DIR* dirp = opendir("TSP_data_files/");
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) 
    {
        fname = dp->d_name;
        if (fname.substr(fname.find('.')) == ".tsp")
        {   
            if(validate_file("TSP_data_files/"+fname))
                valid_fnames.push_back(fname);
        }
    }
    closedir(dirp);
    return valid_fnames;
}

vector<vector<float>> read_data_file(string filename, float *maxima)  //reads nodes from data file 
{
    string curr_line;  //reads the file
    ifstream MyReadFile(filename); 
    unordered_map<string, string> file_specs;  //gets the tsplib file specifics from the first few lines into a dictionary
    while(true)
    {
        getline(MyReadFile, curr_line);
        int idx = curr_line.find(":");
        if (idx<curr_line.length())
        {   
            string s1 = curr_line.substr(0,idx);
            string s2 = curr_line.substr(idx+1);
            file_specs[strip(s1)] = strip(s2);
        }
        else if (curr_line=="NODE_COORD_SECTION")  
            break;
    }
    int dim = stoi(file_specs["DIMENSION"]);  //gets the number of nodes/points in the file
    vector<vector<float> > nodes(dim, vector<float>(2));
    float xmin=INFINITY, xmax=-INFINITY, ymin=INFINITY, ymax=-INFINITY;
    for(int i=0; i<dim; i++)
    {
        getline (MyReadFile, curr_line);
        curr_line = strip(curr_line);
        curr_line = curr_line.substr(curr_line.find(' ')+1);
        curr_line = strip(curr_line);
        nodes[i][0] = stof(curr_line.substr(0, curr_line.find(' ')));  //gets the 1st and 2nd coordinate values
        nodes[i][1] = stof(curr_line.substr(curr_line.find(' ')+1));
        if(nodes[i][0]>xmax)
            xmax=nodes[i][0];
        if(nodes[i][0]<xmin)
            xmin=nodes[i][0];
        if(nodes[i][1]>ymax)
            ymax=nodes[i][1];
        if(nodes[i][1]<ymin)
            ymin=nodes[i][1];
    }
    *maxima=max(abs(ymax-ymin),abs(xmax-xmin));
    MyReadFile.close();  //close the file
    return nodes;
}

vector<int> _2opt_step(vector<int> ham_path, int x1, int x2, int x3, int x4)  //performs a 2 opt exchange (args needs to be in a cyclic order)
{//edges 'path[x1]->path[x2]', 'path[x3]->path[x4]' are modified to 'path[x1]->path[x3]', 'path[x2]->path[x4]' where path denotes hamiltonian path
    int n=ham_path.size();
    int x_1=min(x1,x3), x_3=max(x1,x3);
    int x_2=(x_1+1)%n, x_4=(x_3+1)%n;
    if(x_4==0)
        x_4=n;
    vector<int> new_ham_path(n);  //0->(n-1) directional loop  
    for(int i=0;i<=x_1;i++)
        new_ham_path[i]=ham_path[i];
    for(int i=0;i<=x_3-x_2;i++)
        new_ham_path[x_1+1+i]=ham_path[x_3-i];
    for(int i=0;i<n-x_4;i++)
        new_ham_path[x_4+i]=ham_path[x_4+i];
    return new_ham_path;
}

vector<int> read_opt_tour_file(string filename)  //reads optimal tour
{
    string curr_line;  //reads the file
    ifstream MyReadFile(filename); 
    unordered_map<string, string> file_specs;  //gets the file specifications in the first few lines of the tsplib file into a dictionary
    while(true)
    {
        getline(MyReadFile, curr_line);
        int idx = curr_line.find(":");
        if (idx<curr_line.length())
        {   
            string s1 = curr_line.substr(0,idx);
            string s2 = curr_line.substr(idx+1);
            file_specs[strip(s1)] = strip(s2);
        }
        else if (curr_line=="TOUR_SECTION")  
            break;
    }
    int dim = stoi(file_specs["DIMENSION"]);  //gets the number of nodes/points in the file
    vector<int> path(dim);
    for(int i=0; i<dim; i++)
    {
        getline (MyReadFile, curr_line);
        path[i]=stoi(curr_line)-1;
    }
    MyReadFile.close();  //closes the file
    return (path);
}

vector<vector<float>> make_adjacency_matrix(vector<vector<float>> nodes, float (*function)(float,float,float,float))  
{//creates a square adjacency matrix based on distance function passed in argument
    int dim = nodes.size();
    vector<vector<float> > adj_mat(dim, vector<float>(dim));
    for(int i=0; i<dim; i++)
        for(int j=0; j<dim; j++)
            adj_mat[i][j] = adj_mat[j][i] = function(nodes[i][0], nodes[i][1], nodes[j][0], nodes[j][1]);
    return adj_mat;
}

void print_adj_mat(vector<vector<float> > adj_mat)  //prints square adjacency matrix
{
    cout<<endl<<"adjacency matrix: "<<endl;
    for(int p=0;p<adj_mat.size();p++)
    {
        for(int q=0;q<adj_mat.size();q++)
            cout<<adj_mat[p][q]<<',';
        cout<<endl;
    }       
}

void print_path(vector<int> path)  //prints any path
{
    cout<<endl<<"path: "<<endl;
    for(int p=0;p<path.size();p++)
        cout<<path[p]<<',';
    cout<<endl;
}

vector<int> rand_ham_path(int size)  //creates a random hampath (0->(size-1))! 
{
    srand(time(NULL));
    vector<int> path(size);
    int index;
    for(int i=0;i<size;i++)
        path[i]=i;
    for(int i=0;i<size;i++)
    {
        index=rand()%(size-i);
        swap(path[index],path[size-i-1]);
    }
    return path;
}

int get_index(vector<int> path, int point)  //returns index of a point in a ham-path vector
{
    vector<int>::iterator it=find(path.begin(),path.end(),point);
    if(it!=path.end())
    {
        int index=it-path.begin();
        return (index);
    }
    else
        return (-1);
}

float scatter(vector<vector<float> > adj_mat, vector<int> ham_path, int n, int* min_node)  //finds the scatter (minimum edge) of a path
{
    float min_scatter=INFINITY, temp_dist;
    for(int i=0;i<n;i++)
    {
        temp_dist = adj_mat[ham_path[i]][ham_path[(i+1)%n]];
        if (temp_dist<min_scatter)
        {
            min_scatter = temp_dist;
            *min_node = i;
        }
    }
    return (min_scatter);
}

int main()  //start
{
    vector<string> valid_fnames = find_valid_files();
    vector<int> opt_path;
    float scatter_bound, naive_greedy_result, dirac_result=INFINITY, pure_2opt_randomize_result, pure_2opt_result, naive_weave_result, hoffmann_weave_result;
    vector<vector<float>> adj_mat, new_adj_mat, nodes, new_nodes;
    clock_t start, end;
    int MAX_ITER=1000, counter=100, nodes_size; 
    float maxima, minima, beta=0.5;
    ofstream MyFile("final_results.txt");  //log file
    float mean=0.0, std_dev=1.0;
    float a=0.0, b=1.0, std_dev_e,std_dev_s;
    default_random_engine generator;  //std object for random number generation
    int randomness=1;  //randomness={0: 'no randomness', 1: 'gaussian (mean, std_dev)', 2: 'uniform (a,b)'}
    int emsun=0;  //flag variable to distinguish two different variances used
    for(emsun=0; emsun<=1;emsun+=1)
    {
        for(randomness=1; randomness<=2; randomness+=1)
        {
            if(randomness==0)
                MyFile<<"null"<<endl;
            else if (randomness==1)
                MyFile<<"gaussian"<<endl;
            else if (randomness==2)
                MyFile<<"uniform"<<endl;
            else
            {
                cout<<"unknown randomness! exiting!"<<endl;        
                exit(0);
            }
            for(int i=0; i<valid_fnames.size(); i++) 
            {
                nodes = read_data_file("TSP_data_files/"+valid_fnames[i],&maxima);  //read data-coordinates
                nodes_size=nodes.size();
                if(nodes.size()>2000)
                {
//                  cout<<"Ignored due to large size: "<<valid_fnames[i]<<endl;
                    continue;
                }
                cout<<"file "<<valid_fnames[i]+": "<<endl;  //prints file name
                MyFile << "file: "<<valid_fnames[i]<<endl;
                adj_mat = make_adjacency_matrix(nodes, &manhattan_dist);  //calculate adjacency matrix from coordinates and distance metric
                vector<float> mini=*min_element(adj_mat.begin(),adj_mat.end());
                sort(mini.begin(),mini.end());
                minima=mini[1];
                std_dev_e=beta*minima;
                std_dev_s=0.5*maxima/sqrt(nodes_size*log(nodes_size));
                if(emsun==0)
                {
                    std_dev=std_dev_s;
                    a=-std_dev_s;
                    b=std_dev_s;  //'emsun'=0 => variance = max square width * (1/(2*sqrt(n*logn)))
                }
                else if(emsun==1)
                {
                    std_dev=std_dev_e;
                    a=-std_dev_e;
                    b=std_dev_e;  //'emsun'=1 => variance = beta * min_dist between any pair of nodes
                }                    
                if (randomness==1)
                    MyFile <<"gaussian: "<<std_dev<<endl;
                else if (randomness==2)
                    MyFile <<"uniform: "<<a<<','<<b<<endl;
                for(int k=0;k<counter;k++)  //perturbation iterations over a fixed (mean, sigma)
                {
                    scatter_bound = pow(2,16)-1;  //INFINITE initial scatter bound (minimum among all 2nd max edge of all vertices)
                    if (randomness==1)
                        new_nodes=gaussian_noise(nodes,&generator,mean,std_dev);  //gaussian perturbation
                    else if (randomness==2)
                        new_nodes=uniform_noise(nodes,&generator,a,b);  //uniform perturbation
                    else
                        new_nodes=nodes;
                    if (randomness)
                        new_adj_mat = make_adjacency_matrix(new_nodes, &manhattan_dist);  //calculate adjacency matrix from coordinates and distance metric
                    else
                        new_adj_mat = adj_mat;
                    vector<double> time_taken;
                    start = clock();  //start time to evaluate runtime of 'naive_greedy'
                    naive_greedy_result = naive_greedy(new_adj_mat);
                    end = clock();  //end time to evaluate runtime of 'naive_greedy'
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));  //storing execution time in a vector
                    start=clock();
                    opt_path = dirac(new_adj_mat, &dirac_result, &scatter_bound);
                    end = clock();
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));
                    start=clock();
                    opt_path =_2opt(opt_path, new_adj_mat, &MAX_ITER, &pure_2opt_result,0);  //randomize=0
                    end = clock();
                    cout<<"count1: "<<MAX_ITER<<endl;  //maximum #iterations for "pure-2opt" to converge
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));
                    start=clock();
                    opt_path =_2opt(opt_path, new_adj_mat, &MAX_ITER, &pure_2opt_randomize_result,1);  //randomize=1   
                    end = clock();
                    cout<<"count2: "<<MAX_ITER<<endl;  //maximum #iterations taken for "randomized-2opt" to converge
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));
                    start=clock();
                    naive_weave_result=max(max(naive_weave(new_nodes,1), naive_weave(new_nodes,0.5)), naive_weave(new_nodes,2));
                    end = clock();
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));
                    start=clock();
                    hoffmann_weave_result=max(max(hoffmann_weave(new_nodes,1), hoffmann_weave(new_nodes,0.5)), hoffmann_weave(new_nodes,2));
                    end = clock();
                    time_taken.push_back(double(end-start)/double(CLOCKS_PER_SEC));
                    cout<<"k: "<<k<<endl;  //printing obtained results in terminal, starting with perturbation's iteration count
                    cout<<"0. scatter_bound: "<<scatter_bound<<endl;
                    cout<<"1. naive greedy: "<<naive_greedy_result<<" ;alpha: "<<naive_greedy_result/scatter_bound<<" ;time_taken: "<<time_taken[0]<<endl;
                    cout<<"2. dirac: "<<dirac_result<<" ;alpha: "<<dirac_result/scatter_bound<<" ;time_taken: "<<time_taken[1]<<endl;
                    cout<<"3. pure 2 opt (dirac init): "<<pure_2opt_result<<" ;alpha: "<<pure_2opt_result/scatter_bound<<" ;time_taken: "<<time_taken[2]<<endl;
                    cout<<"4. pure 2 opt (randomized init): "<<pure_2opt_randomize_result<<" ;alpha: "<<pure_2opt_randomize_result/scatter_bound<<" ;time_taken: "<<time_taken[3]<<endl;
                    cout<<"5. naive weave: "<<naive_weave_result<<" ;alpha: "<<naive_weave_result/scatter_bound<<" ;time_taken: "<<time_taken[4]<<endl;
                    cout<<"6. hoffmann weave: "<<hoffmann_weave_result<<" ;alpha: "<<hoffmann_weave_result/scatter_bound<<" ;time_taken: "<<time_taken[5]<<endl;
                    MyFile <<"k: "<<k<<endl;  //storing obtained results in a file in the below order
                    MyFile << "scatter_bound: "<<scatter_bound<<endl;
                    MyFile << "results: "<<naive_greedy_result<<','<<dirac_result<<','<<pure_2opt_result<<','<<pure_2opt_randomize_result<<','<<naive_weave_result<<','<<hoffmann_weave_result<<endl;
                    MyFile << "times: "<<time_taken[0]<<','<<time_taken[1]<<','<<time_taken[2]<<','<<time_taken[3]<<','<<time_taken[4]<<','<<time_taken[5]<<endl;
                }
            }
        }
    }
    MyFile.close();
    return 0;
}
//ran the setup in google-cloud and parsed 'final_results' data log file to obtain results in paper.
