/*
 * MN.cpp
 *
 *
 *      Authors: Sara Rouhani,Tahrima Rahman,Vibhav Gogate. The University of Texas at Dallas
 *      Contacts: {sara.rouhani,tahrima.rahman,vibhav.gogate}@utdallas.edu
 *
 *      MIT License
 *
 *      Permission is hereby granted, free of charge, to any person obtaining a copy
 *      of this software and associated documentation files (the "Software"), to deal
 *      in the Software without restriction, including without limitation the rights
 *      to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *      copies of the Software, and to permit persons to whom the Software is
 *      furnished to do so, subject to the following conditions:
 *
 *      The above copyright notice and this permission notice shall be included in all
 *      copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *      IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *      AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *      LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *      OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *      SOFTWARE.
 */

#include "MN.h"
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <iterator>




extern long double
greedy_solve_MCKP(vector<vector<long double> > &weights, vector<vector<long double> > &profits, long double max_cost,
                 vector<int> &solution);


// Read the Markov network
void MN::readMN(string filename)
{
    ifstream infile(filename);
    int num_variables;
    string tmp_string;
    infile >> tmp_string;
    if (tmp_string.compare("MARKOV") != 0) {
        cerr << "Not a Markov network\n";
        exit(-1);
        return;
    }
    infile >> num_variables;
    // Read domains
    variables = vector<Variable*>(num_variables);
    for (int i = 0; i < num_variables; i++) {
        int domain_size;
        infile >> domain_size;
        variables[i] = new Variable(i, domain_size);
    }
    int num_functions;
    infile >> num_functions;
    vector < vector<Variable*> > scope(num_functions);
    for (int i = 0; i < num_functions; i++) {
        // Read parents of variables
        int num_vars_in_func;
        infile >> num_vars_in_func;
        scope[i] = vector<Variable*>(num_vars_in_func);
        for (int j = 0; j < num_vars_in_func; j++) {
            int temp;
            infile >> temp;
            scope[i][j] = variables[temp];
        }
    }
    potentials = vector<Potential*>(num_functions);
    //srand(100000000L);
    for (int i = 0; i < num_functions; i++) {
        int num_entries;
        infile >> num_entries;
        potentials[i] = new Potential();
        potentials[i]->variables = scope[i];
        int num_values = Variable::getDomainSize(scope[i]);
        potentials[i]->table = vector<long double>(num_values);
        for (int j = 0; j < num_values; j++) {
            Variable::setAddress(scope[i], j);
            long double value;
            infile >> value;
            int entry = Variable::getAddress(potentials[i]->variables);
            if (value > 0.0)
                potentials[i]->table[entry] = log(value);
            else {
                cerr << "Cannot handle zeros: Log-potentials\n";
                exit(-1);
            }
        }
    }
    infile.close();
}

void MN::readMN2(string filename,MN& mn1)
{
    ifstream infile(filename);
    int num_variables;
    string tmp_string;
    infile >> tmp_string;
    if (tmp_string.compare("MARKOV") != 0) {
        cerr << "Not a Markov network\n";
        exit(-1);
        return;
    }
    infile >> num_variables;
    if(num_variables!=mn1.variables.size()){
        cerr << "Markov networks do not match in number of variables\n";
        exit(-1);
        return;
    }
    // Read domains
    variables = mn1.variables;
    for (int i = 0; i < num_variables; i++) {
        int domain_size;
        infile >> domain_size;
        if(variables[i]->domain_size!=domain_size){
            cerr << "Variables in Markov networks do not match; different domains\n";
            exit(-1);
            return;
        }
    }
    int num_functions;
    infile >> num_functions;
    vector < vector<Variable*> > scope(num_functions);
    for (int i = 0; i < num_functions; i++) {
        // Read parents of variables
        int num_vars_in_func;
        infile >> num_vars_in_func;
        scope[i] = vector<Variable*>(num_vars_in_func);
        for (int j = 0; j < num_vars_in_func; j++) {
            int temp;
            infile >> temp;
            scope[i][j] = variables[temp];
        }
    }
    potentials = vector<Potential*>(num_functions);
    //srand(100000000L);
    for (int i = 0; i < num_functions; i++) {
        int num_entries;
        infile >> num_entries;
        potentials[i] = new Potential();
        potentials[i]->variables = scope[i];
        int num_values = Variable::getDomainSize(scope[i]);
        potentials[i]->table = vector<long double>(num_values);
        for (int j = 0; j < num_values; j++) {
            Variable::setAddress(scope[i], j);
            long double value;
            infile >> value;
            int entry = Variable::getAddress(potentials[i]->variables);
            if (value > 0.0)
                potentials[i]->table[entry] = log(value);
            else {
                cerr << "Cannot handle zeros: Log-potentials\n";
                exit(-1);
            }
        }
    }
    infile.close();
}
vector <set<int>> MN::findKseparator(int k, vector<Variable*>& cut_variables)
{
    cut_variables = vector<Variable*>();
    vector < set<int> > graph(variables.size());
    vector<int> degree(variables.size(), 0);
    //Construct the graph
    for (auto & potential : potentials) {
        for (int j = 0; j < potential->variables.size(); j++) {
            int var1_id = potential->variables[j]->id;
            for (int k = j + 1; k < potential->variables.size(); k++) {
                int var2_id = potential->variables[k]->id;
                graph[var1_id].insert(var2_id);
                graph[var2_id].insert(var1_id);
            }
        }
    }
    // Initialize the degrees
    for (int i = 0; i < graph.size(); i++) {
        degree[i] = graph[i].size();
    }
    // remove all vertices with degree >= k
    while (true) {
        //Find the node with the highest degree
        int max_degree_id = distance(degree.begin(),max_element(degree.begin(),degree.end()));
        int max_degree = degree[max_degree_id];
        // there is no vertex with degree >= k
        if (max_degree < k)
            break;
        cut_variables.push_back(variables[max_degree_id]);
        // Remove the variable max_degree_id from the graph and update degree
        for (auto i = graph[max_degree_id].begin(); i != graph[max_degree_id].end(); i++) {
            graph[*i].erase(max_degree_id);
            degree[*i]--;
        }
        degree[max_degree_id] = 0;
        graph[max_degree_id] = set<int>();
    }
    //remove from components with more than k vertices
    while (true){
        int max_component = 0;
        int max_component_id = -1;
        //find the components
        vector <set<int>> components;
        components = this->connectedComponents(graph);
        //cout<<"Number of components = "<<components.size()<<endl;
        //find the component with max size
        for (int i = 0; i < components.size(); i++) {
            if (components[i].size() > max_component) {
                max_component = components[i].size();
                max_component_id = i;
            }

        }
        //cout<<"Max component size = "<<max_component<<endl;
        // there is no vertices with degree >= k
        if (max_component <= k)
            break;
        //find max degree vertices in max size component
        int max_degree = 0;
        int max_degree_id = -1;
        for (std::__1::__tree_const_iterator<int, std::__1::__tree_node<int, void *> *, long>::value_type i : components[max_component_id]) {
            if (degree[i] > max_degree) {
                max_degree = degree[i];
                max_degree_id = i;
            }
        }
        cut_variables.push_back(variables[max_degree_id]);
        // Remove the variable max_degree_id from the graph and update degree
        for (auto i = graph[max_degree_id].begin(); i != graph[max_degree_id].end(); i++) {
            graph[*i].erase(max_degree_id);
            degree[*i]--;
        }
        degree[max_degree_id] = 0;
        graph[max_degree_id] = set<int>();
    }
    //remove cut variables from components
    vector <set<int>> output_components;
    vector <set<int>> final_components;
    final_components = this->connectedComponents(graph);
    set<int> cut_set;
    for (auto & cut_variable : cut_variables){
        cut_set.insert(cut_variable->id);
    }
    for (auto & final_component : final_components) {
        if (final_component.size() ==1)
            for (auto j = final_component.begin(); j != final_component.end(); j++) {
                if (cut_set.find(*j) == cut_set.end())
                    output_components.push_back(final_component);
            }
        else{
            output_components.push_back(final_component);
        }
    }
    return output_components;
    /*
    //Merge some components
    vector<bool> visited(output_components.size(),false);
    vector <set<int>> merged_components;
    for(int i=0;i<output_components.size();i++){
        if(visited[i]) continue;
        visited[i]=true;
        merged_components.push_back(output_components[i]);
        int a=merged_components.size()-1;
        for(int j=i+1;j<output_components.size();j++){
            if(!visited[j]){
                if(output_components[j].size()+merged_components[a].size()<=k){
                    merged_components[a].insert(output_components[j].begin(),output_components[j].end());
                    visited[j]=true;
                    if(merged_components[a].size()==k) {
                        break;
                    }
                }
            }
        }
    }

    //for(int i=0;i<merged_components.size();i++){
      //  cout<<" size of component "<<i+1<<" is "<<merged_components[i].size()<<endl;
    //}
    return merged_components;
     */

}

vector <set<int>> MN::connectedComponents(vector <set<int>> adj)
{
    vector < set<int> > components;
    bool *visited = new bool[adj.size()];
    for(int v = 0; v < adj.size(); v++)
        visited[v] = false;
    for (int v=0; v<adj.size(); v++)
    {
        if (!visited[v])
        {
            components.emplace_back();
            components[components.size()-1].insert(v);
            DFSUtil(v, visited, adj, components);
        }
    }
    delete[] visited;
    return components;
}

void MN::DFSUtil(int v, bool visited[], vector <set<int>> adj, vector <set<int>>& components){
    visited[v] = true;
    for(auto i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
        {   components[components.size()-1].insert(*i);
            DFSUtil(*i, visited, adj, components);}
}

vector <set<int>> MN::generate_buckets(vector <set<int>>& components)
{
    vector <set<int>> buckets (components.size()+1);
    set<int> cut_vars;
    vector <set<int>> potential_vars (potentials.size());
    set<int> potentials_ids;
    //a set of cpts ids
    for (int i = 0; i < potentials.size(); i++){
        potentials_ids.insert(i);
    }
    //a set of cpts var ids
    for (int i = 0; i < potentials.size(); i++) {
        for (int j = 0; j < potentials[i]->variables.size(); j++)
            potential_vars[i].insert(potentials[i]->variables[j]->id);
    }
    //generate buckets
    for (int i = 0; i < potentials.size(); i++){
        for (int j = 0; j < components.size(); j++){
            set<int> intersect;
            std::set_intersection(potential_vars[i].begin(),potential_vars[i].end(),components[j].begin(),components[j].end(),
                                  std::inserter(intersect,intersect.begin()));
            if (!intersect.empty()){
                buckets[j].insert(i);
                potentials_ids.erase(i);
                break;
            }
        }
    }
    //put cut var functions in last bucket
    for (std::set<int>::iterator it=potentials_ids.begin(); it!=potentials_ids.end(); ++it)
        buckets[buckets.size()-1].insert(*it);
    return buckets;
}


bool MN::knapsack_greedy(long double logq, vector<Potential>& functions, vector<Potential>& functions_c, vector<int>& var_assignment, long double& best_prob)
{

    //generating weights
    vector<vector<long double>> weights(functions.size());
    vector<int> assignment;
    for(int i=0;i<functions.size();i++){
        weights[i]=functions[i].table;
    }
    //generating values
    vector<vector<long double>> values(functions_c.size());
    for(int i=0;i<functions_c.size();i++){
        values[i]=functions_c[i].table;
    }
    best_prob=greedy_solve_MCKP(weights,values,logq,assignment);

    for (int t = 0; t < functions.size(); t++) {
        //set the best assignment for weight and value
        Variable::setAddress(functions[t].variables,assignment[t]);
        Variable::setAddress(functions_c[t].variables,assignment[t]);
        for (int l = 0; l < functions[t].variables.size(); l++)
            var_assignment[functions[t].variables[l]->id] = functions[t].variables[l]->value;
    }
    return true;
}

long double MN::run_experiments_neurips(MN& mn_c, long double logq, int k, ostream& out1, int max_time){
    //creating a graphical model from the original one
    //Uncomment the following line if you want to generate knapsack networks randomly
    //this->create_knapsack(mn_c);
    //MN_constructed mn_c(potentials, variables);
    long double best_prob=0;
    long double greedy_output;
    vector<int> assignment(variables.size());
    vector<int> best_assignment(variables.size());
    out1.precision(20);
    //run for all values less than k
    for (int h = 1; h < k+1; h +=2) {
        std::time_t start_time = std::time(nullptr);
        std::time_t write_time = std::time(nullptr);
        //initial with the worst answer
        best_prob = -1 * std::numeric_limits<long double>::max();
        vector<Variable*> cut_variables;
        vector <set<int>> components;
        components = findKseparator(h,cut_variables);
        vector <set<int>> buckets;
        buckets = generate_buckets(components);
        cout<<"Statistics:"<<endl;
        cout<<"K = "<<h<<endl;
        cout<<"Number of Variables in the K-separator = "<<cut_variables.size()<<endl;
        cout<<"Number of components = "<<buckets.size()-1<<endl;
        for (int i = 0; i<buckets.size()-1; i++){
            cout<<"number of variables in component "<<i+1<<" "<<components[i].size()<<endl;
        }
        for(int i=0;i<variables.size();i++){
            variables[i]->value=rand()%variables[i]->domain_size;
        }
        //run for sampling number
        int num_assignments_explored=0;
        while (true) {
            num_assignments_explored++;
            //random values to k sep variables for both MNs
            long double current_value=mn_c.getValue();
            long double current_weight=this->getValue();
            int change_variable=-1;
            int change_value=-1;
            // If the current solution is not feasible move towards a feasible solution
            if (current_weight > logq) {
                for (int j = 0; j < cut_variables.size(); j++) {
                    int index = cut_variables[j]->value;
                    for (int k = 0; k < cut_variables[j]->domain_size; k++) {
                        if (index == k) continue;
                        cut_variables[j]->value = k;
                        long double sol_value = mn_c.getValue();
                        long double sol_weight = getValue();
                        // If you have already found a feasible solution move towards better objective
                        if (current_weight <= logq) {
                            if (sol_weight <= logq && sol_value > current_value) {
                                current_value = sol_value;
                                current_weight = sol_weight;
                                change_variable = j;
                                change_value = k;
                            }
                        } else {
                            if (sol_weight < current_weight) {
                                current_value = sol_value;
                                current_weight = sol_weight;
                                change_variable = j;
                                change_value = k;
                            }
                        }
                    }
                    cut_variables[j]->value = index;
                }
            }
            else {
                for (int j = 0; j < cut_variables.size(); j++) {
                    int index = cut_variables[j]->value;
                    for (int k = 0; k < cut_variables[j]->domain_size; k++) {
                        if (index == k) continue;
                        cut_variables[j]->value = k;
                        long double sol_value = mn_c.getValue();
                        long double sol_weight = getValue();
                        if (sol_weight <= logq && sol_value > current_value) {
                            current_value = sol_value;
                            change_variable = j;
                            change_value = k;
                            if (sol_value > best_prob) {
                                best_prob = sol_value;
                            }
                        }
                    }
                    cut_variables[j]->value = index;
                }
            }
            // Check for Local maxima
            if (change_variable==-1){
                    //Escape the local maxima by making random assignments to cut variables
                    // With 10% probability make a random global move
                    if (rand()%100>=90){
                        for (int j = 0; j < cut_variables.size(); j++) {
                            cut_variables[j]->value = rand() % cut_variables[j]->domain_size;
                            assignment[cut_variables[j]->id] = cut_variables[j]->value;
                        }
                    }
                    else {
                        // With 90% probability make a local random move
                        int j=rand()%cut_variables.size();
                        cut_variables[j]->value = rand() % cut_variables[j]->domain_size;
                        assignment[cut_variables[j]->id] = cut_variables[j]->value;
                    }
            }
            else {
                // No local maxima: Accept the move
                assignment[cut_variables[change_variable]->id] = change_value;
                cut_variables[change_variable]->value=change_value;
            }

            //generating functions for greedy for both MNs
            vector<Potential> greedy_functions;
            vector<Potential> greedy_functions_c;
            for (int i = 0; i<buckets.size()-1; i++){
                vector<Potential*> functions;
                vector<Potential*> functions_c;
                vector<Variable*> component_variables;
                for (std::set<int>::iterator it=buckets[i].begin(); it!=buckets[i].end(); ++it){
                    functions.push_back(this->potentials[*it]);
                    functions_c.push_back(mn_c.potentials[*it]);
                }
                for (std::set<int>::iterator it=components[i].begin(); it!=components[i].end(); ++it)
                    component_variables.push_back(this->variables[*it]);
                greedy_functions.emplace_back(functions, component_variables);
                greedy_functions_c.emplace_back(functions_c, component_variables);
            }
            //conditioning logq on k-sep vars of original MN
            long double q_sep = 0.0;
            for (std::set<int>::iterator it=buckets[buckets.size()-1].begin(); it!=buckets[buckets.size()-1].end(); ++it)
                q_sep += this->potentials[*it]->getValue();
            long double obj_c = 0.0;
            for (std::set<int>::iterator it=buckets[buckets.size()-1].begin(); it!=buckets[buckets.size()-1].end(); ++it)
                obj_c += mn_c.potentials[*it]->getValue();

            long double new_logq = logq - q_sep;
            knapsack_greedy(new_logq, greedy_functions, greedy_functions_c, assignment, greedy_output);
            if ((greedy_output + obj_c) > best_prob){
                best_prob = greedy_output + obj_c;
            }
            std::time_t curr_time = std::time(nullptr);
            if ((curr_time - start_time) % GlobalSearchOptions::print_interval == 0 && write_time != curr_time) {
                write_time = std::time(nullptr);

                out1 << logq << "," << h << "," <<  best_prob << ',' << curr_time - start_time << "\n";
                cerr << std::setprecision(20) <<logq << "," << h << "," <<  best_prob << ',' << num_assignments_explored<<","<<curr_time - start_time << "\n";
            }
            if ((curr_time - start_time) >= max_time) {
                break;
            }
        }
    }
    return best_prob;
};



