/*
 * MN.h
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

#ifndef MN_H_
#define MN_H_

//neurips
#include <set>
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;
struct GlobalSearchOptions{
    static int print_interval;
};

/*
 * struct Variable
 * id:			the variable number (numbering starts from 0)
 * 				Default=-1 means that the variable is invalid
 * domain_size: the number of values in the domain of the variable.
 * 				For example, for binary variables domain_size=2
 * 				Default=-1 means that variable is in invalid state
 * 				and needs to be initialized
 * value:		the value currently assigned to the variable.
 * 				value should be in the set {0,...,domain_size}
 * 				Default=-1 means that variable is not assigned a value
 */
struct Variable {
    int id;
    int domain_size;
    int value;

    /*
     * Functions
     */
    Variable() :
            id(-1),domain_size(-1), value(-1) {
    }
    Variable(int id_,int domain_size_) :
            id(id_),domain_size(domain_size_), value(-1) {
    }
    ~Variable(){}
    /*
     * Useful Static Functions
     */
    inline static int getAddress(const vector<Variable*>& variables){
        int add_ress = 0;
        int multiplier = 1;
        for (auto variable : variables) {
            add_ress += (multiplier * variable->value);
            multiplier *= variable->domain_size;
        }
        return add_ress;
    }

    inline static int getDomainSize(const vector<Variable*>& variables) {
        int domain_size = 1;
        for (auto variable : variables)
            domain_size *= variable->domain_size;
        return domain_size;
    }
    // Get the maximum domain size of the set of variables
    inline static void setAddress(const vector<Variable*>& variables, const int add_ress_) {
        int add_ress = add_ress_;
        for (auto variable : variables) {
            variable->value = add_ress % variable->domain_size;
            add_ress /= variable->domain_size;
        }
    }
};

/*
 * Note that all potentials are log-potentials.
 * A log-potential is a pair <X,table> where
 *       X is a the scope of the potential and
 *       table is a set of weights, one for each possible assignment to X
 */
struct Potential
{
    vector<Variable*> variables;
    vector<long double> table;
    Potential()= default;
    long double getValue() {
        return table[Variable::getAddress(variables)];
    }
    Potential(const vector <Potential*>& potentials, const vector<Variable*>&component_variables) {
        variables = component_variables;
        int num_values = Variable::getDomainSize(variables);
        table = vector<long double>(num_values, 0.0);
        for (int i = 0; i < num_values; i++) {
            Variable::setAddress(component_variables, i);
            for (auto & potential : potentials)
                table[i] += potential->getValue();
        }
    }
};

struct MN{
    vector<Variable*> variables;
    vector<Potential*> potentials;
    MN(){}
    void readMN(string filename_);
    void readMN2(string filename_,MN& mn1);

    inline long double getValue(){long double logp=0.0;for(int i=0;i<potentials.size();i++) logp+=potentials[i]->getValue(); return logp;}

    vector <set<int>> findKseparator(int k, vector<Variable*>& cut_variables);
    static void DFSUtil(int v, bool visited[], vector <set<int>> adj, vector <set<int>>& components);
    vector <set<int>> connectedComponents(vector <set<int>> adj);
    vector <set<int>> generate_buckets(vector <set<int>>& components);
    bool knapsack_greedy(long double logq, vector<Potential>& functions, vector<Potential>& functions_c, vector<int>& var_assignment, long double& best_prob);
    long double run_experiments_neurips(MN& mn_c, long double logq, int k=15, ostream& out1=cout, int max_time=1200);
    void writeMPS(MN& mn_c, long double logp, ostream &out1=cout);
};
#endif /* MN_H_ */
