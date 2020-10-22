/*
 * main_mckp_test.cpp
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

#include "Knapsack.h"
#include <algorithm>
#include <iostream>
#include <ortools/linear_solver/linear_solver.h>


extern long double
greedy_solve_MCKP(vector<vector<long double> > &weights, vector<vector<long double> > &profits, long double max_cost,
                 vector<int> &solution);

void print_help(const string &program_name) {
    cerr << "Usage: " << program_name << " -n <num-bins> -s <size of each bin> -i <max-int-size> -seed <seed> \n";
    cerr << "-------------------------------------------------------------------------\n";
}

int main(int argc, char *argv[]) {

    int num_bins = 5, size_bin = 5, max_int = 1000, seed = 10000L;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            num_bins = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            size_bin = atoi(argv[i + 1]);

        } else if (strcmp(argv[i], "-i") == 0) {
            max_int = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-seed") == 0) {
            seed = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv[0]);
            exit(-1);
        }
    }
    srand(seed);

    //Generate random knapsack
    vector<vector<long double> > values(num_bins);
    vector<vector<long double> > weights(num_bins);

    long double max_cost = 0.0;
    for (int i = 0; i < num_bins; i++) {
        values[i] = vector<long double>(size_bin);
        weights[i] = vector<long double>(size_bin);
        for (int j = 0; j < size_bin; j++) {
            values[i][j] = rand() % max_int + 1;
            //weights[i][j] = rand() % max_int + 1;
            long double rand_val=static_cast <long double> (rand()) /(static_cast <long double> (RAND_MAX)+1);
            rand_val/=100.0;
            if(rand()%2==0)
                weights[i][j]=values[i][j]+rand_val;
            else
                weights[i][j]=values[i][j]-rand_val;
        }
        max_cost += weights[i][rand() % size_bin];
    }
    vector<int> solution;
    cout<<"Greedy algorithm solution  = "<<greedy_solve_MCKP(weights, values, max_cost, solution)<<endl;

    using namespace operations_research;
    MPSolver solver("simple_mip_program",
                    MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    const double infinity = solver.infinity();
    // x[j] is an array of non-negative, integer variables.
    vector<vector<const MPVariable *> > x(num_bins);
    for (int i = 0; i < num_bins; ++i) {
        x[i] = vector<const MPVariable *>(values[i].size());
        for (int j = 0; j < x[i].size(); j++) {
            //x[i][j]=solver.MakeNumVar(0.0,1.0,"");
            x[i][j] = solver.MakeBoolVar("");
        }
    }
    // Add the cost constraint
    MPConstraint *constraint1 = solver.MakeRowConstraint(-infinity, max_cost, "");
    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[i].size(); j++) {
            constraint1->SetCoefficient(x[i][j], weights[i][j]);
        }
    }
    // Add the constraint that exactly one item from each bin must be chosen
    for (int i = 0; i < num_bins; i++) {
        MPConstraint *constraint = solver.MakeRowConstraint(1.0, 1.0, "");
        for (int j = 0; j < values[i].size(); j++) {
            constraint->SetCoefficient(x[i][j], 1.0);
        }
    }

    // Create the objective function.
    MPObjective *const objective = solver.MutableObjective();
    for (int i = 0; i < values.size(); i++) {
        for (int j = 0; j < values[i].size(); j++) {
            objective->SetCoefficient(x[i][j], values[i][j]);
        }
    }
    objective->SetMaximization();
    solver.set_time_limit(30 * 1000);
    //solver.SetPrimalTolerance(1e-100);
    const MPSolver::ResultStatus result_status = solver.Solve();
    // Check that the problem has an optimal solution.
    if (result_status != MPSolver::OPTIMAL) {
        LOG(FATAL) << "The problem does not have an optimal solution.";
    }
    LOG(INFO) << "Solution:";
    LOG(INFO) << "Optimal objective value = " << objective->Value();

    long double cbc_solution_verify = 0.0;
    long double cbc_cost_verify = 0.0;
    for (int i = 0; i < values.size(); i++) {
        for (int j = 0; j < values[i].size(); j++) {
            if (x[i][j]->solution_value() > 0.99) {
                cbc_solution_verify += values[i][j];
                cbc_cost_verify += weights[i][j];
                break;
            }
        }
    }
    cout << "CBC solution =" << cbc_solution_verify << " CBC cost = " << cbc_cost_verify << ", Max cost = " << max_cost
         << endl;
}

