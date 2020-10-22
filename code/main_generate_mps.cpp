/*
 * main_generate_mps.cpp
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

#include <vector>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <ortools/linear_solver/linear_solver.h>

#include "MN.h"

using namespace std;

void print_help(const string &program_name) {
    cerr << "Usage: " << program_name << " -m1 <uaifilename1> -m2 <uaifilename2> -o <outfilename> -q <q-value>\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Required Option\n";
    cerr << "\t\t\t uaifilename1 and uaifilename2: are evidence instantiated Markov networks in UAI format\n";
    cerr << "\t\t\t outfilename: mpsfile will be stored here\n";
    cerr << "\t\t\t q-value: (Real number): constraint on weight of the assignment in CMPE\n";
    cerr << "-------------------------------------------------------------------------\n";

}


int GlobalSearchOptions::print_interval = 1; //print stat every 100 seconds


using namespace operations_research;
void MN::writeMPS(MN& mn_c, long double logp, ostream &out1) {

    MPSolver solver("simple_mip_program",
                    MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    vector<Potential*>& g=potentials;
    vector<Potential*>& h=mn_c.potentials;
    vector<vector<int> > var2funcids(variables.size());
    for(int i=0;i<g.size();i++){
        for(int j=0;j<g[i]->variables.size();j++)
            var2funcids[g[i]->variables[j]->id].push_back(i);
    }
    const double infinity = solver.infinity();
    // x[j] is an array of non-negative, integer variables.
    vector<vector<const MPVariable*> > x(potentials.size());
    for (int i = 0; i < potentials.size(); ++i) {
        x[i]=vector<const MPVariable*> (g[i]->table.size());
        for(int j=0;j<x[i].size();j++){
            //x[i][j]=solver.MakeNumVar(0.0,1.0,"");
            x[i][j]=solver.MakeBoolVar("");
        }
    }
    // Write the constraint that \sum_{i,j} g[i][j]*x[i][j] should be <=logp
    MPConstraint* constraint1 = solver.MakeRowConstraint(-infinity, logp, "");
    for(int i=0;i<g.size();i++){
        int domain_size=Variable::getDomainSize(g[i]->variables);
        for(int j=0;j<domain_size;j++){
            Variable::setAddress(g[i]->variables,j);
            int entry=Variable::getAddress(g[i]->variables);
            constraint1->SetCoefficient(x[i][j],g[i]->table[entry]);
        }
    }
    // Hard constraint to make sure that exactly one value is chosen from each potential
    // Write the constraint that \sum_j x[i][j]=1 for each i
    for(int i=0;i<g.size();i++){
        MPConstraint* constraint = solver.MakeRowConstraint(1.0, 1.0, "");
        for(int j=0;j<g[i]->table.size();j++){
            constraint->SetCoefficient(x[i][j],1.0);
        }
    }
    // Write the following constraints:
    //  for each variable A
    //      for each function i such that i mentions A
    //          for each function j such that j mentions A i not-equal-to j
    //              Constraint for A=0: \sum_{k|A=0} x[i][k] + \sum_{k|A=1} x[j][k] = 1
    //              Constraint for A=1: \sum_{k|A=1} x[i][k] + \sum_{k|A=0} x[j][k] = 1

    for(int a=0;a<variables.size();a++){
        int A=variables[a]->id;
        for(int b=0;b<var2funcids[A].size();b++){
            int i=var2funcids[A][b];
            int dsize_i=Variable::getDomainSize(g[i]->variables);
            for(int c=b+1;c<var2funcids[A].size();c++){
                int j=var2funcids[A][c];
                int dsize_j=Variable::getDomainSize(g[j]->variables);
                for(int d=0;d<2;d++){
                    MPConstraint *constraint = solver.MakeRowConstraint(1.0,1.0 , "");
                    for(int e=0;e<dsize_i;e++){
                        Variable::setAddress(g[i]->variables,e);
                        if (variables[a]->value==d){
                            constraint->SetCoefficient(x[i][e], 1.0);
                        }
                    }
                    for(int e=0;e<dsize_j;e++){
                        Variable::setAddress(g[j]->variables,e);
                        if (variables[a]->value==1-d){
                            constraint->SetCoefficient(x[j][e], 1.0);
                        }
                    }
                }
            }
        }
    }


    // Create the objective function.
    MPObjective* const objective = solver.MutableObjective();
    for(int i=0;i<h.size();i++){
        int domain_size=Variable::getDomainSize(h[i]->variables);
        for(int j=0;j<domain_size;j++){
            Variable::setAddress(h[i]->variables,j);
            int entry=Variable::getAddress(h[i]->variables);
            objective->SetCoefficient(x[i][j],-h[i]->table[entry]);
        }
    }
    objective->SetMinimization();
    string model_str;
    solver.ExportModelAsMpsFormat(false, false, &model_str);
    out1 << model_str;

    /*
     * Uncomment the following if you want to solve the mps using cbc
     *


    solver.set_time_limit(30*1000);
    //solver.SetPrimalTolerance(1e-100);
    const MPSolver::ResultStatus result_status = solver.Solve();
    // Check that the problem has an optimal solution.
    if (result_status != MPSolver::OPTIMAL) {
        LOG(FATAL) << "The problem does not have an optimal solution.";
    }
    LOG(INFO) << "Solution:";
    LOG(INFO) << "Optimal objective value = " << objective->Value();

    // Find the value of the objective function based on the solution
    for(int i=0;i<h.size();i++) {
        int domain_size = Variable::getDomainSize(h[i]->variables);
        for (int j = 0; j < domain_size; j++) {
            if(x[i][j]->solution_value()>0.99){
                Variable::setAddress(h[i]->variables,j);
                break;
            }
        }
    }

    long double constraint=0.0;
    for(int i=0;i<g.size();i++) {
        constraint+=g[i]->getValue();
    }

    cout<<"Value of constraint = "<<constraint<<" and q value = "<<logp<<endl;

    long double objective_value=0.0;
    for(int i=0;i<h.size();i++) {
        objective_value+=h[i]->getValue();
    }
    cout<<std::setprecision(20)<<"Value of objective = "<<objective_value<<endl;
     */
}

/*
 * This program can be run in two modes
 * Mode 1: Generate 20 q values
 * Mode 2: Run a particular algorithm for K minutes for given q values
 */
int main(int argc, char *argv[]) {
    srand(1000000L);
    //Default mode is 2
    int mode = 2;
    string uai_filename1;
    string uai_filename2;
    string out_filename;
    int max_time = 1200;
    int k = 15;
    int sampling_number = 1000;
    long double q;
    bool uaioption1 = false, uaioption2 = false, outoption = false, qoption = false;
    if (argc == 1) {
        print_help(argv[0]);
        exit(-1);
    }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m1") == 0) {
            uai_filename1 = argv[i + 1];
            uaioption1 = true;
        } else if (strcmp(argv[i], "-m2") == 0) {
            uai_filename2 = argv[i + 1];
            uaioption2 = true;
        } else if (strcmp(argv[i], "-o") == 0) {
            out_filename = argv[i + 1];
            outoption = true;
        } else if (strcmp(argv[i], "-q") == 0) {
            q = atof(argv[i + 1]);
            qoption = true;
        } else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv[0]);
            exit(-1);
        }
    }
    if (!uaioption1) {
        cerr << "UAI file1 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!uaioption2) {
        cerr << "UAI file2 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!outoption) {
        cerr << "Output file not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!qoption) {
        cerr << "Q not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    MN mn1, mn2;
    mn1.readMN(uai_filename1);
    mn2.readMN2(uai_filename2,mn1);
    if (mn1.variables.size() != mn2.variables.size()) {
        cerr << "Variable size mismatch\n";
        cerr << "Code requires the two Markov networks be defined over the same set of variables\n";
        exit(-1);
    }
    ofstream out(out_filename, ofstream::out);
    mn1.writeMPS(mn2,q,out);
    out.close();
    return 0;
}
