/*
 * main_generate_q.cpp
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

#include "MN.h"

using namespace std;

void print_help(const string &program_name) {
    cerr << "Usage: " << program_name << " -m <uaifilename> -o <outfilename>\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Required Option\n";
    cerr << "\t\t\t uaifilename: is an evidence instantiated Markov network in UAI format\n";
    cerr << "\t\t\t outfilename: q-values are stored in outfilename\n";
}


int GlobalSearchOptions::print_interval = 1; //print stat every 100 seconds

/*
 * This program can be run in two modes
 * Mode 1: Generate 20 q values
 * Mode 2: Run a particular algorithm for K minutes for given q values
 */
int main(int argc, char *argv[]) {
    string uai_filename,out_filename;
    bool uaioption = false, outoption = false;
    if (argc == 1) {
        print_help(argv[0]);
        exit(-1);
    }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0) {
            uai_filename = argv[i + 1];
            uaioption = true;
        } else if (strcmp(argv[i], "-o") == 0) {
            out_filename = argv[i + 1];
            outoption = true;
        } else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv[0]);
            exit(-1);
        }
    }
    if (!uaioption) {
        cerr << "UAI file1 not specified\n";
        print_help(argv[0]);
        exit(-1);
    }
    if (!outoption) {
        cerr << "Output file not specified\n";
        print_help(argv[0]);
        exit(-1);
    }


    ofstream out(out_filename);
    MN mn;
    mn.readMN(uai_filename);
    int num_samples = 1000000;
    vector<long double> q(num_samples);
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < mn.variables.size(); j++) {
            mn.variables[j]->value = rand() % mn.variables[j]->domain_size;
        }
        q[i] = mn.getValue();
    }
    sort(q.begin(), q.end());
    out.precision(20);
    for (int i = 0; i < 21; i++) {
        if (i == 0) {
            out << q[i] << endl;
        } else {
            out << q[i * 50000 - 1] << endl;
        }
    }
    out.close();
    return 0;
}
