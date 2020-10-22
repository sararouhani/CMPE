/*
 * main.cpp
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
    cerr << "Usage: " << program_name << " -m1 <uaifilename1> -m2 <uaifilename2> -o <outfilename> -q <q-value>\n";
    cerr << "\t Other Options: [-t -k -s -si ]\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Required Option\n";
    cerr << "\t\t\t uaifilename1 and uaifilename2: are evidence instantiated Markov networks in UAI format\n";
    cerr << "\t\t\t outfilename: Results of experiments will be stored here\n";
    cerr << "\t\t\t q-value: (Real number): constraint on weight of the assignment in CMPE\n";
    cerr << "-------------------------------------------------------------------------\n";
    cerr << "\t\t Details on Other Options and Default values\n";
    cerr << "\t\t\t -t     [int]: max time for which each k is run; default 2\n";
    cerr << "\t\t\t -k     [int]: max k-seperator size. The code iterates from 1 to k in increments of 2\n";
    cerr << "\t\t\t -s     [int]: Seed for Repeatability; default 1000000L\n";
    cerr << "\t\t\t -si    [int]: print status every integer seconds; default 1\n";
    //cerr << "\t\t\t -w  [string]: Write file in MPS format and store it in string\n";
}


int GlobalSearchOptions::print_interval = 1; //print stat every 100 seconds

/*
 * This program can be run in two modes
 * Mode 1: Generate 20 q values
 * Mode 2: Run a particular algorithm for K minutes for given q values
 */
int main(int argc, char *argv[]) {
    srand(1000000L);
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
        } else if (strcmp(argv[i], "-t") == 0) {
            max_time = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-k") == 0) {
            k = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            int seed = atoi(argv[i + 1]);
            srand(seed);
        } else if (strcmp(argv[i], "-si") == 0) {
            GlobalSearchOptions::print_interval = atoi(argv[i + 1]);
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
    mn1.run_experiments_neurips(mn2, q, k, out, max_time);
    out.close();
    return 0;
}
