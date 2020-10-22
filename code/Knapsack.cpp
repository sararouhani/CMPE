/*
 * Knapsack.cpp
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


bool Item_sorter(Item const &lhs, Item &rhs) {
    return lhs.cost < rhs.cost;
}

void sortBin(Bin &bin) {
    sort(bin.begin(), bin.end(), &Item_sorter);
}

void RemoveDominatedItems(Bin &bin_) {
    if (bin_.empty()) return;
    sortBin(bin_);
    //return;
    Bin bin = bin_;
    bin_ = vector<Item>();
    long double max_profit_so_far = bin[0].profit;
    bin_.emplace_back(bin[0]);
    for (int i = 1; i < bin.size(); i++) {
        if (bin[i].profit <= max_profit_so_far) continue;
        max_profit_so_far = bin[i].profit;
        bin_.push_back(bin[i]);
    }
}

void print_mckp(MCKP &mckp) {
    cout << "Num bins = " << mckp.size() << endl;
    for (int i = 0; i < mckp.size(); i++) {
        cout << "Bin = " << i + 1 << ": ";
        for (int j = 0; j < mckp[i].size(); j++) {
            cout << "(" << mckp[i][j].cost << "," << mckp[i][j].profit << ") ";
        }
        cout << endl;
    }
}

// Returns the current best value and stores the best solution so far in solution
long double
greedy_solve_MCKP(vector<vector<long double> > &weights, vector<vector<long double> > &profits, long double max_cost,
                 vector<int> &solution) {
    // Begin: Construct the MCKP from profits and weights
    if (weights.size() != profits.size()) {
        cerr << "Mismatch in the number of Bins\n";
        exit(-1);
    }
    MCKP mckp = vector<Bin>(weights.size());
    for (int i = 0; i < weights.size(); i++) {
        mckp[i] = vector<Item>(weights[i].size());
        if (weights[i].size() != profits[i].size()) {
            cerr << "Mismatch in the number of items in Bin " << i << "\n";
            exit(-1);
        }
        for (int j = 0; j < weights[i].size(); j++) {
            mckp[i][j].profit = profits[i][j];
            mckp[i][j].cost = weights[i][j];
            mckp[i][j].pos_in_bin = j;
        }
    }
    //print_mckp(mckp);
    // End: Cosntruct MCKP

    // Construct Greedy solution to MCKP
    int num_bins = weights.size();
    // Step 1. Remove Dominated items in each bin

    vector<int> multi_item_bin_ids;
    for (int i = 0; i < num_bins; i++) {
        RemoveDominatedItems(mckp[i]);
        if ((int) mckp[i].size() > 1) {
            multi_item_bin_ids.emplace_back(i);
        }
    }

    //print_mckp(mckp);
    vector<Item> current_solution(num_bins);
    long double current_total_profit = 0.0;
    long double current_total_cost = 0.0;
    solution = vector<int>(num_bins);
    vector<int> current_solution_index(num_bins, 0);
    for (int i = 0; i < num_bins; i++) {
        current_solution[i] = mckp[i][0];
        current_total_cost += mckp[i][0].cost;
        current_total_profit += mckp[i][0].profit;
        // Check if problem is infeasible
        if (current_total_cost > max_cost) {
            return -1 * std::numeric_limits<long double>::max();
        }
    }

    for(int i=0;i<num_bins;i++){
        for (int j = 1; j < mckp[i].size(); j++) {
            long double new_total_cost = current_total_cost + mckp[i][j].cost - current_solution[i].cost;
            long double new_total_profit = current_total_profit + mckp[i][j].profit - current_solution[i].profit;
            if (new_total_cost <= max_cost && new_total_profit>current_total_profit) {
                //cout<<current_total_profit<<" "<<new_total_profit<<endl;
                current_total_cost = new_total_cost;
                current_solution[i] = mckp[i][j];
                current_total_profit = new_total_profit;
                current_solution_index[i] = j;
            }
        }
    }
    /*
    cout << "First Greedy objective = " << current_total_profit << ", cost = " << current_total_cost
         << ", max-cost = " << max_cost << endl;
         */
    long double best_total_profit = current_total_profit;
    vector<Item> best_solution = current_solution;
    long double best_total_cost = current_total_cost;

    // Step 3. Perform Local Search
    // Step 2. Construct a greedy solution multiple times by replacing low profit items with high profit items
    // Fill up residual capacity with the bin where the iteration terminates
    for(int num_restarts=0;num_restarts<100;num_restarts++) {
        current_total_profit = 0.0;
        current_total_cost = 0.0;
        current_solution_index=vector<int>(num_bins, 0);
        for (int i = 0; i < num_bins; i++) {
            current_solution[i] = mckp[i][0];
            current_total_cost += mckp[i][0].cost;
            current_total_profit += mckp[i][0].profit;
        }
        for (int iter = 0; iter < 1000; iter++) {
            int i = multi_item_bin_ids[rand() % multi_item_bin_ids.size()];
            int mi = current_solution_index[i];
            for (int j = 0; j < mckp[i].size(); j++) {
                if (j == mi) continue;
                long double new_total_cost = current_total_cost + mckp[i][j].cost - current_solution[i].cost;
                long double new_total_profit = current_total_profit + mckp[i][j].profit - current_solution[i].profit;
                if (new_total_cost <= max_cost && new_total_profit > current_total_profit) {
                    //cout << current_total_profit << " " << new_total_profit << endl;
                    current_total_cost = new_total_cost;
                    current_solution[i] = mckp[i][j];
                    current_total_profit = new_total_profit;
                    current_solution_index[i] = j;
                }
            }
            if (current_total_profit > best_total_profit) {
                best_total_profit = current_total_profit;
                best_solution = current_solution;
                best_total_cost = current_total_cost;
            }
        }
    }
    //cout<<endl;

    //cout << "After local search: objective = " << best_total_profit << ", cost = " << best_total_cost << ", max-cost = "
      //   << max_cost << endl;

    for (int i = 0; i < num_bins; i++) {
        solution[i] = best_solution[i].pos_in_bin;
    }
    return best_total_profit;
}

