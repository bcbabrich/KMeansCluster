// KMeansCluster.cpp : B Chase Babrich 2020, personal project

#include "vector"    // gives us vectors
#include <string>    // gives us getline method
#include <iostream>  // gives us cout, cin
#include <fstream>   // gives us file readin (fstream, ofstream)
#include <sstream>   // gives us ss(line), i.e., turns strings into streams :)
#include <math.h>    // gives us sqrt
#include <stdlib.h>  // gives us srand, rand (random numbers for guesing centers)
#include <time.h>    // use in conjunctio with rand
using namespace std; // don't have to write std:: everywhere


// D processing. Turns data file (csv) into datastructure (2-D vector of doubles).
//
// IN: name of a csv file with D in the following format:
//          a,b,...,n,L where a,b,...,n are ints >=0 and L is in [A-Z]
//
// OUT: a 2-d vector of ints representing the D in the csv file
vector<vector<double>> get_D(string filename)
{
    // initialize variables
    vector<vector<double>> D;
    string line;
    ifstream myfile(filename);

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            vector<double> row;

            // Create a stringstream of the current line
            // i.e., we'll parse a line on a per character basis
            std::stringstream ss(line);

            // Extract each comma separated value (string) as a double
            char c;
            string numstr = "";
            while (ss >> c) {
                // read in chars as ints
                numstr += c;

                // push double onto row
                if (ss.peek() == ',') {
                    row.push_back(stod(numstr));
                    numstr = "";
                    ss.ignore();
                }

            }
            
            // lines are not terminated by commas
            // therefore we missed pushing this numstr onto row
            row.push_back(stod(numstr));
            
            // push row into D vector
            D.push_back(row);
        }
        myfile.close();
    }

    else cout << "Unable to open file";

    return D;

}

// Euclidean distance computing for two arbitrarily long vectors of doubles. 
//
// IN: Two vectors **of equal length!** containing doubles x, y
//
// OUT: A double corresponding to the distance between x and y 
double distance(vector<double> &x, vector<double> &y) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) { sum += pow(y[i] - x[i], 2); }
    return sqrt(sum);
}

// Returns the closest center to a single datapoint
// I.e., calcates distance(dp, c) for each c in C,
// and returns the index of the c with the smallest distance
//
// IN: A vector of doubles "D points" and a vector of double vectors "centers"
//
// OUT: A single integers corresponding to an index in C
int assign_center(vector<double> &dp, vector< vector<double>> &C) {
    int I;
    double min_dist;
    for (int i = 0; i < C.size(); i++) {
        double dist = distance(dp, C[i]);
        if ( i == 0 || dist < min_dist ) { min_dist = dist; I = i; }
    }
    return I;
}

// Returns the centroid of a set of datapoints
// Does this by adding together all the datapoints and dividing by the numer of datapoints
//
// IN: A vector of double vectors "datapoints"
//
// OUT: A single double vector "centroid"
vector<double> compute_centroid(vector< vector<double>> &dps) {
    vector<double> sum = dps[0];
    for (int i = 1; i < dps.size(); i++) {
        for (int j = 0; j < dps[i].size(); j++) {
            sum[j] += dps[i][j];
        }
    }
    for (int j = 0; j < sum.size(); j++) { sum[j] /= dps.size(); }
    return sum;
}

// Convert a vector of doubles to a string. Debugging/Testing helper function.
// Use the actual debugger for anything in-depth!
//
// IN: A vector of doubles
//
// OUT: A string of all those doubles joined with ', '
string print_dp(vector<double> dp) {
    string str = "";
    for (int i = 0; i < dp.size(); i++) { str += to_string(dp[i]) + ", "; }
    return str;
}

// Convert a group of double vectors to a string. Debugging/Testing helper function.
// Use the actual debugger for anything in-depth!
//
// IN: A vector of double vectors
//
// OUT: A string of print_dp()'d vectors joined with '\n'
string print_grp(vector<vector<double>>& grp) {
    string str = "";
    for (int i = 0; i < grp.size(); i++) { str += print_dp(grp[i]) + "\n"; }
    return str;
}

// Generate starting centers for K-Means Clustering
//      - Method 1: Generate k random numbers in range of each dimension in the datapoint space
//      - Method 2: Use average of datapoints to make educated guesses..?
//
// IN: A 2D vector of doubles: the data
//
// OUT: A 2D vector of doubles: centers
//      - k vectors of size equal to the num of dims in datapoint space
vector< vector<double>> generate_centers(vector< vector<double>> &D, int K) {

    // -------------- Get range for each dimension in the datapoint space --------------

    // datapoints are vectors in R^n
    // i.e., D = { [a_1, a_2, ..., a_i, ..., a_n], [b_1, b_2, ..., b_i, ..., b_n], ...} 
    // find the max and the min for each dimension i
    vector <double> maxes;
    vector <double> mins;
    for (int d = 0; d < D.size(); d++) {        // iterate over datapoints
        for (int i = 0; i < D[d].size(); i++) { // iterate over elements

            // populate maxes, mins with first datapoint's vals
            if (d == 0) { maxes.push_back(D[d][i]); mins.push_back(D[d][i]); }

            // otherwise, check if new maxes or mins are found
            else {
                if (maxes[i] < D[d][i]) { maxes[i] = D[d][i]; }
                if (mins[i] > D[d][i]) { mins[i] = D[d][i]; }
            }
        }
    }

    // -------------- Generate and return k random centers in ranges found --------------

    vector< vector<double>> C;
    for (int k = 0; k < K; k++) {
        vector<double> c;
        for (int i = 0; i < D[0].size(); i++) {
            srand(time(NULL) + k + i); // re-initialize random seed
            c.push_back(double(rand() % int(maxes[i]) + int(mins[i])));
        }
        C.push_back(c);
    }

    return C;
}

// Generate labels for data via K-Means Clustering
//
// IN: A 2D vector of doubles: the data D
//
// OUT: A vector of ints: the labels for each datapoint in D
//      - same length as first dimension of D
vector<int> KMeansCluster(vector< vector<double>>& D) {

    // -------------- Setting centers up --------------

    // Generate K centers, "C"
    int K = 3;
    vector<vector<double>> C = generate_centers(D,K);

    cout << "centers:\n" << print_grp(C) << "\n";

    // Allocate space for our datapoint centers, "dp_C" so we can immediately index into it
    // We need both C and dp_C because we reassign one using the other
    vector< vector<double>> dp_C;
    for (int i = 0; i < D.size(); i++) { dp_C.push_back({ 0.0 }); }


    // -------------- Cluster until convergence --------------

    // indexes of convergence
    int eps = 1; // distance
    int N = 3; // consecutive count
    int v = 0; // current count

    while (v < N) {
        // assign a center to each datapoint
        for (int i = 0; i < D.size(); i++) {
            dp_C[i] = C[assign_center(D[i], C)];
        }

        // partition datapoints based on their assigned centers
        // i.e., define new clusters based on current centers
        vector<vector<vector<double>>> dp_P;
        for (int c = 0; c < C.size(); c++) {
            vector<vector<double>> interim;
            for (int d = 0; d < D.size(); d++) {
                if (dp_C[d] == C[c]) { interim.push_back(D[d]); }
            }
            dp_P.push_back(interim);
        }

        // recompute centers based on datapoint partitions
        // i.e., define new centers based on current clusters
        vector<vector<double>> C_prev(C.begin(), C.end()); // save before changing for convergence calculation
        for (int p = 0; p < dp_P.size(); p++) {
            // it is possible for a partition to be empty
            // this happens when there are no datapoints for which that partition's center is the closest
            if (!dp_P[p].empty()) { C[p] = compute_centroid(dp_P[p]); }
        }
        cout << "new centers:\n" << print_grp(C) << "\n";

        // check for convergence of centers:
        // if all centers do not change more than epsilon for N iterations in a row, they have converged.
        bool pass = true;
        for (int c = 0; c < C.size(); c++) {
            if (eps < distance(C[c], C_prev[c])) { pass = false; break; }
        }
        if (pass) { v++; }

    }

    cout << "converged\n";

    // -------------- Gather and return labels --------------

    // the label for a datapoint is the index of its assigned center in C
    vector<int> lbls;
    for (int d = 0; d < D.size(); d++) {
        lbls.push_back(find(C.begin(), C.end(), dp_C[d]) - C.begin());
    }

    return lbls;
}

int main()
{
    // -------------- Setting data up --------------

    // point string to file name
    string filename = "data/1.csv";

    // convert file data to local data structure "D"
    vector<vector<double>> D = get_D(filename);

    cout << "data:\n" << print_grp(D) << "\n";


    // -------------- Get and display labels for data using K-Means Clustering --------------

    // Get
    vector<int> lbls = KMeansCluster(D);

    // Display
    cout << "\npredictions:\n";
    for (int d = 0; d < D.size(); d++) { cout << print_dp(D[d]) << ": " << lbls[d] << "\n"; }

}
