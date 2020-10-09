// KMeansCluster.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <string>
#include <iostream>
#include <fstream>
#include "vector"
#include <sstream>
using namespace std;


// IN: name of a csv file with data in the following format:
// X,Y,L
// where X,Y are ints >=0 and L is in [A-Z]
//
// OUT: a 2-d vector of ints representing the data in the csv file
vector<vector<int>> get_data(string filename)
{
    // initialize variables
    vector<vector<int>> data;
    string line;
    ifstream myfile(filename);

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            vector<int> row;

            // Create a stringstream of the current line
            std::stringstream ss(line);

            // Keep track of the current column index
            int colIdx = 0;

            // Extract each integer
            char c;
            int cur_val;
            int prev_val = 0;
            int tot_val = 0;
            int dec = 1;
            while (ss >> c) {
                // read in chars as ints
                cur_val = c - 48;

                // If the next token is a comma, ignore it and move on
                if (ss.peek() == ',') {
                    ss.ignore();

                    // place int val in current row
                    tot_val += prev_val*dec + cur_val;
                    row.push_back(tot_val);
                    tot_val = 0;
                    prev_val = 0;
                    dec = 1;
                }
                else {
                    // otherwise, multiple the previous val by current decimal place and add it to total
                    tot_val += prev_val * dec + cur_val;
                    prev_val = cur_val;
                    dec *= 10;
                }

            }
            
            // push row into data vector
            row.push_back(tot_val);
            data.push_back(row);
        }
        myfile.close();
    }

    else cout << "Unable to open file";

    return data;

}

int main()
{
    cout << "Hello! I want to be a K-Means Clustering Algorithm When I Grow Up!\n";

    string filename = "data/1.csv";

    vector<vector<int>> data = get_data(filename);

    cout << "Here is the data from " << filename << "\n";

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            cout << data[i][j] << ",";
        }
        cout << "\n";
    }

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
