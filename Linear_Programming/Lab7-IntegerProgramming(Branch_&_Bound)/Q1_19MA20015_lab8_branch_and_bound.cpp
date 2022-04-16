/*
    Name: Gaddam Yogesh
    Roll No: 19MA20015
    Lab No: 8
    Q1 branch and bound method
*/
// Menu driven program to solve a LPP using Branch and Bound method for integer solution
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <random>
#include <queue>
#include <algorithm>
#include <functional>
#include <string>
#define nline "\n"
#define INT_MAX 1e9;
#define INT_MIN -1e9;

using namespace std;

typedef vector<int> VI;
typedef vector<float> VF;
typedef vector<VF> VVF;
typedef vector<string> VS;

int n, m; // Global variables to store the dimensions of the simplex table
vector<map<int, float>> X;
VF Z;


// Function to convert an integer to a sstring
string convert_to_string(int num);
// Function to get he count of a variable
int get_var_cnt(char c, VS &basic, VS &nbasic);
// Function to print the simplex table
void print_table(VVF &coeff_mat, VS &basic, VS &nbasic);
// Function to perform an iteration on the simplex table
int simplex_method(VVF &coeff_mat, VS &basic, VS &nbasic, bool p2);
// Function to add surplus variables to the simplex table
void put_surplus_var(VVF &coeff_mat, VS &nbasic, VF &ObjFn, VF &artObjFn, int &surplus, int index);
// Function to add artificial variables to the simplex table
void put_artificial_var(VVF &coeff_mat, VS &basic, VF &artObjFn, VF &RHS, float &SOL, int &art);
// Function to take input fron the user and add it to the table
int take_input(VVF &coeff_mat, VS &basic, VS &nbasic);
// Function to handle all the cases for the menu driven program
int two_phase_simplex(VVF coeff_mat, VS basic, VS nbasic);
// function to solve the LPP using Branch and Bound method
int solve(VVF coeff_mat, VS basic, VS nbasic);
// Function to print the simplex table
void print_table(VVF &coeff_mat, VS &basic, VS &nbasic)
{
    for (int i = 0; i < m - 1; i++)
    {
        cout << setw(10);
        cout << "-" + nbasic[i];
    }
    cout << setw(10) << "1";
    cout << setw(10) << "NBV/BV" << nline;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << setw(10) << coeff_mat[i][j];
        cout << setw(10) << basic[i] << nline;
    }
}
// Function to convert an integer to a sstring
string convert_to_string(int num)
{
    string res = "";
    bool flag = false;
    if (num < 0)
    {
        num *= -1;
        flag = true;
    }
    while (num != 0)
    {
        res = char((num % 10) + '0') + res;
        num /= 10;
    }
    if (flag)
        res = '-' + res;
    return res;
}

// Function to get he count of a variable
int get_var_cnt(char c, VS &basic, VS &nbasic)
{
    int cnt = 0;
    for (auto e : basic)
        if (e[0] == c)
            cnt++;
    for (auto e : nbasic)
        if (e[0] == c)
            cnt++;
    return cnt;
}

// Function to perform an iteration on the simplex table
// returns 1 if iteration is successful
// returns 0 if optimal solution is reached
// returns -1 if unbounded
// returns -2 if infeasible
int simplex_method(VVF &coeff_mat, VS &basic, VS &nbasic, bool p2)
{
    // checking for the largest negative number to get the required column
    float min = INT_MAX;
    int col = -1;
    for (int i = 0; i < m - 1; i++)
    {
        if (coeff_mat[n - 1][i] < min)
        {
            min = coeff_mat[n - 1][i];
            col = i;
        }
    }
    if (coeff_mat[n - 1][col] >= 0)
    {
        // alternate_optimal(coeff_mat, basic, nbasic);
        for (int i = 0; i < n - 1; i++)
        {
            if (basic[i][0] == 'a' && coeff_mat[i][m - 1] != 0)
            {
                return -1;
            }
        }
        return 0;
    }
    else
    {
        // checking for the row with the pivot element
        VF ratios;
        float min = INT_MAX;
        int row = 0;
        for (int i = 0; i < n - 1; i++)
        {
            if (coeff_mat[i][m - 1] >= 0 && coeff_mat[i][col] > 0)
            {
                float rat = coeff_mat[i][m - 1] / coeff_mat[i][col];
                ratios.push_back(rat);
                if (rat == min && basic[row][0] == 'a')
                    row = i;
                else if (rat < min)
                {
                    min = rat;
                    row = i;
                }
            }
        }
        // if no positive ratios found => unbounded
        if (ratios.size() == 0)
            return -1;
        else
        {
            // swapping the non basic and basic variables
            swap(nbasic[col], basic[row]);
            float pivot_val = coeff_mat[row][col];
            coeff_mat[row][col] = 1 / coeff_mat[row][col];
            // s* = (ps-qr)/p
            int len = n;
            if (p2)
                len++;
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (i != row && j != col)
                    {
                        coeff_mat[i][j] = (coeff_mat[i][j] * pivot_val - coeff_mat[i][col] * coeff_mat[row][j]) / pivot_val;
                    }
                }
            }

            // modifying the column elements
            for (int i = 0; i < len; i++)
            {
                if (i != row)
                    coeff_mat[i][col] = -1 * coeff_mat[i][col] / pivot_val;
            }

            // modifying the row elements
            for (int i = 0; i < m; i++)
            {
                if (i != col)
                    coeff_mat[row][i] = coeff_mat[row][i] / pivot_val;
            }
        }
    }
    return 1;
}

// Function to add surplus variables to the simplex table
void put_surplus_var(VVF &coeff_mat, VS &nbasic, VF &ObjFn, VF &artObjFn, int &surplus, int index)
{
    surplus++;
    for (int i = 0; i < n; i++)
    {
        if (i == index)
            coeff_mat[i][m] = -1;
        else
            coeff_mat[i][m] = 0;
    }
    ObjFn.push_back(0);
    artObjFn.push_back(0);
    nbasic.push_back("s" + convert_to_string(surplus));
    m++;
}

// Function to add artificial variables to the simplex table
void put_artificial_var(VVF &coeff_mat, VS &basic, VF &artObjFn, VF &RHS, float &SOL, int &art)
{
    art++;
    int row = basic.size();
    for (int i = 0; i < m; i++)
    {
        artObjFn[i] -= coeff_mat[row][i];
    }
    SOL -= RHS[row];
    basic.push_back("a" + convert_to_string(art));
}

// Function to take input fron the user and add it to the table
int take_input(VVF &coeff_mat, VS &basic, VS &nbasic)
{
    // number of variables and equations
    cout << "Enter the number of variables: ";
    cin >> m;
    cout << "Enter the number of inequalities: ";
    cin >> n;
    coeff_mat.resize(n + 2, VF(m + n + 1));
    // Taking the objective of the problem, minimize or maximize
    string type;
    cout << "Provide the objective (min/max): ";
    cin >> type;
    if (type != "min" && type != "max")
    {
        cout << "INVALID type!!";
        return -1;
    }

    // taking input for objective function
    VF ObjFn(m);
    VF artObjFn(m, 0);
    float SOL = 0;
    cout << "\nEnter the coefficients of the objective function:\n";
    for (int i = 0; i < m; i++)
    {
        cout << "coefficient of x" << i + 1 << " : ";
        cin >> ObjFn[i];
    }

    // taking inputs for the inequalitites
    VF RHS;
    string sym;
    VS ineq;
    float temp;
    int surp = 0;
    cout << "\nProvide the conditions :\n\n";
    for (int i = 0; i < n; i++)
    {
        cout << "Condition #" << i + 1 << ":\n";
        for (int j = 0; j < m; j++)
        {
            cout << "coefficient of x" << j + 1 << " : ";
            cin >> coeff_mat[i][j];
        }
        cout << "inequality (>=, = , <=): ";
        cin >> sym;
        if (sym == ">=")
            surp++;
        cout << "RHS: ";
        cin >> temp;
        if (temp < 0)
        {
            for (int j = 0; j < m; j++)
            {
                coeff_mat[i][j] *= -1;
            }
            temp *= -1;
            if (sym == ">=")
                sym = "<=";
            else if (sym == "<=")
                sym = ">=";
        }
        ineq.push_back(sym);
        RHS.push_back(temp);
        cout << nline;
    }

    // Printing the given input
    // printing the objective function
    cout << "Given Objective function:\n"
         << type << ": Z = ";
    for (int i = 0; i < m - 1; i++)
    {
        cout << ObjFn[i] << "(x" << i + 1 << ") + ";
        if (type == "max")
            ObjFn[i] *= -1;
    }
    cout << ObjFn[m - 1] << "(x" << m << ")" << nline;
    if (type == "max")
        ObjFn[m - 1] *= -1;
    // printing the conditions
    cout << "Given conditions are:\n";
    for (int i = 0; i < n; i++)
    {
        cout << i + 1 << ") ";
        for (int j = 0; j < m - 1; j++)
        {
            cout << coeff_mat[i][j] << "(x" << j + 1 << ") + ";
        }
        cout << coeff_mat[i][m - 1] << "(x" << m << ") " << ineq[i] << " " << RHS[i] << nline;
    }

    // storing the variable names
    for (int i = 0; i < m; i++)
        nbasic.push_back("x" + convert_to_string(i + 1));

    // introducing slack, artificial and surplus variables into the problem
    int art = 0, slack = 0, surplus = 0;
    for (int i = 0; i < ineq.size(); i++)
    {
        // For 'greater than' type, introducing a surplus and an artificial variable
        if (ineq[i] == ">=")
        {
            put_surplus_var(coeff_mat, nbasic, ObjFn, artObjFn, surplus, i);
            put_artificial_var(coeff_mat, basic, artObjFn, RHS, SOL, art);
        }
        // For 'equal to' type, introducing an artificial variable
        else if (ineq[i] == "=")
        {
            put_artificial_var(coeff_mat, basic, artObjFn, RHS, SOL, art);
        }
        // For 'less than' type, introducing a slack variable
        else if (ineq[i] == "<=")
        {
            slack++;
            basic.push_back("z" + convert_to_string(slack));
        }
        else
        {
            cout << "Looks like you have provided an INVALID inequality!!\n";
            return -1;
        }
    }
    basic.push_back("-f");
    if (type == "max")
        basic.push_back("Z");
    else
        basic.push_back("-Z");

    // combining all the components to make the final table
    for (int i = 0; i < n; i++)
    {
        coeff_mat[i][m] = RHS[i];
    }
    artObjFn.push_back(SOL);
    ObjFn.push_back(0);
    m++;
    for (int i = 0; i < m; i++)
    {
        coeff_mat[n][i] = artObjFn[i];
        coeff_mat[n + 1][i] = ObjFn[i];
    }
    n++;
    RHS.clear();
    artObjFn.clear();
    ineq.clear();
    return 0;
}

// function to solve a branch using simplex method
int solve_branch(VVF coeff_mat, VS basic, VS nbasic)
{
    int res = 1;
    while (res != 0)
    {
        cout << "\n--------------------------------------------------*\n";
        print_table(coeff_mat, basic, nbasic);
        // iterating the simplex table untill an optimal solution is reached or the solution is not feasible or unbounded
        res = simplex_method(coeff_mat, basic, nbasic, false);
        if (res == -1 || res == -2)
            return res;
    }
    // storing the obtained solution
    if (basic[n - 1][0] == '-')
        Z.push_back(-1 * coeff_mat[n - 1][m - 1]);
    else
        Z.push_back(coeff_mat[n - 1][m - 1]);

    // storing the corresponding variable values for the solution
    map<int, float> vals;
    for (int i = 0; i < n - 1; i++)
        if (basic[i][0] == 'x')
            vals[stoi(basic[i].substr(1)) - 1] = coeff_mat[i][m - 1];

    for (int i = 0; i < m - 1; i++)
        if (nbasic[i][0] == 'x')
            vals[stoi(nbasic[i].substr(1)) - 1] = 0;

    X.push_back(vals);
    return 0;
}

// Function to apply two phase simplex method for the simplex table with artificial variables
int two_phase_simplex(VVF coeff_mat, VS basic, VS nbasic)
{
    int res = 1;
    // Running phase 1
    while (res != 0)
    {
        cout << "\n--------------------------------------------------\n";
        print_table(coeff_mat, basic, nbasic);
        // iterating the simplex table untill an optimal solution is reached or the solution is not feasible or unbounded
        res = simplex_method(coeff_mat, basic, nbasic, true);
        if (res == -1 || res == -2)
            return res;
    }
    for (int i = 0; i < n - 1; i++)
    {
        // if an artificial variable is non zero then non feasible solution
        if (basic[i][0] == 'a')
            return -2;
    }
    // filling the row with artificial variable with 0s
    for (int j = 0; j < m - 1; j++)
    {
        if (nbasic[j][0] == 'a')
        {
            for (int i = 0; i < n + 1; i++)
                coeff_mat[i][j] = 0;
        }
    }
    for (int i = 0; i < m; i++)
        swap(coeff_mat[n - 1][i], coeff_mat[n][i]);

    swap(basic[n - 1], basic[n]);
    return solve_branch(coeff_mat, basic, nbasic);
}

// function to solve the LPP using Branch and Bound method
int solve(VVF coeff_mat, VS basic, VS nbasic)
{
    // taking the previously calculated optimal solution to proceed to the ext step
    int sol_cnt = Z.size();
    int des_var = X[0].size();
    bool flag = true;
    int var = -1;
    int next, prev;
    for (int j = 0; j < des_var; j++)
    {
        // looking for a basic variable which has a non integer value
        double num = abs(round(X[sol_cnt - 1][j]) - X[sol_cnt - 1][j]);
        if (num > 0.001)
        {
            // storing the value of the integer greater and less than the obtained value
            next = ceil(X[sol_cnt - 1][j]);
            prev = floor(X[sol_cnt - 1][j]);
            var = j;
            flag = false;
            break;
        }
    }
    // all values are integers
    if (var == -1)
        return 0;
    // storing the initial table in temporary variable for further reference
    VVF temp_mat = coeff_mat;

    /*
        right branch i.e integers <= the floor value
    */
    // resizing the matrix to accomodate the new inequality
    coeff_mat.resize(n + 1, VF(m));
    // shifting the objective function values to the bottom row to accomodate the new variable
    for (int i = 0; i < m; i++)
    {
        coeff_mat[n][i] = coeff_mat[n - 1][i];
        coeff_mat[n - 1][i] = 0;
    }
    // finding the number of slack variables already present
    int slack = get_var_cnt('z', basic, nbasic);
    // adding the slack variable to the basic variable list
    basic.resize(n + 1);
    basic[n] = "z" + convert_to_string(slack + 1);
    swap(basic[n - 1], basic[n]);
    // increasing the n value due to the added row
    n++;
    // entering the values for the new row i.e x <= c
    coeff_mat[n - 2][var] = 1;
    coeff_mat[n - 2][m - 1] = prev;
    // solving the newly obtained table to obtain a new solution
    if (solve_branch(coeff_mat, basic, nbasic) == 0)
        solve(coeff_mat, basic, nbasic); // continuing the process with the condition
    // returning to the original table for other branches
    swap(basic[n - 1], basic[n - 2]);
    basic.pop_back();
    // replacing the table with the initial table
    coeff_mat = temp_mat;
    n--;

    /*
        left branch i.e. integers >= the ceil value
    */
    // resizing the matrix to accomodate the new inequality
    coeff_mat.resize(n + 2, VF(m + 1));
    // shifting the objective function values to the bottom row to accomodate the new variable
    for (int i = 0; i < m + 1; i++)
    {
        coeff_mat[n][i] = 0;
        coeff_mat[n + 1][i] = 0;
    }
    // entering the values for the new row i.e x >= c
    coeff_mat[n + 1][var] = 1;
    coeff_mat[n + 1][m - 1] = next;
    // shifting the RHS values to the right most to accomodate the surplus variable
    for (int i = 0; i < n + 2; i++)
    {
        coeff_mat[i][m] = coeff_mat[i][m - 1];
        coeff_mat[i][m - 1] = 0;
    }
    coeff_mat[n + 1][m - 1] = -1;
    for (int i = 0; i < m + 1; i++)
    {
        swap(coeff_mat[n - 1][i], coeff_mat[n + 1][i]);
        coeff_mat[n][i] = -1 * coeff_mat[n - 1][i];
    }
    // getting the count of the variables to add a new one
    int art = get_var_cnt('a', basic, nbasic);
    int surp = get_var_cnt('s', basic, nbasic);
    // adding the artificial variable to the basic variable list
    basic.resize(n + 2);
    // adding the surplus variable to the basic variable list
    nbasic.push_back("s" + convert_to_string(surp + 1));
    basic[n] = "a" + convert_to_string(art + 1);
    // adding the artificial objective function to the table for two phase simplex method
    basic[n + 1] = "-f";
    swap(basic[n - 1], basic[n]);
    swap(basic[n], basic[n + 1]);
    m++;
    n++;
    // solving the newly obtained table to obtain a new solution
    if (next != 0 && two_phase_simplex(coeff_mat, basic, nbasic) == 0)
        solve(coeff_mat, basic, nbasic); // continuing the process with the condition
    // replacing the table with the initial table
    coeff_mat = temp_mat;
    // returning to the original table for other branches
    basic.pop_back();
    swap(basic[n - 1], basic[n - 2]);
    basic.pop_back();
    nbasic.pop_back();
    n--;
    m--;
    return 0;
}

int main()
{
    // taking user input
    VVF coeff_mat;
    VS basic, nbasic;
    if (take_input(coeff_mat, basic, nbasic) == -1)
        return -1;
    // obtaining an initial non integer optimal solution
    if (two_phase_simplex(coeff_mat, basic, nbasic) != 0)
        return 0;
    for (int i = 0; i < m; i++)
        coeff_mat[n - 1][i] = coeff_mat[n][i];
    swap(basic[n - 1], basic[n]);
    basic.resize(n);
    solve(coeff_mat, basic, nbasic);

    // printing all the obtained solutions, integer and non-integer solutions
    int sol_cnt = Z.size();
    int des_var = X[0].size();
    cout << "\n--------------------------------------------------\n";
    cout << "All the obtained solutions:\n\n";
    for (int i = 0; i < des_var; i++)
        cout << setw(10) << "x" + convert_to_string(i + 1);
    cout << setw(10) << "Z\n";
    for (int i = 0; i < sol_cnt; i++)
    {
        for (int j = 0; j < des_var; j++)
            cout << setw(10) << X[i][j];
        cout << setw(10) << Z[i];
        cout << nline;
    }
    cout << "--------------------------------------------------\n";
    bool flag;
    float max = INT_MIN;
    int sol_index = -1;
    for (int i = 0; i < sol_cnt; i++)
    {
        flag = true;
        for (int j = 0; j < des_var; j++)
        {
            double num = abs(round(X[i][j]) - X[i][j]);
            if (num > 0.0001)
                flag = false;
        }
        if (flag && Z[i] >= max)
        {
            max = Z[i];
            sol_index = i;
        }
    }

    if (sol_index == -1)
    {
        cout << "Optimal Solution doesnt exist.\n";
        return 0;
    }
    cout << "Optimal Solution: Z = " << Z[sol_index] << nline;
    for (int j = 0; j < des_var; j++)
    {
        cout << "x" << j + 1 << " = " << X[sol_index][j] << nline;
    }
    cout << "--------------------------------------------------\n";
    return 0;
}
