/*
    Name: Gaddam Yogesh
    Roll No: 19MA20015
    Lab No: 8
    Q1 cutting plane method
*/
// program to perform cutting plane method to find integer solution of an LPP
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
#define nline "\n"

using namespace std;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef vector<string> VS;

int n, m; // Global variables to store the dimensions of the simplex table

// Function to convert an integer to a string
string convert_to_string(int num);
// function to roundoff a number upto 4 places after decimal
double round_off(double num);
// Function to print the simplex table
void print_table(VVD &coeff_mat, VS &basic, VS &nbasic);
// Function to perform an iteration on the simplex table
int simplex_method(VVD &coeff_mat, VS &basic, VS &nbasic);
// Function to perform dual simplex method
int Dual_simplex_method(VVD &coeff_mat, VS &basic, VS &nbasic);
// Function to add surplus variables to the simplex table
void put_surplus_var(VVD &coeff_mat, VS &nbasic, VD &ObjFn, VD &artObjFn, int &surplus, int index);
// Function to add artificial variables to the simplex table
void put_artificial_var(VVD &coeff_mat, VS &basic, VD &artObjFn, VD &RHS, double &SOL, int &art);
// Function to take input fron the user and add it to the table
int take_input(VVD &coeff_mat, VS &basic, VS &nbasic);
// Function to handle all the cases for the menu driven program
int get_optimal(VVD &coeff_mat, VS &basic, VS &nbasic, int &curr_iter);
// Funtion to Implement the cutting plane method to get an integer solution
int cutting_plane_method(VVD coeff_mat, VS basic, VS nbasic, int &curr_iter, bool optimal);
// Function to solve the given LPP
int solve(VVD coeff_mat, VS basic, VS nbasic, int iter, bool optimal);
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
// function to roundoff a number upto 4 places after decimal
double round_off(double num)
{
    return round(num * 10000) / 10000;
}
// Function to print the simplex table
void print_table(VVD &coeff_mat, VS &basic, VS &nbasic)
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

// Function to perform an iteration on the simplex table
// returns 1 if iteration is successful
// returns 0 if optimal solution is reached
// returns -1 if unbounded
// returns -2 if infeasible
int simplex_method(VVD &coeff_mat, VS &basic, VS &nbasic)
{
    // checking for the largest negative number to get the required column
    double min = INT_MAX;
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
        for (int i = 0; i < n - 1; i++)
        {
            if (basic[i][0] == 'a' && coeff_mat[i][m - 1] != 0)
                return -2;
        }
        return 0;
    }
    else
    {
        // checking for the row with the pivot element
        VD ratios;
        double min = INT_MAX;
        int row = 0;
        for (int i = 0; i < n - 1; i++)
        {
            if (coeff_mat[i][m - 1] >= 0 && coeff_mat[i][col] > 0)
            {
                double rat = coeff_mat[i][m - 1] / coeff_mat[i][col];
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
            double pivot_val = coeff_mat[row][col];
            coeff_mat[row][col] = 1 / coeff_mat[row][col];

            // s* = (ps-qr)/p
            for (int i = 0; i < n + 1; i++)
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
            for (int i = 0; i < n + 1; i++)
            {
                if (i != row)
                    coeff_mat[i][col] = (-1 * coeff_mat[i][col] / pivot_val);
            }

            // modifying the row elements
            for (int i = 0; i < m; i++)
            {
                if (i != col)
                    coeff_mat[row][i] = (coeff_mat[row][i] / pivot_val);
            }
        }
    }
    return 1;
}

// Function to perform dual simplex method
int Dual_simplex_method(VVD &coeff_mat, VS &basic, VS &nbasic)
{
    // checking for the largest negative number to get the required row
    double min = INT_MAX;
    int row = 0;
    for (int i = 0; i < n - 1; i++)
    {
        if (coeff_mat[i][m - 1] < min)
        {
            min = coeff_mat[i][m - 1];
            row = i;
        }
    }
    if (coeff_mat[row][m - 1] >= -0)
    {
        return 0;
    }
    else
    {
        // checking for the row with the pivot element
        VD ratios;
        double max = INT_MIN;
        int col = -1;
        for (int i = 0; i < m - 1; i++)
        {
            if (coeff_mat[row][i] < 0)
            {
                double rat = coeff_mat[n - 1][i] / coeff_mat[row][i];
                if (rat > 0)
                    continue;
                ratios.push_back(rat);
                if (rat > max)
                {
                    max = rat;
                    col = i;
                }
            }
        }
        // if no negative ratios found => unbounded
        if (col == -1)
        {
            cout << "Integer solution doesn't exist.\n";
            return -1;
        }
        else
        {
            // printing the pivot element
            if (true)
            {
                cout << "\nPivot element: " << coeff_mat[row][col] << nline;
            }
            // printing the incoming and outgoing variables
            if (true)
            {
                cout << "Incoming variable: " << nbasic[col] << nline;
            }
            if (true)
            {
                cout << "Outgoing variable: " << basic[row] << nline;
            }
            // swapping the non basic and basic variables
            swap(nbasic[col], basic[row]);
            double pivot_val = coeff_mat[row][col];
            coeff_mat[row][col] = (1 / coeff_mat[row][col]);

            // s* = (ps-qr)/p
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (i != row && j != col)
                    {
                        coeff_mat[i][j] = ((coeff_mat[i][j] * pivot_val - coeff_mat[i][col] * coeff_mat[row][j]) / pivot_val);
                    }
                }
            }

            // modifying the column elements
            for (int i = 0; i < n; i++)
            {
                if (i != row)
                    coeff_mat[i][col] = (-1 * coeff_mat[i][col] / pivot_val);
            }

            // modifying the row elements
            for (int i = 0; i < m; i++)
            {
                if (i != col)
                    coeff_mat[row][i] = (coeff_mat[row][i] / pivot_val);
            }
        }
    }
    return 1;
}

// Function to add surplus variables to the simplex table
void put_surplus_var(VVD &coeff_mat, VS &nbasic, VD &ObjFn, VD &artObjFn, int &surplus, int index)
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
void put_artificial_var(VVD &coeff_mat, VS &basic, VD &artObjFn, VD &RHS, double &SOL, int &art)
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
int take_input(VVD &coeff_mat, VS &basic, VS &nbasic)
{
    // number of variables and equations
    cout << "Enter the number of variables: ";
    cin >> m;
    cout << "Enter the number of inequalities: ";
    cin >> n;
    coeff_mat.resize(n + 2, VD(m + n + 1));
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
    VD ObjFn(m);
    VD artObjFn(m, 0);
    double SOL = 0;
    cout << "\nEnter the coefficients of the objective function:\n";
    for (int i = 0; i < m; i++)
    {
        cout << "coefficient of x" << i + 1 << " : ";
        cin >> ObjFn[i];
    }

    // taking inputs for the inequalitites
    VD RHS;
    string sym;
    VS ineq;
    double temp;
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

// Function to handle all the cases for the menu driven program
int get_optimal(VVD &coeff_mat, VS &basic, VS &nbasic, int &curr_iter)
{
    int res = 1;
    // Running phase 1
    while (res != 0)
    {
        cout << "Iteration #" << curr_iter << nline << nline;
        print_table(coeff_mat, basic, nbasic);
        cout << "\n-------------------------------------------------\n";
        curr_iter++;
        res = simplex_method(coeff_mat, basic, nbasic);
        if (res == -1)
        {
            cout << "unbounded\n";
            return -1;
        }
        if (res == -2)
        {
            cout << "The solution is not feasible.\n";
            return -2;
        }
    }
    for (int i = 0; i < n - 1; i++)
    {
        if (basic[i][0] == 'a')
        {
            cout << "The solution is not feasible..\n";
            return -2;
        }
    }
    for (int j = 0; j < m - 1; j++)
    {
        if (nbasic[j][0] == 'a')
        {
            for (int i = 0; i < n + 1; i++)
            {
                coeff_mat[i][j] = 0;
            }
        }
    }
    for (int i = 0; i < m; i++)
    {
        swap(coeff_mat[n - 1][i], coeff_mat[n][i]);
    }
    swap(basic[n - 1], basic[n]);
    res = 1;
    // Running phase 2
    while (res != 0)
    {
        cout << "Iteration #" << curr_iter << nline << nline;
        print_table(coeff_mat, basic, nbasic);
        cout << "\n-------------------------------------------------\n";
        res = simplex_method(coeff_mat, basic, nbasic);
        curr_iter++;
        if (res == -1)
        {
            cout << "unbounded\n";
            return -1;
        }
        if (res == -2)
        {
            cout << "The solution is not feasible...\n";
            return -2;
        }
    }
    return 0;
}

// Funtion to Implement the cutting plane method to get an integer solution
int cutting_plane_method(VVD coeff_mat, VS basic, VS nbasic, int &curr_iter, bool optimal)
{
    double max = INT_MIN;
    int row = -1;
    for (int i = 0; i < n - 1; i++)
    {
        double num;
        double fract = modf(coeff_mat[i][m - 1], &num);
        if (fract > max && fract != 0)
        {
            max = fract;
            row = i;
        }
    }
    // no non integer solution
    if (row == -1)
    {
        if (optimal)
        {
            cout << "\n-------------------------------------------------\n";
            cout << "Optimal Solution: ";
            if (basic[n - 1][0] == '-')
                cout << "Z = " << -1 * round_off(coeff_mat[n - 1][m - 1]) << nline;
            else
                cout << "Z = " << round_off(coeff_mat[n - 1][m - 1]) << nline;

            for (int i = 0; i < n - 1; i++)
                if (basic[i][0] == 'x')
                    cout << basic[i] << "=" << round_off(coeff_mat[i][m - 1]) << " ";

            for (int i = 0; i < m - 1; i++)
                if (nbasic[i][0] == 'x')
                    cout << nbasic[i] << "=" << 0 << " ";
            cout << "\n-------------------------------------------------\n";
        }
        return curr_iter;
    }
    // print the cutting plane equation
    if (true)
    {
        cout << "-------------------------------------------------\n";
        cout << "Equation of cutting plane:\n";
        cout << "(" << nbasic[0] << " " << coeff_mat[row][0] << ")";
        for (int i = 1; i < m - 1; i++)
        {
            cout << " + (" << nbasic[i] << " " << coeff_mat[row][i] << ")";
        }
        cout << " = " << coeff_mat[row][m - 1];
        cout << nline;
    }

    int cnt = 0;
    for (auto e : basic)
        if (e[0] == 'z')
            cnt++;
    for (auto e : nbasic)
        if (e[0] == 'z')
            cnt++;
    coeff_mat.resize(n + 1, VD(m));
    for (int i = 0; i < m; i++)
    {
        double num;
        coeff_mat[n][i] = coeff_mat[n - 1][i];
        coeff_mat[n - 1][i] = -1 * modf(coeff_mat[row][i], &num);
        if (coeff_mat[n - 1][i] > 0)
            coeff_mat[n - 1][i] = coeff_mat[n - 1][i] - 1;
    }
    basic.resize(n + 1);
    basic[n] = "z" + convert_to_string(cnt + 1);
    swap(basic[n - 1], basic[n]);
    n++;
    int res = 1;
    while (res != 0)
    {
        cout << "\nIteration #" << curr_iter << nline << nline;
        print_table(coeff_mat, basic, nbasic);
        res = Dual_simplex_method(coeff_mat, basic, nbasic);
        curr_iter++;
        if (res == -1)
            return curr_iter;
    }
    for (int i = 0; i < n - 1; i++)
    {
        double num = abs(round(coeff_mat[i][m - 1]) - coeff_mat[i][m - 1]);
        if (num > 0.0001 && basic[i][0] == 'x')
        {
            return cutting_plane_method(coeff_mat, basic, nbasic, curr_iter, optimal);
        }
    }
    if (optimal)
    {
        cout << "\n-------------------------------------------------\n";
        cout << "Optimal Solution: ";
        if (basic[n - 1][0] == '-')
            cout << "Z = " << -1 * round_off(coeff_mat[n - 1][m - 1]) << nline;
        else
            cout << "Z = " << round_off(coeff_mat[n - 1][m - 1]) << nline;

        for (int i = 0; i < n - 1; i++)
            if (basic[i][0] == 'x')
                cout << basic[i] << "=" << round_off(coeff_mat[i][m - 1]) << " ";

        for (int i = 0; i < m - 1; i++)
            if (nbasic[i][0] == 'x')
                cout << nbasic[i] << "=" << 0 << " ";
        cout << "\n-------------------------------------------------\n";
    }
    return curr_iter;
}

// Function to solve the given LPP
int solve(VVD coeff_mat, VS basic, VS nbasic, bool optimal)
{
    int curr_iter = 0;
    if (get_optimal(coeff_mat, basic, nbasic, curr_iter) != 0)
        return curr_iter;
    curr_iter = 0;
    cout << "Using cutting plane method:\n";
    return cutting_plane_method(coeff_mat, basic, nbasic, curr_iter, optimal);
}

int main()
{
    // taking user input
    VVD coeff_mat;
    VS basic, nbasic;
    if (take_input(coeff_mat, basic, nbasic) == -1)
        return -1;
    cout << "\n-------------------------------------------------\n";
    int max_iter = solve(coeff_mat, basic, nbasic, true) - 1;
    return 0;
}
