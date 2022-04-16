/*
    Name: Gaddam Yogesh
    Roll No: 19MA20015
    Lab No: 6
*/

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
typedef vector<float> VF;
typedef vector<VF> VVF;
typedef vector<string> VS;

int n, m;

// Function to convert an integer to a sstring
string convert_to_string(int num);
// Function to take input fron the user and add it to the table
int take_input(VVF &P, VF &C, VVF &CBT, VVF &B, VVF &b, VS &basic, VS &nbasic);
// Function to handle all the cases for the menu driven program
void solve(VVF P, VF C, VVF CBT, VVF B, VVF b, VS basic, VS nbasic);
// Function to find the inverse of a matrix using Product Form of Inverse
VVF PFI(VVF B, VVF Binv, VVF Bnew, int col);
// function to find the product of two matrices
VVF matrix_multiplication(VVF A, VVF B);
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
// Function to find the inverse of a matrix using Product Form of Inverse
VVF PFI(VVF B, VVF Binv, VVF Bnew, int col)
{
    int dim = B.size();
    VVF c(dim, VF(1));
    for (int i = 0; i < dim; i++)
    {
        c[i][0] = Bnew[i][col];
    }
    VVF e = matrix_multiplication(Binv, c);
    float pivot = e[col][0];
    for (int i = 0; i < dim; i++)
        e[i][0] /= (-1 * pivot);
    e[col][0] = 1 / pivot;
    VVF E(dim, VF(dim, 0));
    for (int i = 0; i < dim; i++)
        E[i][i] = 1;
    for (int i = 0; i < dim; i++)
        E[i][col] = e[i][0];
    return matrix_multiplication(E, Binv);
}
// function to find the product of two matrices
VVF matrix_multiplication(VVF A, VVF B)
{
    int l = A.size();
    int r = A[0].size();
    if (r != B.size())
    {
        cout << "[ERROR!!] Matrix multiplication not feasible\n";
        return VVF(0);
    }
    int n = B[0].size();
    VVF res(l, VF(n));
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < n; j++)
        {
            float sum = 0;
            for (int k = 0; k < r; k++)
            {
                sum += (A[i][k] * B[k][j]);
            }
            res[i][j] = sum;
        }
    }
    return res;
}
// Function to take input fron the user and add it to the table
int take_input(VVF &P, VF &C, VVF &CBT, VVF &B, VVF &b, VS &basic, VS &nbasic)
{
    cout << "Enter the number of variables: ";
    cin >> n;
    cout << "Enter the number of inequalities: ";
    cin >> m;
    B.resize(m, VF(m, 0));
    for (int i = 0; i < m; i++)
    {
        B[i][i] = 1;
    }
    CBT.resize(1, VF(m, 0));
    P.resize(m, VF(n, 0));
    b.resize(m, VF(1, 0));
    string type;
    cout << "Provide the objective (min/max): ";
    cin >> type;
    if (type != "min" && type != "max")
    {
        cout << "INVALID type!!";
        return -1;
    }
    float temp;
    float SOL = 0;
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++)
    {
        cout << "coefficient of x" << i + 1 << " : ";
        cin >> temp;
        C.push_back(temp);
        nbasic.push_back("x" + convert_to_string(i + 1));
    }

    int slack = 0, art = 0;
    string sym;
    cout << "\nProvide the conditions :\n\n";
    for (int i = 0; i < m; i++)
    {
        cout << "Condition #" << i + 1 << ":\n";
        for (int j = 0; j < n; j++)
        {
            cout << "coefficient of x" << j + 1 << " : ";
            cin >> P[i][j];
        }
        cout << "inequality (>=, = , <=): ";
        cin >> sym;
        cout << "RHS: ";
        cin >> b[i][0];
        if (sym == ">=")
        {
            b[i][0] *= -1;
            for (int j = 0; j < n; j++)
            {
                P[i][j] *= -1;
            }
            basic.push_back("s" + convert_to_string(slack + 1));
            slack++;
        }
        else if (sym == "<=")
        {
            basic.push_back("s" + convert_to_string(slack + 1));
            slack++;
        }
        else if (sym == "=")
        {
            basic.push_back("a" + convert_to_string(art + 1));
            art++;
        }
        else
        {
            cout << "INVALID SYMBOL!!!\n";
            return -1;
        }
    }
    return 0;
}
// Function to implement the revised simplex method
int revised_simplex(VVF &P, VF &C, VVF &CBT, VVF &B, VVF &b, VVF &Binv, VS &basic, VS &nbasic, bool altopt)
{
    VVF Y = matrix_multiplication(CBT, Binv);
    int col = 0;
    float min = INT_MAX;
    for (int j = 0; j < n; j++)
    {
        float sum = 0;
        for (int i = 0; i < m; i++)
        {
            sum += (Y[0][i] * P[i][j]);
        }
        sum -= C[j];
        if (sum < min)
        {
            min = sum;
            col = j;
        }
    }
    if (!altopt && min == 0)
        return 2;
    if (min > 0)
        return 0; // optimal solution is reached
    // Step 2
    VVF XB = matrix_multiplication(Binv, b);
    VVF Pcol(m, VF(1));
    for (int i = 0; i < m; i++)
        Pcol[i][0] = P[i][col];
    VVF alpha = matrix_multiplication(Binv, Pcol);
    min = INT_MAX;
    int row = -1;
    for (int i = 0; i < m; i++)
    {
        if (alpha[i][0] > 0 && XB[i][0] >= 0)
        {
            float ratio = XB[i][0] / alpha[i][0];
            if (ratio < min)
            {
                min = ratio;
                row = i;
            }
        }
    }
    if (row == -1)
    {
        return -1; // unbounded
    }
    swap(basic[row], nbasic[col]);
    VVF Bnew = B;
    for (int i = 0; i < m; i++)
    {
        Bnew[i][row] = P[i][col];
    }
    Binv = PFI(B, Binv, Bnew, row);
    swap(C[col], CBT[0][row]);
    for (int i = 0; i < m; i++)
    {
        swap(B[i][row], P[i][col]);
    }
    return 1;
}
// Function to handle all the cases for the menu driven program
void solve(VVF P, VF C, VVF CBT, VVF B, VVF b, VS basic, VS nbasic)
{
    VVF Binv(m, VF(m, 0));
    for (int i = 0; i < m; i++)
        Binv[i][i] = 1;

    int res = 1;
    while (res != 0 && res != 2)
    {
        res = revised_simplex(P, C, CBT, B, b, Binv, basic, nbasic, false);
        if (res == -1)
            break;
    }
    if (res == 2 || res == 0)
    {
        VVF XB = matrix_multiplication(Binv, b);
        for (int i = 0; i < m; i++)
        {
            if (basic[i][0] == 'a' && XB[i][0] != 0)
            {
                cout << "infeasible solution.\n";
                return;
            }
        }
        cout << "Optimal Solution: " << matrix_multiplication(CBT, XB)[0][0] << nline;
        for (int i = 0; i < n; i++)
        {
            if (nbasic[i][0] == 'x')
                cout << nbasic[i] << "=0 ";
        }
        for (int i = 0; i < m; i++)
        {
            if (basic[i][0] == 'x')
                cout << basic[i] << "=" << XB[i][0] << " ";
        }
        cout << nline;
        if (res != 0)
        {
            res = revised_simplex(P, C, CBT, B, b, Binv, basic, nbasic, true);
            VVF XB = matrix_multiplication(Binv, b);
            for (int i = 0; i < m; i++)
            {
                if (basic[i][0] == 'a' && XB[i][0] != 0)
                {
                    cout << "infeasible solution.\n";
                    return;
                }
            }
            cout << "Alternate Optimal Solution:\n";
            for (int i = 0; i < n; i++)
            {
                if (nbasic[i][0] == 'x')
                    cout << nbasic[i] << "=0 ";
            }
            for (int i = 0; i < m; i++)
            {
                if (basic[i][0] == 'x')
                    cout << basic[i] << "=" << XB[i][0] << " ";
            }
            cout << nline;
        }
    }
    else
        cout << "unbounded solution";
}

int main()
{
    // taking user input
    VVF P;
    VF C;
    VVF b;
    VVF CBT;
    VVF B;
    VS basic, nbasic;
    if (take_input(P, C, CBT, B, b, basic, nbasic) == -1)
        return -1;
    solve(P, C, CBT, B, b, basic, nbasic);
}
