/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
    TRANSPORTATION PROBLEM
    ---------------------------
    Modi Method For Optimizing
*/

/*
    LOGIC :-
    ---------
    Occupied Cell - cells in which initial basic feasible solutions are present
    Non Occupied Cells - cells withot initial basic feasible solution

    Degneracy condition - 
    Degenrate soln are converted to Non Degenerate solution
    Modi method applied to only Non degenrate solutions
    When :-
    (Number of occupied cells = Number of Rows + Number of columns - 1 = n + m -1)
    Then that initial basic feasible solution is NonDegenrate Solution

    To conver NonDegenreate toDegenrate initial BFS
    when Occupied Cells < n + m - 1 -> Degenerate solution
    Choose the lease cost among non occupied cell and assign a small quantity.

    From the Initial BFS we will find the new values of Supply and Demand
    costValue = u[i] + v[j]
    for all the cells in BFS

    Since Number of Variables > Number of Equations
    Assign the variable occuring most as 0
    Find all other values

    Find the Oppurtunity cost Table
    If all values in Oppurtunity Cost Table > 0
    Then We have achiever Optimal Solution 
    else repeat

    Inputs - according to Question given :-
    Input 1 - 
    ---------
    3 4 1
    2 3 5 1
    7 3 4 6
    4 1 7 2
    8 10 20
    6 8 9 15


    Input 2 -
    ----------
    3 4 1
    19 30 50 10
    70 30 40 60
    40 8 70 20
    7 9 18
    5 8 7 14


    Check out input and output of the given question at bottom


*/

#include <bits/stdc++.h>
using namespace std;

// global variables to store the size of the transportation problem table
int n,m;


// function to convert a float into string
string to_String(float num);
// function to print the transportation table
void display(vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<vector<float>> &mul_mat);
// function to balance the transportation by adding a dummy column or row
void balance_and_join(vector<vector<float>> &costArray, vector<float> &stock, vector<float> &demand);
// function to take input from user
void take_input(vector<vector<float>> &costArray);
// function to find the minimum element in the transportation table
pair<int,int> find_min(vector<vector<float>> &costArray, vector<bool> &row, vector<bool> &col);
// north-west corner rule for phase I
void NWCR(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, float &total);
// Least common method for phase I
void MMM(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, float &total);
// function to find penalties for VAM method
void find_penalties(vector<vector<float>> &costArray, vector<bool> &row, vector<bool> &col, vector<float> &Prow, vector<float> &Pcol);
// Vogel's Approximation Method for phase I
void VAM(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, vector<bool> &row, vector<bool> &col, float &total);
// find the values of u recursively
void find_U(int row, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v);
// find the values of v recursively
void find_V(int col, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v);
// finding the entering variable i.e. max element to be minimised
void find_entering_var(int &x, int &y, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v);
// function to find a loop in the transportation table
bool find_loop(bool row_search, vector<pair<int,int>> &loop, vector<vector<bool>> &basis);
// Modified Distribution (MODI) method for finding optimal solution
void MODI(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, vector<vector<bool>> &basis, float &total);
// function to convert a BFS into non-degenerate
void degenerate_handler(vector<vector<float>> &costArray, vector<vector<bool>> &basis);

// function to convert a float into string
string to_String(float num)
{
    string ans = to_string(num);
    bool flag = false;
    for (auto s : ans)
    {
        if (s == '.')
        {
            flag = true;
            break;
        }
    }
    if (flag)
    {
        string new_ans = "";
        int n = ans.length();
        int i = n - 1;
        for (; i >= 0; i--)
        {
            if (ans[i] == '.')
            {
                i -= 1;
                break;
            }
            if (ans[i] != '0')
                break;
        }
        for (; i >= 0; i--)
        {
            new_ans = ans[i] + new_ans;
        }
        return new_ans;
    }
    return ans;
}

// function to print the transportation table
void display(vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<vector<float>> &mul_mat)
{
    cout<<"__________________________________________________________________________________"<<endl;
    cout << left << setw(10) << "Src/Dest"<<" | ";
    float sum = 0;
    for (int i = 0; i < m - 1; i++)
    {
        cout << left<< setw(10) << "D[" + to_string(i + 1) + "]"<<" | ";
        sum += costArray[n - 1][i];
    }
    cout << left<< setw(10) << "" << endl;
    cout<<"__________________________________________________________________________________"<<endl;

    for (int i = 0; i < n - 1; i++)
    {
        cout << left<< setw(10) << "S[" + to_string(i + 1) + "]"<<" | ";
        for (int j = 0; j < m - 1; j++)
        {
            if (basis[i][j])
            {
                cout << left<< setw(10) << to_String(costArray[i][j]) + "(" + to_String(mul_mat[i][j]) + ")"<<" | ";
            }
            else
            {
                cout << left<< setw(10) << costArray[i][j] <<" | ";
            }
        }
        cout << left<< setw(10) << costArray[i][m - 1] << endl;
    }
    cout<<"__________________________________________________________________________________"<<endl;

    cout << left<< setw(10) << ""<<" | ";
    for (int j = 0; j < m - 1; j++)
        cout << left<< setw(10) << costArray[n - 1][j]<<" | ";
    cout << left<< setw(10) << sum << endl;
    cout<<"__________________________________________________________________________________"<<endl;

}

// // function to balance the transportation by adding a dummy column or row
// void balance_and_join(vector<vector<float>> &costArray, vector<float> &stock, vector<float> &demand)
// {
//     float sum1 = 0, sum2 = 0, sum = 0;
//     for (auto e : stock)
//         sum1 += e;
//     for (auto e : demand)
//         sum2 += e;
//     // stock is less than demand
//     if (sum1 < sum2)
//     {
//         for (int j = 0; j < m; j++)
//         {
//             costArray[n][j] = 0;
//         }
//         stock.push_back(sum2 - sum1);
//         n++;
//         sum = sum2;
//     }
//     // stock is more than demand
//     else if (sum1 > sum2)
//     {
//         for (int i = 0; i < n; i++)
//         {
//             costArray[i][m] = 0;
//         }
//         demand.push_back(sum1 - sum2);
//         m++;
//         sum = sum1;
//     }
//     n++;
//     m++;
//     for (int j = 0; j < m - 1; j++)
//     {
//         costArray[n - 1][j] = demand[j];
//     }
//     for (int i = 0; i < n - 1; i++)
//     {
//         costArray[i][m - 1] = stock[i];
//     }
//     costArray[n - 1][m - 1] = sum;
// }

// north-west corner rule for phase I
void NWCR(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, float &total)
{
    int i = 0, j = 0;
    while (i < n && j < m)
    {
        float diff = abs(costArray[n - 1][j] - costArray[i][m - 1]);
        if (costArray[n - 1][j] < costArray[i][m - 1])
        {
            total += (costArray[i][j] * costArray[n - 1][j]);
            mul_mat[i][j] = costArray[n - 1][j];
            costArray[n - 1][j] = 0;
            costArray[i][m - 1] = diff;
            j++;
        }
        else if (costArray[n - 1][j] > costArray[i][m - 1])
        {
            total += (costArray[i][j] * costArray[i][m - 1]);
            mul_mat[i][j] = costArray[i][m - 1];
            costArray[n - 1][j] = diff;
            costArray[i][m - 1] = 0;
            i++;
        }
        else
        {
            total += (costArray[i][j] * costArray[n - 1][j]);
            mul_mat[i][j] = costArray[n - 1][j];
            costArray[n - 1][j] = 0;
            costArray[i][m - 1] = 0;
            i++;
            j++;
        }
    }
}

// function to find the minimum element in the transportation table
pair<int,int> find_min(vector<vector<float>> &costArray, vector<bool> &row, vector<bool> &col)
{
    float min = INT_MAX;
    pair<int,int> res = {-1, -1};
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            if (row[i] && col[j] && costArray[i][j] <= min)
            {
                min = costArray[i][j];
                res = {i, j};
            }
        }
    }
    return res;
}

// Least common method for phase I
void MMM(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, float &total)
{
    vector<bool> row(n - 1, true), col(m - 1, true);
    pair<int,int> ind = find_min(costArray, row, col);
    while (ind.first != -1 && ind.second != -1)
    {
        float diff = abs(costArray[n - 1][ind.second] - costArray[ind.first][m - 1]);
        if (costArray[n - 1][ind.second] < costArray[ind.first][m - 1])
        {
            total += (costArray[ind.first][ind.second] * costArray[n - 1][ind.second]);
            mul_mat[ind.first][ind.second] = costArray[n - 1][ind.second];
            costArray[n - 1][ind.second] = 0;
            costArray[ind.first][m - 1] = diff;
            col[ind.second] = false;
        }
        else if (costArray[n - 1][ind.second] > costArray[ind.first][m - 1])
        {
            total += (costArray[ind.first][ind.second] * costArray[ind.first][m - 1]);
            mul_mat[ind.first][ind.second] = costArray[ind.first][m - 1];
            costArray[n - 1][ind.second] = diff;
            costArray[ind.first][m - 1] = 0;
            row[ind.first] = false;
        }
        else
        {
            total += (costArray[ind.first][ind.second] * costArray[n - 1][ind.second]);
            mul_mat[ind.first][ind.second] = costArray[n - 1][ind.second];
            costArray[n - 1][ind.second] = 0;
            costArray[ind.first][m - 1] = 0;
            row[ind.first] = false;
            col[ind.second] = false;
        }
        ind = find_min(costArray, row, col);
    }
}

// function to find penalties for VAM method
void find_penalties(vector<vector<float>> &costArray, vector<bool> &row, vector<bool> &col, vector<float> &Prow, vector<float> &Pcol)
{
    for (int i = 0; i < n - 1; i++)
    {
        vector<float> temp;
        for (int j = 0; j < m - 1; j++)
        {
            if (row[i] && col[j])
            {
                temp.push_back(costArray[i][j]);
            }
        }
        sort(temp.begin(), temp.end());
        int len = temp.size();
        if (len == 0)
            Prow[i] = INT_MIN;
        else if (len == 1)
            Prow[i] = temp[0];
        else
            Prow[i] = temp[1] - temp[0];
        temp.clear();
    }
    for (int j = 0; j < m - 1; j++)
    {
        vector<float> temp;
        for (int i = 0; i < n - 1; i++)
        {
            if (row[i] && col[j])
            {
                temp.push_back(costArray[i][j]);
            }
        }
        sort(temp.begin(), temp.end());
        int len = temp.size();
        if (len == 0)
            Prow[j] = INT_MIN;
        else if (len == 1)
            Pcol[j] = temp[0];
        else
            Pcol[j] = temp[1] - temp[0];
        temp.clear();
    }
}

// Vogel's Approximation Method for phase I
void VAM(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, vector<bool> &row, vector<bool> &col, float &total)
{
    vector<float> Prow(n - 1, -1), Pcol(m - 1, -1);
    find_penalties(costArray, row, col, Prow, Pcol);
    int row_max = max_element(Prow.begin(), Prow.end()) - Prow.begin();
    int col_max = max_element(Pcol.begin(), Pcol.end()) - Pcol.begin();
    int i = -1, j = -1;
    // finding the pivot element
    if (Prow[row_max] > Pcol[col_max])
    {
        i = row_max;
        float min = INT_MAX;
        for (int tmp = 0; tmp < m - 1; tmp++)
        {
            if (row[row_max] && col[tmp] && costArray[row_max][tmp] <= min)
            {
                min = costArray[row_max][tmp];
                j = tmp;
            }
        }
    }
    else
    {
        j = col_max;
        float min = INT_MAX;
        for (int tmp = 0; tmp < n - 1; tmp++)
        {
            if (row[tmp] && costArray[tmp][col_max] <= min)
            {
                min = costArray[tmp][col_max];
                i = tmp;
            }
        }
    }
    if (i == -1 || j == -1)
        return;
    // delivering from the source to destination
    float diff = abs(costArray[n - 1][j] - costArray[i][m - 1]);
    if (costArray[n - 1][j] < costArray[i][m - 1])
    {
        total += (costArray[i][j] * costArray[n - 1][j]);
        mul_mat[i][j] = costArray[n - 1][j];
        costArray[n - 1][j] = 0;
        costArray[i][m - 1] = diff;
        col[j] = false;
    }
    else if (costArray[n - 1][j] > costArray[i][m - 1])
    {
        total += (costArray[i][j] * costArray[i][m - 1]);
        mul_mat[i][j] = costArray[i][m - 1];
        costArray[n - 1][j] = diff;
        costArray[i][m - 1] = 0;
        row[i] = false;
    }
    else
    {
        total += (costArray[i][j] * costArray[n - 1][j]);
        mul_mat[i][j] = costArray[n - 1][j];
        costArray[n - 1][j] = 0;
        costArray[i][m - 1] = 0;
        col[j] = false;
        row[i] = false;
    }
    VAM(costArray, mul_mat, row, col, total);
}

// find the values of u recursively
void find_U(int col, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v)
{
    for (int i = 0; i < n - 1; i++)
    {
        if (basis[i][col] && u[i] == INT_MIN)
        {
            u[i] = costArray[i][col] - v[col];
            find_V(i, costArray, basis, u, v);
        }
    }
}

// find the values of v recursively
void find_V(int row, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v)
{
    for (int j = 0; j < m - 1; j++)
    {
        if (basis[row][j] && v[j] == INT_MIN)
        {
            v[j] = costArray[row][j] - u[row];
            find_U(j, costArray, basis, u, v);
        }
    }
}

// finding the entering variable i.e. max element to be minimised
void find_entering_var(int &x, int &y, vector<vector<float>> &costArray, vector<vector<bool>> &basis, vector<float> &u, vector<float> &v)
{
    int max = 0;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            float diff = u[i] + v[j] - costArray[i][j];
            if (!basis[i][j] && diff > max)
            {
                max = diff;
                x = i;
                y = j;
            }
        }
    }
}

// function to find a loop in the transportation table
bool find_loop(bool row_search, vector<pair<int,int>> &loop, vector<vector<bool>> &basis)
{
    if (row_search)
    {
        auto pt = loop.back();
        for (int j = 0; j < m - 1; j++)
        {
            if (basis[pt.first][j] && j != pt.second)
            {
                loop.push_back({pt.first, j});
                if (j == loop.front().second)
                    return true;
                if (find_loop(false, loop, basis))
                    return true;
                else
                    loop.pop_back();
            }
        }
        return false;
    }
    else
    {
        auto pt = loop.back();
        for (int i = 0; i < n - 1; i++)
        {
            if (basis[i][pt.second] && i != pt.first)
            {
                loop.push_back({i, pt.second});
                if (i == loop.front().first)
                    return true;
                if (find_loop(true, loop, basis))
                    return true;
                else
                    loop.pop_back();
            }
        }
        return false;
    }
}

// Modified Distribution (MODI) method for finding optimal solution
void MODI(vector<vector<float>> &costArray, vector<vector<float>> &mul_mat, vector<vector<bool>> &basis, float &total)
{
    vector<float> u(n - 1, INT_MIN), v(m - 1, INT_MIN);
    u[0] = 0;
    // finding the values of all u and v
    find_V(0, costArray, basis, u, v);
    // finding the entering variable
    int x = -1, y = -1;
    find_entering_var(x, y, costArray, basis, u, v);
    if (x == -1 || y == -1)
        return;
    // finding if a loop exists
    vector<pair<int,int>> loop;
    loop.push_back({x, y});
    find_loop(true, loop, basis);
    for(pair<int,int> x : loop) {
        cout<<x.first<<" "<<x.second<<endl;
    }
    int len = loop.size();
    float min = INT_MAX;
    for (int i = 1; i < len; i += 2)
    {
        int a = loop[i].first;
        int b = loop[i].second;
        if (mul_mat[a][b] < min)
        {
            min = mul_mat[a][b];
            x = a;
            y = b;
        }
    }
    for (int i = 0; i < len; i++)
    {
        int a = loop[i].first;
        int b = loop[i].second;
        if (i % 2 == 0)
            mul_mat[a][b] += min;
        else
            mul_mat[a][b] -= min;
    }
    basis[loop[0].first][loop[0].second] = true;
    basis[x][y] = false;
    float sum = 0;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            if (basis[i][j])
                sum += (mul_mat[i][j] * costArray[i][j]);
        }
    }
    total = sum;
    MODI(costArray, mul_mat, basis, total);
}

// function to convert a BFS into non-degenerate
void degenerate_handler(vector<vector<float>> &costArray, vector<vector<bool>> &basis)
{
    multimap<float, pair<int,int>, greater<float>> elts;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            if (!basis[i][j])
            {
                elts.insert({costArray[i][j], {i, j}});
            }
        }
    }
    int len = elts.size();
    auto it = elts.rbegin();
    for (; it != elts.rend(); it++)
    {
        vector<pair<int,int>> loop;
        loop.push_back(it->second);
        find_loop(true, loop, basis);
        if (loop.size() == 1)
            break;
    }
    cout << "Degeneracy handled by adding element at (";
    cout << it->second.first << ", " << it->second.second << ") to the basis.\n";
    basis[it->second.first][it->second.second] = true;
}


int main()
{
    vector<vector<float>> costArray;
    float total = 0;
    
    // Taking Input - Supply Sources & Dest
    cout<<"Enter number of Supply Sources : ";
    cin>>n;
    cout<<"Enter number of Destinaton Demands : ";
    cin>>m;

    // Taking Input - Cost coffecient array
    costArray.resize(n + 2, vector<float>(m + 2));
    cout<<"Enter the Cost Values :-"<<endl;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            cout<<"Source : S["<<i<<"] -> Destination : D["<<j<<"] : ";
            cin>>costArray[i][j];
        }
    }

    // Taking Input - stock array
    vector<float> stock(n);
    cout<<"Enter Values of Supply from :- "<<endl;
    for(int i=1; i<=n; i++) {
        cout<<"Source ["<<i<<"] : ";
        cin>>stock[i];
    }

    // Taking Input - demand array
    vector<float> demand(m);
    cout<<"Enter Value of Demand for :- "<<endl;
    for(int j=1; j<=m; j++) {
        cout<<"Destination ["<<j<<"] : ";
        cin>>demand[j];
    }


    vector<vector<float>> mul_mat(n, vector<float>(m, 0));
    vector<vector<bool>> basis(n - 1, vector<bool>(m - 1, false));


    // ------------------------------------------------------------
    // Displaying the Initial Transportation Table
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Initial Transportation Table :-"<<endl;    
    display(costArray, basis, mul_mat);
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    // ------------------------------------------------------------    
    cout << "PHASE - I"<<endl;
    cout << "------------------------------------------------------------";


    // Option to calculate the Method to be used to Find Initial Basic Feasible Solution
    int option;
    cout<<"Enter 1 for Least Cost Method | 2 for North West Corner Methos | 3 for Vogel's Approximation Method"<<endl;
    cin >> option;
    
    
    
    
    // Performing phase 1 for BFS
    cout << "\nBFS of the given transportation problem:\n";
    if (option == 1)
    {
        /// Matrix Minima Method / Least Cost Method to find initial BFS
        MMM(costArray, mul_mat, total);
        cout << " Using Matrix Minima Method => Cost =  ";
        cout << total << endl;
    }
    else if (option == 2)
    {
        // North West Corner Rule to find initial BFS
        NWCR(costArray, mul_mat, total);
        cout << " Using North-West Corner Rule => Cost = ";
        cout << total << endl;
    }
    else if (option == 3)
    {
        vector<bool> row(n - 1, true), col(m - 1, true);
        // Vogel's Approximation Method to find initial BFS
        VAM(costArray, mul_mat, row, col, total);
        cout << " Using Vogel's Approximation Method => Cost = ";
        cout << total << endl;
    }
    else
    {
        cout << "INVALID option.\n";
        return 0;
    }
    cout << "\n------------------------------------------------------------\n";
    // modifying the basis matrix
    int cnt = 0;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            if (mul_mat[i][j] > 0)
            {
                basis[i][j] = true;
                cnt++;
            }
            else
                basis[i][j] = false;
        }
    }
    // printing the obtained degenerate solution
    cout << "\nBasic Feasible Solution:\n\n";
    display(costArray, basis, mul_mat);
    cout << endl;
    // checking if the BFS is degenerate or not
    while (cnt != n + m - 3)
    {
        cout << "------------------------------------------------------------\n";
        cout << "Degeneracy observed in the BFS of the transportation problem.";
        cout << "\n------------------------------------------------------------\n";
        // adding another variable to the basis to make the table non-degenerate
        degenerate_handler(costArray, basis);
        cout << endl;
        display(costArray, basis, mul_mat);
        cnt = 0;
        for (int i = 0; i < n - 1; i++)
            for (int j = 0; j < m - 1; j++)
                if (basis[i][j])
                    cnt++;
    }
    cout << "____________________________________________________________\n\n";
    // Performing phase 2 for optimal solution
    cout << "PHASE - II";
    cout << "\n------------------------------------------------------------\n";
    cout << "Optimal solution by using MODI method:";
    // using MODI method for phase 2
    MODI(costArray, mul_mat, basis, total);
    cout << " Cost = " << total << endl;
    cout << "------------------------------------------------------------\n";
    display(costArray, basis, mul_mat);
    cout << "____________________________________________________________\n";
}