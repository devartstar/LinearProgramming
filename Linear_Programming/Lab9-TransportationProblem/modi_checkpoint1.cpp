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

int n, m;
double costArray[100][100];
double supply[100];
double demand[100];
double costArrayCopy[100][100], supplyCopy[100], demandCopy[100];
bool covered[100][100];
long long INF = 1e9;

// for Modi method
double newSupply[100], newDemand[100];
bool occupied[100][100];


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function To Display the Table for Transportation Problem
void display(double A[][100], double B[], double C[]) {
    cout<<left<<setw(10)<<""<<" | ";
    for(int j=1;j<=m;j++){
        cout<<left<<setw(10)<<"D["+to_string(j)+"]"<<" | ";
    }
    cout<<endl;
    cout<<"__________________________________________________________________"<<endl;
    for(int i=1; i<=n; i++){
        cout<<left<<setw(10)<<"s["+to_string(i)+"]"<<" | ";
        for(int j=1; j<=m; j++){
            cout<<left<<setw(10)<<A[i][j]<<" | ";
        }
        cout<<left<<setw(10)<<B[i]<<endl;
    }
    cout<<"__________________________________________________________________"<<endl;
    
    cout<<left<<setw(10)<<""<<" | ";
    for(int j=1; j<=m; j++){
        cout<<left<<setw(10)<<C[j]<<" | ";
    }
    cout<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// To check if the given Transportation Problem is Balanced or Not
// Balanced Transportation Problem => ( Total Supply = Total Demand )
bool checkIfBalanced(double S[], double D[]){
    double totalSupply = 0;
    double totalDemand = 0;
    for(int i=1; i<=n; i++){
        totalSupply += S[i];
    }
    for(int j=1; j<=m; j++){
        totalDemand += D[j];
    }
    if(abs(totalDemand - totalSupply) <= 0.001){
        return true;
    }else{
        return false;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function To Create Copy costArray, supplyArray, demandArray
void createCopy() {
    // ------------------------------------------------------------
    // creating a copy of CostArray, supply, demand
    // CostArrayCopy -> will contain Number of Unots for each source to Destination
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            costArrayCopy[i][j] = costArray[i][j];
            covered[i][j] = false;
        }
    }
    for(int i=1; i<=n; i++){
        supplyCopy[i] = supply[i];
    }
    for(int j=1; j<=m; j++){
        demandCopy[j] = demand[j];
    }
    // ------------------------------------------------------------
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// After every Iteration we need to mind the cell in Transportation Table will min Cost Values
// Returns the position of the Min Cost cell among the left cells
vector<pair<int, int>> findPosOfMin(double costArray[][100], bool covered[][100]){
    vector<pair<int,int>> allpos;
    int min = INF;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            if(covered[i][j] == false){
                if(costArray[i][j] < min){
                    allpos.clear();
                    allpos.push_back({i,j});
                    min = costArray[i][j];
                }else if(costArray[i][j] == min){
                    allpos.push_back({i,j});
                }
            }
        }
    }
    return allpos;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function To get the Basic Feasible solution for Transportation Problem using Least Cost ethod
void LeastCostMethod() {
    createCopy();
    int iterationCount = 0;
    while (iterationCount <= n+m-1) {
        vector<pair<int, int>> posArr = findPosOfMin(costArray, covered);
        if((int)posArr.size()==0){
            break;
        }
        pair<int, int> pos = {-1,-1};
        double mx = -INF;
        // now in all eleemts in alll pos maximize allocation
        for(pair<int, int> x : posArr) {
            double temp;
            if(supplyCopy[x.first] >= demandCopy[x.second]){
                temp = supplyCopy[x.first];
            }else{
                temp = demandCopy[x.second];
            }

            if(temp > mx){
                pos = x;
                mx = temp;
            }
        }
        if(pos.first == -1 && pos.second == -1){
            break;
        }

        covered[pos.first][pos.second] = true;
        if(supplyCopy[pos.first] <= demandCopy[pos.second]) {
            costArrayCopy[pos.first][pos.second] = supplyCopy[pos.first];
            for(int j=1; j<=m; j++){
                if(j==pos.second){
                    continue;
                }else{
                    if(covered[pos.first][j] == true){
                        continue;
                    }else{
                        covered[pos.first][j] = true;
                        costArrayCopy[pos.first][j] = 0;
                    }
                }
            }
            demandCopy[pos.second] -= supplyCopy[pos.first];
            supplyCopy[pos.first] = 0;
        }else{
            costArrayCopy[pos.first][pos.second] = demandCopy[pos.second];
            for(int i=1; i<=n; i++){
                if(i==pos.first){
                    continue;
                }else{
                    if(covered[i][pos.second] == true){
                        continue;
                    }else{
                        covered[i][pos.second] = true;
                        costArrayCopy[i][pos.second] = 0;
                    }
                }
            }
            supplyCopy[pos.first] -= demandCopy[pos.second];
            demandCopy[pos.second] = 0;
        }
        display(costArrayCopy, supplyCopy, demandCopy);
        iterationCount++;
    }
    display(costArrayCopy, supplyCopy, demandCopy);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function To get the Basic Feasible solution for Transportation Problem using North West Corner Method
void NorthWestCornerMethod() {
    createCopy();
    int iterationCount = 0;
    pair<int, int> pos = {1,1};
    while(iterationCount <= n+m-1){
        covered[pos.first][pos.second]  =  true;
        if( supplyCopy[pos.first] <= demandCopy[pos.second] ){
            costArrayCopy[pos.first][pos.second] = supplyCopy[pos.first];
            for(int j=1; j<=m; j++){
                if(j==pos.second){
                    continue;
                }else{
                    if(covered[pos.first][j] == true){
                        continue;
                    }else{
                        // covered[pos.first][j] = true;
                        costArrayCopy[pos.first][j] = 0;
                    }
                }
            }
            demandCopy[pos.second] -= supplyCopy[pos.first];
            supplyCopy[pos.first] = 0;
            if(pos.first+1 > n){
                break;
            }
            pos = {pos.first+1, pos.second};
        }else{
            costArrayCopy[pos.first][pos.second] = demandCopy[pos.second];
            for(int i=1; i<=n; i++){
                if(i==pos.first){
                    continue;
                }else{
                    if(covered[i][pos.second] == true){
                        continue;
                    }else{
                        // covered[i][pos.second] = true;
                        costArrayCopy[i][pos.second] = 0;
                    }
                }
            }
            supplyCopy[pos.first] -= demandCopy[pos.second];
            demandCopy[pos.second] = 0;
            if(pos.second+1 > m){
                break;
            }
            pos = {pos.first, pos.second+1};            
        }

        display(costArrayCopy, supplyCopy, demandCopy);
        iterationCount++;
    }
    display(costArrayCopy, supplyCopy, demandCopy);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// function to find loop

// function to find a loop in the transportation table
bool find_loop(bool row_search, vector<pair<int,int>> &loop, double costArrayTemp[][100])
{
    if (row_search)
    {
        auto pt = loop.back();
        for (int j = 1; j <= m; j++)
        {
            if (costArrayTemp[pt.first][j] == 0 && j != pt.second)
            {
                loop.push_back({pt.first, j});
                if (j == loop.front().second)
                    return true;
                if (find_loop(false, loop, costArrayTemp))
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
        for (int i = 1; i <= n; i++)
        {
            if (costArrayTemp[i][pt.second]==0 && i != pt.first)
            {
                loop.push_back({i, pt.second});
                if (i == loop.front().first)
                    return true;
                if (find_loop(true, loop, costArrayTemp))
                    return true;
                else
                    loop.pop_back();
            }
        }
        return false;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// After getting the Initial Basic Feasible Solution - Modi Method to Optimize it
void ModiMethod() {
    // We keep on Iterating till 
    while(1){
        // newSupply and newDemand will contain the values of the 
        // supply and demand values calculated in each iteration
        for(int k=0; k<100; k++){
            newSupply[k] = INF;
            newDemand[k] = INF;
        }

        // --------------------------------------------------------------
        // Occupied and un-occupied cells
        for(int i=1; i<=n; i++){
            for(int j=1; j<=m; j++){
                if(costArrayCopy[i][j] > 0){
                    occupied[i][j] = true;
                }else{
                    occupied[i][j] = false;
                }
            }
        }

        // --------------------------------------------------------------
        // since Number of Variables > Number of Equation
        // we will assign some values 0 - corresponding to the variable which occures most
        // finding the row or col with max occupied cells
        int mxRow = 0, mxCol = 0, rowNo, colNo;
        for(int i=1; i<=n; i++){
            int cnt = 0;
            for(int j=1; j<=m; j++){
                if(occupied[i][j]) cnt++;
            }
            if(cnt > mxRow){
                mxRow = cnt;
                rowNo = i;
            }
        }
        for(int j=1; j<=m; j++){
            int cnt = 0;
            for(int i=1; i<=n; i++){
                if(occupied[i][j]) cnt++;
            }
            if(cnt > mxCol){
                mxCol = cnt;
                colNo = j;
            }
        }

        // --------------------------------------------------------------
        // This part calculates the remaning new values of supply and demand
        //since each equation has 2 variables, so if we get value of row(newSupply) cell
        // then pushing it in a vector so new we will calculate all col(newDemand) cell values
        // and similarly vice-versa till we find all values
        vector<int> rows, cols;
        if(mxRow >= mxCol) {
            newSupply[rowNo] = 0;
            rows.push_back(rowNo);
        }else{
            newDemand[colNo] = 0;
            cols.push_back(colNo);
        }

        while((int)rows.size() > 0 || (int)cols.size() > 0) {
            vector<int> tempRows, tempCols;
            for(int x : rows) {
                for(int j=1; j<= m; j++){
                    if(occupied[x][j] && newDemand[j] == INF){
                        newDemand[j] = costArray[x][j] - newSupply[x];
                        tempCols.push_back(j);
                    }
                }
            }
            for(int y : cols) {
                for(int i=1; i<=n; i++) {
                    if(occupied[i][y] && newSupply[i] == INF){
                        newSupply[i] = costArray[i][y] - newDemand[y];
                        tempRows.push_back(i);
                    }
                }
            }
            rows = tempRows;
            cols = tempCols;
        }
        cout<<"Values of U : "<<endl;
        for(int i=1; i<=n; i++) {
            cout<<"U["<<i<<"] = "<<newSupply[i]<<" | ";
        }
        cout<<endl;
        cout<<"Values of V : "<<endl;
        for(int j=1; j<=m; j++) {
            cout<<"V["<<j<<"] = "<<newDemand[j]<<" | ";
        }
        cout<<endl;
        // --------------------------------------------------------------
        // costArrayTeamp will containd values of Oppurtunity Cost
        double costArrayTemp[100][100];
        bool ok = true;
        for(int i=1; i<=n; i++){
            for(int j=1; j<=m; j++) {
                if(occupied[i][j] == false){
                    costArrayTemp[i][j] = costArray[i][j] - newSupply[i] - newDemand[j];
                }else{
                    costArrayTemp[i][j] = 0;
                }
                if(costArrayTemp[i][j] < 0) {
                    ok = false;
                }
            }
        }
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"Oppurtunity Cost Table :-"<<endl;
        display(costArrayTemp, newSupply, newDemand);
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;


        if(ok){
            cout<<"Found Optimal solution"<<endl;
            return;
        }else{
            // Among all the Negative in Oppurtunity Cost Table ->
            // finding the Most negative one
            int mn = INF;
            pair<int, int> pos;
            for(int i=1; i<=n; i++) {
                for(int j=1; j<=m; j++) {
                    if(costArrayTemp[i][j] < mn){
                        mn = costArrayTemp[i][j];
                        pos = {i,j};
                    }
                }
            }

            // LOGIC TO FIND CLOSED PATH :-
            // now we have to make a closed path starting and ending at this pos
            // Cuurently only coded for rectangular closed path
            // But general closed path -> Using Bfs on the Adjacency matrix of Oppurtunity Cost with
            // the cells having value 0 as nodes
            // in that grapg we can search for a cycle (starting node be the most negavive cell)
            // assign alternate nodes to add and sub vectors


            // finding if a loop exists
            vector<pair<int,int>> loop;
            loop.push_back(pos);
            find_loop(true, loop, costArrayTemp);
            for(pair<int,int> x : loop) {
                cout<<x.first<<","<<x.second<<endl;
            }
            int sz = (int)loop.size();
            double mnn = INF;
            for(int i=0;i<sz;i++){
                if(i&1) {
                    int x = loop[i].first;
                    int y = loop[i].second;
                    if(costArrayCopy[x][y] < mnn){
                        mnn = costArrayCopy[x][y];
                    }
                }
            }
            cout<<mnn<<endl;
            for(int i=0;i<sz;i++){
                int x = loop[i].first;
                int y = loop[i].second;
                if(i&1) {
                    costArrayCopy[x][y] -= mnn;
                }else{
                    costArrayCopy[x][y] += mnn;
                }
            }


            // vector<int> xx, yy;
            // for(int j=1; j<=m; j++){
            //     if(costArrayCopy[pos.first][j] > 0)
            //         yy.push_back(j);
            // }
            // for(int i=1; i<=n; i++){
            //     if(costArrayCopy[i][pos.second] > 0)
            //         xx.push_back(i);
            // }
            // vector<pair<int,int>> add, sub;
            // add.push_back(pos);
            // for(int x : xx){
            //     for(int y : yy) {
            //         if(costArrayCopy[x][y] > 0){
            //             sub.push_back({x,pos.second});
            //             sub.push_back({pos.first,y});
            //             add.push_back({x,y});
            //         }
            //     }
            // }
            // if(sub.size() < 2){
            //     cout<<"Cannot find a closed path"<<endl;
            //     break;
            // }
            // double val;
            // if(costArrayCopy[sub[0].first][sub[0].second] < costArrayCopy[sub[1].first][sub[1].second]) {
            //     costArrayCopy[sub[1].first][sub[1].second] -= costArrayCopy[sub[0].first][sub[0].second];
            //     costArrayCopy[add[0].first][add[0].second] += costArrayCopy[sub[0].first][sub[0].second];
            //     costArrayCopy[add[1].first][add[1].second] += costArrayCopy[sub[0].first][sub[0].second];
            //     costArrayCopy[sub[0].first][sub[0].second] = 0;
            // }else{
            //     costArrayCopy[sub[0].first][sub[0].second] -= costArrayCopy[sub[1].first][sub[1].second];
            //     costArrayCopy[add[0].first][add[0].second] += costArrayCopy[sub[1].first][sub[1].second];
            //     costArrayCopy[add[1].first][add[1].second] += costArrayCopy[sub[1].first][sub[1].second];
            //     costArrayCopy[sub[1].first][sub[1].second] = 0;
            // }

            cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
            cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
            cout<<"Transportation Table after Interation"<<endl;
            display(costArrayCopy, supplyCopy, demandCopy);
            cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        }
        
    }

}


int main(){
    // ------------------------------------------------------------
    // Taking the input
    cout<<"Enter number of Supply Sources : ";
    cin>>n;
    cout<<"Enter number of Destinaton Demands : ";
    cin>>m;

    // Option to calculate the Method to be used to Find Initial Basic Feasible Solution
    int option;
    cout<<"Enter 1 for North West Corner Methos | 2 for Least Cost Method"<<endl;
    cin>>option;

    cout<<"Enter the Cost Values :-"<<endl;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            cout<<"Source : S["<<i<<"] -> Destination : D["<<j<<"] : ";
            cin>>costArray[i][j];
        }
    }

    cout<<"Enter Values of Supply from :- "<<endl;
    for(int i=1; i<=n; i++) {
        cout<<"Source ["<<i<<"] : ";
        cin>>supply[i];
    }

    cout<<"Enter Value of Demand for :- "<<endl;
    for(int j=1; j<=m; j++) {
        cout<<"Destination ["<<j<<"] : ";
        cin>>demand[j];
    }

    // ------------------------------------------------------------
    // Displaying the Initial Transportation Table
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Initial Transportation Table :-"<<endl;
    display(costArray, supply, demand);
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
    // ------------------------------------------------------------
    // Chekcking if the problem is Balanced
    bool isBalanced = checkIfBalanced(supply, demand);

    // If not balanced -> display message and Exit
    if(!isBalanced) {
        cout<<"The given Transportation problem is not balanced"<<endl;
        return 0;
    }

    
    if(option == 1){
        NorthWestCornerMethod();
    }else{
        LeastCostMethod();
    }

    
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Initial Basic Feasible Solution  :- "<<endl;
    double totalCost = 0;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            if(costArrayCopy[i][j] > 0){
                cout<<"S["<<i<<"] -> D["<<j<<"] = "<<costArrayCopy[i][j]<<endl;
                totalCost = totalCost + costArrayCopy[i][j] * costArray[i][j];
            }
        }
    }
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Initial Transportation Cost :- "<<endl;
    cout<<totalCost<<endl;
    cout<<"----------------------------------------------------------"<<endl;


    // Calling Modi Method to Optimize the Solution
    ModiMethod();


    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Solution After Using Modi Method :- "<<endl;
    totalCost = 0;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            if(costArrayCopy[i][j] > 0){
                cout<<"S["<<i<<"] -> D["<<j<<"] = "<<costArrayCopy[i][j]<<endl;
                totalCost = totalCost + costArrayCopy[i][j] * costArray[i][j];
            }
        }
    }
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Transportation Cost :- "<<endl;
    cout<<totalCost<<endl;
    cout<<"----------------------------------------------------------"<<endl;
}

/*
Input 1 -
3 4 1
2 3 5 1
7 3 4 6
4 1 7 2
8 10 20
6 8 9 15


Output 1 -
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
The Initial Transportation Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 2          | 3          | 5          | 1          | 8
s[2]       | 7          | 3          | 4          | 6          | 10
s[3]       | 4          | 1          | 7          | 2          | 20
__________________________________________________________________
           | 6          | 8          | 9          | 15         |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Enter 1 for North West Corner Methos | 2 for Least Cost Method
1
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 3          | 5          | 1          | 2
s[2]       | 0          | 3          | 4          | 6          | 10        
s[3]       | 0          | 1          | 7          | 2          | 20        
__________________________________________________________________
           | 0          | 8          | 9          | 15         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 3          | 4          | 6          | 10        
s[3]       | 0          | 1          | 7          | 2          | 20        
__________________________________________________________________
           | 0          | 6          | 9          | 15         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 4          | 6          | 4
s[3]       | 0          | 0          | 7          | 2          | 20        
__________________________________________________________________
           | 0          | 0          | 9          | 15         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 4          | 0          | 0
s[3]       | 0          | 0          | 7          | 2          | 20
__________________________________________________________________
           | 0          | 0          | 5          | 15         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 4          | 0          | 0
s[3]       | 0          | 0          | 5          | 2          | 15        
__________________________________________________________________
           | 0          | 0          | 0          | 15         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 4          | 0          | 0
s[3]       | 0          | 0          | 5          | 15         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          |
----------------------------------------------------------
----------------------------------------------------------
----------------------------------------------------------
Initial Basic Feasible Solution  :-
S[1] -> D[1] = 6
S[1] -> D[2] = 2
S[2] -> D[2] = 6
S[2] -> D[3] = 4
S[3] -> D[3] = 5
S[3] -> D[4] = 15
----------------------------------------------------------
Initial Transportation Cost :-
117
----------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       | 
__________________________________________________________________
s[1]       | 0          | 0          | 1          | 2          | 0
s[2]       | 5          | 0          | 0          | 7          | 0
s[3]       | -1         | -5         | 0          | 0          | 3
__________________________________________________________________
           | 2          | 3          | 4          | -1         |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Transportation Table after Interation
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 2          | 0          | 0          | 0
s[2]       | 0          | 1          | 9          | 0          | 0
s[3]       | 0          | 5          | 0          | 15         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 0          | 0          | 1          | -3         | 3
s[2]       | 5          | 0          | 0          | 2          | 3
s[3]       | 4          | 0          | 5          | 0          | 1
__________________________________________________________________
           | -1         | 0          | 1          | 1          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Transportation Table after Interation
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 6          | 0          | 0          | 2          | 0
s[2]       | 0          | 1          | 9          | 0          | 0
s[3]       | 0          | 7          | 0          | 13         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 0          | 3          | 4          | 0          | 0
s[2]       | 2          | 0          | 0          | 2          | 3
s[3]       | 1          | 0          | 5          | 0          | 1
__________________________________________________________________
           | 2          | 0          | 1          | 1          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Found Optimal solution
----------------------------------------------------------
----------------------------------------------------------
----------------------------------------------------------
Solution After Using Modi Method :-
S[1] -> D[1] = 6
S[1] -> D[4] = 2
S[2] -> D[2] = 1
S[2] -> D[3] = 9
S[3] -> D[2] = 7
S[3] -> D[4] = 13
----------------------------------------------------------
Transportation Cost :-
86
----------------------------------------------------------
*/

/*
Input 2 -
3 4 1
19 30 50 10
70 30 40 60
40 8 70 20
7 9 18
5 8 7 14

Output 2 -
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
The Initial Transportation Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 19         | 30         | 50         | 10         | 7
s[2]       | 70         | 30         | 40         | 60         | 9
s[3]       | 40         | 8          | 70         | 20         | 18
__________________________________________________________________
           | 5          | 8          | 7          | 14         |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Enter 1 for North West Corner Methos | 2 for Least Cost Method
1
           | D[1]       | D[2]       | D[3]       | D[4]       |  
__________________________________________________________________
s[1]       | 5          | 30         | 50         | 10         | 2
s[2]       | 0          | 30         | 40         | 60         | 9
s[3]       | 0          | 8          | 70         | 20         | 18        
__________________________________________________________________
           | 0          | 8          | 7          | 14         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 30         | 40         | 60         | 9
s[3]       | 0          | 8          | 70         | 20         | 18        
__________________________________________________________________
           | 0          | 6          | 7          | 14         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 40         | 60         | 3
s[3]       | 0          | 0          | 70         | 20         | 18        
__________________________________________________________________
           | 0          | 0          | 7          | 14         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 3          | 0          | 0
s[3]       | 0          | 0          | 70         | 20         | 18
__________________________________________________________________
           | 0          | 0          | 4          | 14         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 3          | 0          | 0
s[3]       | 0          | 0          | 4          | 20         | 14
__________________________________________________________________
           | 0          | 0          | 0          | 14         |
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 6          | 3          | 0          | 0
s[3]       | 0          | 0          | 4          | 14         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          | 
----------------------------------------------------------
----------------------------------------------------------
----------------------------------------------------------
Initial Basic Feasible Solution  :-
S[1] -> D[1] = 5
S[1] -> D[2] = 2
S[2] -> D[2] = 6
S[2] -> D[3] = 3
S[3] -> D[3] = 4
S[3] -> D[4] = 14
----------------------------------------------------------
Initial Transportation Cost :-
1015
----------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       | 
__________________________________________________________________
s[1]       | 0          | 0          | 10         | 20         | 0
s[2]       | 51         | 0          | 0          | 70         | 0
s[3]       | -9         | -52        | 0          | 0          | 30
__________________________________________________________________
           | 19         | 30         | 40         | -10        |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Transportation Table after Interation
           | D[1]       | D[2]       | D[3]       | D[4]       | 
__________________________________________________________________
s[1]       | 5          | 2          | 0          | 0          | 0
s[2]       | 0          | 2          | 7          | 0          | 0
s[3]       | 0          | 4          | 0          | 14         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 0          | 0          | 10         | -32        | 30        
s[2]       | 51         | 0          | 0          | 18         | 30
s[3]       | 43         | 0          | 52         | 0          | 8
__________________________________________________________________
           | -11        | 0          | 10         | 12         |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Transportation Table after Interation
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 5          | 0          | 0          | 2          | 0
s[2]       | 0          | 2          | 7          | 0          | 0
s[3]       | 0          | 6          | 0          | 12         | 0
__________________________________________________________________
           | 0          | 0          | 0          | 0          |
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Oppurtunity Cost Table :-
           | D[1]       | D[2]       | D[3]       | D[4]       |
__________________________________________________________________
s[1]       | 0          | 32         | 42         | 0          | 0
s[2]       | 19         | 0          | 0          | 18         | 32
s[3]       | 11         | 0          | 52         | 0          | 10
__________________________________________________________________
           | 19         | -2         | 8          | 10         | 
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Found Optimal solution
----------------------------------------------------------
----------------------------------------------------------
----------------------------------------------------------
Solution After Using Modi Method :-
S[1] -> D[1] = 5
S[1] -> D[4] = 2
S[2] -> D[2] = 2
S[2] -> D[3] = 7
S[3] -> D[2] = 6
S[3] -> D[4] = 12
----------------------------------------------------------
Transportation Cost :-
743
----------------------------------------------------------
*/