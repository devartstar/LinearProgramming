/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
    LAB 7 - INTEGER PROGRAMING PROBLEM
    INITIAL SUBMISSION
*/

/*
    LOGIC :-
    ----------------
    Data Structure used
    1. To represent a tree i will be using an array of structure
    Each structure in the array will represent a node and will contain all information :-
        The Inequalities, Equation, number of unknows, everything
        Also all the simplex tables and their solutions
    2. Inside the array(every node) we have a structure 
        ie. (node is a structure and insde a structure we have a nestedarray of structure)
        each item of this array of structure table will store all the deatils of the i'th iteration of that node/branch

    STEP 1 :-
        Take in all the input given by the user inside tree[1] <--- first node
    STEP 2 :-
        make a recursive call
        parent (index : x) ----> if not integer solution 
                        ------> left child (index : 2x) (add the extra <= inequality)  
                        ------> right child (index : 2x + 1) (add the extra >= inequality)  
        
        inside the recursive call 
            function transformEqn() -> make the inequalities in desired format (ie. all rhs > 0)
            function conInqualityToEquality() -> will add the slack and arificial variable and convert to system of equalities

        Solve this LPP question with big M method

        if infeasible solution return false;
        check if all answer are integers :-
            YES - check if its the best optimal solution and return true
            NO - branch the problem in left and righ node
                and recursively proceed

        we return TRUE from a node if atleast any branch lower than it (in both left and right subtree) 
        has an integer optimal solution
        else we return FALSE
        
        finally if we return TRUE from node 1 -> we got an integer optimal solution
        else Not feasible solution
            


*/

#include <bits/stdc++.h>
using namespace std;

long long M = 1000000;  // a large number ( to be multiples with artificial variable in Objective function)

// A structure to define the table used for every iteration
struct table{
    double coefMarix[100][100];         // n, m
    double rhs[100];                    // n
    double argumentedMatrix[100][100];  // n, m+1
    pair<int, double> CB[100];          // <index, value>
    double ratio[100];                  // n
    double Zj[100];                     // m
    double CjSubZj[100];                // m
    int kRow;
    int kCol;
    int kElement;

};

// Lets define a tree
// [i][j] ---> ith node in a tree has 2*i 2*i+1 as its children
// so in branching we will call the children from the parent node
struct treeNode{
    int n;  // number of Equations
    int m;  // number of Variables
    int a, last_var_ind;  // a will have the count of number of Artificial Variables
    int equality_count = 0;
    int mm;
    double Equation[100];
    double Inequalities[100][100];
    double RHS[100];
    int chc[100];
    bool toMaximize = false;

    int copy_n;
    int copy_m;
    double copy_Equation[100];
    double copy_Inequalities[100][100];
    double copy_RHS[100];
    int copy_chc[100];
    bool copy_toMaximize = false;

    table arr[100];
    int iterationCount = 0;
};

treeNode tree[100];
int node = 1;

double optimalValue;
int optimalSolutionNode;

///////////////////////////////////////////////////////////////////////////////////
// Display the Inequalities
void displayInEqualities(int node){
    if(tree[node].toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= tree[node].m; j++){
            cout<<" + ("<<tree[node].Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= tree[node].m; j++){
            cout<<" + ("<<-1*tree[node].Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }

    cout<<"Inequality Constrains :- "<<endl;
    for(int i = 1; i <= tree[node].n; i++){
        cout<<"Inq "<<i<<" :- ";
        for(int j = 1; j <= tree[node].m; j++){
            cout<<" + ("<<tree[node].Inequalities[i][j]<<"x["<<j<<"])";
        }
        cout<<(tree[node].chc[i]==1?" <= ":" >= ")<<tree[node].RHS[i]<<endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////
// Displaying The Transformed Equations after inserting slack variables
void displayEqualities(int node){
    if(tree[node].toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= tree[node].mm; j++){
            cout<<" + ("<<tree[node].Equation[j]<<"x["<<j<<"])";
        }
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= tree[node].mm; j++){
            cout<<" + ("<<-1*tree[node].Equation[j]<<"x["<<j<<"])";
        }
    }
    for(int j = tree[node].mm+1; j <= tree[node].m; j++){
        cout<<" + ("<<tree[node].Equation[j]<<"x["<<j<<"])";
    }
    cout<<endl;

    cout<<"Simultaneous equations :- "<<endl;
    for(int i = 1; i <= tree[node].n; i++){
        cout<<"Eq "<<i<<" :- ";
        for(int j = 1; j <= tree[node].m; j++){
            cout<<" + ("<<tree[node].Inequalities[i][j]<<"x["<<j<<"])";
        }
        cout<<" = "<<tree[node].RHS[i]<<endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////
// Function to Display the Table
void displayTable(table t, int node){
    cout<<left<<setw(18)<<"Basic-Variable"<<left<<setw(10)<<"Coef";
    for(int j = 1; j <= tree[node].m; j++){
        cout<<left<<setw(10)<<("x["+to_string(j)+"]");

    }
    cout<<left<<setw(10)<<"Value"<<left<<setw(10)<<"Ratio"<<endl;

    for(int i = 1; i <= tree[node].n; i++){
        cout<<left<<setw(18)<<("x["+to_string(t.CB[i].first)+"]")<<left<<setw(10)<<t.CB[i].second;
        for(int j = 1; j <= tree[node].m; j++){
            cout<<left<<setw(10)<<t.argumentedMatrix[i][j];
        }
        cout<<left<<setw(10)<<t.argumentedMatrix[i][tree[node].m+1]<<left<<setw(10)<<(t.ratio[i]!=100005?to_string(t.ratio[i]):"_")<<endl;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"Z :- ";
    for(int j = 1; j <= tree[node].m+1; j++){
        cout<<left<<setw(10)<<t.Zj[j];
    }
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"C - Z :- ";
    for(int j = 1; j <= tree[node].m+1; j++){
        cout<<left<<setw(10)<<t.CjSubZj[j];
    }
    cout<<endl;
}

///////////////////////////////////////////////////////////////////////////////////
// Function to construct the Initial Table
table constructInitialTable(int node){
    table t;

    for(int i = 1; i<= tree[node].n; i++){
        for(int j = 1; j <= tree[node].m; j++){
            t.coefMarix[i][j] = tree[node].Inequalities[i][j];
        }
    }

    // Constructing the SOl column
    for(int i = 1; i<= tree[node].n; i++){
        t.rhs[i] = tree[node].RHS[i];
    }

    // Constructing the Argumented Matrix
    for(int i = 1; i<= tree[node].n; i++){
        for(int j = 1; j <= tree[node].m; j++){
            t.argumentedMatrix[i][j] = t.coefMarix[i][j];
        }
        t.argumentedMatrix[i][tree[node].m+1] = t.rhs[i];
    }

    // initial base variables last n variables (may include artificial + slack variables)
    for(int i = tree[node].n; i>=1; i--){
        t.CB[tree[node].n-i+1].first = tree[node].m + 1 - i;
        t.CB[tree[node].n-i+1].second = tree[node].Equation[tree[node].m + 1 - i];
    }

    // calculating Zj and cjSubZj for all m variables
    for(int j = 1; j<= tree[node].m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= tree[node].n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.CjSubZj[j] = (tree[node].Equation[j] - t.Zj[j]);
    }

    // kCOl is the max of all CjSubZj
    double mx = -100005;
    for(int j = 1; j <= tree[node].m; j++){
        if(t.CjSubZj[j] > mx){
            mx = t.CjSubZj[j];
            t.kCol = j;
        }
    }

    // Calculating the Ratios
    for(int i = 1; i <= tree[node].n; i++){
        if(t.argumentedMatrix[i][t.kCol] > 0){
            t.ratio[i] = t.argumentedMatrix[i][tree[node].m+1] / t.argumentedMatrix[i][t.kCol];
        }else{
            t.ratio[i] = 100005;
        }
    }

    // kRow is the min of all the ratios
    double mn = 100005;
    for(int i = 1; i <= tree[node].n; i++){
        if(t.ratio[i] < mn){
            mn = t.ratio[i];
            t.kRow = i;
        }
    }

    return t;
}


///////////////////////////////////////////////////////////////////////////////////
// Function to construct the simplex table for next iteration from current table
table constructNextSimplexTable(table t_prev, int node){
    table t_current;
    int leavingVariable  = t_prev.CB[t_prev.kRow].first;
    int enteringVariable = t_prev.kCol;
    double pivotElement = t_prev.argumentedMatrix[t_prev.kRow][t_prev.kCol];

    // Constructing the Base Variable and Corresponding ARGUMENTED MATRIX
    int k = 1;
    for(int i = 1; i <= tree[node].n; i++){
        if(t_prev.CB[i].first == leavingVariable)   continue;
        t_current.CB[k] = t_prev.CB[i];
        // row in the argumented matrix according to the base variable
        for(int j = 1; j <= tree[node].m+1; j++){
            t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[i][j] - (t_prev.argumentedMatrix[t_prev.kRow][j]*t_prev.argumentedMatrix[i][t_prev.kCol])/pivotElement;
        }
        k++;
    }
    t_current.CB[k] = {enteringVariable, tree[node].Equation[enteringVariable]};
    for(int j = 1; j <= tree[node].m+1; j++){
        t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[t_prev.kRow][j] / pivotElement;
    }

    // Calculating the Zj and corresponding ZjSubCj
    for(int j = 1; j<= tree[node].m+1; j++){
        t_current.Zj[j] = 0;
        for(int i = 1; i <= tree[node].n; i++){
            t_current.Zj[j] = t_current.Zj[j] + t_current.CB[i].second * t_current.argumentedMatrix[i][j];
        }
        t_current.CjSubZj[j] = (tree[node].Equation[j] - t_current.Zj[j]);
    }

    // kCOl is the max of all CjSubZj
    double mx = -100005;
    for(int j = 1; j <= tree[node].m; j++){
        if(t_current.CjSubZj[j] > mx){
            mx = t_current.CjSubZj[j];
            t_current.kCol = j;
        }
    }

    // Calculating the Ratios
    for(int i = 1; i <= tree[node].n; i++){
        if(t_current.argumentedMatrix[i][t_current.kCol] > 0){
            t_current.ratio[i] = t_current.argumentedMatrix[i][tree[node].m+1] / t_current.argumentedMatrix[i][t_current.kCol];
        }else{
            t_current.ratio[i] = 100005;    // a very large value (bcs to find kRow we have to select min of ratios)
        }
    }

    // kRow is the min of all the ratios
    double mn = 100005;
    for(int i = 1; i <= tree[node].n; i++){
        if(t_current.ratio[i] < mn){
            mn = t_current.ratio[i];
            t_current.kRow = i;
        }
    }

    return t_current;
}



///////////////////////////////////////////////////////////////////////////////////
// We return true if all Elements in Ci - Zi (CjSubZj) are negative (then we have found optimal solution)
bool stopingCondition(table t, int node){
    bool ok = true;
    for(int j = 1; j <= tree[node].m; j++){
        if(t.CjSubZj[j] > 0){
            ok = false;
            break;
        }
    }
    for(int j = 1; j<=tree[node].n; j++){
        if(t.CB[j].first > tree[node].last_var_ind){
            ok = false;
            break;
        }
    }
    return ok;
}

/////////////////////////////////////////////////////////////////////////////////////
// constructing the set of Inequations of the child node from parent node
void constructEqn(int node){
    // number of inequalities inc by 1
    tree[node].n = tree[node/2].copy_n+1;   
    tree[node].m = tree[node/2].copy_m;

    for(int i = 1; i <= tree[node/2].copy_n; i++){
        for(int j = 1; j <= tree[node/2].copy_m; j++){
            tree[node].Inequalities[i][j] = tree[node/2].copy_Inequalities[i][j];
        }
        tree[node].RHS[i] = tree[node/2].copy_RHS[i];
        tree[node].chc[i] = tree[node/2].copy_chc[i];
    }
    // Add that one extra inequality
    // .... ??? node%2 == 0--->
    // .... ??? node%2 == 1--->
    // yet to be coded

    tree[node].toMaximize = tree[node/2].copy_toMaximize;

    for(int j = 1; j <= tree[node/2].copy_m; j++){
        tree[node].Equation[j] = tree[node/2].copy_Equation[j];
    }
}


////////////////////////////////////////////////////////////////////////////////////////
void transformEqn(int node){
    // making sure all inequality is of the form less than equals to
    for(int i=1; i<=tree[node].n; i++){
        if(tree[node].chc[i] == 3){
            tree[node].equality_count ++;
        }
        // since RHS cannot be negative (multiplying LHS and RHS by -1 and changing equality sign)
        // and we alter chc
        if(tree[node].RHS[i] < 0){
            for(int j = 1; j <= tree[node].m; j++){
                tree[node].Inequalities[i][j] = -1*tree[node].Inequalities[i][j];
            }
            tree[node].RHS[i] = -1*tree[node].RHS[i];
            if(tree[node].chc[i] == 1)    tree[node].chc[i] = 2;
            else if(tree[node].chc[i] == 2) tree[node].chc[i] = 1;
        }
    }
    
    // Since we want the eqn to be to maximized
    // so if we have to minimize - we just multiply the equation by -1 to convert to maximization problem
    if(tree[node].toMaximize == false){
        for(int j = 1; j <= tree[node].m; j++){
            tree[node].Equation[j] = -1 * tree[node].Equation[j];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
void conInequalityToEquality(int node){
    tree[node].mm = tree[node].m;
    tree[node].m = tree[node].m + tree[node].n - tree[node].equality_count;
    tree[node].last_var_ind = tree[node].m;

    // NOW WE HAVE TO TRANSFORM THE INEQUALITES TO EQUATIONS
    tree[node].a = 1;
    int k = 1;
    for(int i = 1; i <= tree[node].n; i++){
        if(tree[node].chc[i] == 3){
            // Equality  (add only artificial variables)
            tree[node].Inequalities[i][tree[node].m+tree[node].a] = 1;   // aftificial variables
            tree[node].a++;
        }else if(tree[node].chc[i] == 2){
            // >=  (add slack and artificial variables)
            tree[node].Inequalities[i][tree[node].mm+k] = -1; // slack variable
            tree[node].Inequalities[i][tree[node].m+tree[node].a] = 1;   // aftificial variables
            tree[node].a++;
            k++;
        }else if(tree[node].chc[i] == 1){
            // <=  (add slack variables)
            tree[node].Inequalities[i][tree[node].mm+k] = 1;
            k++;
        }
    }
    // NO also insert artificial variables in Objective Function
    for(int i = tree[node].m+1; i <= tree[node].m+tree[node].a-1; i++){
        tree[node].Equation[i] = -M;
    }
    tree[node].m = tree[node].m + tree[node].a - 1;
}


bool checkIfIntegerSolution(int node){
    // ... yet to code.
    // just iterate over all the unknows in 
    // for i = 1 to tree[node].n
    //      if( tree[node].arr[tree[node].iterationcount].CB[i] is integer or not )
    //             if not integer return FALSE else continue
    // if we havent retun yet means all integer so return TRUE
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// will return true if it found integer solution 
// on solving the simplex for that node
bool recursiveSolve(int node){
    if(node > 100){
        // we keep on getting fractions as answer
        // then we stop
        return false;
    }
    if(node > 1) constructEqn(node);

    // converts to proper desired Inequality format
    transformEqn(node);
    // now we can solely focus on maximization problem
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    displayInEqualities(node);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;

    // converts Inequality to Equality
    conInequalityToEquality(node);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Transformed Equations are :- "<<endl;
    displayEqualities(node);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // constructing the initial simplex table
    tree[node].arr[0] = constructInitialTable(node);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE :-"<<endl;
    displayTable(tree[node].arr[0], node);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // After Constructing Initial Simplex Table 
    // Iterate to construst simplex table of i+1'th iteration using i'th iteration
    // till the stopping condition is achievied
    tree[node].iterationCount = 0; 
    while(stopingCondition(tree[node].arr[tree[node].iterationCount], node) == false){
        tree[node].iterationCount++;
        tree[node].arr[tree[node].iterationCount] = constructNextSimplexTable(tree[node].arr[tree[node].iterationCount-1], node);
        cout<<"--------------------------------------------------------------------"<<endl;
        cout<<"THE SIMPLEX TABLE AFTER "<<tree[node].iterationCount<<"'th ITERATION :-"<<endl;
        displayTable(tree[node].arr[tree[node].iterationCount], node);
        cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(tree[node].iterationCount >= 10){    // 10 for debugging make it to 100 or more later
            return false;
        }
    }

    // check for the answers ---->
    bool ans = false, ans1 = false, ans2 = false;
    if(checkIfIntegerSolution(node)){   // all ans of basic variable integer
        double finval = tree[node].arr[tree[node].iterationCount].Zj[tree[node].m+1];
        if(tree[node].toMaximize){
            if(finval > optimalValue)
                optimalValue = finval;
        }else{
            if(finval < optimalValue)
                optimalValue = finval;
        }
        // no need to branch further
        // we can go back comparing
        return true;
    }else{
        ans1 = recursiveSolve(2*node);
        if(ans1 == true){
            double finval = tree[2*node].arr[tree[2*node].iterationCount].Zj[tree[2*node].m+1];
            if(tree[2*node].toMaximize){
                if(finval > optimalValue){
                    optimalValue = finval;
                    optimalSolutionNode = 2*node;
                }
            }else{
                if(finval < optimalValue){
                    optimalValue = finval;
                    optimalSolutionNode = 2*node;
                }
            }
        }
        ans2 = recursiveSolve(2*node+1);
        if(ans2 == true){
            double finval = tree[2*node+1].arr[tree[2*node+1].iterationCount].Zj[tree[2*node+1].m+1];
            if(tree[2*node+1].toMaximize){
                if(finval > optimalValue){
                    optimalValue = finval;
                    optimalSolutionNode = 2*node+1;
                }
            }else{
                if(finval < optimalValue){
                    optimalValue = finval;
                    optimalSolutionNode = 2*node+1;
                }
            }
        }
        // if any 1 of the children nodes has got integer solution, then we can say from below that node there exists a solution
        // so return true
        // if both dosent get infeasible solution then return false 
        return ans1 | ans2;
    }
}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

int main(){
    ////////////////////////////////////////////////////////////////////
    // I am assuming non negativity to hold & else if x < 0 then we have to replace that
    // variable x by another variable p such that p = (-x), so p > 0     
    cout<<"Enter Number of Inequalities (Apart from Non Negativity): - ";
    cin>>tree[1].n;
    cout<<"Enter Number of Variables : - ";
    cin>>tree[1].m;
    tree[1].copy_n = tree[1].n;
    tree[1].copy_m = tree[1].m;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;

    for(int i = 1; i <= tree[1].n; i++){
        cout<<"Enter details of "<<i<<"'th inequality"<<endl;
        for(int j = 1; j <= tree[1].m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
            cin>>tree[1].Inequalities[i][j];
            tree[1].copy_Inequalities[i][j] = tree[1].Inequalities[i][j];
        }
        cout<<"Enter Value of Rhs of inequation : ";
        cin>>tree[1].RHS[i];
        tree[1].copy_RHS[i] = tree[1].RHS[i];

        cout<<"RHS is less equal(Enter 1) || greater equal(Enter 2) || equal(Enter 3): ";
        cin>>tree[1].chc[i];
        tree[1].copy_chc[i] = tree[1].chc[i];

    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter 1 for Maximization || 2 for Minimization : ";
    int choice; cin>>choice;
    if(choice==1)   tree[1].toMaximize = true;
    tree[1].copy_toMaximize = tree[1].toMaximize;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter the Equation to Maximize / Minimize :- "<<endl;
    for(int j = 1; j <= tree[1].m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
        cin>>tree[1].Equation[j];
        tree[1].copy_Equation[j] = tree[1].Equation[j];
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(tree[1].toMaximize){
        optimalValue = -1e9;
    }else{
        optimalValue = 1e9;
    }
    optimalSolutionNode = -1;

    bool ok = recursiveSolve(1);

    if(ok == true){
        cout<<"Found Integer Solution"<<endl;
        cout<<"The optimal Value is = "<<optimalValue<<endl;
        cout<<"The value of the unkows are = "<<endl;
        // for that optimal node tree[optimalNode] show the values of basic variables

    }else{
        cout<<"Not Fount Integer Solution"<<endl;
    }

    
}