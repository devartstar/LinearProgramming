/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
*/


/*
    LOGIC :- 
    Solve the Linear Programing Problem using Simplex Method

    Consider an example :-
    Max : Z : 12x(1) + 16x(2)
    Cond: 10x(1) + 20x(2) <= 120
          8x(1) + 8x(2) <= 80
          x1 >= 0
          x2 >= 0

    Step 1 :-
        Convert the question(inequality) into standard form(equality) :-
        Lets add slack variables - x(3) & x(4)
        Max : Z : 12x(1) + 16x(2) + x(3) + x(4)
        Cond:   10x(1) + 20x(2) + x(3) = 120
                8x(1) + 8x(2) + x(4) = 80

    Step 2 :-
        Input Simplex Table :-
        _____________________________________________________________________________________________
        CB | (BasicVariable   )    | (x(1),12) | (x(2),16) | (x(3),0) | (x(4),0) | Soln | Ratio |
        _____________________________________________________________________________________________
        0  |      x(3)             |     10    |      20   |    1     |     0    | 120  | 120/20 = 6
        0  |      x(4)             |      8    |      8    |    0     |     1    | 80   | 80/8 = 8
        _____________________________________________________________________________________________
        Zj = Sum(CB(i)*a(i)(j))    |     0     |      0    |    0     |     0    |  0   |
        _____________________________________________________________________________________________
        Cj - Zj                    |     12    |     16    |    0     |     0    |
        _____________________________________________________________________________________________

        The Max vaklue of (Cj - Zj) ROW gives the k col -> The basic Variable in k col (Incomming Element)
        The Min value of (Ratio) COL gives the k row -> The basic Variable in k row (Outgoing Element)
        The Intersection element gives the k element 

    Step 3 :-
        Construst the new Simplex table :- with the k element being 1 and all the elements in k col be 0
        New Value = OLd Value -  Corresponding (k col value * k row value) / k element
        _____________________________________________________________________________________________
        CB | (BasicVariable)       | (x(1),12) | (x(2),16) | (x(3),0) | (x(4),0) | Soln | Ratio |
        _____________________________________________________________________________________________
        16 |      x(2)             |    1/2    |      1    |  1/20    |     0    |   6  |   6/1/2 = 12
        0  |      x(4)             |      4    |      0    |   -2/5   |     1    |  32  | 32/4 = 8
        _____________________________________________________________________________________________
        Zj = Sum(CB(i)*a(i)(j))    |     8     |     16    |    4/5   |     0    |  
        _____________________________________________________________________________________________
        Cj - Zj                    |     4     |     0     |  -4/5    |     0    |
        _____________________________________________________________________________________________

        For Maimization problem :- Stopping condition is all (Cj - Zj) <= 0
        Outgoing Variable :- x(4)
        INcomming Variable x(1)

        _____________________________________________________________________________________________
        CB | (BasicVariable)       | (x(1),12) | (x(2),16) | (x(3),0) | (x(4),0) | Soln | Ratio |
        _____________________________________________________________________________________________
        16 |      x(2)             |      0    |      1    |    0     |   -1/8   |   2  |   
        12 |      x(1)             |      1    |      0    |   -1/10  |     8    |   8  | 
        _____________________________________________________________________________________________
        Zj = Sum(CB(i)*a(i)(j))    |     12    |     16    |   2/5    |     1    |   128
        _____________________________________________________________________________________________
        Cj - Zj                    |     0     |     0     |  -2/5    |     -1    |
        _____________________________________________________________________________________________

    Since all Cj - Zj is <= 0 so we have optimal vlaue
    x(1) = 12
    x(2) = 16
    Max Z : 128

*/


#include <bits/stdc++.h>
using namespace std;


long long M = 1000000;  // a large number ( to be multiples with artificial variable in Objective function)
int a, last_var_ind;  // a will have the count of number of Artificial Variables
int equality_count = 0;
int n, m, mm;
double Equation[100];
double Inequalities[100][100];
double RHS[100];
bool toMaximize = false;

int chc[100];

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

// an array of these tables ( we can also create array size dynamically insead of fixing max Iterations)
table arr[3];

///////////////////////////////////////////////////////////////////////////////////
// Display the Inequalities
void displayInEqualities(){
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<-1*Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }

    cout<<"Inequality Constrains :- "<<endl;
    for(int i = 1; i <= n; i++){
        cout<<"Inq "<<i<<" :- ";
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<Inequalities[i][j]<<"x["<<j<<"])";
        }
        cout<<(chc[i]==1?" <= ":" >= ")<<RHS[i]<<endl;
    }
}


///////////////////////////////////////////////////////////////////////////////////
// Displaying The Transformed Equations after inserting slack variables
void displayEqualities(){
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
        }
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<-1*Equation[j]<<"x["<<j<<"])";
        }
    }
    for(int j = mm+1; j <= m; j++){
        cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
    }
    cout<<endl;

    cout<<"Simultaneous equations :- "<<endl;
    for(int i = 1; i <= n; i++){
        cout<<"Eq "<<i<<" :- ";
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<Inequalities[i][j]<<"x["<<j<<"])";
        }
        cout<<" = "<<RHS[i]<<endl;
    }
}


///////////////////////////////////////////////////////////////////////////////////
// Function to Display the Table
void displayTable(table t){
    cout<<left<<setw(18)<<"Basic-Variable"<<left<<setw(10)<<"Coef";
    for(int j = 1; j <= m; j++){
        cout<<left<<setw(10)<<("x["+to_string(j)+"]");

    }
    cout<<left<<setw(10)<<"Value"<<left<<setw(10)<<"Ratio"<<endl;

    for(int i = 1; i <= n; i++){
        cout<<left<<setw(18)<<("x["+to_string(t.CB[i].first)+"]")<<left<<setw(10)<<t.CB[i].second;
        for(int j = 1; j <= m; j++){
            cout<<left<<setw(10)<<t.argumentedMatrix[i][j];
        }
        cout<<left<<setw(10)<<t.argumentedMatrix[i][m+1]<<left<<setw(10)<<(t.ratio[i]!=100005?to_string(t.ratio[i]):"_")<<endl;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"Z :- ";
    for(int j = 1; j <= m+1; j++){
        cout<<left<<setw(10)<<t.Zj[j];
    }
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"C - Z :- ";
    for(int j = 1; j <= m+1; j++){
        cout<<left<<setw(10)<<t.CjSubZj[j];
    }
    cout<<endl;
}


///////////////////////////////////////////////////////////////////////////////////
// Function to construct the Initial Table
table constructInitialTable(){
    table t;

    for(int i = 1; i<= n; i++){
        for(int j = 1; j <= m; j++){
            t.coefMarix[i][j] = Inequalities[i][j];
        }
    }

    // Constructing the SOl column
    for(int i = 1; i<= n; i++){
        t.rhs[i] = RHS[i];
    }

    // Constructing the Argumented Matrix
    for(int i = 1; i<= n; i++){
        for(int j = 1; j <= m; j++){
            t.argumentedMatrix[i][j] = t.coefMarix[i][j];
        }
        t.argumentedMatrix[i][m+1] = t.rhs[i];
    }

    // initial base variables last n variables (may include artificial + slack variables)
    for(int i = 1; i<=n; i++){
        t.CB[i].first = m + 1 - i;
        t.CB[i].second = Equation[m + 1 - i];
    }

    // calculating Zj and cjSubZj for all m variables
    for(int j = 1; j<= m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.CjSubZj[j] = (Equation[j] - t.Zj[j]);
    }

    // kCOl is the max of all CjSubZj
    double mx = -100005;
    for(int j = 1; j <= m; j++){
        if(t.CjSubZj[j] > mx){
            mx = t.CjSubZj[j];
            t.kCol = j;
        }
    }

    // Calculating the Ratios
    for(int i = 1; i <= n; i++){
        if(t.argumentedMatrix[i][t.kCol] > 0){
            t.ratio[i] = t.argumentedMatrix[i][m+1] / t.argumentedMatrix[i][t.kCol];
        }else{
            t.ratio[i] = 100005;
        }
    }

    // kRow is the min of all the ratios
    double mn = 100005;
    for(int i = 1; i <= n; i++){
        if(t.ratio[i] < mn){
            mn = t.ratio[i];
            t.kRow = i;
        }
    }

    return t;
}

///////////////////////////////////////////////////////////////////////////////////
// Function to construct the simplex table for next iteration from current table
table constructNextSimplexTable(table t_prev){
    table t_current;
    int leavingVariable  = t_prev.CB[t_prev.kRow].first;
    int enteringVariable = t_prev.kCol;
    double pivotElement = t_prev.argumentedMatrix[t_prev.kRow][t_prev.kCol];

    // Constructing the Base Variable and Corresponding ARGUMENTED MATRIX
    int k = 1;
    for(int i = 1; i <= n; i++){
        if(t_prev.CB[i].first == leavingVariable)   continue;
        t_current.CB[k] = t_prev.CB[i];
        // row in the argumented matrix according to the base variable
        for(int j = 1; j <= m+1; j++){
            t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[i][j] - (t_prev.argumentedMatrix[t_prev.kRow][j]*t_prev.argumentedMatrix[i][t_prev.kCol])/pivotElement;
        }
        k++;
    }
    t_current.CB[k] = {enteringVariable, Equation[enteringVariable]};
    for(int j = 1; j <= m+1; j++){
        t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[t_prev.kRow][j] / pivotElement;
    }

    // Calculating the Zj and corresponding ZjSubCj
    for(int j = 1; j<= m+1; j++){
        t_current.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t_current.Zj[j] = t_current.Zj[j] + t_current.CB[i].second * t_current.argumentedMatrix[i][j];
        }
        t_current.CjSubZj[j] = (Equation[j] - t_current.Zj[j]);
    }

    // kCOl is the max of all CjSubZj
    double mx = -100005;
    for(int j = 1; j <= m; j++){
        if(t_current.CjSubZj[j] > mx){
            mx = t_current.CjSubZj[j];
            t_current.kCol = j;
        }
    }

    // Calculating the Ratios
    for(int i = 1; i <= n; i++){
        if(t_current.argumentedMatrix[i][t_current.kCol] > 0){
            t_current.ratio[i] = t_current.argumentedMatrix[i][m+1] / t_current.argumentedMatrix[i][t_current.kCol];
        }else{
            t_current.ratio[i] = 100005;    // a very large value (bcs to find kRow we have to select min of ratios)
        }
    }

    // kRow is the min of all the ratios
    double mn = 100005;
    for(int i = 1; i <= n; i++){
        if(t_current.ratio[i] < mn){
            mn = t_current.ratio[i];
            t_current.kRow = i;
        }
    }

    return t_current;
}

///////////////////////////////////////////////////////////////////////////////////
// We return true if all Elements in Ci - Zi (CjSubZj) are negative (then we have found optimal solution)
bool stopingCondition(table t){
    bool ok = true;
    for(int j = 1; j <= m; j++){
        if(t.CjSubZj[j] > 0){
            ok = false;
            break;
        }
    }
    for(int j = 1; j<=n; j++){
        if(t.CB[j].first > last_var_ind){
            ok = false;
            break;
        }
    }
    return ok;
}

///////////////////////////////////////////////////////////////////////////////////
// check if the current set of sulution is Basic Feasible solution
bool checkIsBFS(vector<double> v){
    bool ok = true;
    for(int j = 1; j <= m; j++){
        if(v[j] < 0){
            ok = false;
            break;
        }
        if(v[j] > 0 && j > last_var_ind){
            ok = false;
            break;
        }
    }
    return ok;
}

int main(){
    ////////////////////////////////////////////////////////////////////
    // I am assuming non negativity to hold & else if x < 0 then we have to replace that
    // variable x by another variable p such that p = (-x), so p > 0     
    cout<<"Enter Number of Inequalities (Apart from Non Negativity): - ";
    cin>>n;
    cout<<"Enter Number of Variables : - ";
    cin>>m;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    
    for(int i = 1; i <= n; i++){
        cout<<"Enter details of "<<i<<"'th inequality"<<endl;
        for(int j = 1; j <= m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
            cin>>Inequalities[i][j];
        }
        cout<<"Enter Value of Rhs of inequation : ";
        cin>>RHS[i];
        cout<<"RHS is less equal(Enter 1) || greater equal(Enter 2) || equal(Enter 3): ";
        cin>>chc[i];
        if(chc[i] == 3){
            equality_count ++;
        }
        // since RHS cannot be negative (multiplying LHS and RHS by -1 and changing equality sign)
        // and we alter chc
        if(RHS[i] < 0){
            for(int j = 1; j <= m; j++){
                Inequalities[i][j] = -1*Inequalities[i][j];
            }
            RHS[i] = -1*RHS[i];
            if(chc[i] == 1)    chc[i] = 2;
            else if(chc[i] == 2) chc[i] = 1;
        }
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter 1 for Maximization || 2 for Minimization : ";
    int choice; cin>>choice;
    if(choice==1)   toMaximize = true;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter the Equation to Maximize / Minimize :- "<<endl;
    for(int j = 1; j <= m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
        cin>>Equation[j];
    }
    // Since we want the eqn to be to maximized
    // so if we have to minimize - we just multiply the equation by -1 to convert to maximization problem
    if(toMaximize == false){
        for(int j = 1; j <= m; j++){
            Equation[j] = -1 * Equation[j];
        }
    }
    // now we can solely focus on maximization problem
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    displayInEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mm = m;
    m = m + n - equality_count;  // adding n slack variables x[m+1], x[m+2], ...., x[m+n]
    last_var_ind = m;
    // NOW WE HAVE TO TRANSFORM THE INEQUALITES TO EQUATIONS
    a = 1;
    int k = 1;
    for(int i = 1; i <= n; i++){
        if(chc[i] == 3){
            // Equality  (add only artificial variables)
            Inequalities[i][m+a] = 1;   // aftificial variables
            a++;
        }else if(chc[i] == 2){
            // >=  (add slack and artificial variables)
            Inequalities[i][mm+k] = -1; // slack variable
            Inequalities[i][m+a] = 1;   // aftificial variables
            a++;
            k++;
        }else if(chc[i] == 1){
            // <=  (add slack variables)
            Inequalities[i][mm+k] = 1;
            k++;
        }
    }
    // NO also insert artificial variables in Objective Function
    for(int i = m+1; i <= m+a-1; i++){
        Equation[i] = -M;
    }
    m = m + a - 1;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Transformed Equations are :- "<<endl;
    displayEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STEP 2 - Construxt the initial simplex table
    arr[0] = constructInitialTable();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE :-"<<endl;
    displayTable(arr[0]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    arr[1] = constructNextSimplexTable(arr[0]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE SIMPLEX TABLE AFTER ITERATION 1 :-"<<endl;
    displayTable(arr[1]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    arr[2] = constructNextSimplexTable(arr[1]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE SIMPLEX TABLE AFTER ITERATION 2 :-"<<endl;
    displayTable(arr[2]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    while(1){

            int menu;
            cout<<"**********************************************************************"<<endl;
            cout<<"MENU : - "<<endl;
            cout<<"Enter 1  :- List of all BFS"<<endl;
            cout<<"Enter 2  :- Print Initial Simplex Table"<<endl;
            cout<<"Enter 3  :- List of Non Basic Variables in Initial Simplex Table"<<endl;
            cout<<"Enter 4  :- List of Basic Variable & Min Ratio in First Iteration"<<endl;
            cout<<"Enter 5  :- Simplex Table Of Second Iteration"<<endl;
            cout<<"Enter -1 :- To Exit"<<endl;
            cout<<"**********************************************************************"<<endl;
            cout<<"Enter Your Choice :- "<<endl;
            cin>>menu;
            cout<<"**********************************************************************"<<endl;            
            if(menu == -1){
                break;
            }

            // Check for invalid choice Input
            if(menu < 1 || menu > 6){
                cout<<"Invalid User Input"<<endl;
            }

            // If choice input = 1
            if(menu == 1){
                for(int i = 0; i <= 2; i++){
                    // storing the solution of the Basic & Non Basic variables in the BFS vector
                    vector<double> BFS(m+1, 0);
                    vector<bool> isBV(m+1, false);
                    for(int k = 1; k <= n; k++){
                        BFS[arr[i].CB[k].first] = arr[i].argumentedMatrix[k][m+1]; 
                        isBV[arr[i].CB[k].first] = true;
                    }

                    // fot the i'th iteration we check if it is BFS (all basic Variables >= 0)
                    bool isBFS = checkIsBFS(BFS);

                    // if solution for i'th iteration is BFS then printing it
                    if(isBFS){
                        // displayTable(arr[i]);
                        cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
                        cout<<"THE BFS OBTAINED IN "<<i<<"'th ITERATION"<<endl;
                        cout<<"BASIC VARIABLES :- "<<endl;
                        for(int k = 1; k <= m; k++){
                            if(isBV[k]){
                                cout<<"x["<<k<<"] = "<<BFS[k]<<endl;
                            }
                        }
                        cout<<"NON - BASIC VARIABLES :- "<<endl;
                        for(int k = 1; k <= m; k++){
                            if(!isBV[k]){
                                cout<<"x["<<k<<"] = "<<BFS[k]<<endl;
                            }
                        }
                        cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
                    }

                }
            }

            // If choice input = 2
            if(menu == 2){
                
                cout<<"---------------------------------------------------------------"<<endl;
                cout<<"Initial Simplex Table"<<endl;
                displayTable(arr[0]);
                cout<<"---------------------------------------------------------------"<<endl;
            }

            // If choice input = 3
            if(menu == 3){
                int i = 0;
                vector<bool> isNonBasic(m+1, true);
                for(int k = 1; k <= n; k++){
                    isNonBasic[arr[i].CB[k].first] = false;
                }
                cout<<"List of Non Basic Variables - "<<endl;
                for(int k = 1; k <= m; k++){
                    if(isNonBasic[k]){
                        cout<<"x["<<k<<"]"<<endl;
                    }
                }
            }

            // If choice input = 4
            if(menu == 4){
                int i = 1;
                cout<<"AFTER "<<i<<"'th ITERATION :-"<<endl;
                cout<<left<<setw(18)<<"Basic-Variable"<<left<<setw(18)<<"Ratio"<<endl;
                cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
                for(int k = 1; k <= n; k++){
                    cout<<left<<setw(18)<<"x["+to_string(arr[i].CB[k].first)+"]"<<left<<setw(18)<<(arr[i].ratio[k]!=100005?to_string(arr[i].ratio[k]):"_")<<endl;
                }
                cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
            }

            // If choice input = 5
            if(menu == 5){
                cout<<"---------------------------------------------------------------"<<endl;
                cout<<"Simplex Table after 2nd Iteration"<<endl;
                displayTable(arr[2]);
                cout<<"---------------------------------------------------------------"<<endl;
            }
    }
}

