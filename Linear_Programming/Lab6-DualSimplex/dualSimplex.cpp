/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
*/

/*
    When we can use dual simplex :-
    The 2 conditions must be true -
    1. For Maximization Problem :-
                    All Inequalities should be less than equall
                    All coffecients in objective Function is Minimum
    2. For Minimization Problem 
                    All Inequalities should be greater than equall
                    All coffecients in objective Function is Maximum

    Primial :-
        In the primial -- convert in proper format
        If it does not satisfy the proper condition then ->
        Convert it in the Dual form
        If it does not satisfy the proper condition

    STEP1 ---> Take Input from the User
    STEP2 ---> Convert the Primial into desired form
    STEP3 ---> Check if the Primial can be solved
                If not Convert the Primial into Dual Form
    STEP4 ---> Check if the Dual can be solved

    STEP5 ---> Solve the system of Equations

*/

/*
    Example :-
    ------------
    Max Z = -2x1 - 0x2 - x3
    Inequalities
        x1 + x2 - x3 >= 5
        x1 - 2x2 + 4x3 >= 8
    ---------------------------------
    Transformed Ineq :-
        -x1 - x2 + x3 <= -5
        -x1 - 2x2 - 4x3 <= -8
    ----------------------------------
    Max Z = -2x1 - 0x2 - x3 + 0x4 + 0x5
    Transformation to Equalities
        -x1 - x2  + x3  + x4  + 0x5 = -5
        -x1 - 2x2 - 4x3 + 0x4 + x5  = -8
    ------------------------------------

    Initial Table :-
    ______________________________________________________________________________________________________________________
    CB  | (BasicVariable   )    | (x(1),-2) | (x(2),-0) | (x(3),-1) | (x(4),0) | (x(5),0)  | solXb| 
    _____________________________________________________________________________________________________
    0   |      x(4)             |     -1    |     -1    |     1    |     1    |      0     | -5  | 
    0   |      x(5)             |     -1    |      2    |    -4    |    0     |      1     | -8  |  
    ___________________________________________________________________________________________________________
    Zj                         |     0      |     0     |    0     |    0    |     0      |  0
    ___________________________________________________________________________________________________________
    Zj - Cj                    |     2     |     0     |    1     |    0     |     0      | 
    ______________________________________________________________________________________________________________________
    Basic Solution :-
        x1 = 0 x2 = 0 x3 = 0 x4 = -5 x5 = -8
    It is infeasible but optimal

    Leaving Vector --> Min(Xb) --> x(5) = Row(2)
    Entering Vector --> max( Ratio (Zj-Cj / Leaving Row) <-- of only the negative --> x(3) = Col(3)

    Same way of Row Tranformation to make key element = 1 and all other eleemnts as 0
    ______________________________________________________________________________________________________________________
    CB  | (BasicVariable   )    | (x(1),-2) | (x(2),-0) | (x(3),-1) | (x(4),0) | (x(5),0)  | solXb| 
    _____________________________________________________________________________________________________
    0   |      x(4)             |     -5/4   |    -1/2   |     0    |     1    |      1/4  | -7  | 
    -1  |      x(3)             |     1/4    |    -1/2   |    1     |    0     |      -1/4 |  2  |  
    ___________________________________________________________________________________________________________
    Zj                         |     -1/4    |     1/2   |    -1     |    0    |     1/4   |  -2
    ___________________________________________________________________________________________________________
    Zj - Cj                    |     7/4     |     1/2     |    0     |    0   |     1/4   | 
    ______________________________________________________________________________________________________________________
    Basic Solution :-
        x1 = 0 x2 = 0 x3 = 2 x4 = -7 x5 = 0
        Moving from Infeasibility to feasibility
    Leaving Vector --> Min(Xb) --> x(5) = Row(1)
    Entering Vector --> max( Ratio (Zj-Cj / Leaving Row) <-- of only the negative --> x(2) = Col(2)

    ______________________________________________________________________________________________________________________
    CB  | (BasicVariable   )    | (x(1),-2) | (x(2),-0) | (x(3),-1) | (x(4),0) | (x(5),0)  | solXb| 
    _____________________________________________________________________________________________________
    0   |      x(2)             |      5/2   |     1    |     0    |     -2    |    -1/2   | 14  | 
    -1  |      x(3)             |     3/4    |    0     |    1    |     -1     |     -1/2  |  9  |  
    ___________________________________________________________________________________________________________
    Zj                         |     -3/4    |     0    |     1     |    1   |     1/2   |  -9
    ___________________________________________________________________________________________________________
    Zj - Cj                    |     1/2     |     0     |    0     |    1   |     1/2   | 
    ______________________________________________________________________________________________________________________
    Basic Solution :-
        x1 = 0 x2 = 14 x3 = 9 x4 = 0 x5 = 0
    all values are positive

*/



#include <bits/stdc++.h>
using namespace std;


int a, last_var_ind;  // a will have the count of number of Artificial Variables
int equality_count = 0;
int n, m, mm;
double Equation[100];
double Inequalities[100][100];
double RHS[100];

double EquationPrimial[100];
double InequalitiesPrimial[100][100];
double RHSPrimial[100];

double EquationDual[100];
double InequalitiesDual[100][100];
double RHSDual[100];

bool solvabilityPrimial = true;
bool solvabilityDual = true;

bool toMaximize = false;
bool actualtoMaximize = false;

int chc[100];
bool isFreeVariable[100];

// A structure to define the table used for every iteration
struct table{
    double coefMarix[100][100];         // n, m
    double rhs[100];                    // n
    double argumentedMatrix[100][100];  // n, m+1
    pair<int, double> CB[100];          // <index, value>
    double ratio[100];                  // n
    double Zj[100];                     // m
    double ZjSubCj[100];                // m
    int kRow;
    int kCol;
    int kElement;

};

// an array of these tables ( we can also create array size dynamically insead of fixing max Iterations)
table arr[100];

///////////////////////////////////////////////////////////////////////////////////
// Display the Inequalities
void displayInEqualities(){
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;        
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
    }
    for(int j = 1; j <= m; j++){
        cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
    }
    cout<<endl;

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
// Display the Inequalities
void displayInEqualitiesPrimial(){
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<EquationPrimial[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<-1*EquationPrimial[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }

    cout<<"Inequality Constrains :- "<<endl;
    for(int i = 1; i <= n; i++){
        cout<<"Inq "<<i<<" :- ";
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<InequalitiesPrimial[i][j]<<"x["<<j<<"])";
        }
        cout<<(chc[i]==1?" <= ":" >= ")<<RHSPrimial[i]<<endl;
    }
}


///////////////////////////////////////////////////////////////////////////////////
// Display the Inequalities
void displayInEqualitiesDual(){
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<EquationDual[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<-1*EquationDual[j]<<"x["<<j<<"])";
        }
        cout<<endl;
    }

    cout<<"Inequality Constrains :- "<<endl;
    for(int i = 1; i <= n; i++){
        cout<<"Inq "<<i<<" :- ";
        for(int j = 1; j <= m; j++){
            cout<<" + ("<<InequalitiesDual[i][j]<<"x["<<j<<"])";
        }
        cout<<(chc[i]==1?" <= ":" >= ")<<RHSDual[i]<<endl;
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
    cout<<left<<setw(10)<<"Value"<<endl;

    for(int i = 1; i <= n; i++){
        cout<<left<<setw(18)<<("x["+to_string(t.CB[i].first)+"]")<<left<<setw(10)<<t.CB[i].second;
        for(int j = 1; j <= m; j++){
            cout<<left<<setw(10)<<t.argumentedMatrix[i][j];
        }
        cout<<left<<setw(10)<<t.argumentedMatrix[i][m+1]<<endl;
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
        cout<<left<<setw(10)<<t.ZjSubCj[j];
    }
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"Ratios :- ";
    for(int j = 1; j <= m; j++){
        cout<<left<<setw(10)<<(t.ratio[j]!=-100005?to_string(t.ratio[j]):"_");
    }
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    
    cout<<"Leaving Variable : "<<t.CB[t.kRow].first<<endl;
    cout<<"Entering Variable : "<<t.kCol<<endl;

    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
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
        t.CB[i].first = mm + i;
        t.CB[i].second = Equation[m + i];
    }

    // calculating Zj and cjSubZj for all m variables
    for(int j = 1; j<= m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.ZjSubCj[j] = (t.Zj[j] - Equation[j]);
    }

    // kCOl is the max of all ZjSubCj
    double mx = -100005;
    for(int j = 1; j <= m; j++){
        if(t.ZjSubCj[j] > mx){
            mx = t.ZjSubCj[j];
            t.kCol = j;
        }
    }

    // leaving vector --> MIN of all Xb
    double mn = 100005;
    for(int j=1; j<=n; j++){
        if(t.argumentedMatrix[j][m+1] < mn){
            t.kRow = j;
            mn = t.argumentedMatrix[j][m+1];
        }
    }
    int tt = t.kRow;

    // Calculating the Ratios
    for(int j = 1; j <= m; j++){
        if(t.argumentedMatrix[tt][j] < 0 && t.ZjSubCj[j] > 0){
            t.ratio[j] = t.ZjSubCj[j] / t.argumentedMatrix[tt][j];
        }else{
            t.ratio[j] = -100005;
        }
    }


    // kCol is the max of all the ratios
    mx = -100005; 
    for(int j = 1; j <= m; j++){
        if(t.ratio[j] > mx){
            mx = t.ratio[j];
            t.kCol = j;
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

    // calculating Zj and cjSubZj for all m variables
    for(int j = 1; j<= m+1; j++){
        t_current.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t_current.Zj[j] = t_current.Zj[j] + t_current.CB[i].second * t_current.argumentedMatrix[i][j];
        }
        t_current.ZjSubCj[j] = (t_current.Zj[j] - Equation[j]);
    }

    // kCOl is the max of all ZjSubCj
    double mx = -100005;
    for(int j = 1; j <= m; j++){
        if(t_current.ZjSubCj[j] > mx){
            mx = t_current.ZjSubCj[j];
            t_current.kCol = j;
        }
    }

    // leaving vector --> MIN of all Xb
    double mn = 0;
    for(int j=1; j<=n; j++){
        cout<<t_current.argumentedMatrix[j][m+1]<<endl;
        if(t_current.argumentedMatrix[j][m+1] < mn){
            mn = t_current.argumentedMatrix[j][m+1];
            t_current.kRow = j;
        }
    }
    int tt = t_current.kRow;

    // Calculating the Ratios
    for(int j = 1; j <= m; j++){
        if(t_current.argumentedMatrix[tt][j] < 0 && t_current.ZjSubCj[j] > 0){
            t_current.ratio[j] = t_current.ZjSubCj[j] / t_current.argumentedMatrix[tt][j];
        }else{
            t_current.ratio[j] = -100005;
        }
    }


    // kCol is the max of all the ratios
    mx = -100005; 
    for(int j = 1; j <= m; j++){
        if(t_current.ratio[j] > mx){
            mx = t_current.ratio[j];
            t_current.kCol = j;
        }
    }

    return t_current;
}

///////////////////////////////////////////////////////////////////////////////////
// We return true if all Elements in Ci - Zi (CjSubZj) are negative (then we have found optimal solution)
bool stopingCondition(table t){
    bool ok = true;
    for(int i=1; i<=n; i++){
        if(t.argumentedMatrix[i][m+1] < 0 && isFreeVariable[t.CB[i].first] == false){
            // add extra condition for free variables
            ok = false;
            break;
        }
    }
    return ok;
}


///////////////////////////////////////////////////////////////////////////////////
// Function to make Primial in desired form and check if it can be solved
bool solvePrimial(){
    // modifying the equations given as input to Primial form
    for(int i=1; i<=m; i++){
        EquationPrimial[i] = Equation[i];
    }

    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= m; j++){
            InequalitiesPrimial[i][j] = Inequalities[i][j];
        }
        RHSPrimial[i] = RHS[i];
    }

    if(toMaximize){
        for(int i=1; i<=n; i++){
            if(chc[i] == 2){
                for(int j = 1; j <= m; j++){
                    InequalitiesPrimial[i][j] = -1 * InequalitiesPrimial[i][j];
                }
                chc[i] = 1;
                RHSPrimial[i] = -1 * RHSPrimial[i];
            }
        }
    }else{
        for(int i=1; i<=n; i++){
            if(chc[i] == 1){
                for(int j = 1; j <= m; j++){
                    InequalitiesPrimial[i][j] = -1 * InequalitiesPrimial[i][j];
                }
                chc[i] = 2;
                RHSPrimial[i] = -1 * RHSPrimial[i];
            }
        }
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Equation in Primial Form"<<endl;
    displayInEqualitiesPrimial();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // checkSolvability
    // If maximization ---> then all Negative
    // If minimization ---> then all Positive
    
    for(int i=1; i<=m; i++){
        if(toMaximize && EquationPrimial[i] > 0){
            solvabilityPrimial = false;
            return solvabilityPrimial;
        }
        if(!toMaximize && EquationPrimial[i] < 0){
            solvabilityPrimial = false;
            return solvabilityPrimial;
        }
    }


    // Now we can solve the primial set of equations
    if(!toMaximize){
        for(int i = 1; i<=m; i++){
            Equation[i] = -1 * EquationPrimial[i];
        }
        for(int i = 1; i <= n; i++){
            if(chc[i] == 2){
                for(int j = 1; j <= m; j++){
                    Inequalities[i][j] = -1*InequalitiesPrimial[i][j];
                }
                RHS[i] = -1*RHSPrimial[i];
                chc[i] = 1;
            }
        }
        toMaximize = !toMaximize;
    }else{
        for(int i = 1; i<=m; i++){
            Equation[i] = EquationPrimial[i];
        }
        for(int i = 1; i <= n; i++){
            for(int j = 1; j <= m; j++){
                Inequalities[i][j] = InequalitiesPrimial[i][j];
            }
            RHS[i] = RHSPrimial[i];
        }
    }
    

    return solvabilityPrimial;
}

///////////////////////////////////////////////////////////////////////////////////
// Function that will convert Primial to Dual (in desired form) and check if it can be solved
bool solveDual(){
    // now we will use the primial equation to get the new set of equations for dual

    for(int i=1; i<=n; i++){
        EquationDual[i] = RHSPrimial[i];
    }
    
    for(int i=1; i<=m; i++){
        for(int j=1; j<=n; j++){
            InequalitiesDual[i][j] = InequalitiesPrimial[j][i];
        }
        if(toMaximize)  chc[i] = 2;
        else            chc[i] = 1;
        RHSDual[i] = EquationPrimial[i];
    }
    swap(n,m);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Equation in Dual Form"<<endl;
    displayInEqualitiesDual();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i=1; i<=m; i++){
        if(toMaximize && EquationPrimial[i] > 0){
            solvabilityDual = false;
            return solvabilityDual;
        }
        if(!toMaximize && EquationPrimial[i] < 0){
            solvabilityDual = false;
            return solvabilityDual;
        }
    }

    // Now we can solve the Dual set of equations
    if(!toMaximize){
        for(int i = 1; i<=m; i++){
            Equation[i] = -1 * EquationDual[i];
        }
        for(int i = 1; i <= n; i++){
            if(chc[i] == 2){
                for(int j = 1; j <= m; j++){
                    Inequalities[i][j] = -1*InequalitiesDual[i][j];
                }
                RHS[i] = -1*RHSDual[i];
                chc[i] = 1;
            }
        }
    }else{
        for(int i = 1; i<=m; i++){
            Equation[i] = EquationDual[i];
        }
        for(int i = 1; i <= n; i++){
            for(int j = 1; j <= m; j++){
                Inequalities[i][j] = InequalitiesDual[i][j];
            }
            RHS[i] = RHSDual[i];
        }
    }

    return solvabilityDual;
}


int main(){

    //////////////   STEP 1

    ////////////////////////////////////////////////////////////////////  
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
    if(choice==1){
        toMaximize = true;
        actualtoMaximize = true;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter the Equation to Maximize / Minimize :- "<<endl;
    for(int j = 1; j <= m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
        cin>>Equation[j];
    }


    // asking user which variables are free variables
    cout<<"For each variable \nEnter 1 - Free Variable | Enter 2 - Non negativity constrain holds "<<endl;
    for(int i=1; i<=m ;i++){
        int ttt;
        cout<<"For x["<<i<<"] : ";
        cin>>ttt;
        if(ttt==1){
            isFreeVariable[i] = true;
        }else{
            isFreeVariable[i] = false;
        }
    }

    // /////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Input we got :- "<<endl;
    displayInEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    // /////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////   STEP 2  and  STEP 3 and STEP 4
    bool solved = solvePrimial();
    if(!solved){
        cout<<"The Primial Equation is not in solvable form"<<endl;
        solved = solveDual();
        if(!solved){
            cout<<"The Dual Equation is not in solvable form"<<endl;
            // given system of equation cant be solved
            cout<<"Hence given system of equation cant be solved"<<endl;
            return 0;
        }
    }


    //////////////   STEP 5
    // /////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Transformed Inequality to Solve :- "<<endl;
    displayInEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    // /////////////////////////////////////////////////////////////////////////////////////////////////

    
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
    m = m + a - 1;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Transformed Equations are :- "<<endl;
    displayEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STEP 2 - Construxt the initial simplex table
    arr[0] = constructInitialTable();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE :-"<<endl;
    displayTable(arr[0]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // After Constructing Initial Simplex Table 
    // Iterate to construst simplex table of i+1'th iteration using i'th iteration
    // till the stopping condition is achievied
    int iterationCount = 0; 
    while(stopingCondition(arr[iterationCount]) == false){
        iterationCount++;
        arr[iterationCount] = constructNextSimplexTable(arr[iterationCount-1]);
        cout<<"--------------------------------------------------------------------"<<endl;
        cout<<"THE SIMPLEX TABLE AFTER "<<iterationCount<<"'th ITERATION :-"<<endl;
        displayTable(arr[iterationCount]);
        cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(iterationCount >= 100){
            break;
        }
    }


    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The optimal solution of Objective Function "<<endl;   

    // Minimization Problem (-Ans) else Ans (because we had earlier transformed the equation)
    if(actualtoMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
        cout<<"Maximum Value : "<<arr[iterationCount].Zj[m+1]<<endl;
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<-1*Equation[j]<<"x["<<j<<"])";
        }
        cout<<endl;
        cout<<"Minimum Value : "<<-1*arr[iterationCount].Zj[m+1]<<endl;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;


    while(1){
        int menu;
        cout<<"****************************************************************************************************"<<endl;
        cout<<"MENU :- "<<endl;
        cout<<"Enter 1  :- Initial Table"<<endl;
        cout<<"Enter 2  :- Incomming Variable | Outgoing Variable | Objective Value"<<endl;
        cout<<"Enter 3  :- Optimal Value of Problem"<<endl;
        cout<<"Enter -1 :- To EXIT"<<endl;
        cout<<"****************************************************************************************************"<<endl;
        cout<<"Enter Your Choice :- "<<endl;
        cin>>menu;
        cout<<"****************************************************************************************************"<<endl;

        if(menu == -1){
            break;
        }
        
        // Check for invalid choice Input
        if(menu < 1 || menu > 3){
            cout<<"Invalid User Input"<<endl;
        }

        // if choice input = 1
        if(menu == 1){
            cout<<"----------------------------------------------------------------------------------------------"<<endl;
            cout<<"Initial Table :- "<<endl;
            displayTable(arr[0]);
            cout<<"----------------------------------------------------------------------------------------------"<<endl;
        }

        // if choice input = 2
        if(menu == 2){
            cout<<"----------------------------------------------------------------------------------------------"<<endl;
            for(int k = 1; k<= iterationCount; k++){
                cout<<"----------------------------------------------------------------------------------------------"<<endl;
                cout<<"For the "<<k<<"'th iteration"<<endl;
                cout<<"Leaving Variable :- x["<<arr[k].CB[arr[k].kRow].first<<"]"<<endl;
                cout<<"Leaving Variable :- x["<<arr[k].kCol<<"]"<<endl;
                cout<<"Optimal Value"<<endl;
                for(int j = 1; j <= mm; j++){
                    cout<<" + ("<<-1*Equation[j]<<"x["<<j<<"])";
                }
                cout<<" = "<<arr[k].Zj[m+1]<<endl;
                cout<<"----------------------------------------------------------------------------------------------"<<endl;
            }
            cout<<"----------------------------------------------------------------------------------------------"<<endl;
        }

        // if choice input = 3
        if(menu == 3){
            cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
            cout<<"The optimal solution of Objective Function "<<endl;   

            // Minimization Problem (-Ans) else Ans (because we had earlier transformed the equation)
            if(actualtoMaximize){
                cout<<"The Equation To Maximize :- "<<endl;
                for(int j = 1; j <= mm; j++){
                    cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
                }
                cout<<endl;
                cout<<"Maximum Value : "<<arr[iterationCount].Zj[m+1]<<endl;
            }else{
                cout<<"The Equation To Minimize :- "<<endl;
                for(int j = 1; j <= mm; j++){
                    cout<<" + ("<<-1*Equation[j]<<"x["<<j<<"])";
                }
                cout<<endl;
                cout<<"Minimum Value : "<<-1*arr[iterationCount].Zj[m+1]<<endl;
            }
            cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
        }
    }

}