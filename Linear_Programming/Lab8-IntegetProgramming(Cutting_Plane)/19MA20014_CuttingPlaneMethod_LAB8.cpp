/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
    INTEGER PROGRAMING PROBLEM
    ---------------------------
    CUTTING PLANE METHOD
*/


/*
    LOGIC :-
    ----------------

    STEP 1 :-
     Simplex Method :-
    Changle the inequalities to equalities using slack variables
    Make initial simplex Table and keep iterating till we get a solution
    --- Stopping Cond z-c >= 0 for all 
        else select most minimum in Z-C
        Make the Ratios - min of ratios


    STEP 2 :-
    Now the solution acieved mignt not be integer
    --- So choose the variable whose fraction value is Maximum = k
    Consider the row with hiighest fraction value,
        Conside the positive fractions only
        k = - Summation(f * x) + g1
    Add the newly formed equation as a row

    STEP 3 :-
    Now apply Dual Simplex method to the new set of equation

    Iterate over STEP 2 and STEP 3 until we get all the solutions as integer

*/

#include <bits/stdc++.h>
using namespace std;

long long M = 1000000; 

int a, last_var_ind;    
//last_var_ind is the last index for basic +slack varibles
int equality_count = 0;
int n, m, mm;
double Equation[100];
double Inequalities[100][100];
double RHS[100];
bool toMaximize = false;
bool actualtoMaximize = false;
bool isSlackArtificial[100];
int chc[100];
int iterationCount;
bool solvable = true;

// A structure to define the table used for every iteration
struct table{
    double coefMarix[100][100];
    double rhs[100];
    double argumentedMatrix[100][100];
    pair<int, double> CB[100];
    bool isCB[100];
    double Zj[100];
    double ZjSubCj[100];
    double ratio[100];
    int kRow;
    int kCol;
    int kElement;
};

//////////////////////////////////////////////////////////////////////////////////////////
// an arrSimplexay of these tables ( we can also create arrSimplexay size dynamically insead of fixing max Iterations)
table arrSimplex[100];

//////////////////////////////////////////////////////////////////////////////////////////
// the tables formed after Adding Gomorion Constrains
table arrDual[100];


// START -------------------------------------------------------------------------
////////////////////////  FUNCTIONS FOR SIMPLE SIMPLEX METHODS ////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Display the Inequalities of the Simplex Equation Form
void displayInEqualitiesSimplex(){
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
void displayEqualitiesSimplex(){
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
void displayTableDualSimplex(table t){
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(18)<<"Basic-Variable"<<left<<setw(10)<<"Coef";
    for(int j = 1; j <= m; j++){        
        cout<<left<<setw(10)<<("x["+to_string(j)+"]");
    }
    cout<<left<<setw(10)<<"Value"<<left<<setw(10)<<"Ratio"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    for(int i = 1; i <= n; i++){
        cout<<left<<setw(18)<<("x["+to_string(t.CB[i].first)+"]")<<left<<setw(10)<<t.CB[i].second;
        for(int j = 1; j <= m; j++){
            if(t.isCB[j] == true){
                cout<<left<<setw(10)<<t.argumentedMatrix[i][j];
            }else{
                cout<<left<<setw(10)<<"-";
            }
        }
        cout<<left<<setw(10)<<t.argumentedMatrix[i][m+1]<<left<<setw(10)<<(t.ratio[i]!=100005?to_string(t.ratio[i]):"_")<<endl;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"Z :- ";
    for(int j = 1; j <= m; j++){
        if(t.isCB[j] == true){
            cout<<left<<setw(10)<<t.Zj[j];
        }else{
            cout<<left<<setw(10)<<"-";
        }
    }
    cout<<left<<setw(10)<<t.Zj[m+1];
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<left<<setw(28)<<"Z - C :- ";
    for(int j = 1; j <= m; j++){
        if(t.isCB[j] == true){
            cout<<left<<setw(10)<<t.ZjSubCj[j];
        }else{
            cout<<left<<setw(10)<<"-";
        }
    }
    cout<<endl;
}



///////////////////////////////////////////////////////////////////////////////////
// Function to construct the Initial Table
table constructInitialTableSimplex(){
    table t;
    for(int i = 0; i<= m; i++){
        t.isCB[i] = true;
    }
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
    int k = 1;
    for(int i = 1; i<=m; i++){
        if(isSlackArtificial[i]){
            t.CB[k].first = i;
            t.CB[k].second = Equation[i];
            k++;
        }
    }

    // calculating Zj and ZjSubCj for all m variables
    for(int j = 1; j<= m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.ZjSubCj[j] = (t.Zj[j] - Equation[j]);
    }

    // kCOl is the min of all ZjSubCj
    double mn = 100005;
    for(int j = 1; j <= m; j++){
        if(t.isCB[j] == true && t.ZjSubCj[j] < mn){
            mn = t.ZjSubCj[j];
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


    // kRow is the min of all the Ratios RHS
    mn = 100005;
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
table constructNextTableSimplex(table t_prev){
    table t_current;
    int leavingVariable  = t_prev.CB[t_prev.kRow].first;
    int enteringVariable = t_prev.kCol;
    double pivotElement = t_prev.argumentedMatrix[t_prev.kRow][t_prev.kCol];
    for(int i = 1; i <= m; i++){
        if(i > last_var_ind){
            if(i==leavingVariable || t_prev.isCB[i] == false){
                t_current.isCB[i] = false;
            }
        }else{
            t_current.isCB[i] = true;
        }
    }

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
        t_current.ZjSubCj[j] = (t_current.Zj[j] - Equation[j]);
    }

    // kCOl is the min of all ZjSubCj
    double mn = 100005;
    for(int j = 1; j <= m; j++){
        if(t_current.isCB[j] && t_current.ZjSubCj[j] < mn){
            mn = t_current.ZjSubCj[j];
            t_current.kCol = j;
        }
    }

    // Calculating the Ratios
    for(int i = 1; i <= n; i++){
        if(t_current.argumentedMatrix[i][t_current.kCol] > 0){
            t_current.ratio[i] = t_current.argumentedMatrix[i][m+1] / t_current.argumentedMatrix[i][t_current.kCol];
        }else{
            t_current.ratio[i] = 100005;
        }
    }


    // kRow is the min of all the Ratios RHS
    mn = 100005;
    for(int i = 1; i <= n; i++){
        if(t_current.ratio[i] < mn){
            mn = t_current.ratio[i];
            t_current.kRow = i;
        }
    }

    return t_current;
}



///////////////////////////////////////////////////////////////////////////////////
// We return true if all Elements in Ci - Zi (ZjSubCj) are negative (then we have found optimal solution)
bool stoppingConditionSimplex(table t){
    bool ok = true;
    for(int j = 1; j <= m; j++){
        if(t.isCB[j]){
            if(t.ZjSubCj[j] < 0){
                ok = false;
                break;
            }
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

// END ------------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------------



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
// Function to check if all solutions are integers
bool checkIfIntegerSoln(table t){
    for(int i = 1; i <= n; i++){
        int index = t.CB[i].first;
        double val = t.argumentedMatrix[i][m+1];
        if(index < mm && (int)ceil(val) != (int)floor(val)){
            return false;
        }
    }
    return true;
}




// START -------------------------------------------------------------------------
////////////////////////  FUNCTIONS FOR DUAL SIMPLEX METHODS ////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Function to Display the Table
void displayTableDual(table t){
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
    
    if(t.CB[t.kRow].first > 0 && t.kCol > 0){
        cout<<"Leaving Variable : "<<t.CB[t.kRow].first<<endl;
        cout<<"Entering Variable : "<<t.kCol<<endl;
    }

    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
}



///////////////////////////////////////////////////////////////////////////////////
// Construct the New table at every step of new variable added
table constructNewTable(table t_prev){

    pair<int, double> pp;
    double mnn = -1;
    int n1 = n-1, m1 = m-1;
    for(int i = 1; i <= n1; i++){
        int index = t_prev.CB[i].first;
        double val = t_prev.argumentedMatrix[i][m1+1];
        // fractional part
        double temp = val - (int)floor(val); 
        if(temp < 0) {
            temp = 1 + temp;
        }
        if(temp > mnn){
            mnn = temp;
            pp = {i, mnn};
        }
    }

    table t;

    for(int i=1; i<=n-1; i++) {
        for(int j=1; j<=m-1; j++){
            t.argumentedMatrix[i][j] = t_prev.argumentedMatrix[i][j];
        }
        t.argumentedMatrix[i][m] = 0;
        t.argumentedMatrix[i][m+1] = t_prev.argumentedMatrix[i][m];
    }

    // now we have to construct for the n'th (new) row
    int ii = pp.first;
    for(int j=1; j<=m-1; j++){
        if((int)ceil(t.argumentedMatrix[ii][j]) == (int)floor(t.argumentedMatrix[ii][j])){
            t.argumentedMatrix[n][j] = 0;
            continue;
        }
        double temp = t.argumentedMatrix[ii][j] - (int)floor(t.argumentedMatrix[ii][j]);
        t.argumentedMatrix[n][j] = -1*temp;
    }
    t.argumentedMatrix[n][m] = 1;
    t.argumentedMatrix[n][m+1] = -1*pp.second;
    

    for(int i = 1; i<=n-1; i++){
        t.CB[i].first = t_prev.CB[i].first;
        t.CB[i].second = t_prev.CB[i].second;
    }
    t.CB[n].first = m;
    t.CB[n].second = 0;

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
// We return true if all Elements in Ci - Zi (CjSubZj) are negative (then we have found optimal solution)
bool stopingConditionDual(table t){
    bool ok = true;
    for(int i=1; i<=n; i++){
        if(t.argumentedMatrix[i][m+1] < 0){
            ok = false;
            break;
        }
    }
    return ok;
}



///////////////////////////////////////////////////////////////////////////////////
// Function to construct the simplex table for next iteration from current table
table constructNextTableDual(table t_prev){
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

    double mx = -100005;
    for(int j = 1; j <= m; j++){
        if(t_current.ZjSubCj[j] > mx){
            mx = t_current.ZjSubCj[j];
            t_current.kCol = j;
        }
    }

    // leaving vector --> MIN of all Xb
    double mn = 0;
    t_current.kRow = 0;
    for(int j=1; j<=n; j++){
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

// END ------------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------------
//  ------------------------------------------------------------------------------------



int main(){
    ////////////////////////////////////////////////////////////////////
    // I am assuming non negativity to hold & else if x < 0 then we have to replace that
    // variable x by another variable p such that p = (-x), so p > 0     

    ///////////////////////////////////////////////////////////////////////
    // STEP1
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
    if(choice==1){
        toMaximize = true;
        actualtoMaximize = true;
    }
    else{
        toMaximize = false;
        actualtoMaximize = false;
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Enter the Equation to ";
    if(toMaximize) cout<<"Maximize :- "<<endl;
    else           cout<<"Minimize :- "<<endl;
    for(int j = 1; j <= m; j++){
        cout<<"enter coffecient of x["<<j<<"] : ";
        cin>>Equation[j];
    }

    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    displayInEqualitiesSimplex();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // STEP2
    ///////////////////////////////////////////////////////////////////////
    mm = m;
    m = m + n;  // adding n slack variables x[m+1], x[m+2], ...., x[m+n]
    last_var_ind = m;
    // NOW WE HAVE TO TRANSFORM THE INEQUALITES TO EQUATIONS
    for(int i=0; i<=m; i++){
        isSlackArtificial[i] = false;
    }
    a = 1;
    int k = 1;
    for(int i = 1; i <= n; i++){
        if(chc[i] == 2 || chc[i] == 3){
            // >=  (add slack and artificial variables)
            Inequalities[i][mm+k] = -1; // surplus variable
            Inequalities[i][m+a] = 1;   // aftificial variables
            isSlackArtificial[m+a] = true;
            Equation[mm+k] = 0;
            Equation[m+a] = 1;
            a++;
            k++;
        }else if(chc[i] == 1){
            // <=  (add slack variables)
            Inequalities[i][mm+k] = 1;
            Equation[mm+k] = 0;
            isSlackArtificial[mm+k] = true;
            k++;
        }
    }
    m = m + a - 1;
    // Since we want the eqn to be to maximized
    // so if we have to minimize - we just multiply the equation by -1 to convert to maximization problem
    if(toMaximize == false){
        for(int j = 1; j <= m; j++){
            Equation[j] = -1 * Equation[j];
        }
        toMaximize = true;
    }
    // now we can solely focus on maximization problem

    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"The Transformed Equations are :- "<<endl;
    displayEqualitiesSimplex();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STEP 3 - Construxt the initial simplex table
    arrSimplex[0] = constructInitialTableSimplex();
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE :-"<<endl;
    displayTableDualSimplex(arrSimplex[0]);
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // After Constructing Initial Simplex Table 
    // Iterate to construst simplex table of i+1'th iteration using i'th iteration
    // till the stopping condition is achievied
    int iterationCount = 0; 
    while(stoppingConditionSimplex(arrSimplex[iterationCount]) == false){
        iterationCount++;
        arrSimplex[iterationCount] = constructNextTableSimplex(arrSimplex[iterationCount-1]);
        cout<<"--------------------------------------------------------------------"<<endl;
        cout<<"THE SIMPLEX TABLE AFTER "<<iterationCount<<"'th ITERATION :-"<<endl;
        displayTableDualSimplex(arrSimplex[iterationCount]);
        cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(iterationCount >= 99){
            solvable = false;
            break;
        }
    }  

    // If we Reach to Infeasible or diverging solution
    if(solvable == false) {
        cout<<"No solution can be found"<<endl;
        return 0;
    }  


    bool gotIntegerSolution = checkIfIntegerSoln(arrSimplex[iterationCount]);


    // If we directly found integer solution by Simple Simplex Methos
    if(gotIntegerSolution) {
        cout<<" ----------------------------------------------------------- "<<endl;
        cout<<" ----------------------------------------------------------- "<<endl;
        cout<<" ----------------------------------------------------------- "<<endl;
        cout<<"Optimal Solution Obtained at :-"<<endl;
        for(int i=1; i<=n; i++){
            if(arrSimplex[iterationCount].CB[i].first <= mm){
                cout<<"x["<<arrSimplex[iterationCount].CB[i].first<<"] = "<<arrSimplex[iterationCount].argumentedMatrix[i][m+1]<<endl;
            }
        }
        cout<<"Optimal Value = "<<arrSimplex[iterationCount].Zj[m+1]<<endl;
        cout<<" ----------------------------------------------------------- "<<endl;
        cout<<" ----------------------------------------------------------- "<<endl;
        cout<<" ----------------------------------------------------------- "<<endl;
        return 0;
    }


    table prevTable = arrSimplex[iterationCount];
    int GomorialIneqCount = 0;
    // Now will keep on iterating and adding Gomorian inequalities at each step till we find solution
    while(gotIntegerSolution == false){
        bool solvable = true;
        n++;
        m++;
        GomorialIneqCount ++;
        arrDual[0] = constructNewTable(prevTable);

        cout<<"---------------------------------------------------------"<<endl;
        cout<<"---------------------------------------------------------"<<endl;
        cout<<"---------------------------------------------------------"<<endl;
        cout<<"After Introducing "<<GomorialIneqCount<<" new Inequalities :"<<endl;
        cout<<"---------------------------------------------------------"<<endl;
        cout<<"THE DUAL SIMPLEX TABLE FROM PREVIOUS SET OF EQUATIONS :-"<<endl;

        displayTableDual(arrDual[0]);

        int itercount = 0;
        while(stopingConditionDual(arrDual[itercount]) == false){
            itercount++;
            arrDual[itercount] = constructNextTableDual(arrDual[itercount - 1]);
            
            cout<<"--------------------------------------------------------------------"<<endl;
            cout<<"THE SIMPLEX TABLE AFTER "<<itercount<<"'th ITERATION :-"<<endl;
            displayTableDual(arrDual[itercount]);
            cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
            if(itercount > 99){
                solvable = false;
                break;
            }
        }
        if(solvable == false){
            cout<<"No solution can be found"<<endl;
            return 0;
        }else{
            gotIntegerSolution = checkIfIntegerSoln(arrDual[itercount]);
            if(gotIntegerSolution){
                cout<<" ----------------------------------------------------------- "<<endl;
                cout<<" ----------------------------------------------------------- "<<endl;
                cout<<" ----------------------------------------------------------- "<<endl;
                cout<<"Optimal Solution Obtained at :-"<<endl;
                for(int i=1; i<=n; i++){
                    if(arrDual[itercount].CB[i].first <= mm){
                        cout<<"x["<<arrDual[itercount].CB[i].first<<"] = "<<arrDual[itercount].argumentedMatrix[i][m+1]<<endl;
                    }
                }
                cout<<"Optimal Value = "<<arrDual[itercount].Zj[m+1]<<endl;
                cout<<" ----------------------------------------------------------- "<<endl;
                cout<<" ----------------------------------------------------------- "<<endl;
                cout<<" ----------------------------------------------------------- "<<endl;
                return 0;            
            } else{
                prevTable = arrDual[itercount];
            }
        }

        if(GomorialIneqCount > 50){
            cout<<"No Solution Found"<<endl;
            break;
        }

    }
}



/*********************************************
 * Example Input (Copy and Paste in Terminal)
 
2
2
-1
3
6
1
7
1
35
1
1
7
9

 * Output :-
 -----------------------------------------------------------
 -----------------------------------------------------------
 -----------------------------------------------------------
    Optimal Solution Obtained at :-
    x[2] = 3
    x[1] = 4
    Optimal Value = 55
 -----------------------------------------------------------
 -----------------------------------------------------------
 -----------------------------------------------------------
**************************************************/