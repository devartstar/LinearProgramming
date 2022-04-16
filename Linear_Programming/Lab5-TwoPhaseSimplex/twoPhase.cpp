/*
    ---------------------------
    Name    - Devjit Choudhury
    Roll No - 19MA20014
    ---------------------------
*/

/*
    LOGIC :-
    
    STEP 1 :-
    Take Input :-
        1.  If Objective Function is Maximization -> conver to Minimization
        2.  If the RHS of inequality is  < 0 --> Multiply -1 both sides

    STEP 2:-
        1. Now will convert the inequality to equality constrains
        2. for (>=) OR (=) constrains - subtract slack variables | add artificial variables
        3. for (<=) add slack variable
    
    STEP 3 :-  << PHASE 1 >>
        1. For the revised objective function 
        2. Make initial simplex table for PHASE 1
        3. Iterate over the simplex table till optimality condition is reached
    
        For Phase 1 :-
            Optimality condition is reached when all Cj -Zj >= 0
            All vbalues of Objective function should be 0 
                            -> then proceed to PHASE 2
                            -> else show not optimal solution

    STEP 4:- <<PHASE 2>>
        1. New Objective function with no artificial variables
        2. Make the initial simplex table for PHASE 2
            using the final simplex table PHASE 1
        3. Iterate over the simplex table till optimality is reached

*/

#include <bits/stdc++.h>
using namespace std;

int a, last_var_ind;    
//last_var_ind is the last index for basic +slack varibles
int equality_count = 0;
int n, m, mm;
double Equation[100];
double OptimalityEquation[100];    // for PHASE 1
double OptimalityEquation2[100];    // for PHASE 2
double Inequalities[100][100];
double RHS[100];
bool toMaximize = false;
bool actualtoMaximize = false;
int chc[100];
int iterationCount;
bool isSlackArtificial[100];
bool isInFeasibleInPahse1 = false;
bool isInFeasibleInPahse2 = false;
string reason = "";

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

table phase1Arr[100];
table phase2Arr[100];

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
            cout<<" + ("<<Equation[j]<<"x["<<j<<"])";
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
            cout<<" + ("<<OptimalityEquation[j]<<"x["<<j<<"])";
        }
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<-1*OptimalityEquation[j]<<"x["<<j<<"])";
        }
    }
    for(int j = mm+1; j <= m; j++){
        cout<<" + ("<<OptimalityEquation[j]<<"x["<<j<<"])";
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
table constructInitialTable(){
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
            t.CB[k].second = OptimalityEquation[i];
            k++;
        }
    }

    // calculating Zj and ZjSubCj for all m variables
    for(int j = 1; j<= m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.ZjSubCj[j] = (t.Zj[j] - OptimalityEquation[j]);
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
table constructNextSimplexTablePhase1(table t_prev){
    table t_current;
    int leavingVariable  = t_prev.CB[t_prev.kRow].first;
    int enteringVariable = t_prev.kCol;
    double pivotElement = t_prev.argumentedMatrix[t_prev.kRow][t_prev.kCol];
    cout<<leavingVariable<<endl;
    cout<<enteringVariable<<endl;
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
    t_current.CB[k] = {enteringVariable, OptimalityEquation[enteringVariable]};
    for(int j = 1; j <= m+1; j++){
        t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[t_prev.kRow][j] / pivotElement;
    }

    // Calculating the Zj and corresponding ZjSubCj
    for(int j = 1; j<= m+1; j++){
        t_current.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t_current.Zj[j] = t_current.Zj[j] + t_current.CB[i].second * t_current.argumentedMatrix[i][j];
        }
        t_current.ZjSubCj[j] = (t_current.Zj[j] - OptimalityEquation[j]);
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
bool stoppingConditionPhase1(table t){
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


///////////////////////////////////////////////////////////////////////////////////
// Function to construct the Initial Table
table constructInitialTablePhase2(){
    table t;
    for(int i=1; i<= m; i++){
        t.isCB[i] = phase1Arr[iterationCount].isCB[i];
    }
    // for(int i = last_var_ind+1; i<= m; i++){
    //     t.isCB[i] = false;
    // }

    for(int i = 1; i<= n; i++){
        for(int j = 1; j <= m; j++){
            t.coefMarix[i][j] = phase1Arr[iterationCount].argumentedMatrix[i][j];
        }
    }

    // Constructing the SOl column
    for(int i = 1; i<= n; i++){
        t.rhs[i] = phase1Arr[iterationCount].rhs[i];
    }

    // Constructing the Argumented Matrix
    for(int i = 1; i<= n; i++){
        for(int j = 1; j <= m; j++){
            t.argumentedMatrix[i][j] = t.coefMarix[i][j];
        }
        t.argumentedMatrix[i][m+1] = phase1Arr[iterationCount].argumentedMatrix[i][m+1];
    }

    // initial base variables last n variables (may include artificial + slack variables)
    for(int i = 1; i<=n; i++){
        t.CB[i].first = phase1Arr[iterationCount].CB[i].first;
        t.CB[i].second = OptimalityEquation[phase1Arr[iterationCount].CB[i].first];
    }

    // calculating Zj and ZjSubCj for all m variables
    for(int j = 1; j<= m+1; j++){
        t.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t.Zj[j] = t.Zj[j] + t.CB[i].second * t.argumentedMatrix[i][j];
        }
        t.ZjSubCj[j] = (t.Zj[j] - OptimalityEquation[j]);
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
table constructNextSimplexTablePhase2(table t_prev){
    table t_current;
    int leavingVariable  = t_prev.CB[t_prev.kRow].first;
    int enteringVariable = t_prev.kCol;
    double pivotElement = t_prev.argumentedMatrix[t_prev.kRow][t_prev.kCol];
    cout<<leavingVariable<<endl;
    cout<<enteringVariable<<endl;
    for(int i = 1; i <= m; i++){
        t_current.isCB[i] = t_prev.isCB[i];
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
    t_current.CB[k] = {enteringVariable, OptimalityEquation[enteringVariable]};
    for(int j = 1; j <= m+1; j++){
        t_current.argumentedMatrix[k][j] = t_prev.argumentedMatrix[t_prev.kRow][j] / pivotElement;
    }

    // Calculating the Zj and corresponding ZjSubCj
    for(int j = 1; j<= m+1; j++){
        t_current.Zj[j] = 0;
        for(int i = 1; i <= n; i++){
            t_current.Zj[j] = t_current.Zj[j] + t_current.CB[i].second * t_current.argumentedMatrix[i][j];
        }
        t_current.ZjSubCj[j] = (t_current.Zj[j] - OptimalityEquation[j]);
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
// Checking condition for Infeasibility - No valid ratio
bool checkInfeasibilityRatio(table t){
    bool isInfeas = true;
    for(int i=1; i<=n; i++){
        if(t.ratio[i] > 0 && t.ratio[i] != 100005){
            isInfeas = false;
            break;
        }
    }
    return isInfeas;
}

//////////////////////////////////////////////////////////////////////////////////
bool checkInfeasibilityPhase1(table t){
    bool isInFeas = false;
    for(int i=1; i<=n; i++){
        if(t.CB[i].first > last_var_ind){
            isInFeas = true;
        }
    }
    return isInFeas;
}


//////////////////////////////////////////////////////////////////////////////////
// Function to call the series of Operations to solve for Phase 1
void solvePhase1(){
     // THe Objective FUnction for Phase 1
    for(int i = 0; i <= last_var_ind; i++){
        OptimalityEquation[i] = 0;
    }
    for(int  i = last_var_ind+1; i <= m; i++){
        OptimalityEquation[i] = -1;
    }

    cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"===================================="<<endl;
    cout<<"PHASE 1"<<endl;
    cout<<"===================================="<<endl;
    cout<<"The Transformed Equations are :- "<<endl;
    displayEqualities();
    cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    
    phase1Arr[0] = constructInitialTable();
    cout<<"********************************************************************************************************************************************"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE FOR PHASE 1 :-"<<endl;
    displayTable(phase1Arr[0]);
    cout<<"********************************************************************************************************************************************"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // cout<<phase1Arr[0].<<endl;
    // After Constructing Initial Simplex Table 
    // Iterate to construst simplex table of i+1'th iteration using i'th iteration
    // till the stopping condition is achievied
    iterationCount = 0;
    while(stoppingConditionPhase1(phase1Arr[iterationCount]) == false){
        iterationCount++;
        phase1Arr[iterationCount] = constructNextSimplexTablePhase1(phase1Arr[iterationCount-1]);
        if(checkInfeasibilityRatio(phase1Arr[iterationCount])){
            isInFeasibleInPahse1 = true;
            reason += "All infeasible rations in Phase 1 at "+to_string(iterationCount)+"'th iteration";
            break;
        }
        cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"THE SIMPLEX TABLE AFTER "<<iterationCount<<"'th ITERATION :-"<<endl;
        displayTable(phase1Arr[iterationCount]);
        cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(iterationCount >= 100){
            break;
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////
// Function to call the series of Operations to solve for Phase 1
void solvePhase2(){
    // for(int i=1; i <= mm; i++){
    //     OptimalityEquation[i] = Equation[i];
    // }
    // for(int i = 1; i <= n; i++){
    //     if(chc[i] == 2 || chc[i] == 3){
    //         OptimalityEquation[mm+i] = -1;
    //     }else{
    //         OptimalityEquation[mm+i] = 1;
    //     }
    // }
    // for(int  i = last_var_ind+1; i <= m; i++){
    //     OptimalityEquation[i] = 0;
    // }
    for(int i=1; i<=m; i++){
        if(i>last_var_ind)  OptimalityEquation[i] = 0;
        else                OptimalityEquation[i] = Equation[i];
    }
    cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"===================================="<<endl;
    cout<<"PHASE 2"<<endl;
    cout<<"===================================="<<endl;
    cout<<"New Optimality Equation is :- "<<endl;
    if(toMaximize){
        cout<<"The Equation To Maximize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<OptimalityEquation[j]<<"x["<<j<<"])";
        }
    }else{
        cout<<"The Equation To Minimize :- "<<endl;
        for(int j = 1; j <= mm; j++){
            cout<<" + ("<<-1*OptimalityEquation[j]<<"x["<<j<<"])";
        }
    }
    for(int j = mm+1; j <= m; j++){
        cout<<" + ("<<OptimalityEquation[j]<<"x["<<j<<"])";
    }
    cout<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    phase2Arr[0] = constructInitialTablePhase2();
    cout<<"********************************************************************************************************************************************"<<endl;
    cout<<"THE INITIAL SIMPLEX TABLE FOR PHASE 2 :-"<<endl;
    displayTable(phase2Arr[0]);
    cout<<"********************************************************************************************************************************************"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    iterationCount = 0; 
    while(stoppingConditionPhase1(phase2Arr[iterationCount]) == false){
        iterationCount++;
        phase2Arr[iterationCount] = constructNextSimplexTablePhase2(phase2Arr[iterationCount-1]);
        if(checkInfeasibilityRatio(phase2Arr[iterationCount])){
            isInFeasibleInPahse2 = true;
            reason += "All infeasible rations in Phase 2 at "+to_string(iterationCount)+"'th iteration";
            break;
        }
        cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"THE SIMPLEX TABLE AFTER "<<iterationCount<<"'th ITERATION :-"<<endl;
        displayTable(phase2Arr[iterationCount]);
        cout<<"------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(iterationCount >= 100){
            break;
        }
    }
}



int main(){
    ////////////////////////////////////////////////////////////////////
    // I am assuming non negativity to hold & else if x < 0 then we have to replace that
    // variable x by another variable p such that p = (-x), so p > 0     

    ///////////////////////////////////////////////////////////////////////
    // STEP1
    ///////////////////////////////////////////////////////////////////////
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
    displayInEqualities();
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

    // cout<<"-------------------------------------------";
    // cout<<"List of slack and artificial"<<endl;
    // for(int i=1; i<=m; i++){
    //     if(isSlackArtificial[i]){
    //         cout<<"x["<<i<<"] ";
    //     }
    // }
    // cout<<endl;
    // cout<<"-------------------------------------------";
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // STEP 3 PHASE 1
    ///////////////////////////////////////////////////////////////////////
    solvePhase1();
    
    // function to check if we got feasible solution for phase one
    if(iterationCount < 100 && isInFeasibleInPahse1 == false){
        ///////////////////////////////////////////////////////////////////////
        // STEP 4 ---> PHASE 2
        ///////////////////////////////////////////////////////////////////////
        solvePhase2();
        if(iterationCount >= 100){
            isInFeasibleInPahse2 = true;
            reason += " Solution in Phase 2 does not converge";
        }
    }else{
        isInFeasibleInPahse1 = true;
        reason += " Solution in Phase 1 does not converge";
    }
    
    if(isInFeasibleInPahse1 || isInFeasibleInPahse2){
        cout<<"---------------------------------------------------"<<endl;
        cout<<"Not Feasable solution"<<endl;
        cout<<"---------------------------------------------------"<<endl;
    }else{
        cout<<"---------------------------------------------------"<<endl;
        if(actualtoMaximize){
            for(int i=1; i<=mm-1; i++){
                cout<<"("<<Equation[i]<<"x["<<i<<"]) + ";
            }
            cout<<"("<<Equation[mm]<<"x["<<mm<<"])"<<endl;
            cout<<"The Value of the Objective Function is :- "<<endl;
            cout<<phase2Arr[iterationCount].Zj[m+1]<<endl;
        }else{
            for(int i=1; i<=mm-1; i++){
                cout<<"("<<-1*Equation[i]<<"x["<<i<<"]) + ";
            }
            cout<<"("<<-1*Equation[mm]<<"x["<<mm<<"])"<<endl;
            cout<<"The Value of the Objective Function is :- "<<endl;
            cout<<-1*phase2Arr[iterationCount].Zj[m+1]<<endl;
        }
    }



        while(1){

            int menu;
            cout<<"****************************************************************************************************"<<endl;
            cout<<"MENU :- "<<endl;
            cout<<"Enter 1  :- Initial table for Phase I"<<endl;
            cout<<"Enter 2  :- Initial table for Phase II"<<endl;
            cout<<"Enter 3  :- Print optimal solution if exists"<<endl;
            cout<<"Enter -1 :- To EXIT"<<endl;
            cout<<"****************************************************************************************************"<<endl;
            cout<<"Enter Your Choice :- "<<endl;
            cin>>menu;
            cout<<"****************************************************************************************************"<<endl;
            if(menu == -1){
                break;
            }
            // Check for invalid choice Input
            if(menu < 1 || menu > 6){
                cout<<"Invalid User Input"<<endl;
            }
            // If choice input = 1 Displaying the initial simplex table for Phase 1
            if(menu == 1){
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"The Initial Simplex Table for Phase 1 is :-"<<endl;
                displayTable(phase1Arr[0]);
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
            }
            // If choice input = 2 Displaying the initial simplex table for Phase 2
            if(menu == 2){
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"The Initial Simplex Table for Phase 2 is :-"<<endl;
                displayTable(phase2Arr[0]);
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
            }
            // If choice input = 2 Displaying the initial simplex table for Phase 2
            if(menu == 3){
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"The Solution using Simplex method :-"<<endl;
                if(isInFeasibleInPahse1 || isInFeasibleInPahse2){
                    cout<<"---------------------------------------------------"<<endl;
                    cout<<"Not Feasable solution"<<endl;
                    cout<<reason<<endl;
                    cout<<"---------------------------------------------------"<<endl;
                }else{
                    cout<<"---------------------------------------------------"<<endl;
                    if(actualtoMaximize){
                        for(int i=1; i<=mm-1; i++){
                            cout<<"("<<Equation[i]<<"x["<<i<<"]) + ";
                        }
                        cout<<"("<<Equation[mm]<<"x["<<mm<<"])"<<endl;
                        cout<<"The Value of the Objective Function is :- "<<endl;
                        cout<<phase2Arr[iterationCount].Zj[m+1]<<endl;
                        cout<<endl;
                        for(int i=1; i<=mm; i++){
                            cout<<"x["<<phase2Arr[iterationCount].CB[i].first<<"] = "<<phase2Arr[iterationCount].argumentedMatrix[i][m+1]<<endl;
                        }
                    }else{
                        for(int i=1; i<=mm-1; i++){
                            cout<<"("<<-1*Equation[i]<<"x["<<i<<"]) + ";
                        }
                        cout<<"("<<-1*Equation[mm]<<"x["<<mm<<"])"<<endl;
                        cout<<"The Value of the Objective Function is :- ";
                        cout<<-1*phase2Arr[iterationCount].Zj[m+1]<<endl;
                        cout<<endl;
                        for(int i=1; i<=mm; i++){
                            cout<<"x["<<phase2Arr[iterationCount].CB[i].first<<"] = "<<phase2Arr[iterationCount].argumentedMatrix[i][m+1]<<endl;
                        }
                    }
                }
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
                cout<<"------------------------------------------------------------------------------------------------"<<endl;
            }

        }
}
