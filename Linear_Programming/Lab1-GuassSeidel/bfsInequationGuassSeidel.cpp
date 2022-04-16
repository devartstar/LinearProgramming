/*
    __________________________
    Name :- Devjit Choudhury
    ROll :- 19MA20014
    __________________________

    Q4 Solve the Linear Programming Problem.
*/

/*
    Consider system of inequalities :-
        x[1] + 4x[2] <= 24
        3x[1] + x[2] <= 21
        x[1] + x[2] <= 9
    subject to confition :-{ x[1], x[2] >=0 }
    Objective Function :- Maximize(2x[1] + 5x[2])

        lets convert the inequality into equality by adding slack variables
         x[1] + 4x[2] + x[3] = 24
        3x[1] +  x[2] + x[4] = 21
         x[1] +  x[2] + x[5] = 9
    subject to confition :-{ x[1], x[2], x[3], x[4], x[5] >=0 }

    Now we have a system of linear equations with <n equations> and <m variables>
    n = number of equations = 3
    m = number of variables = 5
    .... so we will have 5C3 => 10 cases of solutions

    Among them identify the basic feasible solution and apply it in objective function
     And the solution set for which we get maximum value

     <nan> - means the system cant be solved using Guass - Seidel method
*/

#include <bits/stdc++.h>
using namespace std;

int n, m, m1;
// n = number of equations
// m = number of variable (m1 is compy of m)
int iteration_count = 100;

//*******************************************************************************
// Function to take input from the user
void takeInput(double A[][100], double B[]){
    for(int i = 1; i <= n; i++){
        cout<<"---------------------------------"<<endl;
        cout<<"Enter the "<<i<<"th equation"<<endl;
        for(int j = 1; j <= m; j++){
            cout<<"Enter the coffecient for x["<<j<<"] : ";
            cin>>A[i][j];
        }
        cout<<"Enter the value in RHS :- ";
        cin>>B[i];
    }
}
//*******************************************************************************

//*******************************************************************************
// function to conver inequality to equality in LPP questions
void convertInequalityToEquality(double A[][100], double B[]){
    int jj = m+1;
    for(int i = 1; i<=n; i++){
        A[i][jj] = 1;
        jj++;
    }
}
//*******************************************************************************

//*******************************************************************************
// Function to display the system of linear equations
void displayEquations(double arr[][100], double B[]){
    // Lets set the precision to 5 digits after decimal
    cout<<setprecision(5);

    for(int i = 1; i<=n; i++){
        cout<<"Eq "<<i<<" ";
        for(int j = 1; j<=m; j++){
            if(j < m)
                cout<<"("<<arr[i][j]<<"x["<<j<<"]) + ";
            else
                cout<<"("<<arr[i][j]<<"x["<<j<<"])";
        }
        cout<<" = "<<B[i]<<endl;
    }
}
//*******************************************************************************


//*******************************************************************************
// Check if coffecient matrix is diagonally dominant
bool checkDiagonalDominance(double A[][100]){
    for(int i =1; i<=n; i++){
        int sum = 0;
        for(int j = 1; j<=n; j++){
            if(i!=j){
                sum += abs(A[i][j]);
            }
        }
        if(abs(A[i][i]) < sum){
            return false;
        }
    }
    return true;
}
//*******************************************************************************



//*******************************************************************************
// check if the ans calculated is valid answer
bool checkIfAnswerSatisfies(double A[][100], double B[], double ans[]){
    for(int j=1; j<=n; j++){
        if(ans[j] > 10000){
            return false;
        }
    }
    for(int i = 1; i<=n; i++){
        double sum = 0;
        for(int j=1; j<=n; j++){
            sum = sum + (A[i][j] * ans[j]);
        }
        if(abs(sum - B[i]) > 0.1){
            return false;
        }
    }
    return true;

}
//*******************************************************************************


//*******************************************************************************
// Function to solve system of linear equations using Guass Seidel
// it should have equal number of equations and variables
bool solveGuassSeidel(double A[][100], double B[], double X_CURR[]){
    // Lets set the precision to 5 digits after decimal
    cout<<setprecision(5);

    double X_PREV[100];
    //X_CURR -> stores the answer from the Current iterations
    //X_PREV -> stores the answer from the Previous iterations

    // giving the inital values to 1
    for(int i=1; i<=n ;i++){
        X_CURR[i] = 1;
        X_PREV[i] = 1;
    }

    // number of iterations (cal also have error allowed as parameter for number of itrerations)
    iteration_count = 30;   // !!! Increase it to bigger value
    int temp = iteration_count;
    while(iteration_count--){
        for(int i = 1; i <= n; i++){
            double val = B[i];
            for(int j = 1; j<=n; j++){
                if(j <  i){
                    val = val - (A[i][j] * X_CURR[j]);
                }
                if(j > i){
                    val = val - (A[i][j] * X_PREV[j]);
                }
            }
            X_CURR[i] = (val / A[i][i]);
        }

        // X_prev = x_curr for the next iterations
        for(int i=1; i<=n; i++){
            X_PREV[i] = X_CURR[i];
        }
        /*
        cout<<"-----------------------------------------------------"<<endl;
        cout<<"Interation Number : "<<temp-iteration_count<<endl;
        cout<<"[ ";
        for(int i=1; i<=n; i++){
            cout<<"X["<<i<<"]= "<<X_CURR[i]<<" || ";
        }
        cout<<" ]"<<endl;
        cout<<"-----------------------------------------------------"<<endl;
        */
    }
    bool convergentAns = checkIfAnswerSatisfies(A,B,X_PREV);
    return convergentAns;
}
//*******************************************************************************


//*******************************************************************************
// Function to make coffecient matrix diagonally dominant
// SO WILL LOOP THROUGH ALL PERMUTATIONS and if solution converges in any one will return true
bool getAnswer(double A[][100], double B[], double fin_ans[]){
    /*
                A   X   =  B                    AAA
        [a1 b1 c1] [x1]   [y1]              [a1 b1 c1 y1]
        [a2 b2 c2] [x2] = [y2] ------>      [a2 b2 c2 y2]
        [a1 b2 c3] [x3]   [y3]              [a3 b3 c3 y3]
        so when we permute to have diagonal dominance of A we automically have the corresponding B
    */

    // Creating the coeficient + rhs combined maxtix
    double AAA[n+1][n+2];
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            AAA[i][j] = A[i][j];
        }
        AAA[i][n+1] = B[i];
    }


    // maping the index with the correspoing equation
    vector<int> temp;
    map<int, vector<double>> m;
    for(int i=1; i<=n; i++){
        temp.push_back(i);
        for(int j=1; j<=n+1; j++){
            m[i].push_back(AAA[i][j]);
        }
    }

    double coef[100][100];
    double rhs[100];
    do{
        // create the coffecient matrix --- according to the permutation
        for(int i=0;i<n;i++){
            for(int j=0; j<n; j++){
                coef[i+1][j+1] = m[temp[i]][j];
            }
            rhs[i+1] = m[temp[i]].back();
        }

        // apply guass siedel to this equation
        // coef * X = rhs
        // Now we solve for this matrix using guass seidel
        double ans[100];
        bool is_convergent_system = solveGuassSeidel(coef, rhs, ans);
        if(is_convergent_system){
            for(int i=1; i<=n; i++){
                fin_ans[i] = ans[i];
            }
            // cout<<"-----------------------------------------------"<<endl;
            // cout<<"The system of Linear Equations AFTER REARRANGEMENT TO APPLY GUASS SEIDEL:- "<<endl;
            // displayEquations(coef, rhs);
            // cout<<"-----------------------------------------------"<<endl;
            return true;
        }
    }while(next_permutation(temp.begin(), temp.end()));

    return false;
}
//*******************************************************************************

//*******************************************************************************
void copy2DArray(double A[][100], double copy_A[][100]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m1; j++){
            copy_A[i][j] = A[i][j];
        }
    }
}
//*******************************************************************************

//*******************************************************************************
double calculateObjectiveFunction(double objective[], double Z[]){
    double ans = 0.0;
    for(int j = 1; j<=m1; j++){
        ans = ans + (objective[j]*Z[j]);
    }
    return ans;
}
//*******************************************************************************

int main(){

    // Lets set the precision to 5 digits after decimal
    cout<<setprecision(5);

    cout<<"Enter number of Equations :- ";  cin>>n;
    cout<<"Enter number of Variables :- ";  cin>>m;

    m1 = m;

    //---------------------------------------------------------------
    //Defining and initiazing the array
    double A[100][100], B[100];
    memset(A, sizeof(A), 0.0);
    memset(B, sizeof(B), 0.0);
    //--------------------------------------------------------------


    //--------------------------------------------------------------
    // Now lets input the equations from the user
    takeInput(A, B);
    convertInequalityToEquality(A,B);
    //--------------------------------------------------------------


    //--------------------------------------------------------------
    //Objective Function :-
    cout<<"Enter the objective function "<<endl;
    double objective[m1+1];
    for(int j = 1; j <= m1; j++){
        cout<<"Enter the coffecient for x["<<j<<"] : ";
        cin>>objective[j];
    }
    // MAXIMIZE - 1 || MINIMIZE - 2
    cout<<"Enter 1 for Maximization | 2 for Minimizatoin"<<endl;
    int cond;   cin>>cond;
    //--------------------------------------------------------------

    m = m + n;
    m1 = m;

    //--------------------------------------------------------------
    // Displaying the syatem of LInear Equations
    cout<<"---------------------------------------------------------------"<<endl;
    cout<<"The system inequalities after converting to Linear Equations:- "<<endl;
    displayEquations(A,B);
    //--------------------------------------------------------------



    if(n > m){
        cout<<"Solution dont exist as Number of Equations greater than number of unknowns"<<endl;
    }else if(n <= m){
        // we make the extra variables = 0
        m = n;

        int temp_arr[m1+1];
        for(int i=1; i<=m1; i++){
            if(i<=n)    temp_arr[i] = 0;
            else        temp_arr[i] = 1;
        }
        //Now we have an permutation of n 'zeros' and m-n 'ones'
        // the index corresponding to 'zeros' are basic variables   (to be calculated)
        // the index corresponding to 'ones' are non-basic variables (to be assigned 0)
        double objectiveValue;
        if(cond == 1)    objectiveValue = -10000000;
        else if(cond == 2)  objectiveValue = 10000000;

        do{
            double copy_A[100][100];
            int indx[n+1];
            int j_pos = 1;
            for(int k = 1; k<=m1; k++){
                if(temp_arr[k] == 0){
                    indx[j_pos] = k;
                    // x[k] is non basic variable
                    // x[k] = 0 --> so all copy_A[i=1...n][k] = 0
                    // cout<<"k "<<k<<endl;
                    for(int i = 1; i<=n; i++){
                        copy_A[i][j_pos] = A[i][k];
                    }
                    j_pos++;
                }
            }
            cout<<"-----------------------------------------"<<endl;

            // !!!!!!!! CHECK IF EVERYTHING IS RESET BEFORE THE NEXT CALL <-----
            double ans[100];
            bool gotAns = getAnswer(copy_A,B,ans);
            // Non basic Variables
            cout<<"Non Basic Variables :  ";
            for(int i=1; i<=m1; i++){
                if(temp_arr[i] == 1){
                    cout<<"x["<<i<<"] = 0 | ";
                }
            }
            cout<<endl;
            if(gotAns){
                bool isBasicFeasibleSolution = true;
                // Basic Variables --- cant be solved by guass siedel method
                cout<<"Basic Variables : ";
                for(int i=1; i<=n; i++){
                    cout<<"x["<<indx[i]<<"] = "<<ans[i]<<" | ";
                    if(ans[i] < 0){
                        isBasicFeasibleSolution = false;
                    }
                }
                cout<<endl;
                if(isBasicFeasibleSolution == true){
                    double temp_ans = calculateObjectiveFunction(objective,ans);
                    cout<<"Value of Objective function :- "<<temp_ans<<endl;
                    cout<<"THIS IS BASIC FEASIBLE SOLUTION"<<endl;
                    if(cond==1) objectiveValue = max(objectiveValue, temp_ans);
                    else if(cond==2)    objectiveValue = min(objectiveValue, temp_ans);
                }else{
                    cout<<"THIS IS NOT BASIC FEASIBLE SOLUTION"<<endl;
                }
            }else{
                // Basic Variables --- cant be solved by guass siedel method
                cout<<"Basic Variables : Values can be solved using guass seidel method"<<endl;
            }
            cout<<"-----------------------------------------"<<endl;
        }while(next_permutation(temp_arr+1, temp_arr+m1+1));
        cout<<"-------------------------------------------------------------------------"<<endl;
        cout<<"-------------------------------------------------------------------------"<<endl;
        cout<<"FINAL ANSWER FOR OBJECTIVE FUNCTION = "<<objectiveValue<<endl;
        cout<<"-------------------------------------------------------------------------"<<endl;
        cout<<"-------------------------------------------------------------------------"<<endl;
    }


}