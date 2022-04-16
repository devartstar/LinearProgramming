/*
    __________________________
    Name :- Devjit Choudhury
    ROll :- 19MA20014
    __________________________

    Q1. Genreal Program for Guass Seidel Algorithm
*/

/*
    WHEN N < M

        1.  We make number of equations equal to number of variables (n=m) by
            adding extra qeautaions (m-m) of Nom basic variables = 0
            SO WE HAVE M EQUATIONS AND M UNKNOWNS
        ---> OR <---
        2.  We decrease number of variables = n by substituting the non basic variables 
            in system of linear equations to 0
            SO WE HAVE N EQUATIONS AND N UNKNOWNS

        We go through all the arrangements and         
        then we apply Guass Sridel to this set of m*m equations
        
*/
/*
    WHEN N == M
    
    Guass Seidel Algorithm 
        Conditions :-
            1.  The coffecient matrix should be diagonally dominant
            2.  symmetric
            3.  positive define
    Consider system of linear equation :-
        a[1]x + b[1]y + c[1]z = d[1]    cond - a[1] > max(b[1], c[1])   ->  x = (d[1] - b[1]y - c[1]z) / a[1]     --- eq 1
        a[2]x + b[2]y + c[2]z = d[2]    cond - b[2] > max(a[2], c[2])   ->  y = (d[2] - a[2]x - c[2]z) / b[2]     --- eq 2
        a[3]x + b[3]y + c[3]z = d[3]    cond - c[3] > max(a[3], b[3])   ->  z = (d[3] - a[3]x - b[3]y) / c[3]      -- eq 3

        Working Rule :-
        In guass seidel we keep using the calculated value in the nex equation
        and not syop for the next set of iterations to start over to use it (as in Jacobi method)
            FIRST APPROXIMATION :- 
                in eq 1 -> put y=z=0       get x = x1
                in eq 1 -> put x=x1 z=0    get y = y1
                in eq 1 -> put x=x1 y=y1   get z = z1

            SECOND APPROXIMATION :-
                in eq 1 -> put y=y1 z=z1    get x = x2
                in eq 1 -> put x=x2 z=z1    get y = y2
                in eq 1 -> put x=x2 y=y2    get z = z2
                .
                .
                .
                Keep on approximating till we get the desired accuracy (or fix a large number of iterations limit)

*/

#include <bits/stdc++.h>
using namespace std;

int n, m;
// n = number of equations
// m = number of variable
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
// Function to display the system of linear equations
void displayEquations(double arr[][100], double B[]){
    // Lets set the precision to 3 digits after decimal
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
    // Lets set the precision to 3 digits after decimal
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
    iteration_count = 30;   // Can Increase it to bigger value
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
        // !!! UN COMMENT THIS IF TO PRINT ANSWER AFTER EVERY ITERATIONS <-----
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
    // so when we take permutation of arrangements of linear equations 
    // to check if they can be solved by Guass Seidel...
    // we can retrace back easily after rearrangements
    vector<int> temp;
    map<int, vector<double>> m;
    for(int i=1; i<=n; i++){
        temp.push_back(i);
        for(int j=1; j<=n+1; j++){
            m[i].push_back(AAA[i][j]);
        }
    }

    double coef[100][100];  // the coefeciant matrix A after rearrangement
    double rhs[100];        // the right hand side matrix B after similar rearrangement

    // next we loop through all the possible rearrangements (USING next_permutation()) until we get an ans
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
            // cout<<"&&&&&&&&"<<endl;
            // displayEquations(coef, rhs);
            // for(int k =1; k<=n;k++) cout<<fin_ans[k]<<" ";
            // cout<<endl;
            // cout<<"&&&&&&&&"<<endl;
            // cout<<"-----------------------------------------------"<<endl;
            return true;
        }
    }while(next_permutation(temp.begin(), temp.end()));

    return false;
}
//*******************************************************************************


int main(){

    // Lets set the precision to 3 digits after decimal
    cout<<setprecision(5);

    cout<<"Enter number of Equations :- ";  cin>>n;
    cout<<"Enter number of Variables :- ";  cin>>m;

    int m1 = m;

    //---------------------------------------------------------------
    //Defining and initiazing the array
    double A[100][100], B[100];
    memset(A, sizeof(A), 0.0);
    memset(B, sizeof(B), 0.0);
    //--------------------------------------------------------------


    //--------------------------------------------------------------
    // Now lets input the equations from the user
    takeInput(A, B);
    //--------------------------------------------------------------


    //--------------------------------------------------------------
    // Displaying the syatem of LInear Equations
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"The system of Linear Equations GIVEN AS INPUT:- "<<endl;
    displayEquations(A,B);
    cout<<"-------------------------------------------------------"<<endl;
    //--------------------------------------------------------------

    double ans[100];

    if(n > m){
        cout<<"Solution dont exist as Number of Equations greater than number of unknowns"<<endl;
    }else if(n < m){
        // we make the extra variables = 0
        // any extra variable for this q1
        // in q2 and q3 I will make all combinations of extra variables 0
        for(int i=n+1; i<=m; i++){
            ans[i] = 0;
        }
        m = n;
        for(int j = n+1; j<=m; j++){
            for(int i=1; i<=n; i++){
                A[i][j] = 0;
            }
        }
    }

    //--------------------------------------------------------------
    // To get the answer if system of equations can be solved
    
    bool gotAns = getAnswer(A,B,ans);
    //--------------------------------------------------------------
    if(gotAns){
        cout<<"----------------------------------"<<endl;
        cout<<"----------------------------------"<<endl;
        cout<<"FINAL ANSWER :-"<<endl;
        for(int i=1; i<=m1; i++){
            cout<<"X["<<i<<"] = "<<ans[i]<<endl;
        }
        cout<<"----------------------------------"<<endl;
        cout<<"----------------------------------"<<endl;
    }else{
        cout<<"Result didnt Converge"<<endl;
        cout<<"Probably did not find a Diagonally dominant matrix"<<endl;
    }

}