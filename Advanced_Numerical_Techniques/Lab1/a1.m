%{
    NAME : Devjit Choudhury
    ROLL : 19MA20014
%}

function a1()
    a = input("Enter row vector a\n");    % Ex : [2 5 3 7]
    b = input("Enter row vector b\n");    % Ex : [2 4 -9 9 8]
    c = input("Enter row vector c\n");    % Ex : [-1 1 2 4]

    B = input("Enter row vector d\n");    % Ex : [1 2 -3 4 5]

    
    % ------------------------------------------------------
    % Constructing the Tridiagonal matrix
    n = length(b);
    A = zeros(n,n);
    for i = 1 : n
        A(i,i) = b(i);
    end
    for i = 1 : n-1
        A(i, i+1) = c(i);
    end
    for i = 2 : n
        A(i, i-1) = a(i-1);
    end
    % Now we have the tridiagonal matrix
    % AX = B
    % ------------------------------------------------------

    C = zeros(1,n);
    D = zeros(1,n);
    % C stores the values of new ci's after ai's are turned to 0 and bi's
    % are turned to 1 similarly, di's stores these changes accordingly
    
    % ------------------------------------------------------------
    % Thomas Algorithm
    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1,1);

    for i = 2:n
        if (i~=n)
            C(i) = A(i,i+1)/(A(i,i)-A(i,i-1)*C(i-1));
        end
        D(i)=(B(i)-A(i,i-1)*D(i-1))/(A(i,i)-A(i,i-1)*C(i-1));
    end

    y=zeros(n,1);
    y(n)=D(n);

    % Backtracking to reconstruct the solution 
    for i = n-1:-1:1
        y(i)=D(i)-C(i)*y(i+1);
    end
    % ------------------------------------------------------------

    % ------------------------------------------------------------
    % Displaying the output
    fprintf("Answer :- \n");
    for i = 1:n
        fprintf("x %d = %f\n",i,y(i));
    end
    fprintf("\n");
    % ------------------------------------------------------------
end


%{
Input  :-
-------------
[2 5 3 7]
[2 4 -9 9 8]
[-1 1 2 4]
[1 2 -3 4 5]

Output :-
----------
Answer :- 
x 1 = 0.559016
x 2 = 0.118033
x 3 = 0.409836
x 4 = 0.049180
x 5 = 0.581967
%}