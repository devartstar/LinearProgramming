%{
function thomasAlgorithm()
    a = input("Enter row vector a \n");
    b = input("Enter row vector b \n");
    c = input("Enter row vector c \n");

    B = input("Enter RHS vector \n");
    
    %% --------------------------------------
    % The Tridiagonal Matrix
    n = length(b);
    A = zeros(n,n);
    for i = 2 : n
        A(i, i-1) = a(i-1);
    end
    for i = 1 : n
        A(i, i) = b(i);
    end
    for i = 1 : n-1
        A(i, i+1) = c(i);
    end
    % ---------------------------------------

    %% --------------------------------------
    % Upper Triangular form
    C = zeros(1, n);
    D = zeros(1, n);
    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1, 1);

    for i = 2 : n
        if i ~= n
            C(i) = A(i,i+1) / (A(i,i) - A(i,i-1)*C(i-1));
        end
        D(i) = (B(i) - A(i,i-1)*D(i-1)) / (A(i,i)-A(i,i-1)*C(i-1));
    end
    % ---------------------------------------

    %% --------------------------------------
    % Back substitution to get the answer
    y = zeros(1, n);
    y(n) = D(n);
    for i = n-1 : -1 : 1
        y(i) = D(i) - y(i+1)*C(i);
    end

    %% ---------------------------------------
    % Displaying the answer
    fprintf("Answer :- \n");
    for i = 1 : n
        fprintf("x [%d] = %f\n",i, y(i));
    end
    % ----------------------------------------
end
%}

function thomasAlgo()
    %% Given information
    x0 = -1; xn = 1;
    y0 = 0; yn = 0;
    h = 0.25;
    x = [x0+h : h : xn-h];
    x = x';
    n = (xn - x0) / h - 1;

    A = zeros(n, n);
    B = zeros(1, n); 
    for i = 2 : n-1
        A(i, i-1) = 1/h^2;
        A(i, i) = 1 + x(i)^2 - 2/h^2;
        A(i, i+1) = 1/h^2;
        B(i) = -1;
    end
    B(1) = -1 - y0;
    A(1,1) = 1 + x(1)^2 - 2/h^2;
    A(1,2) = 1/h^2;
    B(n) = -1 - yn;
    A(n,n) = 1 + x(n)^2 - 2/h^2;
    A(n,n-1) = 1/h^2;

    y = thomasAlgorithmSolver(A,B);
    display(y);
end

%------------------------------------------------------------
% Thomas Algorithm to solve the tri diagonal system
function y = thomasAlgorithmSolver(A,B)
    [rows,~]=size(A);   
    C = zeros(1,rows);
    D = zeros(1,rows);

    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1,1);
    
    % modifying the matrix to the desired form
    for i = 2:rows
        if i~=rows
            C(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*C(i-1));
        end
        D(i)=(B(i)-A(i,i-1)*D(i-1))/(A(i,i)-A(i,i-1)*C(i-1));
    end
    y = zeros(rows,1);
    y(rows) = D(rows);

    % Backtracking to reconstruct the solution 
    for i = rows-1:-1:1
        y(i) = D(i)-C(i)*y(i+1);
    end
end
%------------------------------------------------------------

