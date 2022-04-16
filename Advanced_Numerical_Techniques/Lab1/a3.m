%{
    NAME : Devjit Choudhury
    ROLL : 19MA20014
%}


%------------------------------------------------------------
function y = a2()
    x0 = -1;      % left limit of interval
    xn = 1;       % right limit of interval
    y0 = 0;       % Boundary Value
    yn = 0;       % Boundary Value

    h = 0.25;               % interval length
    x = [x0+h : h : xn-h];  % values of x to iterate over
    x = x';                 % Transpose to make it column matrix
    n = length(x);

    %{
        ... the ques and simplification
        y'' + (2 + x^2)y + 2 = 0
------------------------------------------------------
    In general :-
        A(i)y(i)'' + B(i)y(i)' + C(i)y(i) = D(i)
        after substitution :
            y'' = y(i+1) - 2y(i) + y(i-1) / h^2
            y' = y(i+1) - y(i-1) / h
        we get --->
-------------------------------------------------------
        a(i)y(i-1) + b(i)y(i) + c(i)y(i+1) = d(i)
        where :-
            a(i) = A(i)/h^2 - B(i)/(2*h)
            b(i) = C(i) - 2*A(i)/h^2
            c(i) = A(i)/h^2 + B(i)/(2*h)
            d(i) = D(i)
        using the substitution
    %}

    % Accorfing to the differential equation given
    A = zeros(length(x),1) + 1;
    B = zeros(length(x),1);
    C = 2 + x.^2;
    D = zeros(length(x),1) - 2;

    % after converting it into tridiagonal system
    a = A/h^2 - B/(2*h);
    b = C - 2*A/h^2;
    c = A/h^2+B/(2*h);
    d = D;

    d(1) = d(1) - (A(1)/h^2 - B(1)/(2*h)) * y0;
    d(n) = d(n) - (A(n)/h^2 + B(n)/(2*h)) * yn;

    % Generating the tri diagonal matrix
    tri_diagonal_matrix = eye(n);       % To get the identity matrix of n*n
    for i = 1 : n
        tri_diagonal_matrix(i,i) = b(i);
    end
    for i = 1 : n-1
        tri_diagonal_matrix(i, i+1) = c(i);
    end
    for i = 2 : n
        tri_diagonal_matrix(i, i-1) = a(i-1);
    end 
    y = thomasAlgorithm(tri_diagonal_matrix,d);
end
%------------------------------------------------------------

%------------------------------------------------------------
% Thomas Algorithm to solve the tri diagonal system
function y = thomasAlgorithm(A,B)
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

%{
Output :-
---------------

    3.2904
    5.9289
    7.6086
    8.1825
    7.6086
    5.9289
    3.2904
%}

