%{
    Solving a BVP with deriviative Boundary Condition
    using finite differencr method --- results in fictitous points
    
    IN GENERAL :- consider the BVP
    -------------------------------
    y'' + p(x)y' + q(x)y = r(x) -------------<1>
    BC  -                       -------------<2>
            a0.y(a) - a1.y'(a) = G1
            b0.y(b) + b1.y'(b) = G2
    Discretizing <1> we get :-
    -------------------------------
    A(i).Y(i-1) + B(i).Y(i) + C(i).Y(i+1) = h^2.r(i)    ----------<3>
    A(i) = 1 - h.p(i)/2
    b(i) = -2 + h^2.q(i)
    c(i) = 1 + h.p(i)/2
        for i = 1,2,3...n (we get unknows y(0) y(1) ... y(n) y(n+1))
    
    -------------------------------
    Discretizing <2> we get :-  
    a0.y(0) - a1.y'(0) = G1 ---> a0.y(0) - a1.(y(1) - y(-1))/2.h = G1
    y(-1) = 2.h.G1/a1 + y(1) - 2.h.a0.y(0)/a1       --------------<4>
    b0.y(n+1) + b1.y'(n+1) = G2 --> b0.y(n+1) + b1.(y(n+2) - y(n))/2.h = G2
    y(n+2) = y(n) - 2.h.b0.y(n+1)/b1 + 2.h.G2/b1    --------------<5>
    
    -------------------------------
    i = 0 in equation <3> and substitute equation <4>
    A(0).Y(-1) + B(0).Y(0) + C(0).Y(1) = h^2.r(0)
    F0 -->
    (B(0) - 2.h.a0.A(0)/a1).Y(0) + (A(0) + C(0)).Y(1) = h^2.r(0) - 2.h.G1.A(0)/a1
    
    i = n+1 in equation <3> and substitute equation <5>
    A(n+1).Y(n) + B(n+1).Y(n+1) + C(n+1).Y(n+2) = h^2.r(n+1)
    FN+1 -->
    (A(n+1) + C(n+1)).Y(n) + (B(n+1) - 2.h.b0.C(n+1)/b1).Y(n+1) = h^2.r(n+1) - 2.h.G2.C(n+1)/b1
    
    F1, F2, .... FN -> using i=1,2,..,n in equation <3>
    ---------------------------------
    WE GET THE COFFECIENT MATRIX :-
    [(B(0) - 2.h.a0.A(0)/a1) (A(0) + C(0))  0    0 ..... ......  0
    [        A1                   B1
    [        0                    A2        B2   C2 .....
    . ....
    . ....
    [        0                    0         0    0 ... (A(n+1) + C(n+1))(B(n+1) - 2.h.b0.C(n+1)/b1) ]
%}

function BVP1()
    % -----------------------------------------------------------------
    % TAKING INPUT
    x0 = 0;
    xn = 1;
    h = 0.25;
    %{ 
        % BC -> 
            a0y(x(0)) - a1y'(x(0)) = G1
            
    %}
    a0 = 1; a1 = -1; G1 = 2;
    b0 = 1; G2 = 2;
    % -----------------------------------------------------------------

    x = [x0 + h : h : xn - h];   % x1, x2, x3, ...., xn
    n = length(x);

    A = zeros(n,1);
    B = zeros(n,1);
    C = zeros(n,1);
    D = zeros(n,1);

    % -----------------------------------------------------------------
    % DISCRETIZATION THE DIFFERENTIAL EQUATION
    % A(i)Y(i-1) + B(i)Y(i) + C(i)Y(i+1) = D(i)
    for i = 1 : n
        A(i) = calcA(h, x(i));
        B(i) = calcB(h, x(i));
        C(i) = calcC(h, x(i));
        D(i) = calcD(h, x(i));
    end
    %-------------------------------------------------------------------

    %-------------------------------------------------------------------
    % HANDELING THE FICTITOUS VARIABLES & BUILDING MATRIX SYSTEM

    coefMatrix = zeros(n+2, n+2);
    coefMatrix(1, 1) = calcB(h, x0) - (2 * h * a0 * calcA(h, x0))/a1;
    coefMatrix(1, 2) = calcA(h, x0) + calcC(h, x0);
    coefMatrix(n+2, n+2) = b0;

    for i = 2 : n+1
        coefMatrix(i, i-1) = A(i-1);
        coefMatrix(i,i) = B(i-1);
        coefMatrix(i,i+1) = C(i-1);
    end
    %-------------------------------------------------------------------

    % Building the Value Matrix in RHS
    rhsMatrix = zeros(n+2, 1);
    rhsMatrix(1) = calcD(h, x0) - (2 * h * G1 * calcA(h,x0))/a1;
    rhsMatrix(n+2) = G2;
    for i = 2 : n+1
        rhsMatrix(i) = D(i-1);
    end
    %-----------------------------------------------------------------
    fprintf("Coef Matrix : \n");
    disp(coefMatrix);
    fprintf("Value Matrix : \n");
    disp(rhsMatrix);
    sol= thomasAlgorithm(coefMatrix, rhsMatrix);
    % Displaying the output
    fprintf("Answer :- \n");
    for i = 1:n+2
        fprintf("x %d = %f\n",i,sol(i));
    end
    fprintf("\n");
    % ------------------------------------------------------------
end

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


function y = coefOfFirstDeriviative(x)
   % y = p(x)
   y = 0;
end

function y = coefOfY(x)
   % y = q(x)
   y = -x;
end

function y = rhs(x)
    % y = r(x)
    y = 2;
end

function y = calcA(h, x)
    y = 1 - (h * coefOfFirstDeriviative(x)) / 2;
end

function y = calcB(h, x)
    y = -2 + (h * h * coefOfY(x));
end

function y = calcC(h,x)
    y = 1 + (h * coefOfFirstDeriviative(x)) / 2;
end

function y = calcD(h, x)
    y = (h * h * rhs(x));
end
