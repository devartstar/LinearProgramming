%{
    NAME : Devjit Choudhury
    ROLL : 19MA20014
%}

%{
    Solving a BVP with deriviative Boundary Condition
    using finite differencr method --- results in fictitous points
    
    IN GENERAL :- consider the BVP
    -------------------------------
    y'' + p(x)y' + q(x)y = r(x) -------------<1>
    BC  -                       -------------<2>
            a0.y(a) - a1.y'(a) = G1
            b0.y(b) = G2
    Discretizing <1> we get :-
    -------------------------------
    A(i).Y(i-1) + B(i).Y(i) + C(i).Y(i+1) = h^2.r(i)    ----------<3>
    A(i) = 1 - h.p(i)/2
    b(i) = -2 + h^2.q(i)
    c(i) = 1 + h.p(i)/2
        for i = 1,2,3...n-1 (we get unknows y(0) y(1) ... y(n) y(n+1))
    
    -------------------------------
    Discretizing <2> we get :-  
    a0.y(0) - a1.y'(0) = G1 ---> a0.y(0) - a1.(y(1) - y(-1))/2.h = G1
    y(-1) = 2.h.G1/a1 + y(1) - 2.h.a0.y(0)/a1       --------------<4>
    
    -------------------------------
    i = 0 in equation <3> and substitute equation <4>
    A(0).Y(-1) + B(0).Y(0) + C(0).Y(1) = h^2.r(0)
    F0 -->
    (B(0) - 2.h.a0.A(0)/a1).Y(0) + (A(0) + C(0)).Y(1) = h^2.r(0) - 2.h.G1.A(0)/a1
    
    i = n in equation <3> and substitute BC :-  b0.Y(n+1) = G2
    A(n).Y(n-1) + B(n).Y(n) + C(n).Y(n+1) = h^2.r(n)
    FN -->
    A(n).Y(n-1) + B(n).Y(n) = h^2.r(n) - C(n).G2/b0 
    
%}

%% Solve the BVP 
function assignment1()
    a = 0; b = 1;       % a and b are boundary limits

    %----------------------------------------------------------------
    h = 0.25;           % h is step size
    sol1 = solve(a,b,h);   % %calling solve function with step size 0.25
    n = (b-a) / h;
    sol1(n+1) = 2;

    % Displaying the output
    fprintf("Answer for 0.25:- \n");
    for i = 1:n+1
        fprintf("x %d (%d) = %f\n",i,a + (i-1)*h, sol1(i));
    end
    fprintf("\n");

    %----------------------------------------------------------------
    h = 0.10;
    sol2 = solve(a, b, h);  %calling solve function with step size 0.1
    n = (b-a) / h;
    sol2(n+1) = 2;
    % Displaying the output
    fprintf("Answer for 0.10:- \n");
    for i = 1:n+1
        fprintf("x %d (%d) = %f\n",i,a + (i-1)*h, sol2(i));
    end
    fprintf("\n");

    %----------------------------------------------------------------
    

    % Plotting the result
    xs = [0:0.25:1];
    xs1 = [0:0.10:1];
    plot(xs, sol1 , '--o');
    hold on;
    plot(xs1, sol2 , '-.*');
    legend('h = 0.25', 'h = 0.10');
    title('Plots for different h');
    hold off
   
end

%% Function to convert the equations to Matrix System
function y = solve(a, b, dx)
    % -----------------------------------------------------------------
    % TAKING INPUT
    x0 = a;
    xn = b;
    h = dx;
    %{ 
        % BC -> 
            a0y(x(0)) - a1y'(x(0)) = G1
            b0y(x(n)) = G2
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
    % HANDELING THE FICTITOUS VARIABLES(derived the eqs) & 
    % BUILDING MATRIX SYSTEM

    coefMatrix = zeros(n+1, n+1);
    coefMatrix(1, 1) = calcB(h, x0) - (2 * h * a0 * calcA(h, x0))/a1;
    coefMatrix(1, 2) = calcA(h, x0) + calcC(h, x0);
    for i = 2 : n
        coefMatrix(i, i-1) = A(i-1);
        coefMatrix(i,i) = B(i-1);
        coefMatrix(i,i+1) = C(i-1);
    end
    coefMatrix(n+1, n) = A(n);
    coefMatrix(n+1, n+1) = B(n);
    %-------------------------------------------------------------------

    % Building the Value Matrix in RHS
    rhsMatrix = zeros(n+1, 1);
    rhsMatrix(1) = calcD(h, x0) - (2 * h * G1 * calcA(h,x0))/a1;
    for i = 2 : n
        rhsMatrix(i) = D(i-1);
    end
    rhsMatrix(n+1) = calcD(h, xn) - C(n)*G2/b0;
    %-----------------------------------------------------------------
    % fprintf("Coef Matrix : \n");
    % disp(coefMatrix);
    % fprintf("Value Matrix : \n");
    % disp(rhsMatrix);
    % CALLING FUNCTION TO SOLVE MATRIX SYSTEM USING GUASS ELIMINATION
    y = gaussElimination(coefMatrix, rhsMatrix);
    % ------------------------------------------------------------
end

%% ------------------------------------------------------------
% the function corresponding to p(x) as writen at top
function y = coefOfFirstDeriviative(x)
   % y = p(x)
   y = 0;
end

%% ------------------------------------------------------------
% the function corresponding to q(x) as writen at top
function y = coefOfY(x)
   % y = q(x)
   y = -x;
end

%% ------------------------------------------------------------
% the function corresponding to r(x) as writen at top
function y = rhs(x)
    % y = r(x)
    y = 2;
end

%% ------------------------------------------------------------
% the function corresponding to A(x) as writen at top after DISCRETIZATION
function y = calcA(h, x)
    y = 1 - (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
% the function corresponding to B(x) as writen at top after DISCRETIZATION
function y = calcB(h, x)
    y = -2 + (h * h * coefOfY(x));
end

%% ------------------------------------------------------------
% the function corresponding to C(x) as writen at top after DISCRETIZATION
function y = calcC(h,x)
    y = 1 + (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
% the function corresponding to RHS of the DISCRETIZED EQUATION
% Calculating value of RHS in discretized equation
function y = calcD(h, x)
    y = (h * h * rhs(x));
end



%% Guass Elimination
% we make the lover triangular half 0 and then calc y(n+1)
% then backtrack to previous solutions
function y = gaussElimination(A,B)
    [n,~]=size(A);
    for i = 2 : n
        const = A(i,i-1) / A(i-1,i-1);
        A(i,i-1)=0;
        A(i,i) = A(i,i) - const*A(i-1,i);
        if (i+1 ~= n+1)
            A(i,i+1)=A(i,i+1)-const*A(i-1,i+1);
        end
        B(i)=B(i)-const*B(i-1);
    end
    y=zeros(n,1);

    y(n)=B(n)/A(n,n);

    for i=n-1:-1:1
        y(i)=(B(i)-A(i,i+1)*y(i+1))/A(i,i);
    end  
end

%{
Answer for 0.25:- 
x 1 (0) = -14.923233
x 2 (2.500000e-01) = -10.629925
x 3 (5.000000e-01) = -6.377709
x 4 (7.500000e-01) = -2.199797
x 5 (1) = 2.000000

Answer for 0.10:- 
x 1 (0) = -14.179587
x 2 (1.000000e-01) = -12.551629
x 3 (2.000000e-01) = -10.916221
x 4 (3.000000e-01) = -9.282647
x 5 (4.000000e-01) = -7.656920
x 6 (5.000000e-01) = -6.041821
x 7 (6.000000e-01) = -4.436931
x 8 (7.000000e-01) = -2.838663
x 9 (8.000000e-01) = -1.240265
x 10 (9.000000e-01) = 0.368211
x 11 (1) = 2.000000
%}

