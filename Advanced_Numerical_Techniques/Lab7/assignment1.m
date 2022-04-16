%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment1 
%}

%{
    Logic :-
    ---------
    I am attaching the pdf where i derived the formulae.
    Assignment1_working.png


    1. Using the given data points make the tridiagonal system in Mi-1 Mi Mi+1
    2. Solve using Thomas Algorithm to get all the values of Mi
    3. Using these values we can make the differen cubic polynomials in each
    interval
    4. For given x to calculate find in which interval it lies
    or which cubic polynomial to be used to calculate its value
%}

function assignment1()
    %% Data provided in the question :-
    % Interval boundaries
    x0 = 1; xn = 4;
    % Values at different data points
    X = [1;2;3;4];
    Y = [1;2;5;11];
    % Step Size
    h = 1;

    % The number of steps
    n = (xn - x0)/h - 1;
    
    % Natural Spline :- Values of M at end interval = 0
    M0 = 0;
    Mn = 0;

    %% Converting it to Matrix Form
    % A - The tri-diagonal matrix
    % B - The RHS matrix
    A = zeros(n,n);
    B = zeros(n,1);
 
    % Let us construct the Tridiagonal System
    %{
        Using the given data points we find value of M1, M2
        Considering Natural spline M0 = M3 = 0
        Using the formulae
        M(i-1) + 4 * M(i) + M(i+1) = 6/h^2 (y(i-1) - 2y(i) + y(i+1))
        i=1 Eq1 -> 4M1 + M2 = 12
        i=2 Eq2 -> M1 + 4M2 = 18
    %}
    %% Writing it in Matrix form -
    for i=1:n
        %ai is the coefficient of delta Mi-1
        %bi is the coefficient of delta Mi
        %ci is the coefficient of delta Mi+1
        %Bi is the rhs side of the equation
        ai = 1;
        bi = 4;
        ci = 1;
        B(i) = 6 * (Y(i) - 2*Y(i+1) + Y(i+2));
        
        % filling up the entries of the matrix
        if i ~= 1
            A(i,i-1) = ai;
        else 
            B(i) = B(i) - ai*M0;
        end

        A(i,i) = bi;

        if i~=n
            A(i,i+1) = ci;
        else
            B(i) = B(i) - ci*Mn;
        end
    end
    
    %% Solving the Matrix System using Thomas Algorithm
    % y = [M1, M2, M3, M4.... Mn]
    y = thomasAlgorithm(A,B);
    
    M = [M0;y;Mn];

    % To print the values asked y(1.5) and y'(3)
    fprintf("Solution :- \n");
    fprintf("y (%f) = %f \n",1.5,calcValue(1.5,X,Y,M,h));
    fprintf("y' (%f) = %f \n",3,calcDiffValue(3,X,Y,M,h));

    xs = linspace(1,4,30000);
    pp = spline(X,Y);
    ys = ppval(pp,xs);
    plot(xs,ys);
    hold on;
    % ploting the values calculated as discrete points
    plot(xs,calcValue(xs,X,Y,M,h));
    title("Cubic Spline Interpolation");
    legend("Actual values","Calculated values");
    hold off;
end


%% -----------------------------------------------------------------
% function to solve tri-diagonal system using Thomas Algorithm
function y=thomasAlgorithm(A,B)
    [rows,~] = size(A);
    
    C = zeros(1,rows);
    D = zeros(1,rows);
    
    % constructing the first row
    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1,1);
    
    % from row 2 to last
    for i = 2 : rows
        if i ~= rows
            C(i) = A(i,i+1) / (A(i,i) - A(i,i-1)*C(i-1));
        end
        D(i) = (B(i) - A(i,i-1)*D(i-1)) / (A(i,i) - A(i,i-1)*C(i-1));
    end

    y = zeros(rows,1);
    y(rows) = D(rows);
    
    %{
        [1 C1  0   0....     ][M1]      [D1]
        [0  1  C2  0....     ][M2]      [D2]
        [0  0  1   C3..      ][M3]  =   [D3]
        [..                  ]      .
        [...................1][Mn]      [Dn]
        this n corresponds to n-1 (but we already subtracted 1 to n
        Mn = Dn
        My(n-1) = D(n-1) + C(n-1)*Mn
        using back substitution... to get the values of Mi
    %}

    for i = rows-1:-1:1
        y(i) = D(i) - C(i)*y(i+1);
    end
end
%---------------------------------------------------------------------

%% -------------------------------------------------------------------
% calculating in which interval x belongs to
% accordingly we will calculate the value in that cubic polynomial
function y = clacIntervalIndex(x,X,M)
    [n,~] = size(M);
    ind = 1;
    for i = 1:n
        if x < X(i)
            ind = i-1;
            break;
        end
    end
    y = ind;
end   
%-----------------------------------------------------------------

%% -----------------------------------------------------------------
% Function to caluculate values of y at discrete points
function y = calcValue(x,X,Y,M,h)
    i = clacIntervalIndex(x,X,M);
    y = (M(i)*((X(i+1)-x).^3))/(6*h) + (M(i+1)*((x-X(i)).^3))/(6*h) + ((X(i+1) - x)*(Y(i) - ((h*h)*M(i))/6))/h + ((x - X(i))*(Y(i+1) - ((h*h)*M(i+1))/6))/h;
end
%-----------------------------------------------------------------


%% -----------------------------------------------------------------
function y = calcDiffValue(x,X,Y,M,h)
    i = clacIntervalIndex(x,X,M);
    y = M(i)*(-3*(X(i+1)-x).^2+1)/(6*h) - Y(i)/h + M(i+1)*(3*(x-X(i))^2-1)/(6*h) + Y(i+1)/h;
end
%-----------------------------------------------------------------

%{
Solution :- 
y (1.500000) = 1.375000 
y' (3.000000) = 4.666667 
%}