%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment1
%}

%{
    Logic :-
    ---------
    I am attaching the pdf where i derived the formulae.
    Assignment1_working.jpg
%}

function assignment1()
    % the Boundary condition for BVP
    x0 = 0; xn = 1;
    y0 = 0; yn = 2;
    % step size
    h = 0.05;

    % The number of steps
    n = (xn - x0)/h - 1;
    %{
        % we do `-1` since we have n-1 unknowns
         -----------------------------------------------------
        x0=0   x1=0.05   x2=0.10                  x(n-1)=0.95  x(n)=1.00
        The unkowns are for y1, y2, .... y(n-1)
        as values of y0 and yn are given
    %}

    x = x0 + h;
    %% Converting it to Matrix Form
    % A - The tri-diagonal matrix
    % B - The RHS matrix
    A = zeros(n,n);
    B = zeros(n,1);
 
    % Let us construct the Tridiagonal System
    %{
        y'' - y = x
        y(0) = 0
        y(1) = 2
        
        Eq 1 --> M(k) = x(k) + y(k)          
        Eq 2 --> M(k-1) + 4*M(k) + M(k+1) = (6/h^2)(y(k-1) - 2*y(k) + y(k+1))
        
        Substitube values from Eq 1 to Eq 2
        [x(k-1) + y(k-1)] + 4*[x(k) + y(k)] + [x(k+1) + y(k+1)] = (6/h^2)(y(k-1) - 2*y(k) + y(k+1))
        Now we have an equation only in the terms of unknows :-
            y(k-1), y(k), y(k+1)
        On rearranging we get :-
        Tri Diagonal System
        -----------------------------------------------------------------------------------------
        y(k-1)(6/h^2 - 1) + y(k)(12/h^2 - 4) + y(k+1)(6/h^2 - 1) = [x(k-1) + 4*x(k) + x(k+1)]
        y(k-1)*ai + y(k)*bi + y(k+1)*c1 = RHS[k]
        -----------------------------------------------------------------------------------------
        [ b1 c1  0  0 0 ...         ][y(1)]     [RHS(1) - a1*y0]
        [ a2 b2 c2  0 0 ...         ][y(2)]     [RHS(2)]
        [ 0  a3 b3 c3 0 ...         ][y(3)]  =  [RHS(3)]
        .
        .
        [ 0  0  0  0  0 ...... bn cn][y(n)]     [RHS(n) - cn*yn]
    %}

    %% Writing it in Matrix form -
    for i = 1:n
        % Functions to calculate values of a, b, c
        % as shown in the above comments
        ai = calcValueOf_a(h);
        bi = calcValueOf_b(h);
        ci = calcValueOf_c(h);

        B(i) = calcValueOf_RHS(x, h);

        %filling up the matrix A and B at each iteration and  
        if i ~= 1
            A(i,i-1) = ai;
        else
            B(i) = B(i) - ai*y0;
        end

        A(i,i) = bi;

        if i ~= n
            A(i,i+1) = ci;
        else
            B(i) = B(i) - ci*yn;
        end
        
        x=x+h;
    end


    %% Solving the Matrix System using Thomas Algorithm
    y=thomasAlgorithm(A,B);
    
    % Printing the Answer
    fprintf("Answer :- \n");
    fprintf("x %d (%f) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%f) = %f\n",i,x0 + i*h, y(i));
    end
    fprintf("x %d (%f) = %f\n",n+1,xn, yn);
    fprintf("\n");

    % calculating the exact solutions
    syms Y;
    Y = dsolve('D2Y-Y=t','Y(0)=0','Y(1)=2');
    ezplot(Y,[0,1])
    
    hold on;
    
    plot([x0:h:xn], [y0;y;yn], '*--');
    hold on;
    %y=2x solves the gicen differential equation and hence we plot those
    %exact values against the calculated values
    title('Solution for h=0.05');
    legend('Calculated values','Actual Values');
    hold off;
end


%%-----------------------------------------------------------------
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
        [1 C1  0   0....     ][del_y1]      [D1]
        [0  1  C2  0....     ][del_y2]      [D2]
        [0  0  1   C3..      ][del_y3]  =   [D3]
        [..                  ]      .
        [...................1][del_yn]      [Dn]
        
        del_yn = Dn
        del_y(n-1) = D(n-1) + C(n-1)*del_yn
        using back substitution... to get the values of del_y
    %}

    for i = rows-1:-1:1
        y(i) = D(i) - C(i)*y(i+1);
    end
end
%---------------------------------------------------------------------



%% functions to calculate the entries of the Tridiagonal Matrix
% --------------------------------------------------------------
function y = calcValueOf_a(h)
    y = 6/h^2 - 1;
end
% --------------------------------------------------------------
function y = calcValueOf_b(h)
    y = - 12/h^2 - 4;
end
% --------------------------------------------------------------
function y = calcValueOf_c(h)
    y = 6/h^2 - 1;
end
% --------------------------------------------------------------
function y = calcValueOf_RHS(x, h)
    y = (x-h) + 4*x + (x+h);
end
% --------------------------------------------------------------


%{
Answer :- 
x 0 (0.000000) = 0.000000
x 1 (0.050000) = 0.077687
x 2 (0.100000) = 0.155693
x 3 (0.150000) = 0.234338
x 4 (0.200000) = 0.313945
x 5 (0.250000) = 0.394838
x 6 (0.300000) = 0.477343
x 7 (0.350000) = 0.561792
x 8 (0.400000) = 0.648521
x 9 (0.450000) = 0.737873
x 10 (0.500000) = 0.830196
x 11 (0.550000) = 0.925846
x 12 (0.600000) = 1.025187
x 13 (0.650000) = 1.128593
x 14 (0.700000) = 1.236447
x 15 (0.750000) = 1.349144
x 16 (0.800000) = 1.467091
x 17 (0.850000) = 1.590708
x 18 (0.900000) = 1.720430
x 19 (0.950000) = 1.856705
x 20 (1.000000) = 2.000000
%}