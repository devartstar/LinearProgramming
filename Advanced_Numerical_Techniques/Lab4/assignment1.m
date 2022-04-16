%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
%}

%{
    Given  f(x) = 0

1.  Find initial root x0 such that f(x0) aprroximately 0
    ie. x0 is near to root of f(x)
    
2.  find f(x0) and f'(x0)
    approximate the root by Newton Raphson method
    using x1 = x0 - f(x0)/f'(x0)

2.  find f(x1) and f'(x1)
    approximate the second root by Newton Raphson method
    using x2 = x1 - f(x1)/f'(x1)

Keep on approximating till we reach a certain accuracy.
General Rule -
    x_n = x_n-1 - f(x_n-1)/f'(x_n-1)

%}

%{
    Given BVP to solve is is :---
    y'' = 4y^3 + 4
    y(1) = 1 
    y(2) = 0.5
    Using 2nd order finite difference method
    -----------------------------------------
     (y(i-1) - 2y(i) + y(i+1))/h^2 = 4*y(i)^3 + 4
    F(y(i-1), y(i), y(i+1)) = (y(i-1) - 2y(i) + y(i+1))/h^2 - 4*y(i)^3 - 4

    a(i) = coef of del F/ del y(i-1) = 1/h^2
    b(i) = coef of del F / del y(i) = -2/h^2 - 12*y(i)^2
    c(i) = coef of del F / del y(i+1) = 1/h^2

    d(i) = - ((y(i-1) - 2y(i) + y(i+1))/h^2 - 4*y(i)^3 - 4)
%}


function assignment1()
    % the Boundary condition for BVP
    x0 = 1; xn = 2;
    y0 = 1; yn = 0.5;
    % step size
    h=0.1;

    % The number of steps
    n = ((xn - x0)/h) - 1;  % -1 since we have n-1 unknowns
    % the unknows are y(1) y(2) y(3) .... y(n)
    yi = zeros(n,1) + 1;
   

    A = zeros(n,n);
    B = zeros(n,1);

    % k = iteration count more iteration cont to decrease error
    for k = 1 : 1000
        % now we will be forming the matrix form using 2nd order Finite
        % difference method ---> AX = B
        for i = 1 : n
            ai = calcValueOf_a(yi,h,i);
            bi = calcValueOf_b(yi,h,i);
            ci = calcValueOf_c(yi,h,i);
            
            % the tridiagonal matrix
            if i ~= 1
                A(i,i-1) = ai;
            end
            A(i,i) = bi;
            if i ~= n
                A(i,i+1)=ci;
            end
            
            % the RHS of matrix eq
            B(i) = calcValueOf_d(yi,h,i);
        end
        % Now we have the tridiagonal matrix system
        %{
            [b(1) c(2)  0   0 ...      ][y(1)]    [B(1)]
            [a(2) b(2) c(2) 0 ...      ][y(2)] =  [B(2)]
            ....    
            [0  0 0 ....      b(n) c(n)][y(n)]    [B(n)]
            solving this tri-diagonal system using thomas algorithm
        %}
        y = thomasAlgorithm(A,B);

        yi = yi + y;
    end


    fprintf("Answer :- \n");
    fprintf("x %d (%d) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%d) = %f\n",i,x0 + i*h, yi(i));
    end
    fprintf("x %d (%d) = %f\n",n+1,xn, yn);
    fprintf("\n");

    plot([x0+h:h:xn-h],yi,'.-');
    title('Solution for h=0.1');
end


%% -------------------------------------------------------------------
% Thomas algorithm to solve the differential equation
%{
    [1 C1  0   0....     ][y1]      [D1]
    [0  1  C2  0....     ][y2]      [D2]
    [0  0  1   C3..      ][y3]  =   [D3]
    [..                  ]      .
    [...................1][yn]      [Dn]
%}
function y = thomasAlgorithm(A,B)
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
    
    % Back substitution to calculate all values of y
    for i = rows-1:-1:1
        y(i) = D(i) - C(i)*y(i+1);
    end
end
%---------------------------------------------------------------------

%% -------------------------------------------------------------------
% function that calc values of a
function y = calcValueOf_a(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
   [n,~]=size(yi);

    if i ~= 1
       ym = yi(i-1);
    else
       ym = 1;
    end
    
    if i ~= n
        yp = yi(i+1);
    else
        yp = 0.5;
    end

    yj = yi(i);
    y = 1/h^2;
end
% ---------------------------------------------------------------------

%% --------------------------------------------------------------------
% function that calculates value of b
function y=calcValueOf_b(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    [n,~]=size(yi);

    if i ~= 1
       ym = yi(i-1);
    else
       ym = 1;
    end
    
    if i ~= n
        yp = yi(i+1);
    else
        yp = 0.5;
    end

    yj=yi(i);
    y = -2/h^2 - 12*yj^2;

end
% -------------------------------------------------------------------

%% --------------------------------------------------------------------
% function that calculates value of c
function y=calcValueOf_c(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    [n,~]=size(yi);

    if i ~= 1
       ym = yi(i-1);
    else
       ym = 1;
    end
    
    if i ~= n
        yp = yi(i+1);
    else
        yp = 0.5;
    end

    yj=yi(i);
    y=1/h^2;
end
% -------------------------------------------------------------------


%% --------------------------------------------------------------------
% function that calculates value of d, rhs
function y = calcValueOf_d(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    [n,~]=size(yi);

    if i ~= 1
       ym = yi(i-1);
    else
       ym = 1;
    end
    
    if i ~= n
        yp = yi(i+1);
    else
        yp = 0.5;
    end


    yj = yi(i);
    y = -1 * ((ym - 2*yj + yp) / h^2 - 4*yj^3 - 4);
end
% -------------------------------------------------------------------


%{
    Answer :- 
    x 0 (1) = 1.000000
    x 1 (1.100000e+00) = 0.746513
    x 2 (1.200000e+00) = 0.549667
    x 3 (1.300000e+00) = 0.399464
    x 4 (1.400000e+00) = 0.291811
    x 5 (1.500000e+00) = 0.225151
    x 6 (1.600000e+00) = 0.198948
    x 7 (1.700000e+00) = 0.213060
    x 8 (1.800000e+00) = 0.267559
    x 9 (1.900000e+00) = 0.362824
    x 10 (2) = 0.500000
%}