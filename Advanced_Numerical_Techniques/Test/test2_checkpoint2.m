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
    y''(y - y') = e^x - 1
    y(0) = 1 
    y(1) = 0
    Using 2nd order finite difference method
    -----------------------------------------
     F = ((y(i-1) - 2y(i) + y(i+1))/h^2) (y(i) - (y(i+1)-y(i-1))/2h) =
     exp(x(i)) - 1

    a(i) = coef of del F/ del y(i-1) 
    b(i) = coef of del F / del y(i) 
    c(i) = coef of del F / del y(i+1)

    d(i) = - F
%}


function assignment1()
    % the Boundary condition for BVP
    x0 = 0; xn = 1;
    y0 = 1; yn = 0;
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
            B(i) = calcValueOf_d(yi,h,i,x0 + i*h);
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
    fprintf("x %d (%f) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%f) = %f\n",i,x0 + i*h, yi(i));
    end
    fprintf("x %d (%f) = %f\n",n+1,xn, yn);
    fprintf("\n");

    
    title('Solution for h=0.1');

    %% using bvp4c to get exact solution
    f=@(x,y)[y(2);(exp(x)-1)/(y(1)-y(2))];
    bc=@(yf,yl)[yf(1)-1;yl(1)-0];
    xmesh=linspace(0,1,1000);
    solinit=bvpinit(xmesh,[0 1]);
    sol=bvp4c(f,bc,solinit);
    plot(sol.x,sol.y(1,:));
    hold on;
    plot([x0:h:xn],[y0;yi;yn],'-.o');
    hold off
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
        yp = 0;
    end

    yj = yi(i);
    y = 1000*ym - 900*yj;
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
        yp = 0;
    end

    yj=yi(i);
    y = -900*ym - 400*yj + 1100*yp;

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
        yp = 0;
    end

    yj=yi(i);
    y=1100*yj - 1000*yp;
end
% -------------------------------------------------------------------


%% --------------------------------------------------------------------
% function that calculates value of d, rhs
function y = calcValueOf_d(yi,h,i,xi)
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
        yp = 0;
    end

    
    yj = yi(i);
    y = ((ym - 2*yj + yp)/(h*h))*(yj - (yp-ym)/(2*h)) - exp(xi) + 1;
    y = -1 * y;
end
% -------------------------------------------------------------------


%{
Answer :- 
x 0 (0.000000) = 1.000000
x 1 (0.100000) = 0.883805
x 2 (0.200000) = 0.768125
x 3 (0.300000) = 0.653599
x 4 (0.400000) = 0.541028
x 5 (0.500000) = 0.431434
x 6 (0.600000) = 0.326149
x 7 (0.700000) = 0.226960
x 8 (0.800000) = 0.136393
x 9 (0.900000) = 0.058339
x 10 (1.000000) = 0.000000
%}

