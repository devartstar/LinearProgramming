%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
%}

%{
    Given BVP to solve is is :---
    y'' = 4yy'
    y(0) = 0.5 
    y(1) = 1
    Using 2nd order finite difference method
    -----------------------------------------
     (y(i-1) - 2y(i) + y(i+1))/h^2 = 4*y(i)*(y(i+1) - y(i-1))/2
    F(y(i-1), y(i), y(i+1)) = (y(i-1) - 2y(i) + y(i+1))/h^2 - 4*y(i)*(y(i+1) - y(i-1))/2

    a(i) = coef of del F/ del y(i-1) = 1/h^2 + 2*y(i)
    b(i) = coef of del F / del y(i) = -2/h^2 - 2*(y(i+1) - y(i-1))
    c(i) = coef of del F / del y(i+1) = 1/h^2 - 2*y(i)

    d(i) = - ((y(i-1) - 2y(i) + y(i+1))/h^2 - 4*y(i)*(y(i+1) - y(i-1))/2)
%}

function assignment2()
    % the Boundary condition for BVP
    x0 = 0; xn = 1;
    y0 = 0.5; yn = 1;
    % step size
    h = 0.1;
    
    % The number of steps
    n = (xn-x0)/h - 1; % -1 since we have n-1 unknowns
    % the unknows are y(1) y(2) y(3) .... y(n)
    yi = zeros(n,1) + 1;
    
    A=zeros(n,n);
    B=zeros(n,1);

    %% using bvp4c to get exact solution
    f = @(x,y)[y(2);4*y(1)*y(2)];
    bc = @(ya,yb)[ya(1)-0.5;yb(1)-1];
    xmesh = linspace(0,1,1000);
    solinit = bvpinit(xmesh,[1 0]);
    sol = bvp4c(f,bc,solinit);

    %% Plotting the graph for the exact solution
    % ploting the exact solution
    plot(sol.x,sol.y(1,:));
    hold on;

    % k = iteration count more iteration cont to decrease error
    for k = 1:1000
        % now we will be forming the matrix form using 2nd order Finite
        % difference method ---> AX = B
        for i = 1:n
            ai = calcValueOf_a(yi,h,i);
            bi = calcValueOf_b(yi,h,i);
            ci = calcValueOf_c(yi,h,i);

            % the tridiagonal matrix
            if i~=1
                A(i,i-1)=ai;
            end
            A(i,i)=bi;
            if i~=n
                A(i,i+1)=ci;
            end

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
        y=thomasAlgorithm(A,B);

        %using the delta values obtained to update the approximate
        %solutions
        yi=yi+y;        
    end

    %% Displaying the answer
    fprintf("Answer :- \n");
    fprintf("x %d (%d) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%d) = %f\n",i,x0 + i*h, yi(i));
    end
    fprintf("x %d (%d) = %f\n",n+1,xn, yn);
    fprintf("\n");


    %% Plotting the graph for the approximated solution

    % ploting the calculated value
    plot([x0:h:xn], [y0;yi;yn],'*--');
    title('Solution for h=0.1');
    legend('Actual Values','Calculated Values');
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
% function that calc values of a which is coef of delta_y(j-1)
function y = calcValueOf_a(yi,h,i)

    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
   if i ~= 1
       ym = yi(i-1);
   else
       ym = 0.5;
   end
   [n,~]=size(yi);
    if i ~= n
        yp = yi(i+1);
    else
        yp = 1;
    end

    yj=yi(i);
    y=1/h^2+2*yj/h;
end


%% -------------------------------------------------------------------
% function that calc values of b which is coef of delta_y(j)
function y = calcValueOf_b(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym=yi(i-1);
   else
       ym=0.5;
    end
   [n,~]=size(yi);
    if i~=n
        yp=yi(i+1);
    else
        yp=1;
    end

    yj=yi(i);
    y=-2/h^2-2*(yp-ym)/h;
end

%% -------------------------------------------------------------------
% function that calc values of c which is coef of delta_y(j+1)
function y = calcValueOf_c(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym=yi(i-1);
   else
       ym=0.5;
    end
   [n,~]=size(yi);
    if i~=n
        yp=yi(i+1);
    else
        yp=1;
    end

    yj=yi(i);
    y=1/h^2-2*yj/h;
end


%% -------------------------------------------------------------------
% function that calc values of d
function y = calcValueOf_d(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym=yi(i-1);
   else
       ym=0.5;
    end
   [n,~]=size(yi);
    if i~=n
        yp=yi(i+1);
    else
        yp=1;
    end

    yj=yi(i);
    y=-1*((ym-2*yj+yp)/h^2-2*yj*(yp-ym)/h);
end

%{
    Answer :- 
    x 0 (0) = 0.500000
    x 1 (1.000000e-01) = 0.513307
    x 2 (2.000000e-01) = 0.529658
    x 3 (3.000000e-01) = 0.549885
    x 4 (4.000000e-01) = 0.575110
    x 5 (5.000000e-01) = 0.606892
    x 6 (6.000000e-01) = 0.647455
    x 7 (7.000000e-01) = 0.700086
    x 8 (8.000000e-01) = 0.769856
    x 9 (9.000000e-01) = 0.865020
    x 10 (1) = 1.000000
%}