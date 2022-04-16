%{
    Name : Devjit Choudhury
    Roll : 19MA20014
%}
%{
    LOGIC :-
    ------------------
    We are given a coupled system of Differnetial Equation
    Eq1 - y'' + (x-1)y' - 6z = x^2
    Eq2 - z'' - 2z' + xy = 4x + 1
    Boundary Conditions :-
    y(0) = y'(2) = 0
    z(0) = z'(2) = 0
    -------------------
    We discretize both the Equations :-
    ---------------------------
    Discretizing Eq 1
    (1/h^2 - (x-1)/2h)y(i-1) + (-2/h^2)y(i) + (1/h^2 + (x-1)/2h)y(i+1) -6z(i) =
    x(i)^2
    --------------------------
    Discretizing Eq 2
    x(i)y(i) + (1/h^2 + 1/h)z(i-1) + (-2/h^2)z(i) + (1/h^2 - 1/h)z(i+1) =
    4x(i) + 1
    ---------------------------
    Now we will construst the block tridiagonal
    Assume Z(i) =  [y(i)]
                   [z(i)]
    Now we will construct the block tridiagonal system 
    and solve it using Thomas Algorithm
    
    %%----------------------------------------------
    Will be attaching a pdf/img of the detailed working
%}

function test1()
    % the Boundary condition for BVP
    x0 = 0; xn = 2;
    y0 = 0; u0 = 0;
    ypn = 0; upn = 0;
    % step size
    h=0.2;
    
    % The number of steps
    n=(xn-x0)/h;

    z0=[y0;u0];

    A=zeros(2*(n),2*(n));
    B=zeros(2*(n),1);

    x=x0+h;

    for i=1:2:2*(n)
        % filling in the block tridiagonal system
        coeffAi = calcValueOf_a(x,h);
        coeffBi = calcValueOf_b(x, h);
        coeffCi = calcValueOf_c(x, h);
        coeffDi = calcValueOf_RHS(x, h);
        if i~=1
            A(i:i+1,i-2:i-1)=coeffAi;
        end
        A(i:i+1,i:i+1)=coeffBi;
        if i<2*(n)-1
            A(i:i+1,i+2:i+3)=coeffCi;
        end
        if i==1
            B(i:i+1,1)=coeffDi-coeffAi*z0;
        elseif i==2*n-1
            A(i:i+1,i-2:i-1)=coeffAi+coeffCi;
            B(i:i+1,1)=coeffDi;
        else
            B(i:i+1,1)=coeffDi;
        end
        x=x+h;
    end

    % calling the function to solve Block Tridiagonal System using Thomas
    % Algorithm
    zs = thomasBlockTriDiagonal(A,B);
    ys = zs(1:2:2*(n));
    us = zs(2:2:2*n);

    for i=1:n
        fprintf("value of y %d (%f) = %f\n",i,x0+i*h,ys(i));
    end
    fprintf("\n\n");

    for i=1:n
        fprintf("value of z %d (%f) = %f\n",i,x0+i*h,us(i));
    end
    
end




%% Instead of element swe have 2*2 matrices in tridiagonal system
function y = thomasBlockTriDiagonal(A,B)
    [n,~] = size(A);

    for i = 1:2:n
       if i == 1
           const = A(i:i+1,i:i+1);
       else
           const = A(i:i+1,i:i+1) - A(i:i+1,i-2:i-1)*A(i-2:i-1,i:i+1);
       end

       if i == 1
           A(i:i+1,i+2:i+3) = inv(const)*A(i:i+1,i+2:i+3);
           B(i:i+1) = inv(const)*B(i:i+1);
       else
           B(i:i+1) = inv(const)*(B(i:i+1)-A(i:i+1,i-2:i-1)*B(i-2:i-1));
           A(i:i+1,i-2:i-1) = zeros(2,2);
           A(i:i+1,i:i+1) = eye(2,2);
           if i~=n-1
               A(i:i+1,i+2:i+3) = inv(const)*A(i:i+1,i+2:i+3);
           end
       end
       A(i:i+1,i:i+1) = eye(2,2);
    end
    z = zeros(n,1);

    % Back substitution to get the values
    for i = n:-2:2
        if i == n
            z(i-1:i,1) = B(i-1:i);
        else
            z(i-1:i,1) = B(i-1:i) - A(i-1:i,i+1:i+2)*z(i+1:i+2);
        end
    end
    y = z;
   
end

%% functions to calculate the entries of the Tridiagonal Matrix
% --------------------------------------------------------------
function y = calcValueOf_a(x, h)
    y = [1/h^2-(x-1)/(2*h) 0;0 1/h^2+1/h];
end
% --------------------------------------------------------------
function y = calcValueOf_b(x, h)
    y = [-2/h^2 -6;x -2/h^2];
end
% --------------------------------------------------------------
function y = calcValueOf_c(x, h)
    y = [1/h^2+(x-1)/(2*h) 0;0 1/h^2-1/h];
end
% --------------------------------------------------------------
function y = calcValueOf_RHS(x, h)
    y = [x^2;4*x+1];
end
% --------------------------------------------------------------

%{
value of y 1 (0.200000) = 0.659073
value of y 2 (0.400000) = 1.390345
value of y 3 (0.600000) = 2.135020
value of y 4 (0.800000) = 2.833403
value of y 5 (1.000000) = 3.435107
value of y 6 (1.200000) = 3.906936
value of y 7 (1.400000) = 4.236838
value of y 8 (1.600000) = 4.433501
value of y 9 (1.800000) = 4.522257
value of y 10 (2.000000) = 4.538777


value of z 1 (0.200000) = -0.169287
value of z 2 (0.400000) = -0.339808
value of z 3 (0.600000) = -0.493396
value of z 4 (0.800000) = -0.617830
value of z 5 (1.000000) = -0.707816
value of z 6 (1.200000) = -0.764551
value of z 7 (1.400000) = -0.794069
value of z 8 (1.600000) = -0.804925
value of z 9 (1.800000) = -0.805889
value of z 10 (2.000000) = -0.804337
%}
