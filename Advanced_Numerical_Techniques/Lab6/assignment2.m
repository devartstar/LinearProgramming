%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment2
%}

%{
    LOGIC :-
    ----------
    S(x) = {p(x), x(k) < x < x(k+1)}

        Eq 1 -->
        P(k) = (M(k)/6)[((x(k+1) - x)^3)/h - h(x(k+1)-x)]
                + (M(k+1)/6)[((x - x(k))^3)/h - h(x-x(k))]
                + y(k)/h (x(k+1) - x) 
        M(k) = S''(x(k))
        
        Eq 2 -->
        M(k-1) + 4M(k) + M(k+1) = 6/h^2 (y(k+1) - 2y(k) + y(k+1))
        k = 1,2,3....,n-1
    
%}

function assignment2()
    % the Boundary condition for BVP
    x0 = 0; xn = 1;
    y0 = 1; yn = 0;
    % step size
    h = 0.2;
   
    % The number of steps
    n = (xn - x0)/h;

    % Let us construct the Block Tridiagonal System
    %{
    
        We get Block Tridiagonal System
-----------------------------------------------------------------------------------------
[0    0  ][M(k-1)]  + [ 1/6-h/3  -1/h][M(k)] + [-h/6  1/h][M(k+1)] = [-5/6]
[h/6 -1/h][y(k-1)]    [ 1/6+h/3   1/h][y(k)]   [ 0    0  ][y(k+1)]   [-5/6]
-----------------------------------------------------------------------------------------
       
    %}

    % using Boundary condition
    z0 = [0;y0];
    zn = [0;yn];

    %% Converting it to Matrix Form
    % A - The Block tri-diagonal matrix
    % B - The RHS matrix
    % A Z = B
    A = zeros(2*(n-1), 2*(n-1));
    B = zeros(2*(n-1),1);

    x = x0 + h;
    for i=1:2:2*(n-1)
        
        % check out the comment above
        Ai = calcValueOf_a(h);
        Bi = calcValueOf_b(h);
        Ci = calcValueOf_c(h);
        Di = calcValueOf_RHS();

        if i~=1
            A(i:i+1,i-2:i-1) = Ai;
        end
        
        A(i:i+1,i:i+1) = Bi;

        if i < (2*(n-1)-1)
            A(i:i+1,i+2:i+3) = Ci;
        end

        if i==1
            B(i:i+1,1) = Di - Ai*z0;
        elseif i==2*n-3
            B(i:i+1,1) = Di - Ci*zn;
        else
            B(i:i+1,1) = Di;
        end

        x = x + h;
    end
    
    fprintf("The block tridiagonal system\n");
    disp(A);

    %storing the answer
    zs = thomasBlockTriDiagonal(A,B);
    ys = [y0;zs(2:2:2*(n-1));yn];
    xs = [x0:h:xn];

    % displaying the answers
    for i = 1 : n+1
        fprintf("x %d = %f || y %d = %f\n",i,xs(i),i,ys(i));
    end
    
    syms Y;
    Y=dsolve('D2Y+6*DY+5=0','Y(0)=1','Y(1)=0');
    ezplot(Y,[0 1]);
    hold on
    plot(xs,ys,'-.*');
    legend('actual values','computed values');
    title('Solution for h=0.2');
    hold off
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
function y = calcValueOf_a(h)
    y = [0 0; h/6 -1/h];
end
% --------------------------------------------------------------
function y = calcValueOf_b(h)
    y = [1/6-h/3 -1/h;1/6+h/3 1/h];
end
% --------------------------------------------------------------
function y = calcValueOf_c(h)
    y = [-h/6 1/h;0 0];
end
% --------------------------------------------------------------
function y = calcValueOf_RHS()
    y = [-5/6;-5/6];
end
% --------------------------------------------------------------

%{
The block tridiagonal system
    0.1000   -5.0000   -0.0333    5.0000         0         0         0         0
    0.2333    5.0000         0         0         0         0         0         0
         0         0    0.1000   -5.0000   -0.0333    5.0000         0         0
    0.0333   -5.0000    0.2333    5.0000         0         0         0         0
         0         0         0         0    0.1000   -5.0000   -0.0333    5.0000
         0         0    0.0333   -5.0000    0.2333    5.0000         0         0
         0         0         0         0         0         0    0.1000   -5.0000
         0         0         0         0    0.0333   -5.0000    0.2333    5.0000

x 1 = 0.000000 || y 1 = 1.000000
x 2 = 0.200000 || y 2 = 0.723851
x 3 = 0.400000 || y 3 = 0.514174
x 4 = 0.600000 || y 4 = 0.336755
x 5 = 0.800000 || y 5 = 0.167400
x 6 = 1.000000 || y 6 = 0.000000
%}

