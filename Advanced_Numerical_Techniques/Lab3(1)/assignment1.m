%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
%}


%% Function to construct the Block Tridiagonal System
function assignment1()

    %------------------------------------------------
    % assigning the values given in question
    h = 0.1;
    x0 = 0; xn = 1;
    y0 = 0; yn = 0;
    ydf0 = 0;   ydfn = 0;
    %------------------------------------------------

    n = (xn - x0) / h;
    z0 = [y0; ydf0];
    zn = [yn; ydfn];

    %------------------------------------------------
    % I have stores matrix inside matrix in this form
    %{
        [[a1 b1]    [a2 b2]]    [a1 b1 a2 b2]
        [[c1 d1]    [c2 d2]]    [c1 d1 c2 d2]
        [                  ] -> [a3 b3 a4 b4]
        [[a3 b3]    [a4 b4]]    [c3 d4 c4 d4]
        [[c3 d3]    [c4 d4]]    
    %}
    % so for block tridiagonal system of n-1 * n-1
    % we need a size 2*(n-1)

    A=zeros(2*(n-1), 2*(n-1));
    B=zeros(2*(n-1), 1);

    x = x0 + h;

    % Constructing the Block Tridiagonal System
    % A(i)X(i-1) + B(i)X(i) + C(i)X(i+1) = D(i)

    for i = 1 : 2 : 2*(n-1)
        coefA_i = calcMatrixA(h);
        coefB_i = calcMatrixB(h);
        coefC_i = calcMatrixC(h);
        coefD_i = [
            0;
            16*(x^2+2);
        ];

        if i ~= 1
            A(i:i+1,i-2:i-1) = coefA_i;
        end
        A(i: i+1, i: i+1) = coefB_i;

        if i < 2*(n-1)-1
            A(i:i+1,i+2:i+3)=coefC_i;
        end
        if i == 1
            B(i:i+1,1)=coefD_i-coefA_i*z0;
        elseif i == 2*n-3
            B(i:i+1,1)=coefD_i-coefC_i*zn;
        else
            B(i:i+1,1)=coefD_i;
        end

        x=x+h;
    end

    zs = blockTDS(A, B);
    sol = zs(1:2:2*(n-1));
    xs = [x0+h:h:xn-h];
    
    % Displaying the output
    fprintf("Answer for 0.2:- \n");
    for i = 1:n-1
        fprintf("x %d (%d) = %f\n",i,x0 + (i)*h, sol(i));
    end
    fprintf("\n");

    %------------------------------------------------------------------
    % Actual Solution
    % solution equation generation taking 5 to 10 sec
    syms Y;
    Y=dsolve('D4Y-4*Y=16*(t^2+2)','Y(0)=0','Y(1)=0','D2Y(0)=0','D2Y(1)=0');
    fprintf("Solution equation:- \n");
    disp(Y);
    %------------------------------------------------------------------

    
    ezplot(Y,[0 1]);
    hold on
    plot(xs,sol,'-.*');
    legend('actual values','computed values');
    title('Numeric values ( h=0.02 )');
    hold off
end

%% Function for solving Block Tridiagonal system using Thomas Algorithm
function y = blockTDS(A, B)
    [n,~] = size(A);
    for i = 1:2:n
       if i==1
           const = A(i:i+1,i:i+1);
       else 
           const = A(i:i+1,i:i+1)-A(i:i+1,i-2:i-1)*A(i-2:i-1,i:i+1);
       end

       if i==1
           A(i:i+1,i+2:i+3) = inv(const)*A(i:i+1,i+2:i+3);
           B(i:i+1) = inv(const)*B(i:i+1);
       else
           B(i:i+1) = inv(const)*(B(i:i+1)-A(i:i+1,i-2:i-1)*B(i-2:i-1));
           A(i:i+1,i-2:i-1) = zeros(2,2);
           A(i:i+1,i:i+1) = eye(2,2);
           if i ~= n-1
               A(i:i+1,i+2:i+3) = inv(const)*A(i:i+1,i+2:i+3);
           end
       end
       A(i:i+1,i:i+1) = eye(2,2);
    end
    z=zeros(n,1);


    for i = n:-2:2
        if i == n
            z(i-1:i,1)=B(i-1:i);
        else
            z(i-1:i,1) = B(i-1:i)-A(i-1:i,i+1:i+2)*z(i+1:i+2);
        end
    end
    y=z;
end

%% Function to return the 2"*2 matrix
% A(i)X(i-1) + B(i)X(i) + C(i)X(i+1) = D(i)
function y = calcMatrixA(h)
    A_i = calcA();
    % B_i = calcB(x);
    C_i = calcC();
    % D_i = calcD(x);
    y = [
            1/(h^2), 0;
            -(C_i/(2*h)), ((1/h^2)-(A_i/(2*h)))
        ];
end

function y = calcMatrixB(h)
    % A_i = calcA(x);
    B_i = calcB();
    % C_i = calcC(x);
    D_i = calcD();
    y = [
            -2/(h^2), -1;
            D_i, (-2/(h^2))+B_i
        ];
end

function y = calcMatrixC(h)
    A_i = calcA();
    % B_i = calcB(x);
    C_i = calcC();
    % D_i = calcD(x);
    y = [
            1/(h^2), 0;
            (C_i/(2*h)), ((1/h^2)+(A_i/(2*h)))
        ];
end

%{
function y = calcMatrixD(h)
    % A_i = calcA(x);
    % B_i = calcB(x);
    % C_i = calcC(x);
    % D_i = calcD(x);
    y = [
            0;
            16*(x^2+2);
        ];
end
%}

function y = calcA()
    % A(x) given in question
    y = 0;
end

function y = calcB()
    % A(x) given in question
    y = 0;
end

function y = calcC()
    % A(x) given in question
    y = 0;
end

function y = calcD()
    % A(x) given in question
    y = -4;
end














