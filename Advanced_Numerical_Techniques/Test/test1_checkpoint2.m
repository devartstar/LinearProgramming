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
        coeffAi = [1/h^2-(x-1)/(2*h) 0;0 1/h^2+1/h];
        coeffBi = [-2/h^2 -6;x -2/h^2];
        coeffCi = [1/h^2+(x-1)/(2*h) 0;0 1/h^2-1/h];
        coeffDi = [x^2;4*x+1];
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

    %solving using thomas algorithm and getting the values of z and y
    %zs stores all those values
    zs=blockTriDiagonal(A,B);
    ys=zs(1:2:2*(n));
    us=zs(2:2:2*n);
    %ys(n)=ys(n-1)+(zs(2*n-2)+zn(2))*h/2
    for i=1:n
        fprintf("value of y at %f is %f\n",x0+i*h,ys(i));
    end
    fprintf("\n\n");

    for i=1:n
        fprintf("value of z at %f is %f\n",x0+i*h,us(i));
    end
    
end




function y=blockTriDiagonal(A,B)
    [n,~]=size(A);
    %standard thomas algorithm but instead of elements we pick us matrices
    %to modify
    for i = 1:2:n
       if i==1
           const=A(i:i+1,i:i+1);
       else 
           const=A(i:i+1,i:i+1)-A(i:i+1,i-2:i-1)*A(i-2:i-1,i:i+1);
       end

       if i==1
           A(i:i+1,i+2:i+3)=inv(const)*A(i:i+1,i+2:i+3);
           B(i:i+1)=inv(const)*B(i:i+1);
       else
           B(i:i+1)=inv(const)*(B(i:i+1)-A(i:i+1,i-2:i-1)*B(i-2:i-1));
           A(i:i+1,i-2:i-1)=zeros(2,2);
           A(i:i+1,i:i+1)=eye(2,2);
           if i~=n-1
               A(i:i+1,i+2:i+3)=inv(const)*A(i:i+1,i+2:i+3);
           end
       end
       A(i:i+1,i:i+1)=eye(2,2);
    end
    z=zeros(n,1);


    for i=n:-2:2
        if i==n
            z(i-1:i,1)=B(i-1:i);
        else
            z(i-1:i,1)=B(i-1:i)-A(i-1:i,i+1:i+2)*z(i+1:i+2);
        end
    end
    y=z;
    
end

