%Solve the following linear boundary value problems by using spline
%interpolation method with the help of Matlab code

% y'' + A(x)y' + B(x)y = C(x)

% y(x0)=alpha

% y(xn) =beta

% Step Size: h

% Unknowns: y1,y2....yn-1


% Defining the parameters
% Definition of the Functions A,B,C

A=@(x)(5);
B=@(x)(6);
C=@(x)(0);

h=0.2;
x0=0;                         % Start Point 
xn=1;                         % End Point
n=(xn-x0)/h;                  % Value of n based on above definition 
x_with_endpoints=x0:h:xn;
x=x0+h:h:xn-h;                % Vector x

alpha=1;                      % Function Values at endpoints
beta=1/exp(3);


% Defining the matrix of the system of Equations in Spline Interpolation
% The Unknowns M0,M1,...Mn,y1,y2....yn-1  (2n unknowns)
% Let L be the LHS matrix. R be RHS matrix.

n_vars=2*n;
L=zeros(n_vars,n_vars);
R=zeros(n_vars,1);

% Case for the 1st row (k=0)

L(1,1)=1/A(x_with_endpoints(1)) - (h/3);
L(1,2)=h/6;
L(1,n+2)=1/h;
R(1,1)=C(x_with_endpoints(1))/A(x_with_endpoints(1)) - alpha*(B(x_with_endpoints(1))/(A(x_with_endpoints(1))) - 1/h);


% Case for last row (k=n)

L(n_vars,n)=h/6;
L(n_vars,n+1)=(h/3 + 1/A(x_with_endpoints(n+1)));
L(n_vars,n_vars)=-1/h;
R(n_vars,1)= C(x_with_endpoints(n+1))/A(x_with_endpoints(n+1)) - beta*(1/h + B(x_with_endpoints(n+1))/A(x_with_endpoints(n+1)));


%disp(L)
%disp(R)

% The Eqns for k=1 

 L(2,2)=1/A(x_with_endpoints(2)) - (h/3);
 L(2,3)=h/6;
 L(2,n+2)=B(x_with_endpoints(2))/A(x_with_endpoints(2)) - 1/h;
 L(2,n+3)=1/h;
 R(2,1)=C(x_with_endpoints(2))/A(x_with_endpoints(2));
 
 L(3,2)=h/3 + 1/A(x_with_endpoints(2));
 L(3,1)=h/6;
 L(3,n+2)=1/h + B(x_with_endpoints(2))/A(x_with_endpoints(2));
 R(3,1)=C(x_with_endpoints(2))/A(x_with_endpoints(2)) - (alpha)*(-1/h);
 
% Case for k=2...n-2  (Two Equations for each k)

for k=2:n-2
    
    % 1st eqn 
    L(2*k,k+1)=1/A(x_with_endpoints(k+1)) - (h/3);
    L(2*k,k+2)=h/6;
    L(2*k,n+k+1)=B(x_with_endpoints(k+1))/A(x_with_endpoints(k+1)) - 1/h;
    L(2*k,n+k+2)=1/h;
    R(2*k,1)=C(x_with_endpoints(k+1))/A(x_with_endpoints(k+1));
    
    % 2nd Eqn
    L(2*k+1,k+1)=h/3 + 1/A(x_with_endpoints(k+1));
    L(2*k+1,k)=h/6;
    L(2*k+1,n+k+1)=1/h + B(x_with_endpoints(k+1))/A(x_with_endpoints(k+1));
    L(2*k+1,n+k)= (-1/h);
    R(2*k+1,1)=C(x_with_endpoints(k+1))/A(x_with_endpoints(k+1));
    
end

% Question>> Write the equations for k=n-1
 
 k=n-1;
 L(k-1,k-1)=1/A(x_with_endpoints(k-1)) - (h/3);
 L(k-1,k-2)=h/6;
 L(k-1,n+2)=B(x_with_endpoints(k-1))/A(x_with_endpoints(k-1)) - 1/h;
 L(k-1,n+3)=1/h;
 R(k-1,k)=C(x_with_endpoints(k-1))/A(x_with_endpoints(k-1));
 
 L(k-2,k-1)=h/3 + 1/A(x_with_endpoints(k-1));
 L(k-2,k)=h/6;
 L(k-2,n+2)=1/h + B(x_with_endpoints(k-1))/A(x_with_endpoints(k-1));
 R(k-2,k)=C(x_with_endpoints(k-1))/A(x_with_endpoints(k-1)) - (alpha)*(-1/h);
 
 
 
 
 
 
 
 

% Display the Resultant LHS and RHS matrix and the Solution  
disp("L=")
disp(L);
disp("R=")
disp(R);

y=L\R;

disp("The numerical Solution obtained:");
disp(y)

