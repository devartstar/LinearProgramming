%Solving non-linear BVP using second order finite difference method and Newton - Raphson method

%BVP :  y'' - y - 2y^3 = x^2
%        y(0) = 0      y(1) = 0
%       step size h = 0.5

h = 0.5;
x0 = 0;
x1 = 1;
n = (x1 - x0)/h ;

x_val = x0 + (0:n)*h ;

y0 = 0;                   %boundary value y at x=0
y2 = 0;                   %boundary value y at x=1

syms y1 z
eqn1 = -2*y1 - 1*(h^2)*y1 + 2*(h^2)*(y1^3) + (h^2)*(x_val(2)^2) == 0;   %Discretisation using second order FDM gives this
eqn = -2*z - 1*(h^2)*z + 2*(h^2)*(z^3) + (h^2)*(x_val(2)^2);

y10 = 0 ;   %Intial guess for Newton Raphson Method

%Q1. Write a command here to call the NRM function and store its value in the variable y1_val
y1 = NRM(y10, eqn);
fprintf('\n The value of y at %f is : %f \n',x_val(2),y1_val);


function sol = NRM(y_initial, func)         %Function for Newton Raphson Method

syms z;
f = func;        
g = diff(f,z);     
epsilon = 10^-(10);      %tolerance

for i=1:100
     %Q2. Write a three line code here for Newton Raphson Method.
     for j = 1 : n
        y1 = y_initial - f(y_initial) / g(y_initial);
     end
     
     %Q3. Calculate error here and save it to err variable.
     err = mod(y_initial - y1);
     
     if err < epsilon 
        break
     end
     y_initial = y1(i);
end
     %Q4. Write a print statement to display the number of iterations executed. 
     fprintf("No of Iterations : %d", i);
     
sol = y1(i);
end


