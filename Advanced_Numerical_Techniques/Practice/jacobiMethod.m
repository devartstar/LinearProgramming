%{
    Most basic iterative method to solve linear systems
    a(1,1)x(1) + a(1,2)x(2) + .... + a(1,n)x(n) = b(1) ---> x(1) = ...
    a(1,1)x(1) + a(2,2)x(2) + .... + a(2,n)x(n) = b(2) ---> x(2) = ...
    .
    .
    a(n,1)x(1) + a(n,2)x(2) + .... + a(n,n)x(n) = b(n) ---> x(n) = ...

    General form :-
    ---------------
    x*(i) = - 1/a(i,i) ( sum(a(i,j) - b(i)) ) i = 1,2...n && i!=j
%}
%{
    Question to solve :-
    4x(1) +  x(2) + 2x(3) -  x(4) = 2
    3x(1) + 6x(2) -  x(3) + 2x(4) = -1
    2x(1) -  x(2) + 5x(3) - 3x(4) = 3
    4x(1) +  x(2) - 3x(3) - 8x(4) = 2
%}
a = [4 1 2 -1;
     3 6 -1 2;
     2 -1 5 -3;
     4 1 -3 -8];
b = [2;-1;3;2];

n = length(b);
x = zeros(n, 1);           % vector of n rows 1 col (all 0)
xnew = zeros(n, 1);

x(:) = 0;
iterlimit = 100;
tol = 0.000001;

for iterations = 1 : iterlimit  % iteration number
    convergence = true;             
    for i = 1 : n           % loop of equations
        s = 0;
        for j = 1 : n
            if j ~= i
                s = s + a(i,j)*x(j);
            end
        end
        xnew = -1/a(i,i) * (s - b(i));
        % even if one of the values > tolerance then convergence is false
        if abs(xnew(i) - x(i)) > tol        
            convergence = false;
        end
    end
    if convergence
        break;
    end
    x = xnew;
end

disp('iterations')
iter;
disp('solutions')
xnew;




















