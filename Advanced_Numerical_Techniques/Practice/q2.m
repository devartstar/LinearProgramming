%{
    GUASS ELIMINATION
    Ax = B
%}

A = [1 1 1; 2 1 3; 3 4 -2];
b = [4; 7; 9];

% Get Argumented Matrix
Ab = [A, b];
n = 3;

%% Guass Elimination
% Row_j = Row_j - alpha(i,j)Row_i
% alpha(i,j) = A(j,i) / A(i,i)

% A(1,1) as pivot element
alpha = Ab(2,1) / Ab(1,1);
Ab(2, :) = Ab(2, :) - alpha*Ab(1, :); % R2 = R2 - alpha*R1

alpha = Ab(3,1) / Ab(1,1);
Ab(3, :) = Ab(3, :) - alpha*Ab(1, :);

% A(2,2) as pivot element
alpha = Ab(3,2) / Ab(2, 2);
Ab(3, :) = Ab(3, :) - alpha*Ab(2, :);

% A(3,3) as pivot element


%% Back Substitution
x = zeros(3, 1);
x(3) = Ab(3,end) / Ab(3,3);
x(2) = (Ab(2,end) - Ab(2,2+1:n)*x(2+1:n)) / Ab(2,2);
% x(1) = (Ab(1,end) - (Ab(1,3)*x(3) + Ab(1,2)*x(2))) / Ab(2,2);
x(1) = (Ab(1,end) - Ab(1,1+1:n)*x(1+1:n)) / Ab(1,1);

x(3) = Ab()
for i = 2 : -1 : 1
    x(i) = (Ab(i,end) - Ab(i,i+1:n)*x(i+1:n)) / Ab(i,i);
end