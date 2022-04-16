%{
% Simple foor loop
for i = 1 : 10
 disp(i*5);
end
%}

% ---------------------------------------------
% Mean Square Error

E = randn(1000, 1);  %random data points
% disp(E)
%{
% ---> USING LOOPS
SSE = 0;
for i = 1 : 1000
    SSE = SSE + E(i) * E(i);
end
%}

% ----> USING DOT PRODUCT
SSE = dot(E, E);
MSE = SSE / 1000;
disp("Mean Square Error ");
disp(MSE);

%------------------------------------------------
% ---> Elements divisible by 2
X = 1 : 10;
Y = zeros(1,10);
for i = 1 : 10
    if mod(X(i),2) == 0
        Y(i) = 1;
    else
        Y(i) = 0;
    end
end
disp(Y);

%------------------------------------------------
% ---> Sum of elements divisible by 3
X = 1 : 10;
S = 0;
%{
for i = 1 : 10
    if mod(X(i), 3) == 0
        S = S + X(i);
    end
end
%}

for x = X
    if mod(x,3) == 0
        S = S + x;
    end
end
disp(S);

% --------------------------------------------
% ---> Find an element within an matrix
X = [2,5,1,3,9,7,2,5];
found = 0; i=0;
while ((~found) && i < length(X))
    i = i + 1;
    if X(i) == 8
        disp("Found It")
        found = 1;
    end
end
if found == 0
    disp("Not Found It");
end

%---------------------------------------------

