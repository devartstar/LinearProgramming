%{
LAGRANGE INTERPOLATION
-----------------------
    y(x) = y_1 * l_1(x) + y_2 * l_2(x) + .... + y_n+1 * l_n+1(x)
    y(x) = sum ( y_i * l_i(x) ) i = 1, 2, ..., n+1

    l_1(x) = (x - x_2)(x - x_3)...(x-x_n+1) / (x_1 - x_2)(x_1 - x_3)....(x_1 - x_n+1)
    l_2(x) = (x - x_1)(x - x_3)...(x-x_n+1) / (x_2 - x_1)(x_2 - x_3)....(x_2 - x_n+1)
    .
    .
    l_n+1(x) = (x - x_1)(x - x_2)...(x-x_n) / (x_n+1 - x_1)(x_n+1 - x_2)....(x_n+1 - x_n)

    l_i(X) = product((x - x_j) / (x_i - x_j))       i!=j
%}

x = [0 20 40 60 80 100];
y = [26.0 48.6 61.6 71.2 74.8 75.2];
n = length(x) - 1;

xp = 40;
sm = 0;
for i = 1 : n+1
    pr = 1;
    for j = 1 : n+1
        if j ~= i
            pr = pr * (xp - x(j)) / (x(i) - x(j));
        end
    end
    sm = sm + y(i) * pr;
end

yp = sm