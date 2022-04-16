%{
EULER METHOD
----------------
    To solve ODE (IVP)
    y(x+h) = y(x) + y(x)'*h + y(x)''*h^(2)/2! + .....
    y(x+h) = y(x) + y(x)'*h
%}
dy = @(x, y) x*y;
f = @(x) exp(x^2 / 2);

x0 = 0;
xn = 2;
y = 1;
h = 0.001;

fprintf("x \t\t y (Euler)\t y (Analytical)\n");
fprintf('%f \t %f\t %f\n', x0, y, f(x0));

for x = x0 : h : xn - h     % every step we calcl y for next x (or x+h)
    y = y + dy(x, y) * h;
    fprintf('%f \t %f\t %f\n', x + h, y, f(x + h));
end












