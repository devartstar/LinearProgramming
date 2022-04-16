%{
SIMPSON 1/3 RD RULE
---------------------
    This method calculates the integral by summing every two division
    at a time, so 3 values of x are taken into account
    to cover the whole domain exactly, the number of strips must be even
    Formulae for the first 2 stips
    A = 1/3*h*(f(x(0)) + 4*f(x(1)) + f(x(2)))
    A(i) = 1/3*h*(f(x(i)) + 4*f(x(i+1)) + f(x(i+2)))
    then sum (A)    i = 0, 1, 2, ... n-2
    sum(A) = 1/3*h*{
                [f(x(0)) + f(x(n))] +
                4*[f(x(1)) + f(x(3)) + ... + f(x(n-1))] +
                2*[f(x(2)) + f(x(4)) + ... + f(x(n-2))]
                }
%}


f = @(x) x * sin(x);
a = 0;
b = pi / 2;
n = 18;
h = (b - a) / n;
s = f(a) + f(b);
for i = 1 : 2 : n-1
    s = s + 4 * f(a + i * h);
end
for i = 2 : 2 : n-2
    s = s + 2 * f(a + i * h);
end

I = h/3 * s