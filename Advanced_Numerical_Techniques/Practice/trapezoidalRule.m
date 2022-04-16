%{
TRAPEZOIDAL RULE
------------------
    Numerical integration method
    A = h * [ f(x(i) + f(x(i + 1)) ] / 2
    
    I = sum(A) i = 0, 1, 2, .. n-1
    I = h * [ 1/2 * (f(x(0) + f((n))) + sum(f(x(i))) ] i =1, 2, 3, ... n-1
%}

f = @(x) x * sin(x);
a = 0;
b = pi / 2;
n = 100;
h = (b - a) / n;
s = 0.5 * (f(a) + f(b));
for i = 1 : n-1
    s = s + f(a + i*h);
end

I = h * s