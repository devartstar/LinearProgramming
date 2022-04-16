y = @(x) 2*x^2 - 5*x + 3;
x1 = input("Enter Value of x1 : ");
x2 = input("Enter Value of x2 : ");
if y(x1) * y(x2) > 0
    fprintf('No roots exist within the given interval\n');
    return;
end

if y(x1) == 0
    fprintf('x1 is one of the roots\n');
    return;
elseif y(x2) == 0
    fprintf('x2 is one of the roots\n');
    return;
end

for i = 1 : 100             % 100 bisections
    x_mid = (x1 + x2) / 2;
    if y(x1) * y(x_mid) < 0
        x2 = x_mid;
    else
        x1 = x_mid;
    end
    if abs(y(x1)) < 1.0E-6
        break;
    end
end
fprintf("The root : %f\nThe number of bisections: %d\n",x1,i)