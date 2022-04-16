%{
    Name : Devjit Choudhury
    Roll : 19MA20014
%}

function assignment2()

    %------------------------------------------------
    % assigning the values given in question
    h = 0.25;
    x0 = 0;
    xn = 2;
    n = (xn - x0)/h;
    %------------------------------------------------

    % Using finite difference to 
    % initially initializing the coffecient matrix
    A = zeros(n-1,n-1); 
    % initially initializing the rhs matrix
    B = zeros(n-1,1);

    x = x0 + h;
    %{
            We will have fictitoue values y-1 and yn+1
            Which is replaced by boundary condition
        The coffecient matrix will be of the form
        [Y1 Y2 Y3 0 0 ....
        [Y1 Y2 Y3 Y4 ......
        [Y1 Y2 Y3 Y4 Y5 ....
        [0 Y1 Y2 Y3 Y4 Y5 ....
                so on we have pentadiagonal system
        [
        [
        [
        [                   0   Yn-3 Yn-2 Yn-1]
    %}
    %-------------------------------------------------------------------------
    % Using the boundary condition to RHS for i = 1, 2
    coef_A = (1/h^4) + (1/h^3);
    coef_B = (-4/h^4) - (2/h^3);
    coef_C = (6/h^4) + 1;
    coef_D = (-4/h^4) + (2/h^3);
    coef_E = (1/h^4) - (1/h^3);
    
    x=x0+h;
    A(1,1) = coef_A+coef_C;
    A(1,2) = coef_D;
    A(1,3) = coef_E;
    B(1) = 2-x^4;

    x=x+h;
    A(2,1) = coef_B;
    A(2,2) = coef_C;
    A(2,3) = coef_D;
    A(2,4) = coef_E;
    B(2) = 2-x^4;

    x=x+h;
    for i=3:n-3
        A(i,i-2) = coef_A;
        A(i,i-1) = coef_B;
        A(i,i) = coef_C;
        A(i,i+1) = coef_D;
        A(i,i+2) = coef_E;
        B(i) = 2-x^4;
        x=x+h;
    end

    % Using the boundary condition to RHS for i = n-1, n-2
    A(n-2,n-4) = coef_A;
    A(n-2,n-3) = coef_B;
    A(n-2,n-2) = coef_C;
    A(n-2,n-1) = coef_D;
    B(n-2) =  2-x^4 -coef_E;

    x=x+h;
    A(n-1,n-3) = coef_A;
    A(n-1,n-2) = coef_B;
    A(n-1,n-1) = coef_C+coef_E;
    B(n-1) =  2-x^4 -coef_D-2*h*coef_E;
    
    %-------------------------------------------------------------------------
    
    ys = guassElimination(A,B);
    xs = [x0:h:xn];
    ys = [0;ys];
    ys = [ys;1];
    disp('The required solution is:');
    ys
   
    %{
    syms Y;
    Y=dsolve('D4Y-2*D3Y+Y=2-t^4','Y(0)=0','Y(2)=1','DY(0)=0','DY(2)=1');
    ezplot(Y,[0 2]);
    hold on;
    plot(xs,ys,'-.*');
    legend('actual values','computed values');
    title('Numeric values h=0.02');
    hold off;
    %}
end



%% Function to solve a matrix equation using guall elimination
function y = guassElimination(A,B)
P=[A B];
[row,col] = size(P);
for n = 1:row-1
    a=P(n,n);  
     P(n,:)= P(n,:)/a;
    for i=n+1:row     
      P(i,:)= P(i,:)- P(i,n)* P(n,:);
    end
end
 a=P(row,row);  
 P(row,:)= P(row,:)/a;
 s=0;
for m=row:-1:2 
    for k=m+1:col
        s=s+P(m-1,k-1)* P(k-1,col);
        P(m-1,col)= P(m-1,col)- s; 
        s=0;
    end
end
  y=P(:,col);
end

