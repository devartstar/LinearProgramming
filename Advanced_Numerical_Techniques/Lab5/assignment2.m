%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment2
%}

%{
    LOGIC :-
    ----------
    I am attaching the pdf where i derived the formulae.
    Assignment2_working.pdf
    
%}

function assignment2()
    % the Boundary condition for BVP
    x0 = 1; xn = 2;
    y0 = 1; yn = 4;
    % step size
    h = 0.1;

    % The number of steps
    n = (xn - x0)/h - 1;
    %{
        % -1 since we have n-1 unknowns
         -----------------------------------------------------
        x0=1   x1=1.1   x2=1.2                  x(n-1)=1.9  x(n)=2.0
        The unkowns are for y1, y2, .... y(n-1)
    %}
    
    % initial approximations
    yi = zeros(n,1) + 3;

    %% Converting it to Matrix Form
    % A - The tri-diagonal matrix
    % B - The RHS matrix
    A = zeros(n,n);
    B = zeros(n,1);

    % k = iteration count more iteration cont to decrease error
    for k = 1:1000
        % now we will be forming the matrix form using 2nd order Finite
        % difference method and Newtons Linearization Method --> AX = B
        for i=1:n
            % Refer to Logic on top for deriviation of ai, bi, ci, Bi
            ai = calcValueOf_a(yi,h,i);
            bi = calcValueOf_b(yi,h,i);
            ci = calcValueOf_c(yi,h,i);

            B(i) = calcValueOf_d(yi,h,i);
            
            % Making the tri-diagonal Matrix and RHS
            if i~=1
                A(i,i-1)=ai;
            else
                B(i) = B(i) - y0*ai;
            end

            A(i,i) = bi;

            if i~=n
                A(i,i+1) = ci;
            else
                B(i) = B(i) - ci*yn;
            end
        end
        
        % Now we have the tridiagonal matrix system
        %{
            [b(1) c(2)  0   0 ...      ][y(1)]    [B(1)]
            [a(2) b(2) c(2) 0 ...      ][y(2)] =  [B(2)]
            ....    
            [0  0 0 ....      b(n) c(n)][y(n)]    [B(n)]
            solving this tri-diagonal system using thomas algorithm
            To get the delta values
        %}
        y = thomasAlgorithm(A,B);
        
        % updating the approximate solutions using the obtained del_y
        yi = y;
    end
    
    % Printing the Answer
    fprintf("Answer :- \n");
    fprintf("x %d (%d) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%d) = %f\n",i,x0 + i*h, yi(i));
    end
    fprintf("x %d (%d) = %f\n",n+1,xn, yn);
    fprintf("\n");

   
    plot([x0:h:xn],[y0;yi;yn],'*--');
    hold on;
    plot([x0:h:xn],[x0:h:xn].^2);
    title('Solution for h=0.1 taking 1000 iterations');
    legend('Calculated values','Actual Values');
    hold off;
end

%%-----------------------------------------------------------------
% function to solve tri-diagonal system using Thomas Algorithm
function y=thomasAlgorithm(A,B)
    [rows,~] = size(A);
    
    C = zeros(1,rows);
    D = zeros(1,rows);
    
    % constructing the first row
    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1,1);
    
    % from row 2 to last
    for i = 2 : rows
        if i ~= rows
            C(i) = A(i,i+1) / (A(i,i) - A(i,i-1)*C(i-1));
        end
        D(i) = (B(i) - A(i,i-1)*D(i-1)) / (A(i,i) - A(i,i-1)*C(i-1));
    end

    y = zeros(rows,1);
    y(rows) = D(rows);
    
    %{
        [1 C1  0   0....     ][del_y1]      [D1]
        [0  1  C2  0....     ][del_y2]      [D2]
        [0  0  1   C3..      ][del_y3]  =   [D3]
        [..                  ]      .
        [...................1][del_yn]      [Dn]
        
        del_yn = Dn
        del_y(n-1) = D(n-1) + C(n-1)*del_yn
        using back substitution... to get the values of del_y
    %}

    for i = rows-1:-1:1
        y(i) = D(i) - C(i)*y(i+1);
    end
end
%---------------------------------------------------------------------


%% -------------------------------------------------------------------
% function that calc values of a
function y=calcValueOf_a(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
   if i~=1
       ym = yi(i-1);
   else
       ym = 1; % using BC
   end

   [n,~]=size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 4; % using BC
    end


    yj = yi(i);
    y = ((ym-yp)/(2*h^2)) - (2*yj/h^2);
end
% ---------------------------------------------------------------------


%% -------------------------------------------------------------------
% function that calc values of b
function y=calcValueOf_b(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym = yi(i-1);
   else
       ym = 1; % using BC
   end
    
   [n,~]=size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 4; % using BC
    end


    yj = yi(i);
    y = (-2*ym/h^2) + (8*yj/h^2) - (2*yp/h^2);

end
% ---------------------------------------------------------------------

%% -------------------------------------------------------------------
% function that calc values of c
function y=calcValueOf_c(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym = yi(i-1);
   else
       ym = 1; % using BC
   end
    
   [n,~] = size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 4; % using BC
    end

    yj = yi(i);
    y= ((yp-ym)/(2*h*h)) - (2*yj/h^2);
end
% ---------------------------------------------------------------------

%% -------------------------------------------------------------------
% function that calc values of d
function y=calcValueOf_d(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym=yi(i-1);
   else
       ym=1; % using BC
   end
    
   [n,~]=size(yi);

    if i~=n
        yp=yi(i+1);
    else
        yp=4; % using BC
    end


    yj = yi(i);
    y= ((yp-ym)/(2*h))^2 - (2*yj*(ym-2*yj+yp)/h^2);
end
% ---------------------------------------------------------------------

%{
Answer :- 
x 0 (1) = 1.000000
x 1 (1.100000e+00) = 1.210000
x 2 (1.200000e+00) = 1.440000
x 3 (1.300000e+00) = 1.690000
x 4 (1.400000e+00) = 1.960000
x 5 (1.500000e+00) = 2.250000
x 6 (1.600000e+00) = 2.560000
x 7 (1.700000e+00) = 2.890000
x 8 (1.800000e+00) = 3.240000
x 9 (1.900000e+00) = 3.610000
x 10 (2) = 4.000000
%}