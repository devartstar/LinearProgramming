%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment1
%}

%{
    Logic :-
    ---------
    I am attaching the pdf where i derived the formulae.
    Assignment1_working.pdf
%}

function assignment1()
    % the Boundary condition for BVP
    x0 = 0; xn = 1;
    y0 = 0; yn = 2;
    % step size
    h = 0.1;

    % The number of steps
    n = (xn - x0)/h - 1;
    %{
        % -1 since we have n-1 unknowns
         -----------------------------------------------------
        x0=0   x1=0.1   x2=0.2                  x(n-1)=0.9  x(n)=1.0
        The unkowns are for y1, y2, .... y(n-1)
    %}
    
    % initial approximations
    yi = zeros(n,1) + 1;

    %% Converting it to Matrix Form
    % A - The tri-diagonal matrix
    % B - The RHS matrix
    A = zeros(n,n);
    B = zeros(n,1);

    % k = iteration count more iteration cont to decrease error
    for k = 1:1000
        % now we will be forming the matrix form using 2nd order Finite
        % difference method and Newtons Linearization Method --> AX = B
        for i = 1:n
            % Refer to Logic on top for deriviation of ai, bi, ci, Bi
            ai = calcValueOf_a(yi,h,i);
            bi = calcValueOf_b(yi,h,i);
            ci = calcValueOf_c(yi,h,i);
            
            % Making the tri-diagonal Matrix
            % for first row there wont be ai
            if i~=1
                A(i,i-1)=ai;
            end
            A(i,i)=bi;
            % for last row(n) there wont be ci
            if i~=n
                A(i,i+1)=ci;
            end

            % the RHS of matrix eq
            B(i) = calcValueOf_d(yi,h,i);
        end

        % Now we have the tridiagonal matrix system
        %{
            [b(1) c(2)  0   0 ...      ][del_y(1)]    [B(1)]
            [a(2) b(2) c(2) 0 ...      ][del_y(2)] =  [B(2)]
            ....    
            [0  0 0 ....      b(n) c(n)][del_y(n)]    [B(n)]
            solving this tri-diagonal system using thomas algorithm
            To get the delta values
        %}
        del_y = thomasAlgorithm(A,B);

        % updating the approximate solutions using the obtained del_y
        yi = yi + del_y;
    end
    
    % Printing the Answer
    fprintf("Answer :- \n");
    fprintf("x %d (%d) = %f\n",0,x0, y0);
    for i = 1:n
        fprintf("x %d (%d) = %f\n",i,x0 + i*h, yi(i));
    end
    fprintf("x %d (%d) = %f\n",n+1,xn, yn);
    fprintf("\n");

    plot([x0:h:xn],[y0;yi;yn],'.-');
    hold on;
    plot([x0:h:xn],2*[x0:h:xn]);
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
function y = calcValueOf_a(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
   if i~=1
       ym = yi(i-1);
   else
       ym = 0; % using BC
   end
    
   [n,~] = size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 2; % using BC
    end

    yj = yi(i);
    y = yj/h^2 - 1/(2*h);
end
% ---------------------------------------------------------------------


%% ---------------------------------------------------------------------
% function that calculates value of b
function y = calcValueOf_b(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym = yi(i-1);
   else
       ym = 0; % using BC
   end
    
   [n,~] = size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 2 ; % using BC
    end


    yj = yi(i);
    y = (yp - 4*yj + ym)/h^2;
end
% ---------------------------------------------------------------------

%% ---------------------------------------------------------------------
function y = calcValueOf_c(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym = yi(i-1);
   else
       ym = 0; % using BC
   end
    
   [n,~]=size(yi);

    if i~=n
        yp = yi(i+1);
    else
        yp = 2; % using BC
    end


    yj = yi(i);
    y = yj/h^2 + 1/(2*h);
end
% ---------------------------------------------------------------------

%% ---------------------------------------------------------------------
function y=calcValueOf_d(yi,h,i)
    % we may need value of y(j-1) & y(j+1) 
    % whode vaued depend ij j-1 and j+1 touch boundary
    if i~=1
       ym=yi(i-1);
   else
       ym=0; % using BV
   end
    
   [n,~]=size(yi);

    if i~=n
        yp=yi(i+1);
    else
        yp=2; % using BV
    end

    yj=yi(i);
    y= 2 - (yj*(ym - 2*yj + yp)/h^2) - ((yp - ym)/(2*h));
end
% ---------------------------------------------------------------------

%{
Answer :- 
x 0 (0) = 0.000000
x 1 (1.000000e-01) = 0.200000
x 2 (2.000000e-01) = 0.400000
x 3 (3.000000e-01) = 0.600000
x 4 (4.000000e-01) = 0.800000
x 5 (5.000000e-01) = 1.000000
x 6 (6.000000e-01) = 1.200000
x 7 (7.000000e-01) = 1.400000
x 8 (8.000000e-01) = 1.600000
x 9 (9.000000e-01) = 1.800000
x 10 (1) = 2.000000
%}