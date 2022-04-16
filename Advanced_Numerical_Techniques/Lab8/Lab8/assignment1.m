%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment1
%}

function assignment1()
    %% Data provided in the question :-
    % Interval boundaries
    x0 = 0; xn = 2; 
    t0 = 0; tn = 1;

    % space step size
    xh = 2/100;

    % time step size
    th = 1/2500;

    % U_t = C * U_xx
    C = 0.21;
    
    % calculating the number of steps
    n = (xn - x0)/xh - 1;

    p = th * C/(2*xh^2);
    
    % previous_u will store the value of u calculated in last step
    previous_u = ft0(x0,xn,xh);
    mesh = [previous_u];
    ts = [4/25 6/25 11/25 19/25];
    xs = [];

    for t = t0+th : th : tn
        %% Converting it to Matrix Form
        % A - The tri-diagonal matrix
        % RHS - right hand side of tridiagonal system
        A = zeros(n,n);
        RHS = zeros(n,1);
        i=1;

        % Let us construct the Tridiagonal System
        %{
            Using the discretized equation
            (uj_n+1 - uj_n) / th = C * (uj+1_n - 2*uj_n + uj-1_n) / xh^2
        %}
        %% Writing it in Matrix form -
        for x = x0+xh : xh : xn-xh
            
            if i~=1
                A(i,i-1) = -p;
            else 
                RHS(i) = RHS(i) + p*f0(x0,t);
            end

            A(i,i) = 1 + 2*p;
        
            if i~=n
                A(i,i+1) = -p;
            else
                RHS(i) = RHS(i) + p*fn(xn,t);
            end

            RHS(i) = previous_u(i+1) + p*(previous_u(i+2) - 2*previous_u(i+1) + previous_u(i));
            i = i+1;
        end
        y = thomasAlgorithm(A, RHS);
        previous_u = [f0(x0,t) y' fn(xn,t)];
        mesh = [mesh;previous_u];
    end
    xs = [mesh(1+4/(25*th),:);mesh(1+6/(25*th),:); mesh(1+11/(25*th),:); mesh(1+19/(25*th),:)];

    %% plotting the graph for given time points
    subplot(2,1,1);
    for i=1:4
        plot([x0:xh:xn],xs(i,:));
        hold on
    end
    title('x values for given time stamps')
    legend('t = 4/25','t = 6/25','t = 11/25','t = 19/25');

    %% plotting the graph for temperature distribution with change in x and t
    subplot(2,1,2);
    [X,Y] = meshgrid(0:xh:2,0:th:1);
    surf(X,Y,mesh);
    colormap("default")
    shading interp
    title('Temperature distribution with change in x and t')
end

%% -----------------------------------------------------------------
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
        [1 C1  0   0....     ][M1]      [D1]
        [0  1  C2  0....     ][M2]      [D2]
        [0  0  1   C3..      ][M3]  =   [D3]
        [..                  ]      .
        [...................1][Mn]      [Dn]
        this n corresponds to n-1 (but we already subtracted 1 to n
        Mn = Dn
        My(n-1) = D(n-1) + C(n-1)*Mn
        using back substitution... to get the values of Mi
    %}

    for i = rows-1:-1:1
        y(i) = D(i) - C(i)*y(i+1);
    end
end
%---------------------------------------------------------------------

%% -----------------------------------------------------------------
% Boundary Condition
function x0 = f0(x,t)
    x0 = 0;
end

function xn = fn(x,t)
    xn = 0;
end
% -------------------------------------------------------------------

%% -----------------------------------------------------------------
% Intermediate Condition
function t0 = ft0(x0,xn,xh)
    xs = [x0:xh:xn];
    t0 = 4*xs.^2-xs.^4;
end
% -------------------------------------------------------------------
