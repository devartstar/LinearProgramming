function CrankNicolsonScheme()
    % Interval Boundary
    x0 = 0; xn = 2;
    t0 = 0; tn = 1;

    % spacial step size
    xh = 2/100;
    % time step size
    th = 1/2500;
    
    % Diffusion constant
    c = 0.21;

    % calculating no. of steps
    n = (xn - x0) / xh - 1;
    p = th * c / xh^2;

    previous_u = ft0(x0,xn,xh);
    mesh = [previous_u];
    ts = [4/25 6/25 11/25 19/25];
    xs = [];

    for t = t0+th : th : tn-th
        A = zeros(n,n);
        RHS = zeros(n,1);
        i = 1;
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
        display(mesh);
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
% Boundary Condition
function x0 = f0(x,t)
    x0 = 0;
end

function xn = fn(x,t)
    xn = 0;
end
% -------------------------------------------------------------------

function y = ft0(x0, xn, xh)
    xs = [x0 : xh : xn];
    y = 4*xs.^2 - xs.^4;
end