%{
    Name : Devjit Choudhury
    Roll No. : 19MA20014
    Assignment1
%}

function assignment1 ()
    %% Data provided in the question :-

    % Spatial Boundary
    x0 = 0; xN = pi; 
    % Spatial Step Size
    xh = pi/128;

    % Time Boundaries
    t0 = 0; tn = 1;
    % Time Step Size
    th = 0.001;

    % calculating number of spatial intervals
    n = (xN - x0)/xh + 1;

    xs = [x0:xh:xN];

    velocity = input("Enter Velocity : ");

    % we find the new value of u using the prev values
    % computing the values at all discrete points of x 
    % using Initial condition
    prev_u = compute(xs, velocity, t0);
    new_u = zeros(1, n);

    p = velocity*th/xh;
    
    i = 0;
    if velocity > 0
        for t = t0+th : th : tn
            for j = n : -1 : 2
                new_u(j) = prev_u(j) - p*(prev_u(j) - prev_u(j-1));
            end
            i = i + 1;
            % using Boundary condition
            new_u(1) = new_u(n);
            
            % assigning the calcualted values to prev values 
            % for next iteration
            prev_u = new_u;

             if i == 250
                 plot(xs, compute(xs,velocity,t));
                 hold on;
                 plot(xs, new_u);
                 hold off;
                 legend('Actual Values', 'Calculated Values');
             end
        end
    end

    if velocity < 0
        for t = t0+th : th : tn
            for j = 1 : n-1
                new_u(j) = prev_u(j) - p*(prev_u(j+1) - prev_u(j));
            end
            i = i+1;
            % using Boundary condition
            new_u(n) = new_u(1);

            % assigning the calcualted values to prev values 
            % for next iteration
            prev_u = new_u;

            if i == 250
                 plot(xs, compute(xs,velocity,t));
                 hold on;
                 plot(xs, new_u);
                 hold off;
                 legend('Actual Values', 'Calculated Values');
            end
        end
    end
end

%% ------------------------------------------------------------------------
% Initial Condition (at t = t0 = 0)
function y = compute(x,velocity,t)
    x = x - velocity*t;
    y = sin(2*x);
end


%% ------------------------------------------------------------------------
%{
    For Stability :-
        |velocity*th/xh| <= 1
    So for Unstability
        |velocity*th/xh| > 1
        -1 > velocity*th/xh > 1
        -xh/xt > velocity > xh/xt
        -24.5436 > velocity > 24.5436
%}