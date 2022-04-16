function practice()
    x0 = 0; xN = pi;
    xh = pi/128;

    t0 = 0; tn = 1;
    th = 0.001;

    n = (xN - x0) / xh + 1;
    xs = [x0 : xh : xN];

    velocity = input("Enter Velocity : ");

    previous_u = compute(xs, velocity, t0);
    new_u = zeros(1,n);

    p = velocity * th / xh;

    i=0;
    if velocity > 0
        for t = t0 + th : th : tn
            for j = n : -1 : 2
                new_u(j) = previous_u(j) + p*(previous_u(j) - previous_u(j-1));
            end
            i = i + 1;
            new_u(1) = new_u(n);
            previous_u = new_u;
            if t == 0.25
                plot(xs, compute(xs, velocity, t));
                hold on;
                plot(xs, new_u);
                hold off;
            end
        end
    end

    if velocity < 0
        for t = t0 + th : th : tn
            for j = 1 : n-1
                new_u(j) = previous_u(j) - p * (previous_u(j+1) - previous_u(j));
            end
            i =  i + 1
            new_u(n) = new_u(1);
            previous_u = new_u;
            if t == 0.25
                plot(xs, compute(xs, velocity, t));
                hold on;
                plot(xs, new_u);
                hold off;
            end
        end
    end
end

function y = compute(x, velocity, t)
    x = x - velocity*t;
    y = sin(2 * x);
end