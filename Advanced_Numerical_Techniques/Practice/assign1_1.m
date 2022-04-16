%{
    NAME : Devjit Choudhury
    ROLL : 19MA20014
%}

%{
    Solving a BVP with deriviative Boundary Condition
    using finite differencr method --- results in fictitous points
    
    IN GENERAL :- consider the BVP
    -------------------------------
    y'' + p(x)y' + q(x)y = r(x) -------------<1>
    BC  -                       -------------<2>
            a0.y(a) - a1.y'(a) = G1
            b0.y(b) + b1.y'(b) = G2
    Discretizing <1> we get :-
    -------------------------------
    A(i).Y(i-1) + B(i).Y(i) + C(i).Y(i+1) = h^2.r(i)    ----------<3>
    A(i) = 1 - h.p(i)/2
    b(i) = -2 + h^2.q(i)
    c(i) = 1 + h.p(i)/2
        for i = 1,2,3...n (we get unknows y(0) y(1) ... y(n) y(n+1))
    
    -------------------------------
    Discretizing <2> we get :-  
    a0.y(0) - a1.y'(0) = G1 ---> a0.y(0) - a1.(y(1) - y(-1))/2.h = G1
    y(-1) = 2.h.G1/a1 + y(1) - 2.h.a0.y(0)/a1       --------------<4>
    b0.y(n+1) + b1.y'(n+1) = G2 --> b0.y(n+1) + b1.(y(n+2) - y(n))/2.h = G2
    y(n+2) = y(n) - 2.h.b0.y(n+1)/b1 + 2.h.G2/b1    --------------<5>
    
    -------------------------------
    i = 0 in equation <3> and substitute equation <4>
    A(0).Y(-1) + B(0).Y(0) + C(0).Y(1) = h^2.r(0)
    F0 -->
    (B(0) - 2.h.a0.A(0)/a1).Y(0) + (A(0) + C(0)).Y(1) = h^2.r(0) - 2.h.G1.A(0)/a1
    
    i = n+1 in equation <3> and substitute equation <5>
    A(n+1).Y(n) + B(n+1).Y(n+1) + C(n+1).Y(n+2) = h^2.r(n+1)
    FN+1 -->
    (A(n+1) + C(n+1)).Y(n) + (B(n+1) - 2.h.b0.C(n+1)/b1).Y(n+1) = h^2.r(n+1) - 2.h.G2.C(n+1)/b1
    
    F1, F2, .... FN -> using i=1,2,..,n in equation <3>
    ---------------------------------
    WE GET THE COFFECIENT MATRIX :-
    [(B(0) - 2.h.a0.A(0)/a1) (A(0) + C(0))  0    0 ..... ......  0
    [        A1                   B1
    [        0                    A2        B2   C2 .....
    . ....
    . ....
    [        0                    0         0    0 ... (A(n+1) + C(n+1))(B(n+1) - 2.h.b0.C(n+1)/b1) ]
%}

function assignment2()
    % -----------------------------------------------------------------
    x0 = 0;
    xn = 2;
    h = 0.02;
    %{ 
        % BC -> 
            a0y(x(0)) - a1y'(x(0)) = G1
            b0y(x(n)) + b1y'(x(n)) = G2
    %}
    a0 = -1; a1 = -1; G1 = 0;
    b0 = 1; b1 = 1; G2 = 28 * exp(8);
    % -----------------------------------------------------------------

    x = [x0 + h : h : xn - h];   % x1, x2, x3, ...., xn
    n = length(x);

    A = zeros(n,1);
    B = zeros(n,1);
    C = zeros(n,1);
    D = zeros(n,1);

    % -----------------------------------------------------------------
    % DISCRETIZATION THE DIFFERENTIAL EQUATION
    % A(i)Y(i-1) + B(i)Y(i) + C(i)Y(i+1) = D(i)
    for i = 1 : n
        A(i) = calcA(h, x(i));
        B(i) = calcB(h, x(i));
        C(i) = calcC(h, x(i));
        D(i) = calcD(h, x(i));
    end
    %-------------------------------------------------------------------

    %-------------------------------------------------------------------
    % HANDELING THE FICTITOUS VARIABLES & BUILDING MATRIX SYSTEM

    coefMatrix = zeros(n+2, n+2);
    coefMatrix(1, 1) = calcB(h, x0) - (2 * h * a0 * calcA(h, x0))/a1;
    coefMatrix(1, 2) = calcA(h, x0) + calcC(h, x0);
    coefMatrix(n+2, n+1) = calcA(h, xn) + calcC(h, xn);
    coefMatrix(n+2, n+2) = calcB(h, xn) - (2 * h * b0 * calcC(h, xn))/b1;

    for i = 2 : n+1
        coefMatrix(i, i-1) = A(i-1);
        coefMatrix(i,i) = B(i-1);
        coefMatrix(i,i+1) = C(i-1);
    end
    %-------------------------------------------------------------------

    % Building the Value Matrix in RHS
    rhsMatrix = zeros(n+2, 1);
    rhsMatrix(1) = calcD(h, x0) - (2 * h * G1 * calcA(h,x0))/a1;
    rhsMatrix(n+2) = calcD(h, xn) - (2 * h * G2 * calcC(h,xn))/b1;
    for i = 2 : n+1
        rhsMatrix(i) = D(i-1);
    end
    %-----------------------------------------------------------------
    %{
    fprintf("Coef Matrix : \n");
    disp(coefMatrix);
    fprintf("Value Matrix : \n");
    disp(rhsMatrix);
    %}
    sol= thomasAlgorithm(coefMatrix, rhsMatrix);
    % Displaying the output
    fprintf("Answer :- \n");
    for i = 1:n+1   
        fprintf("x %d (%d) = %f\n",i,x0 + i*h, sol(i));
    end
    fprintf("\n");
    % ------------------------------------------------------------

    % Trying D solve
    %{
    syms y(t) a
    eqn = diff(y,x,2) - 8*diff(y,x,1) + 16*y == 0;
    ySol(x) = dsolve(eqn)
    %}
    
    % Plotting the result
    xs = [0:h:2];
    plot(xs, (3*xs-1).*exp(4*xs), '--o');
    hold on;
    plot(xs, sol, '-.*');
    legend('Actual Values', 'Computed Values');
    title('Numerical solutions For h = 0.02');
    hold off
end

%% ------------------------------------------------------------
% Thomas Algorithm to solve the tri diagonal system
function y = thomasAlgorithm(A,B)
    [rows,~]=size(A);   
    C = zeros(1,rows);
    D = zeros(1,rows);

    C(1) = A(1,2) / A(1,1);
    D(1) = B(1) / A(1,1);
    
    % modifying the matrix to the desired form
    for i = 2:rows
        if i~=rows
            C(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*C(i-1));
        end
        D(i)=(B(i)-A(i,i-1)*D(i-1))/(A(i,i)-A(i,i-1)*C(i-1));
    end
    y = zeros(rows,1);
    y(rows) = D(rows);

    % Backtracking to reconstruct the solution 
    for i = rows-1:-1:1
        y(i) = D(i)-C(i)*y(i+1);
    end
end


%% ------------------------------------------------------------
function y = coefOfFirstDeriviative(x)
   % y = p(x)
   y = -8;
end

%% ------------------------------------------------------------
function y = coefOfY(x)
   % y = q(x)
   y = 16;
end

%% ------------------------------------------------------------
function y = rhs(x)
    % y = r(x)
    y = 0;
end

%% ------------------------------------------------------------
function y = calcA(h, x)
    y = 1 - (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
function y = calcB(h, x)
    y = -2 + (h * h * coefOfY(x));
end

%% ------------------------------------------------------------
function y = calcC(h,x)
    y = 1 + (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
function y = calcD(h, x)
    y = (h * h * rhs(x));
end


%{
Answer :- 
x 1 = -0.966687
x 2 = -0.984474
x 3 = -0.998506
x 4 = -1.008032
x 5 = -1.012202
x 6 = -1.010057
x 7 = -1.000511
x 8 = -0.982346
x 9 = -0.954188
x 10 = -0.914494
x 11 = -0.861536
x 12 = -0.793375
x 13 = -0.707840
x 14 = -0.602505
x 15 = -0.474660
x 16 = -0.321279
x 17 = -0.138988
x 18 = 0.075972
x 19 = 0.327789
x 20 = 0.621119
x 21 = 0.961142
x 22 = 1.353614
x 23 = 1.804925
x 24 = 2.322170
x 25 = 2.913215
x 26 = 3.586785
x 27 = 4.352546
x 28 = 5.221205
x 29 = 6.204613
x 30 = 7.315886
x 31 = 8.569531
x 32 = 9.981587
x 33 = 11.569780
x 34 = 13.353696
x 35 = 15.354963
x 36 = 17.597459
x 37 = 20.107538
x 38 = 22.914273
x 39 = 26.049732
x 40 = 29.549273
x 41 = 33.451870
x 42 = 37.800471
x 43 = 42.642390
x 44 = 48.029740
x 45 = 54.019900
x 46 = 60.676036
x 47 = 68.067667
x 48 = 76.271285
x 49 = 85.371036
x 50 = 95.459467
x 51 = 106.638342
x 52 = 119.019537
x 53 = 132.726022
x 54 = 147.892931
x 55 = 164.668744
x 56 = 183.216567
x 57 = 203.715549
x 58 = 226.362419
x 59 = 251.373181
x 60 = 278.984957
x 61 = 309.458016
x 62 = 343.077987
x 63 = 380.158279
x 64 = 421.042738
x 65 = 466.108546
x 66 = 515.769391
x 67 = 570.478943
x 68 = 630.734651
x 69 = 697.081893
x 70 = 770.118521
x 71 = 850.499825
x 72 = 938.943965
x 73 = 1036.237911
x 74 = 1143.243932
x 75 = 1260.906694
x 76 = 1390.261020
x 77 = 1532.440370
x 78 = 1688.686109
x 79 = 1860.357638
x 80 = 2048.943466
x 81 = 2256.073311
x 82 = 2483.531313
x 83 = 2733.270490
x 84 = 3007.428512
x 85 = 3308.344948
x 86 = 3638.580103
x 87 = 4000.935598
x 88 = 4398.476844
x 89 = 4834.557598
x 90 = 5312.846778
x 91 = 5837.357752
x 92 = 6412.480319
x 93 = 7043.015643
x 94 = 7734.214393
x 95 = 8491.818391
x 96 = 9322.106087
x 97 = 10231.942209
x 98 = 11228.831972
x 99 = 12320.980253
x 100 = 13517.356200
x 101 = 14827.763745
%}