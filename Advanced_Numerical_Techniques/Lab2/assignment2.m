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
    x0 = 0; xn = 2;     % Boundary points
    h = 0.02;           % step size
    %{ 
        % BC -> 
            a0y(x(0)) - a1y'(x(0)) = G1
            b0y(x(n)) + b1y'(x(n)) = G2
    %}
    a0 = -1; a1 = -1; G1 = 0;           % Variables describing Boundary condition 1
    b0 = 1; b1 = 1; G2 = 28 * exp(8);   % Variables describing Boundary condition 2
    % -----------------------------------------------------------------

    x = [x0 + h : h : xn - h];   % x(1), x(2), x(3), ...., x(n)
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

    % If we discretize the boundary conditions at x(0) = 0 & x(n+1) = 2
    % we get 2 fictitous points x(-1) & x(n+2)
    % using the formulae at top -> we replace x(-1) with x(0) & x(1) &
    % x(n+1) with x(n) & x(n-1)
    
    %-------------------------------------------------------------------
    % BUILDING TRIDIAGONAL MATRIX SYSTEM

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
    %-----------------------------------------------------------------
    % Calling THOMAS FUNCTION to solve the Tri-diagonal System
    sol= thomasAlgorithm(coefMatrix, rhsMatrix);
    % Displaying the output
    fprintf("Answer :- \n");
    for i = 1:n+2
        fprintf("x %d (%d) = %f\n",i,x0 + (i-1)*h, sol(i));
    end
    fprintf("\n");
    % ------------------------------------------------------------


    %-------------------------------------------------------------
    % Plotting the result
    % Using D-solve for plotting the Actual Values
    syms y;
    y = dsolve('D2y - 8*Dy + 16*y = 0', 'Dy(0) - y(0) = 0', 'Dy(2) + y(2) = 28*exp(8)'); 
    fprintf("The Actual Solution equation : \n");
    disp(y);
    fplot(y, [0 2]);
    hold on;    
    xs = [0:h:2];
    plot(xs, sol, '-.*');
    legend('Actual Values', 'Computed Values');
    title('Numerical solutions For h = 0.02');
    hold off
    %-------------------------------------------------------------
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
% the function corresponding to p(x) as writen at top
function y = coefOfFirstDeriviative(x)
   % y = p(x)
   y = -8;
end

%% ------------------------------------------------------------
% the function corresponding to q(x) as writen at top
function y = coefOfY(x)
   % y = q(x)
   y = 16;
end

%% ------------------------------------------------------------
% the function corresponding to r(x) as writen at top
function y = rhs(x)
    % y = r(x)
    y = 0;
end

%% ------------------------------------------------------------
% the function corresponding to A(x) as writen at top after DISCRETIZATION
function y = calcA(h, x)
    y = 1 - (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
% the function corresponding to B(x) as writen at top after DISCRETIZATION
function y = calcB(h, x)
    y = -2 + (h * h * coefOfY(x));
end

%% ------------------------------------------------------------
% the function corresponding to C(x) as writen at top after DISCRETIZATION
function y = calcC(h,x)
    y = 1 + (h * coefOfFirstDeriviative(x)) / 2;
end

%% ------------------------------------------------------------
% the function corresponding to RHS of the DISCRETIZED EQUATION
function y = calcD(h, x)
    y = (h * h * rhs(x));
end


%{
Answer :- 
x 1 (0) = -0.966687
x 2 (2.000000e-02) = -0.984474
x 3 (4.000000e-02) = -0.998506
x 4 (6.000000e-02) = -1.008032
x 5 (8.000000e-02) = -1.012202
x 6 (1.000000e-01) = -1.010057
x 7 (1.200000e-01) = -1.000511
x 8 (1.400000e-01) = -0.982346
x 9 (1.600000e-01) = -0.954188
x 10 (1.800000e-01) = -0.914494
x 11 (2.000000e-01) = -0.861536
x 12 (2.200000e-01) = -0.793375
x 13 (2.400000e-01) = -0.707840
x 14 (2.600000e-01) = -0.602505
x 15 (2.800000e-01) = -0.474660
x 16 (3.000000e-01) = -0.321279
x 17 (3.200000e-01) = -0.138988
x 18 (3.400000e-01) = 0.075972
x 19 (3.600000e-01) = 0.327789
x 20 (3.800000e-01) = 0.621119
x 21 (4.000000e-01) = 0.961142
x 22 (4.200000e-01) = 1.353614
x 23 (4.400000e-01) = 1.804925
x 24 (4.600000e-01) = 2.322170
x 25 (4.800000e-01) = 2.913215
x 26 (5.000000e-01) = 3.586785
x 27 (5.200000e-01) = 4.352546
x 28 (5.400000e-01) = 5.221205
x 29 (5.600000e-01) = 6.204613
x 30 (5.800000e-01) = 7.315886
x 31 (6.000000e-01) = 8.569531
x 32 (6.200000e-01) = 9.981587
x 33 (6.400000e-01) = 11.569780
x 34 (6.600000e-01) = 13.353696
x 35 (6.800000e-01) = 15.354963
x 36 (7.000000e-01) = 17.597459
x 37 (7.200000e-01) = 20.107538
x 38 (7.400000e-01) = 22.914273
x 39 (7.600000e-01) = 26.049732
x 40 (7.800000e-01) = 29.549273
x 41 (8.000000e-01) = 33.451870
x 42 (8.200000e-01) = 37.800471
x 43 (8.400000e-01) = 42.642390
x 44 (8.600000e-01) = 48.029740
x 45 (8.800000e-01) = 54.019900
x 46 (9.000000e-01) = 60.676036
x 47 (9.200000e-01) = 68.067667
x 48 (9.400000e-01) = 76.271285
x 49 (9.600000e-01) = 85.371036
x 50 (9.800000e-01) = 95.459467
x 51 (1) = 106.638342
x 52 (1.020000e+00) = 119.019537
x 53 (1.040000e+00) = 132.726022
x 54 (1.060000e+00) = 147.892931
x 55 (1.080000e+00) = 164.668744
x 56 (1.100000e+00) = 183.216567
x 57 (1.120000e+00) = 203.715549
x 58 (1.140000e+00) = 226.362419
x 59 (1.160000e+00) = 251.373181
x 60 (1.180000e+00) = 278.984957
x 61 (1.200000e+00) = 309.458016
x 62 (1.220000e+00) = 343.077987
x 63 (1.240000e+00) = 380.158279
x 64 (1.260000e+00) = 421.042738
x 65 (1.280000e+00) = 466.108546
x 66 (1.300000e+00) = 515.769391
x 67 (1.320000e+00) = 570.478943
x 68 (1.340000e+00) = 630.734651
x 69 (1.360000e+00) = 697.081893
x 70 (1.380000e+00) = 770.118521
x 71 (1.400000e+00) = 850.499825
x 72 (1.420000e+00) = 938.943965
x 73 (1.440000e+00) = 1036.237911
x 74 (1.460000e+00) = 1143.243932
x 75 (1.480000e+00) = 1260.906694
x 76 (1.500000e+00) = 1390.261020
x 77 (1.520000e+00) = 1532.440370
x 78 (1.540000e+00) = 1688.686109
x 79 (1.560000e+00) = 1860.357638
x 80 (1.580000e+00) = 2048.943466
x 81 (1.600000e+00) = 2256.073311
x 82 (1.620000e+00) = 2483.531313
x 83 (1.640000e+00) = 2733.270490
x 84 (1.660000e+00) = 3007.428512
x 85 (1.680000e+00) = 3308.344948
x 86 (1.700000e+00) = 3638.580103
x 87 (1.720000e+00) = 4000.935598
x 88 (1.740000e+00) = 4398.476844
x 89 (1.760000e+00) = 4834.557598
x 90 (1.780000e+00) = 5312.846778
x 91 (1.800000e+00) = 5837.357752
x 92 (1.820000e+00) = 6412.480319
x 93 (1.840000e+00) = 7043.015643
x 94 (1.860000e+00) = 7734.214393
x 95 (1.880000e+00) = 8491.818391
x 96 (1.900000e+00) = 9322.106087
x 97 (1.920000e+00) = 10231.942209
x 98 (1.940000e+00) = 11228.831972
x 99 (1.960000e+00) = 12320.980253
x 100 (1.980000e+00) = 13517.356200
x 101 (2) = 14827.763745
 
The Actual Solution equation : 
3*t*exp(4*t) - exp(4*t)
%}
