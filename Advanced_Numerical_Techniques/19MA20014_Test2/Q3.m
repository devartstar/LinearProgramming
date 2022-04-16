 %Matlab Program: Heat Diffusion in one dimensional wire using the
% Crank-Nicholson Method

%Consider a thin metal wire of length one meter with no heat exchange with
%its surroundings such that the heat conductivity of wire is 0.5. Initially
%the wire is heated with a temperature distribution sin pi x at time 
%t=0 while both ends of the wire are kept at a fixed temperature of 0 c for
%all times t. Then solve the one-dimensional heat equation to obtain the
%temperature distribution at a later time t using Crank Nicolson method by
%solving a triadiagonal system of linear equations.


clear;
% Parameters to define the heat equation and the range in space and time
L = 1.; % Lenth of the wire
T =1.; % Final time
% Parameters needed to solve the equation within the Crank-Nicholson method
maxk = 2500; % Number of time steps
dt = T/maxk; %Time step
n = 50.; % Number of space steps
dx = L/n; %Space step
cond = 0.25; % Conductivity
r = cond*dt/(dx*dx); % Parameter of the method


% Q.1 Write down a three line code for entering initial temperature
% distribution. 
x0 = 0;
xn = x0+L;
i = 1;
for x = x0 : dx : x0+L
    u(i,0) = sin(x*pi)
    i = i + 1;
end







% Temperature at the boundary (T=0)
for k=1:maxk+1
u(1,k) = 0.;
u(n+1,k) = 0.;
time(k) = (k-1)*dt;
end

% Q.2 Write down a three line code which includes a time loop and call the
%function Tridiag(r, RHS) (See at the end of code) to write the Crank 
%Nicolson iteration scheme. 
for t = t0+th : th : T
    Tridiag(r, u);
end










%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(x,u(:,100),'-',x,u(:,300),'-',x,u(:,500),'-',x,u(:,1500),'-')
title('Temperature within the Crank-Nicholson method')
xlabel('X')
ylabel('T')


% Q.3 Write down a small code to plot a 3-D graph of u vs x and t.

%% plotting the graph for temperature distribution with change in x and t
subplot(2,1,2);
[X,Y] = meshgrid(0:dx:L,0:dt:T); % making the grids for x and t
surf(X,Y,u);
colormap("default")
shading interp
title('Temperature distribution with change in x and t')






%TDMA function
function u = Tridiag(r,RHS)
    n=length(RHS)-1;
    a=zeros(n-1,1);
    b=zeros(n-1,1);
    c=zeros(n-1,1);
    d=zeros(n-1,1);
    nc=zeros(n-1,1);
    nd=zeros(n-1,1);
    u=zeros(n+1,1);    
    for i=1:n-1
        a(i)=r/2;
        b(i)=-(r+1);
        c(i)=r/2;
        d(i)=-r/2 * RHS(i+2) + (r-1) * RHS(i+1) - r/2 * RHS(i);
    end
    
    nc(1)=c(1)/b(1);
    nd(1)=d(1)/b(1);
    for i=2:n-1
        nc(i)=c(i)/(b(i)-a(i)*(nc(i-1)));
        nd(i)=(d(i)-a(i)*nd(i-1))/(b(i)-a(i)*(nc(i-1)));
    end
    
% Q.4 Write down a small code using a for loop (Back substitution loop) for
%returning the value of u to the function Tridiag. 
    u(n+1) = nd(n+1);
    for i = n:-1:1
        u(i) = nd(i) - nc(i)*u(i+1);
    end

end
