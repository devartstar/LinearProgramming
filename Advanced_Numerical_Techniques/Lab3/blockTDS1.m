function calculate()
    h=0.1;
    x0=0;xn=1;
    y0=0;yp0=0;ypn=1;
    
    n=(xn-x0)/h;
    z0=[y0;yp0];
    zn=[0;ypn];
    A=zeros(2*(n-1),2*(n-1));
    B=zeros(2*(n-1),1);
    x=x0+h;
    for i=1:2:2*(n-1)
        coeffAi=[-1 -h/2;0 1/h^2-fA(x)/(2*h)];
        coeffBi=[1 -h/2;fC(x) -2/h^2+fB(x)];
        coeffCi=[0 0;0 1/h^2+fA(x)/(2*h)];
        coeffDi=[0;fD(x)];
        if i~=1
            A(i:i+1,i-2:i-1)=coeffAi;
        end
        A(i:i+1,i:i+1)=coeffBi;
        if i<2*(n-1)-1
            A(i:i+1,i+2:i+3)=coeffCi;
        end
        if i==1
            B(i:i+1,1)=coeffDi-coeffAi*z0;
        elseif i==2*n-3
            B(i:i+1,1)=coeffDi-coeffCi*zn;
        else
            B(i:i+1,1)=coeffDi;
        end
        x=x+h;
    end


    zs=blockTriDiagonal(A,B);
    ys=zs(1:2:2*(n-1));
    ys(n)=ys(n-1)+(zs(2*n-2)+zn(2))*h/2
    xs=[x0+h:h:xn];
end




function y=blockTriDiagonal(A,B)
    [n,~]=size(A);
    for i = 1:2:n
       if i==1
           const=A(i:i+1,i:i+1);
       else 
           const=A(i:i+1,i:i+1)-A(i:i+1,i-2:i-1)*A(i-2:i-1,i:i+1);
       end

       if i==1
           A(i:i+1,i+2:i+3)=inv(const)*A(i:i+1,i+2:i+3);
           B(i:i+1)=inv(const)*B(i:i+1);
       else
           B(i:i+1)=inv(const)*(B(i:i+1)-A(i:i+1,i-2:i-1)*B(i-2:i-1));
           A(i:i+1,i-2:i-1)=zeros(2,2);
           A(i:i+1,i:i+1)=eye(2,2);
           if i~=n-1
               A(i:i+1,i+2:i+3)=inv(const)*A(i:i+1,i+2:i+3);
           end
       end
       A(i:i+1,i:i+1)=eye(2,2);
    end
    z=zeros(n,1);


    for i=n:-2:2
        if i==n
            z(i-1:i,1)=B(i-1:i);
        else
            z(i-1:i,1)=B(i-1:i)-A(i-1:i,i+1:i+2)*z(i+1:i+2);
        end
    end
    y=z;
    
end




function y=fA(x)
    y=4;
end

%to calculate values of B(xk)
function y=fB(x)
    y=1;
end

%to calculate values of C(xk)
function y=fC(x)
    y=-6;
end

%to calculate values of D(xk)
function y=fD(x)
    y=1;
end