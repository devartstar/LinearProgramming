function calculate1()
    h=0.2;
    x0=0;xn=2;
    y0=0;u0=0;upn=0;ypn=0;
    
    n=(xn-x0)/h;
    z0=[y0;u0];
    A=zeros(2*(n),2*(n));
    B=zeros(2*(n),1);
    x=x0+h;
    for i=1:2:2*(n)
        coeffAi=[1/h^2-(x-3)/(2*h) 0;0 1/h^2+1/h];
        coeffBi=[-2/h^2 -6;x -2/h^2];
        coeffCi=[1/h^2+(x-3)/(2*h) 0;0 1/h^2-1/h];
        coeffDi=[x^2;4*x+3];
        if i~=1
            A(i:i+1,i-2:i-1)=coeffAi;
        end
        A(i:i+1,i:i+1)=coeffBi;
        if i<2*(n)-1
            A(i:i+1,i+2:i+3)=coeffCi;
        end
        if i==1
            B(i:i+1,1)=coeffDi-coeffAi*z0;
        elseif i==2*n-1
            A(i:i+1,i-2:i-1)=coeffAi+coeffCi;
            B(i:i+1,1)=coeffDi;
        else
            B(i:i+1,1)=coeffDi;
        end
        x=x+h;
    end

    %solving using thomas algorithm and getting the values of z and y
    %zs stores all those values
    zs=blockTriDiagonal(A,B);
    ys=zs(1:2:2*(n));
    us=zs(2:2:2*n);
    %ys(n)=ys(n-1)+(zs(2*n-2)+zn(2))*h/2
    for i=1:n
        fprintf("value of y at %f is %f\n",x0+i*h,ys(i));
    end
    fprintf("\n\n");

    for i=1:n
        fprintf("value of z at %f is %f\n",x0+i*h,us(i));
    end
    
end




function y=blockTriDiagonal(A,B)
    [n,~]=size(A);
    %standard thomas algorithm but instead of elements we pick us matrices
    %to modify
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