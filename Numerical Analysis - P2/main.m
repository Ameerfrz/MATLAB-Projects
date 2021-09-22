clc;clear;
[num,txt,raw] = xlsread('secondproject.xlsx');
nn = size(num);
n = nn(1,1);
a = zeros(n,n);
for i=1:n
    b(i,1)=num(i,n+1);
    for j=1:n
        a(i,j)=num(i,j);
    end
    xgg(i,1) = num(i,n+2);
end
x = zeros(n,4);
      


%% Gaussian elimination
A = zeros(n,n*n);
B = zeros(n,n);

for i=1:n
    B(i,1)=num(i,n+1);
    for j=1:n
        A(i,j)=num(i,j);
    end
end
for k=1:(n-1)
    for i=(k+1):n
        for j=(k+1):n
            A(i,k*n+j) = A(i,(k-1)*n+j) - (A(i,(k-1)*n+k)/A(k,(k-1)*n+k))*A(k,(k-1)*n+j);
            B(i,k+1) = B(i,k) - (A(i,(k-1)*n+k)/A(k,(k-1)*n+k))*B(k,k);
        end
    end
end

x(n,1) = B(n,n)/A(n,n*n);
for i=(n-1):(-1):1
    s=0;
    for j=i+1:n
        s = s+ A(i,(i-1)*n+j)*x(j,1);
    end
    x(i,1) = (B(i,i) - s)/A(i,(i-1)*n+i);
end
    
%% LU Decomp
l = zeros(n,n);
u = zeros(n,n);

for i=1:n
    l(i,1) = a(i,1);
    u(i,i) = 1;
    u(1,i) = a(1,i)/l(1,1);
end

for j = 2:n
    for i = 2:n
        if i<j
            s=0;
            for k=1:(j-1)
                 s = s + l(i,k)*u(k,j);
            end
            u(i,j) = (a(i,j)-s)/l(i,i);
        end
        
        if j<=i
            s=0;
            for k=1:(j-1)
                s = s + (l(i,k)*u(k,j));
            end
            l(i,j) = a(i,j)- s;
        end
            
        
    end
end
c(1,1) = b(1,1)/l(1,1);
for i=2:n
    s=0;
    for k=1:(i-1)
        s=s+l(i,k)*c(k,1);
    end
    c(i,1) = (b(i,1) - s) / l(i,i);
end

x(n,2) = c(n,1);
for i=(n-1):(-1):1
    s=0;
    for k=(i+1):n
        s=s+u(i,k)*x(k,2);
    end
    x(i,2) = c(i,1) - s;
end

%%
data3 = zeros(20,n);
data4 = zeros(20,n);

%% Jacobi
jj=0;
v = zeros(n,1);
v(1,1)=1;
alpha = 1;
xg = xgg;
while max(v)>0.0001
    jj=jj+1;
    if jj>1000
        disp('the 3th method isnt converging');
        break;
    end

    ss = zeros(n,1);
    for i=1:n
        for j=1:n
            if j~=i
               ss(i,1) = ss(i,1)+ a(i,j)*xg(j,1);
               
            end
        end    
        
        x(i,3) = xg(i,1) + (alpha/a(i,i)) * ( b(i,1)-ss(i,1)- a(i,i)*xg(i,1) );
        
    end
    for i=1:n
            v(i,1) = abs(x(i,3)-xg(i,1));
            xg(i,1) = x(i,3);
    end
end


%% Gauss sidel method
jjj=0;
v = zeros(n,1);
v(1,1)=1;
% alpha = 1;

while max(v)>0.0001
    jjj=jjj+1;
    if jjj>1000
        disp('the 4th method isnt converging');
        break;
    end

    ss = zeros(n,1);
    for i=1:n
        for j=1:n
            if j~=i
               ss(i,1) = ss(i,1)+ a(i,j)*xgg(j,1);
               
            end
        end    
        
        x(i,4) = xgg(i,1) + (alpha/a(i,i)) * ( b(i,1)-ss(i,1)- a(i,i)*xgg(i,1) );
        v(i,1) = abs(x(i,4)-xgg(i,1));
        xgg(i,1) = x(i,4);
        
    end
    
end
x
