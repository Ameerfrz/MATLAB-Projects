clc;clear;
close all;
[num,txt,raw] = xlsread('fp.xlsx');
L = num(1,1); %width
H = num(1,2); %height
alpha = num(1,3);
dx = num(1,4);
dt = num(1,5);
T0 = num(1,6); %T(t=0)
TL = num(1,7); %B.C
TR = num(1,8); %B.C
TU = num(1,9); %B.C
TD = num(1,10); %B.C
t_final = num(1,11);  % zaman nahaii
gdot  = num(1,12); % meghdar shar  
roci  = num(1,13);  %chegali va garmaie vizhe
r = alpha*dt/(dx)^2;

%%
n = floor(H/dx)+1;  %tedad noghat har soton
m = floor(L/dx)+1;  %tedad noghat har satr
q = floor(t_final/dt);  %marahel hal(az nazar zaman)

%% Explicit
if r>0.25
    fprintf('error in explicit');
end

TEX = zeros(n,(q+1)*m);
for i=1:n
    for j=1:m
        TEX(i,j)= T0; %initial value
    end
end

TEX(1,:)=TU;            %boundry condition
TEX(n,:)=TD;            %boundry condition
for nt=0:q
    TEX(:,nt*m+1)=TL;   %boundry condition
    TEX(:,nt*m+m)=TR;   %boundry condition
end


for nt=1:q
    for i=2:m-1
        for j=2:n-1
            TEX(j,i+(nt*m)) = r*( TEX(j,i+((nt-1)*m)+1) + TEX(j,i+((nt-1)*m)-1) + TEX(j+1,i+((nt-1)*m)) + TEX(j-1,i+((nt-1)*m)) )+(1-4*r)*TEX(j,i+((nt-1)*m))+ (gdot*dt/roci);
        end
    end
end
   


%% Implicit
TIM = zeros(n,(2*q+1)*m);
for i=1:n
    for j=1:m
        TIM(i,j)= T0; 
    end
end

TIM(1,:)=TU;
TIM(n,:)=TD;
for nt=0:(2*q)
    TIM(:,nt*m+1)=TL;
    TIM(:,nt*m+m)=TR;
end


majhool = (m-2)*(n-2); %tedad majhoolat

for nt=1:(2*q)
    
    a = zeros(majhool,majhool);
    b = zeros(majhool,1);
    xgg = zeros(majhool,1);
    for i=1:majhool
        xgg(i,1) = T0;
    end
    
    for i=0:(n-3)
        for j=2:(m-1)
            % marhale alal ADI
            if floor(nt/2)<(nt/2)
                % noghat vasat
                if i>0&i<(n-3)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i+1)*(m-2)))) = -r;
                    a(((j-1)+(i*(m-2))),((j-1)+((i-1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+2,j+1+((nt-1)*m)) + (r)*TIM(i+2,j-1+((nt-1)*m)) ;
                    
                % noghat konj
                elseif i==0
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i+1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+2,j+1+((nt-1)*m)) + (r)*TIM(i+2,j-1+((nt-1)*m)) + (r)*TIM(i+1,j+((nt-1)*m));
                elseif i==(n-3)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i-1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+2,j+1+((nt-1)*m)) + (r)*TIM(i+2,j-1+((nt-1)*m)) + (r)*TIM(i+3,j+((nt-1)*m));
                end
                
%%                 
            % marhale dovom ADI    
            elseif floor(nt/2)==(nt/2)
                % noghat vasat
                if j>2&j<(m-1)
                    a(((j-1)+(i*(m-2))),((j-1)+((i)*(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j  )+((i)*(m-2)))) = -r;
                    a(((j-1)+(i*(m-2))),((j-2)+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+3,j+((nt-1)*m)) + (r)*TIM(i+1,j+((nt-1)*m));
                % noghat konj
                elseif j==2
                    a(((j-1)+(i*(m-2))),((j-1)+((i)*(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j  )+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+3,j+((nt-1)*m)) + (r)*TIM(i+1,j+((nt-1)*m)) + (r)*TIM(i+2,j-1+((nt-1)*m));
                elseif j==(m-1)
                    a(((j-1)+(i*(m-2))),((j-1)+((i)*(m-2)))) = 2*r+1;
                    a(((j-1)+(i*(m-2))),((j-2)+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-2*r)*TIM(i+2,j+((nt-1)*m)) + (r)*TIM(i+3,j+((nt-1)*m)) + (r)*TIM(i+1,j+((nt-1)*m)) + (r)*TIM(i+2,j+1+((nt-1)*m));    
                end
                
            end
        end
        
    end
    
       
    
    %% solving with Gauss sidel method
            nn = (n-2)*(m-2);
            jjj=0;
            v = zeros(nn,1);
            v(1,1)=1;
            alpha = 0.8;


            while max(v)>0.0001
                jjj=jjj+1;
                if jjj>30
                    disp('the implicit method isnt converge');
                    break;
                end

                ss = zeros(nn,1);
                for i=1:nn
                    for j=1:nn
                        if j~=i
                           ss(i,1) = ss(i,1)+ a(i,j)*xgg(j,1);
               
                        end
                    end    
        
                    x(i,1) = xgg(i,1) + (alpha/a(i,i)) * ( b(i,1)-ss(i,1)- a(i,i)*xgg(i,1) );
                    v(i,1) = abs(x(i,1)-xgg(i,1));
                    xgg(i,1) = x(i,1);
                end
            end
            

            for i=2:n-1
                for j=2:m-1
                    TIM(i,j+((nt)*m)) = x(((j-1)+((i-2)*(m-2))),1);
                end
            end     
end


 















% 
%% Crank
TCR = zeros(n,(2*q+1)*m);
for i=1:n
    for j=1:m
        TCR(i,j)= T0; 
    end
end

TCR(1,:)=TU;
TCR(n,:)=TD;
for nt=0:(2*q)
    TCR(:,nt*m+1)=TL;
    TCR(:,nt*m+m)=TR;
end


majhool = (m-2)*(n-2); %tedad majhoolat

for nt=1:(2*q)
    
    a = zeros(majhool,majhool);
    b = zeros(majhool,1);
    xgg = zeros(majhool,1);
    for i=1:majhool
        xgg(i,1) = T0;
    end
    
    for i=0:(n-3)
        for j=2:(m-1)
            % marhale alal ADI
            if floor(nt/2)<(nt/2)
                % noghat vasat
                if i>0&i<(n-3)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i+1)*(m-2)))) = -r;
                    a(((j-1)+(i*(m-2))),((j-1)+((i-1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (2*r)*TCR(i+2,j+1+((nt-1)*m)) + (2*r)*TCR(i+2,j-1+((nt-1)*m)) + (r)*TCR(i+3,j+((nt-1)*m)) + (r)*TCR(i+1,j+((nt-1)*m)) ;
                    
                % noghat konj
                elseif i==0
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i+1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (2*r)*TCR(i+2,j+1+((nt-1)*m)) + (2*r)*TCR(i+2,j-1+((nt-1)*m)) + (r)*TCR(i+3,j+((nt-1)*m)) + (2*r)*TCR(i+1,j+((nt-1)*m)) ;
                elseif i==(n-3)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j-1)+((i-1)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (2*r)*TCR(i+2,j+1+((nt-1)*m)) + (2*r)*TCR(i+2,j-1+((nt-1)*m)) + (2*r)*TCR(i+3,j+((nt-1)*m)) + (r)*TCR(i+1,j+((nt-1)*m)) ;
                end



            % marhale dovom ADI    
            elseif floor(nt/2)==(nt/2)
                % noghat vasat
                if j>2&j<(m-1)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j-2)+((i)*(m-2)))) = -r;
                    a(((j-1)+(i*(m-2))),((j  )+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (r)*TCR(i+2,j+1+((nt-1)*m)) + (r)*TCR(i+2,j-1+((nt-1)*m)) + (2*r)*TCR(i+3,j+((nt-1)*m)) + (2*r)*TCR(i+1,j+((nt-1)*m)) ;
                % noghat konj
                elseif j==2
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j  )+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (r)*TCR(i+2,j+1+((nt-1)*m)) + (2*r)*TCR(i+2,j-1+((nt-1)*m)) + (2*r)*TCR(i+3,j+((nt-1)*m)) + (2*r)*TCR(i+1,j+((nt-1)*m)) ;
                elseif j==(m-1)
                    a(((j-1)+(i*(m-2))),((j-1)+( i   *(m-2)))) = 4*r+1;
                    a(((j-1)+(i*(m-2))),((j-2)+((i)*(m-2)))) = -r;
                    b(((j-1)+(i*(m-2))),1)  = (1-4*r)*TCR(i+2,j+((nt-1)*m)) + (2*r)*TCR(i+2,j+1+((nt-1)*m)) + (r)*TCR(i+2,j-1+((nt-1)*m)) + (2*r)*TCR(i+3,j+((nt-1)*m)) + (2*r)*TCR(i+1,j+((nt-1)*m)) ;
                end
                
            end
        end
        
    end
    
    
    
    
    
    
    
    
    %% solving with Gauss sidel method
            nn = (n-2)*(m-2);
            jjj=0;
            v = zeros(nn,1);
            v(1,1)=1;
            alpha = 0.8;


            while max(v)>0.0001
                jjj=jjj+1;
                if jjj>30
                    disp('the crank method isnt converge');
                    break;
                end

                ss = zeros(nn,1);
                for i=1:nn
                    for j=1:nn
                        if j~=i
                           ss(i,1) = ss(i,1)+ a(i,j)*xgg(j,1);
               
                        end
                    end    
        
                    x(i,1) = xgg(i,1) + (alpha/a(i,i)) * ( b(i,1)-ss(i,1)- a(i,i)*xgg(i,1) );
                    v(i,1) = abs(x(i,1)-xgg(i,1));
                    xgg(i,1) = x(i,1);
                end
            end
            

            for i=2:n-1
                for j=2:m-1
                    TCR(i,j+((nt)*m)) = x(((j-1)+((i-2)*(m-2))),1);
                end
            end     
 end


point_result_EX = zeros(1,q+1);
point_result_IM = zeros(1,q+1);
point_result_CR = zeros(1,q+1);
for i=1:(q+1)
    point_result_EX(1,i) = TEX(ceil(n/2),(ceil(m/2)+(i-1)*m));
    point_result_IM(1,i) = TIM(ceil(n/2),(ceil(m/2)+(i-1)*2*m));
    point_result_CR(1,i) = TCR(ceil(n/2),(ceil(m/2)+(i-1)*2*m));
    Time(1,i)            = (i-1)*dt; 
end

plot(Time,point_result_EX,'b','linewidth',6);
hold on
plot(Time,point_result_IM,'r','linewidth',6);
hold on
plot(Time,point_result_CR,'k','linewidth',6);
grid on
xlabel('Time(sec)','fontsize',25);
ylabel('Temperature(^c)','fontsize',25);
title('Center of Steel Plate','fontsize',25);
legend('EX','IM','CR')

%
% ttime=2;
% for i=1:n
%     for j=1:m
%         T_surf(i,j) = TEX(i,j+(m)*ttime);
%     end
% end
% [x,y]=meshgrid(0:dx:L,0:dx:H);
% surf(x,y,T_surf)








