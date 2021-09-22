clc;clear;
[num,txt,raw] = xlsread('tp.xlsx');
m = num(1);     %mass
k = num(2);     %spring constant
c = num(3);     %damper constant
x0 = num(4);    %initial value of x
xx0 = num(5);   %initial value of speed(x dot)
ftime = num(6); %final time(from t=0 up to ftime)
h = num(7);     %delta t
t(1,1) = 0;     %time
%% analytical solution
wn = (k/m)^0.5;
keisi = c/2/((k*m)^0.5);
fi0 = atan(x0*wn*(1-keisi^2)^0.5/(xx0+keisi*wn*x0));
x0A = (x0^2*wn^2+xx0^2+2*x0*xx0*keisi*wn)^0.5/wn/(1-keisi^2)^0.5;
%xA = x0A*exp(-1*keisi*wn*t)*sin((1-keisi^2)^0.5*wn*t+fi0)
analytical_x = @(t) x0A*exp((-1*keisi*wn)*t)*sin((1-keisi^2)^0.5*wn*t+fi0);


%% Modified Euler
xME(1,1) = x0;
xxME(1,1) = xx0;
fME(1,1) = xxME(1,1);                         %x dot
gME(1,1) = (-1*c/m)*xxME(1,1) - (k/m)*xME(1,1);   %x double dot

final_step = ftime/h;
for i=1:final_step
    %predictor
    xpME(i+1,1) = xME(i,1) + h*fME(i,1);   %x
    xxpME(i+1,1) = xxME(i,1) + h*gME(i,1); %x dot
    %corrector
    fME(i+1,1) = xxpME(i+1,1);
    gME(i+1,1) = (-1*c/m)*xxpME(i+1,1) - (k/m)*xpME(i+1,1);
    xME(i+1,1) = xME(i,1) + h/2*(fME(i,1)+fME(i+1,1));
    xxME(i+1,1) = xxME(i,1) + h/2*(gME(i,1)+gME(i+1,1));
    t(i+1,1) = t(i,1)+h;
end


%% Runge kutta
xRK(1,1)  = x0;
yRK(1,1) = xx0;                  %x double dot
f = @ (y) y;                         %x dot
g = @ (x,y) (-1*c/m)*y-(k/m)*x;      %x double dot

for i=1:final_step
    kx1 = h*f(yRK(i,1));
    ky1 = h*g(xRK(i,1),yRK(i,1));
    
    kx2 = h*f(yRK(i,1)+ky1/2);
    ky2 = h*g(xRK(i,1)+kx1/2,yRK(i,1)+ky1/2);
    
    kx3 = h*f(yRK(i,1)+ky2/2);
    ky3 = h*g(xRK(i,1)+kx2/2,yRK(i,1)+ky2/2);
    
    kx4 = h*f(yRK(i,1)+ky3);
    ky4= h*g(xRK(i,1)+kx3,yRK(i,1)+ky3);
    
    xRK(i+1,1) = xRK(i,1)+(kx1+2*kx2+2*kx3+kx4)/6;
    yRK(i+1,1) = yRK(i,1)+(ky1+2*ky2+2*ky3+ky4)/6;
    
end
    

%% Adams-Moulton
for i=1:4                  %initial value + initial 3 step
    xAD(i,1) = xME(i,1);
    yAD(i,1) = xxME(i,1);  %gam haie avalie az javab modified euler mohasabe shode
end

f = @ (y) y;                         %x dot
g = @ (x,y) (-1*c/m)*y-(k/m)*x;      %x double dot

for i=4:final_step
    %predictor
    xAD(i+1,1) = xAD(i,1) + h/24*(55*f(yAD(i,1)) - 59*f(yAD(i-1,1)) + 37*f(yAD(i-2)) - 9*f(yAD(i-3)));
    yAD(i+1,1) = yAD(i,1) + h/24*(55*g(xAD(i,1),yAD(i,1)) - 59*g(xAD(i-1,1),yAD(i-1,1)) + 37*g(xAD(i-2,1),yAD(i-2)) - 9*g(xAD(i-3,1),yAD(i-3)));
    %correcor
    xAD(i+1,1) = xAD(i,1) + h/24*(9*f(yAD(i+1,1)) + 19*f(yAD(i,1)) - 5*f(yAD(i-1)) + f(yAD(i-2)));
    yAD(i+1,1) = yAD(i,1) + h/24*(9*g(xAD(i+1,1),yAD(i+1,1)) + 19*g(xAD(i,1),yAD(i,1)) - 5*g(xAD(i-1,1),yAD(i-1)) + g(xAD(i-2,1),yAD(i-2,1)));
end
 
%% Milen
for i=1:4                  %initial value + initial 3 step
    xMI(i,1) = xME(i,1);
    yMI(i,1) = xxME(i,1);  %gam haie avalie az javab modified euler mohasabe shode
end

f = @ (y) y;                         %x dot
g = @ (x,y) (-1*c/m)*y-(k/m)*x;      %x double dot

for i=4:final_step
    %predictor
    xMI(i+1,1) = xMI(i-3,1) + h*4/3*(2*f(yMI(i,1)) - f(yMI(i-1,1)) + 2*f(yMI(i-2)));
    yMI(i+1,1) = yMI(i,1) + h*4/3*(2*g(xMI(i,1),yMI(i,1)) - g(xMI(i-1,1),yMI(i-1,1)) + 2*g(xMI(i-2,1),yMI(i-2)));
    %correcor
    xMI(i+1,1) = xMI(i-1,1) + h/3*(f(yMI(i+1,1)) + 4*f(yMI(i,1)) + f(yMI(i-1,1)));
    yMI(i+1,1) = yMI(i-1,1) + h/3*(g(xMI(i+1,1),yMI(i+1,1)) + 4*g(xMI(i,1),yMI(i,1)) + g(xMI(i-1,1),yMI(i-1)));
end
    
    

%% plot
fplot(analytical_x,[0,ftime]);
hold on
plot(t,xME,'r');
hold on
plot(t,xRK,'g');
hold on
plot(t,xAD,'c');
hold on
plot(t,xMI,'k');
legend('analytical solution','Modified Euler','Runge kutta#4','Adams-Moulton','Milen')
grid on
xlabel('time(sec)');
ylabel('x(m)');   
    
    
    
    


