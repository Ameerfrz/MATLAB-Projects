clc;clear;
syms f x;
f = -0.0981*x^2+x;
ff = -0.0981*x^2;

result = zeros(30,6);
for i =1:30
    result(i,1) = i;
end

a = 10;  %primary guess
b = 11;

%% Bisection
aa = a;
bb = b;
cc = 0; 
dd = 1;
i=0;

while abs(cc - dd)>0.00001 
    i = i+1;
    dd = cc;
    f1 = subs(f,aa);
    f2 = subs(f,bb);
    cc = (aa+bb)/2;
    f3 = subs(f,cc);
    if f3*f1 < 0
        bb = cc;
    elseif f3*f1 > 0
        aa = cc;
    else
        break;
    end
    result(i,2) = cc;
end
        

%% linear interpolotion
% aa = a;
% bb = b;
% cc = 0; 
% dd = 1;
% i=0;
% 
% while abs(cc - dd)>0.00001 
%     i = i+1;
%     dd = cc;
%     f1 = double(subs(f,aa));
%     f2 = double(subs(f,bb));
%     cc = aa - f1*(bb-aa)/(f2-f1);
%     f3 = subs(f,cc);
%     if f3*f1 < 0
%         bb = cc;
%     elseif f3*f1 > 0
%         aa = cc;
%     else
%         break;
%     end
%     result(i,3) = cc;
% end
% 
% %% Newton Raphson
% aa = a;
% bb = b;
% i=0;
% diff_f = diff(f);
% 
% while abs(bb - aa)>0.00001 
%     i = i+1;
%     bb = aa;
%     f1 = double(subs(f,aa));
%     diff_f1 = double(subs(diff_f,aa));
%     aa = aa - (f1/diff_f1);  
%     result(i,4) = aa;
%     if i>30
%         break
%     end
% end
% 
%  
%  %% modified newton
% % aa = b;
% % bb = b;
% % cc = a;
% % i=0;
% % 
% % while abs(bb - cc)>0.01 
% %     i = i+1;
% %     aa = bb;
% %     bb = cc;
% %     f1 = double(subs(f,aa));
% %     f2 = double(subs(f,bb));
% %     cc = bb - (bb-aa)/(f2-f1);
% %     if i > 29
% %         break;
% %     end
% %     result(i,5) = cc;
% % end
% %  %javab cc ast
% %   
% 
% %% g(x)=x
% x2 = (a+b)/2;
% x1 = x2+1;
% i = 0;
% while abs(x2-x1)>0.01
%     i= i+1;
%     x1 = x2;
%     x2 = double(subs(ff,x2));
%     if i>30
%         break
%     end
%     
%  result(i,5) = x2;
% end
% 
%  %% Muller
% 
%  x1 = max(a,b);
%  x2 = min(a,b);
%  x0 = (x1+x2)/2;
%  xr = x0 ;
%  s  = 1;
%  i = 0;
%  
% while abs(s)>0.00001
% 
%      i=i+1;
%  f0 = double(subs(f,x,x0));
%  f1 = double(subs(f,x,x1));
%  f2 = double(subs(f,x,x2));
%  
%  h1 = abs(x1-x0);
%  h2 = abs(x2-x0);
%  gama = h2/h1;
%  
%  aa = (gama*f1-f0*(1+gama)+f2)/(gama*h1^2*(1+gama));
%  bb = (f1-f0-aa*h1^2)/h1;
%  cc = f0;
%  
%  xr1 = (-2*cc)/(bb+(bb^2-4*aa*cc)^0.5);
%  xr2 = (-2*cc)/(bb-(bb^2-4*aa*cc)^0.5);
%  s = min(abs(xr1),abs(xr2));
%  if abs(xr1)==s
%      ss = sign(xr1);
%  else
%      ss = sign(xr2);
%  end
%  xr = x0 + s*ss;
%  
%  if xr>x0
%      x2 = x0;
%      x0 = xr;
%  else
%      x0 = xr;
%  end
%  if i>30
%      break;
%  end
%  result(i,6) = xr;
%  
%  end
%         
