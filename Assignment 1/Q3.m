%explanation : Use the following method to find two sequences uk and vk such 
%that all the values f(uk),k = 0,1,2,.... have one sign and
%all the values f(vk) have the opposite sign.
%a function f can be assumned anything 
%i used f as f = x^4 -16*x^3 + 89*x^2 - 194*x + 120 (as used in q2)
f = @(x) x^4 -16*x^3 + 89*x^2 - 194*x + 120 ;
%lets assumne uk is the sequence such that f(uk)>0 for all values of uk 
%so we start with u0 = -1 since f(-1) = 420 > 0
%uk is positive sequence
u0 = -1; 
%lets assumne vk is the sequence such that f(vk)>0 for all values of vk 
%so we start with v0 = 1.01 since f(1.01)<0
%vk is negative sequence
v0 = 1.01;
%xlist is the array which stores all the values of sequence uk
xlist(1) = u0; 
%xlist is the array which stores all the values of sequence vk
ylist(1) = v0;
%stopping criteria is k = 10
%since k  starts from 0 to 10 so intotal there will be 11 terms i.e k+1
%terms
k=10;
for i =2:k+1
    % wk = (uk*f(vk) - vk*f(uk))/(f(vk) - f(uk))
    w0 = (u0*f(v0) - v0*f(u0))/(f(v0) - f(u0));
    % if f(wk) and f(uk) of same sign then uk+1 = wk and vk+1 = vk 
    if f(w0)*f(u0)>0
       u0 = w0;
       xlist(i) = u0 ;
       ylist(i) = v0 ;
    else
    % if f(wk) and f(uk) of opposite sign then uk+1 = uk and vk+1 = wk 
       v0 = w0;
       xlist(i) = u0 ;
       ylist(i) = v0;
    end
end
xlist
ylist

