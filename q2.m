%function f = x^4 -16*x^3 + 89*x^2 - 194*x + 120
f = @(x) x^4 -16*x^3 + 89*x^2 - 194*x + 120 ;
%a is initial guess
a = 1.5;
tol = 1e-10;
%a function is called newton rapson method 
[root,itr,xlist] = Newton1(f,a,tol);
%it displays the root of the function f
fprintf('the root = %d \n',root);
% the error en = abs(xn - xn-1)
e(1) = abs(xlist(1)-a);
for i = 2:length(xlist)-1
    e(i) = abs(xlist(i) - xlist(i-1));
end
%q is array to store rate of convergence
q = zeros(length(e)-2 , 1); 
%for each q rate of convergence is stored in tha arry q
%order of convergence is numerically calculated 
%using the formula abs(log(e(n+1)/e(n)) / log(e(n)/e(n-1)))
for n=2:length(e)-1
   q(n-1) = abs((log(e(n+1)/e(n))) / (log(e(n)/e(n-1))));
end
for i=1:length(q)
  fprintf('At itr = %d ,the order of convergence is = %4.2f \n',i+2,q(i));
end
%this prof is for given problem method
%inputf : f(function),a(initial guess), and the tolerance
%output:  root of f ,  number of iterations, 
%          and xlist(it stores value of
%roots at each iteration)
function [sol,no_itr,xlist] = Newton1(f,a,tol)
%xnew = xold - f(xold)/((f(aold+f(aold)) - f(aold))/f(aold))
aold =a;
itr=0;
%first element of xlist will be a
xlist(1) =a;
% it displays initial guess
fprintf(' \n Initial app., a = %7.6f  \n',a);
%stopping condition are used in while loop 
%I have used condition on iter basically it must be less than max_iter 
%else if function diverges then we will get into infinite loop
%so just to not get into infinite loop i used this condition  
max_iter = 200;
while abs(f(aold)) > tol && itr<max_iter
    anew = aold - f(aold)/((f(aold+f(aold)) - f(aold))/f(aold)) ;
    aold = anew;
    itr = itr+1;
    xlist(itr+1) = anew;
    % it displays the iteration number , aold values and f(aold)
    fprintf(' At itr = %d , a = %7.6f , f(a)= %7.6f \n',itr,aold,f(aold));
    
end
sol = anew;
no_itr = itr;
end