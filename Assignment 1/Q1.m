+%explanation : consider f(x) = a -1/x 
%Let x = 1/a be an approximate solution of the equation. At the point
%(x0,f(x0)), we draw a tangent line to the graph of y = f(x). Let
%x1 be the point at which the tangent line intersects the x-axis. It should be an
%improved approximation of the root a.
%now we connect (x1,0) and (x0,f(x0)) and find equation of line we get
%                   x1 = x0(2-a*x0)
%generalising it we get iteration formula x(n+1) = x(n)*(2-a*x(n))
% now let we assume rn = 1 - a*xn 
% next we have to see when will this method converges
% we have error e(n) = (1/a) - xn = rn/a
% the error en converges to zero as n ~ oo if and only if rn converges
% to zero. but rn converges to zero if and only if lr0l < 1
% so we get       -1 < 1- ax0 < 1
%                  0 < x0 < 2/a
% In order that xn converge to 1/a, it is necessary and sufficient that x0 be chosen
% to satisfy  such that 0<x0<2/a
% x0 = initial guess thus 0<x0<2/a
clc
clear all
%Assume function to be f(x) = a-1/x for some a not equal to 0 
%let a to be 10 

f = @(x) 10 - 1/x;
%initial guess
a = 0.01;
tol = 10^(-10);
%the function is called here 
[root,itr] = iterative(f,a,tol);
fprintf('the root is =');
disp(root);


%this prog is for iterative method
%inputf : f(function),a(initial guess), and the tolerance
%output:  root of f , and number of iterations
function [sol,no_itr] = iterative(f,a,tol)
%xnew = xold*(2 - a*xold) i have taken a to be 10
aold =a;
itr=0;
% it displays initial guess
fprintf(' \n Initial app., a = %7.6f  \n',a);

%stopping condition are used in while loop 
%I have used condition on iter basically it must be less than max_iter 
%else if function diverges then we will get into infinite loop
%so just to not get into infinite loop i used this condition
max_iter = 200;
while abs(10*aold-1) > abs(tol*aold) && itr<max_iter
    anew = aold*(2-10*aold);
    aold = anew;
    itr = itr+1;
     % it displays the iteration number , aold values and f(aold)
    fprintf(' At itr = %d , a = %7.6f , f(a)= %7.6f \n',itr,aold,f(aold));
end
sol = anew;
no_itr = itr;
end