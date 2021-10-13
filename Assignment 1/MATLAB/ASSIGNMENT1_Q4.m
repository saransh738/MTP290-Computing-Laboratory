clear all
clc
%function f = tanh(x)
f = @(x) tanh(x);
% df is the differentiation of the function f(x)
df = matlabFunction(diff(sym(f)));
tol = 1e-10;
a = -10; 
b = 15;
s = 0.1;
%a function is called bisection_Newton
[sol, no_iter] =bisection_Newton(f, df, a, b, tol, s);
%it displays the root of the function f
disp(sol);

%this prog is for new(bisection + newton) method
%Input : f(function),df(derivative of f) , [a,b], s,  and  the tolerance
%Output : root of f , and number of iterations

function [sol, no_iter] = bisection_Newton(f, df, a, b, tol, s)
%feval calulates value of function f at a and b
fa = feval(f,a);
fb = feval(f,b);

%if both fa and fb has same sign then no root exists
    if fa*fb > 0
        fprintf("No root lies in the given interval [a,b]");
    end
    
    %if both fa and fb has opposite sign then we use bisection method till
    %the length of interval is greater than s*(b-a) 
    c = (a + b)/2;
    fc = feval(f,c);
    iter = 0;
    % condition to switch from bisection to newton
    swi = s*(b - a);  
    while (b - a) > swi
        % if fa and fc have same sign
        if fa*fc > 0   
            a = c;
            fa = fc;
        else
            b = c;
        end
        c = (a + b)/2;
        fc = feval(f,c);
        iter = iter + 1;
    end
    %when the length of interval becomes less than s*(b-a)
    [sol, no_iter] = Newton(f, df, c, tol);
    no_iter = iter + no_iter;
end

%this prog is Newton Rapson method
%inputf : f(function),df(derivative of f),a(initial guess), and the tolerance
%output:  root of f , and number of iterations
function [sol,no_itr] = Newton(f,df,a,tol)
%xnew = xold - f(xold)/f'(xold)
aold =a;
itr=0;
%stopping condition are used in while loop 
%I have used condition on iter basically it must be less than max_iter 
%else if function diverges then we will get into infinite loop
%so just to not get into infinite loop i used this condition  
max_iter = 200;
while abs(f(aold)) > tol && itr<max_iter
    anew = aold - f(aold)/df(aold);
    aold = anew;
    itr = itr+1;
end
sol = anew;
no_itr = itr;
end
