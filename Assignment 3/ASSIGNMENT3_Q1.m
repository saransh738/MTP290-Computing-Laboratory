clc
clear all
% function is given f(x) = (y/x)- (y/x)^2;
f = @(x,y) (y/x)- (y/x)^2;
%step size
h = 0.01;
%No of steps it is not given in question but we can calculate it
%We know that x increase by h each time ,we are given max value of x i.e. 2
% so x + N*h = 2
%puting all values, we get 1 + N*0.01 = 2
%By solving above equation we get N = 100
N = 100;
%we are given initial condition as y(1) = 1
xinitial = 1;
yinitial = 1;
%function of runge kutta of order 4
[a,b]= Rungekutta(f,h,N,xinitial,yinitial);
%ploting the function 
plot(a,b);
%code for rungekutta
function [a,b] = Rungekutta(f,h,N,xinitial,yinitial)
%input : function f , h stepsize , N , initialconditions
%Output : two vectors a,b which store value of x and y after each iteration
a(1) = xinitial;
b(1) = yinitial;
for j = 1:N
    %vector a stores value of x after each iteration 
    %x will increase by h every time
    a(j+1) = a(j) + h;
    %using 4 eqautions of rungekutta
    %k1 = h*f(x(n), y(n));
    %k2 = h*f(x(n) + h/2 , y(n) + K1/2);
    %k3 = h*f(x(n) + h/2 , y(n) + K2/2);
    %k4 = h*f(x(n) + h, y(n) + K3);
    ama1 = h*f(a(j), b(j));
    ama2 = h*f(a(j) + h/2 , b(j) + ama1/2);
    ama3 = h*f(a(j) + h/2 , b(j) + ama2/2);
    ama4 = h*f(a(j) + h, b(j) + ama3);
    % y(n+1) = y(n) +(1*(k1 + 2*k2 + 2*k3 + k4))/6;
    b(j+1) = b(j) +(1*(ama1 + 2*ama2 + 2*ama3 + ama4))/6;
end
a
b
end
