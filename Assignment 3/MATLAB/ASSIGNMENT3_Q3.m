clc
clear all
% Exaplanation : The given function can be integrated easily
% After integrating we get 2*tan(inv(4)) which is approximately 2.65
% Now we can compare results in part A we get 0.4706 which is very very less than 2.65
% which implies that by using trapezoidal rule we get ans less than actual
% ans in this example
% percentage error ((2.65-0.47)/2.65) = 82% less than actual value 

% In part B we get ans as 5.4902 which is almost double of actual ans
% which implies that by using simpson rule we get ans more than actual
% ans in this example
% percentage error ((2.65-5.49)/2.65) = 107% more than actual value 

% In part C we get ans as 2.6511 by composite trapezoidal rule and 2.6953 by using composite simpson's rule
% which implies that by using composite trapezoidal rule or simpson rule give ans very close to actual
% ans in this example
% percentage error by  composite trapezoidal rule ((2.65-2.65)/2.65) = 0% i.e. equal to actual value 
% percentage error by  composite simpson's rule ((2.65-2.69)/2.65) = 1.5% more than actual value 


% In part D we get ans as 1.2632 which is almost half of actual ans
% which implies that by using two point Gauss-Legendre quadrature we get ans less than actual
% ans in this example 
% percentage error by  composite trapezoidal rule ((2.65-1.26)/2.65) = 52% less than actual value 


% Thus we can conlcude that composite trapezoidal and simpson rule give ans
% almost equal to actual ans while two point Gauss-Legendre quadrature give
% ans which is almost half of actual ans
% Thus composite trapezoidal and simpson rule are more accurate and precise
% than two point Gauss-Legendre quadrature 

%Numerical intergaration
%function f is given f(x) = 1/(1+x^2)
f =@(x) 1/(1+x^2);
%a is lower limit of intergration
a = -4;
%b is upper limit of intergration
b = 4;

%PART A CALCULATION USING TRAPEZOIDAL RULE
h =(b-a);
%using the formula Area = (h*(f(a)+f(b)))/2;
%trapez represents the value of intergral using Trapezoidal rule 
trapez = (h*(f(a)+f(b)))/2;
disp('The value of integral using trapezoidal rule');
disp(trapez);

%PART B CALCULATION USING SIMPSON's RULE
h1 = (b-a)/2;
c = (a+b)/2;
%using the formula Area = (h1*(f(a)+4*f(c)+f(b)))/3;
%simpson represents the value of intergral using Simpson's rule 
simpson = (h1*(f(a)+4*f(c)+f(b)))/3;
disp('The value of integral using simpson rule');
disp(simpson);

%PART C CALCULATION USING COMPOSITE TRAPEZOIDAL RULE
%N = no. of subinterval
N = 10;
compositetrapezoidal = trapezoidal(f,a,b,N);
disp('The value of integral using composite trapezoidal rule');
disp(compositetrapezoidal);
%PART C CALCULATION USING COMPOSITE Simpson's RULE
compositesimpson = compositsimpson(f,a,b,N);
disp('The value of integral using composite simpson rule');
disp(compositesimpson);

%PART D CALCULATION USING Gausian two point quadrature
w1=1;
w2=1;
x1 = ((b-a)/2)*(-1/sqrt(3)) + (b+a)/2;
x2 = ((b-a)/2)*(1/sqrt(3)) + (b+a)/2;
Twopointgausianquadrature = (b-a)/2.0*(w1*f(x1)+w2*f(x2));
disp('The value of integral using two point gausian quadrature ');
disp(Twopointgausianquadrature);

%code for composite trapezoidal rule
function A = trapezoidal(f,a,b,N)
%input : function f , lower limit a and upper limit b
%N = no. of subinterval
%Output : Area A
h2 = (b-a)/N;
A = 0;
for i=1:N-1
    x = a + i*h2;
    A = A + feval(f,x);
end
A = h2*(feval(f,a)+feval(f,b))/2+h2*A;
end

%code for composite simpson's rule
function A1 = compositsimpson(f,a,b,N)
%input : function f , lower limit a and upper limit b
%N = no. of subinterval
%For simpson's 1/3 rd rule N=2k
%Output : Area A
h3 = (b-a)/(N);
s1=0;
s2=0;
A1 =0;
for i=1:N/2
    x1 = a+h3*(2*i-1);
    s1 = s1+f(x1);
end
for j=1:(N/2)-1
    x2 = a+h3*2*j;
    s2 = s2+f(x2);
end
A1 = (h3/3) * (f(a) + f(b) + 4*s1 + 2*s2);
end




