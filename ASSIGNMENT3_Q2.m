clc
clear all
%code to solve -u"+ p(x)u' + q(x)u = f(x);
%this problem is solved using the method of finite difference method
f = @(x) 0;
p = @(x) -3;
q = @(x) -8;
N = 8;
%initial conditions
x0 = 0;
g0 = 1;
x1 = 1;
g1 = 2;
h = (x1-x0)/N;

X = zeros(1,N);
for i=1:N
    X(1,i) = x0+h*(i);
end
    
[A,b] = finitedifferencemethod(f,p,q,N,h,x0,g0,g1);
sol = Gauss(A,b);
disp('The solution of given boundary value problem');
c = sol';
disp(c);
c(8) = 2;
plot(X,c);

%this returns A,b which is solved using gausian elimination 
function [A,b] =  finitedifferencemethod(f,p,q,N,h,x0,g0,g1)
 A = zeros(N-1,N-1);
 b = zeros(N-1,1);
 P = zeros(N-1,1);
 Q = zeros(N-1,1);
 F= zeros(N-1,1);
 for i=1:N-1
     P(i) = p(x0+i*h);
     Q(i) = q(x0+i*h);
     F(i) = f(x0+i*h); 
 end
 for i=1:N-1
     A(i,i) = 2+Q(i)*h*h;
 end
 for i=1:N-2
     A(i,i+1) = -(1-(h*P(i))/2);
 end
 for i=2:N-1
     A(i,i-1) = -(1+(h*P(i))/2);
 end
 b(1,1) = (1+(P(1)*h)/2)*g0+F(1)*h*h;
 b(N-1,1) = (1-(P(N-1)*h)/2)*g1+F(N-1)*h*h;
 for i=2:N-2
     b(i,1) =  F(i)*h^2;
 end
 end

%function for gauss elimination
function  X = Gauss(A,b)
%To solve system AX =b
%Input: A and b
%Output: X
[m,n] = size(A);
if m~=n 
    disp('A is not a square matrix');
end
% X is a matrix of nx1
X = zeros(n,1);
% L be indentity matrix
for i=1:n
    for j =1:n
        if i==j
            L(i,j)=1;
        else
            L(i,j) =0;
        end
    end
end
%Augmented matrix
Aug = [A b];
%code for forward elimination
for k =1:n-1
    for i = k+1:n
        m = Aug(i,k)/Aug(k,k);
        L(i,k) = m;
        Aug(i,k:n+1) = Aug(i,k:n+1) - m*Aug(k,k:n+1);
    end
end
%using backsubsitution
X = backsub(Aug(1:n,1:n),Aug(1:n,1+n));
end

%code for backsubsitution
function X = backsub(U,b)
n = length(b);
% X is a matrix of nx1
X = zeros(n,1);
%U is the upper triangular matrix
X(n,1) = b(n,1)/U(n,n);
for i = n-1: -1 :1
     X(i,1) = (b(i,1) - U(i,i+1:n)*X(i+1:n))/U(i,i); 
end
end
