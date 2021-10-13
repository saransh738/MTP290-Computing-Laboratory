%This is the program to return solution of a triangular system
%take A as triangular matrix A = [2,3,0,0;6,3,9,0;0,2,5,2;0,0,4,3]
clear all
clc
A = [2,3,0,0;6,3,9,0;0,2,5,2;0,0,4,3];
%take D as D = [21,69,34,22]'
D = [21,69,34,22]'; 
n = length(A);

%C matrix is of n rows and 1 columns
%it will stores all values of ci's
%given cn=0
C = zeros(n,1); 
C(n,1) =0;
%By observing Triangular matrix we get ci's
for i=1:n-1
    C(i,1) = A(i,i+1);
end

%A1 matrix is of n rows and 1 columns
%it will stores all values of ai's
%given a1=0
A1 = zeros(n,1);
A1(1,1) = 0;
%By observing Triangular matrix we get ai's
for i=2:n
    A1(i,1) = A(i,i-1);
end

%B matrix is of n rows and 1 columns
%it will stores all values of bi's
B = zeros(n,1);    
%By observing Triangular matrix we get to know that bi's are diagonal
%entries of tridiagonal matrix
for i=1:n
    B(i,1) = A(i,i);
end

%C1 matrix is of n rows and 1 columns
%it will stores all values of ci'
C1 = zeros(n,1);
%given c1' = c1/b1;
C1(1,1) = C(1,1)/B(1,1);
%using the given expression ci' = ci/(bi-ci-1'*ai) we get ci's
for i =2:n-1
     C1(i,1)= C(i,1)/(B(i,1) - A1(i,1)*C1(i-1,1));
end

%D1 matrix is of n rows and 1 columns
%it will stores all values of di'
D1 = zeros(n,1);
%given d1' = d1/b1;
D1(1,1) = D(1,1)/B(1,1);
%using the given expression di' = (di-di-1'*ai)/(bi-ci-1'*ai) we get di's
for i =2:n
     D1(i,1)= (D(i,1)-A1(i,1)*D1(i-1,1))/(B(i,1) - A1(i,1)*C1(i-1,1));
end

% X stores solution of tridiagonal system with n unknowns
% X matrix is of n rows and 1 columns
X = zeros(n,1);
%given xn = dn'
X(n,1) = D1(n,1);
%using the recursion given we can find X 
% xi = di' - ci'*xi+1 
for i = n-1 :-1 :1
     X(i,1)= -X(i+1,1)*C1(i,1) + D1(i,1);
end
X


