% This is for part a and part b
% for part a it will returns L and U
% for part b this code returns Solution set of Ax = b
clear all
clc
%taken matrix  A as A = [4,1,-1;5,1,2;6,1,1]
A = [4,1,-1;5,1,2;6,1,1];
%And take b as b = [-2;4;6]
b = [-2,4,6]';
%U is the upper triangular matrix
%L is the lower triangular matrix
[m,n] = size(A);
%if A is not square matrix

if m~=n 
    disp('A is not a square matrix');
end

%x is matrix with n rows and 1 columns and all zeros
X = zeros(n,1);
%assign all diagonal entry of L as 1
for i=1:n
    for j =1:n
        if i==j
          L(i,j)=1;
        else
          L(i,j) =0;
        end
    end
end
%assign U as A (initially)
U = A;
%after forward elimination we get U as upper triangular matrix
% code for forward elimination
for k =1:n-1
    for i = k+1:n
        % m is the multiplier factor
        m = U(i,k)/U(k,k);
        L(i,k) = m; % for lower triangular matrix
        U(i,k:n) = U(i,k:n) - m*U(k,k:n);
        b(i,1) = b(i,1) - m*b(k,1);
    end
end
U
L
X = backsub(U,b);
X

% code for backward subsitution
function X = backsub(U,b)
n= length(b);
%x is matrix with n rows and 1 columns and all zeros
X = zeros(n,1);
X(n,1) = b(n,1)/U(n,n);
for i = n-1: -1 :1
     X(i,1) = (b(i,1) - U(i,i+1:n)*X(i+1:n))/U(i,i); 
end
end


