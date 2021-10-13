clear all
clc
% Part A
% Solving system of equation using gauss seidel method
% matric A is given in the problem
A = [10,1;1,10];
%b is also given in the problem
b = [11,11];
%initial guess 
X0 = [0.5,0.5]';
%tolerance 
tol=1e-4;
X1 = Seidel(A,b,X0,tol);
disp('The solution is:')
X1

% PART B 
% Solving system of equation using gauss jacobi method
% matric A1 is given in the problem
A1=[4,1,-1;2,7,1;1,-3,12];
%b1 is also given in the problem
b1=[3,19,31];
%initial guess 
X01=[0,0,0]';
%tolerance 
tol=1e-4;
X = Jacobi(A1,b1,X01,tol);
disp('The solution is:')
X

%function for gauss Seidel method
function X = Seidel(A,b,X0,tol)
%input : A(nxn) non singular matrix
%AX = b
%b is a nx1 right hand side vector
%X0 is the inital guess -- Column vector
%tol -- tolerance
%Output -X
n = length(b);
Xnew = zeros(n,1);
itr = 0;
%to start the process
err = 1;
Xold = X0;
%N_maxiter set it to 1000
N_maxiter = 1000;
while err>tol && itr<N_maxiter
    for j=1:n
        if j==1
            Xnew(j) = (b(j)-A(j,2:n)*Xold(2:n))/A(j,j);
        elseif j==n
            Xnew(j) = (b(j)-A(j,1:n-1)*Xnew(1:n-1))/A(j,j);
        else
            Xnew(j) = (b(j)-A(j,1:j-1)*Xnew(1:j-1)-A(j,j+1:n)*Xold(j+1:n))/A(j,j);
        end
    end
    err = max(abs(Xnew-Xold));
    Xold = Xnew;
    itr = itr+1;
end
X =Xnew;
end    

%function for gauss Jacobi method
function X = Jacobi(A,b,X0,tol)
%input : A(nxn) non singular matrix
%AX = b
%b is a nx1 right hand side vector
%X0 is the inital guess -- Column vector
%tol -- tolerance
%Output -X
n = length(b);
Xnew = zeros(n,1);
itr = 0;
%to start the process
error = 1;
Xold = X0;
%N_maxiter set it to 1000
N_maxiter = 1000;
while error>tol && itr<N_maxiter
    for j=1:n
        Xnew(j) = (b(j)-A(j,[1:j-1,j+1:n])*Xold([1:j-1,j+1:n]))/A(j,j);
    end
    error = max(abs(Xnew-Xold));
    Xold = Xnew;
    itr = itr+1;
end
X =Xnew; 
end    
