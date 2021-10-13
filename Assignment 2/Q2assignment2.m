% observations : .
% 1 : The given matrix Hn will be symmetric matrix i.e A' = A 
% (where A' is the transpose of the matrix)
% And if matrix is symmetric then row sum and column sum will be same i.e
% norm1 or column sum and norm_infinity or row sum will be same 
% which is evident form the results we obtained we can see in the output
% that norm1 and norm_infinity are same for n=3,4,5,6

% 2 : We know by the definition of conditional number 
% A matrix with large conditional number is said to be ill-conditioned and
% A matrix with small conditional number is said to be well-conditioned
% In the output we can see that the conditional number is too too much 
% That Means Hn is ill conditioned
% and  also as we increase n the conditional number increses ( can seen in the ouput)

% 3 : as conditional number is product of norm(A)*norm(inv(A)) as if
% conditional number is large that means (A) tends to singular matrix
% AS we know if matrix is singular then conditional number is infinity 
% In our case also the conditional number increases as n increases so as n
% increases we can say that Hn approaches to singular matrix

% 4 : We know if conditional number is high that means 
% for a small change in the inputs there is a large change in the answer 
% and in such cases it becomes difficult to find correct solution of
% equations. In  our case also conditional increases as n increases

% 5 : We can also observe that Euclidean norm is smaller than row or column
% norm 

clear all
clc
%to calculate conditional number
%We know conditional number is norm(A)*norm(inv(A))
%conditional number is displayed for the matrix by taking n=3,4,5,6
H3 = [1,1/2,1/3;1/2,1/3,1/4;1/3,1/4,1/5];
B3 = inv(H3);
%conditional number using norm_infinity or  row norm
conditional_numberinf_3 = norm_inf(H3)*norm_inf(B3);
disp(conditional_numberinf_3);
%conditional number using norm1 or column norm
conditional_number1_3 = norm1(H3)*norm1(B3);
disp(conditional_number1_3);
%conditional number using norm2 or Euclidean norm
conditional_number2_3 = norm2(H3)*norm2(B3);
disp(conditional_number2_3);

H4 = [1,1/2,1/3,1/4;1/2,1/3,1/4,1/5;1/3,1/4,1/5,1/6;1/4,1/5,1/6,1/7];
B4 = inv(H4);
%conditional number using norm_infinity or  row norm
conditional_numberinf_4 = norm_inf(H4)*norm_inf(B4);
disp(conditional_numberinf_4);
%conditional number using norm1 or column norm
conditional_number1_4 = norm1(H4)*norm1(B4);
disp(conditional_number1_4);
%conditional number using norm2 or Euclidean norm
conditional_number2_4 = norm2(H4)*norm2(B4);
disp(conditional_number2_4);

H5 = [1,1/2,1/3,1/4,1/5;1/2,1/3,1/4,1/5,1/6;1/3,1/4,1/5,1/6,1/7;1/4,1/5,1/6,1/7,1/8;1/5,1/6,1/7,1/8,1/9];
B5 = inv(H5);
%conditional number using norm_infinity or  row norm
conditional_numberinf_5 = norm_inf(H5)*norm_inf(B5);
disp(conditional_numberinf_5);
%conditional number using norm1 or column norm
conditional_number1_5 = norm1(H5)*norm1(B5);
disp(conditional_number1_5);
%conditional number using norm2 or Euclidean norm
conditional_number2_5 = norm2(H5)*norm2(B5);
disp(conditional_number2_5);

H6 = [1,1/2,1/3,1/4,1/5,1/6;1/2,1/3,1/4,1/5,1/6,1/7;1/3,1/4,1/5,1/6,1/7,1/8;1/4,1/5,1/6,1/7,1/8,1/9;1/5,1/6,1/7,1/8,1/9,1/10;1/6,1/7,1/8,1/9,1/10,1/11];
B6 = inv(H6);
%conditional number using norm_infinity or  row norm
conditional_numberinf_6 = norm_inf(H6)*norm_inf(B6);
disp(conditional_numberinf_6);
%conditional number using norm1 or column norm
conditional_number1_6 = norm1(H6)*norm1(B6);
disp(conditional_number1_6);
%conditional number using norm2 or Euclidean norm
conditional_number2_6 = norm2(H6)*norm2(B6);
disp(conditional_number2_6);


%this prog is to calculate norm1 or column norm
%Input : A (the matrix)
%Output : a = norm1
function a = norm1(A)
[n,m] = size(A); 
%A matrix having 1 row and m columns
%it stores the column sum in matrix X
X = zeros(1,m);
for i=1:m  
    sum = 0;
    for j=1:n
        sum = sum + abs(A(j,i));
    end
    X(1,i) = sum;
end
%to find maximum column sum
a = X(1,1);
for i = 2:m
    if X(1,i)>a
        a = X(1,i);
    end
end
end

%this prog is to calculate norm2 or Euclidean norm
%Input : A (the matrix)
%Output : a = norm2
function a = norm2(A)
[n,m] = size(A);
sum=0;
for i=1:n
    for j=1:m
        sum = sum+(A(i,j))^2;
    end
end
a = sum^(1/2);
end

%this prog is to calculate norm_infinity or row norm
%Input : A (the matrix)
%Output : a = norm_infinity
function a = norm_inf(A)
[n,m] = size(A);
%A matrix having n row and 1 columns
%it stores the row sum in matrix X
X = zeros(n,1);
for i=1:n
    sum = 0;
    for j=1:m
        sum =sum + abs(A(i,j));
    end
    X(i,1) = sum;
end
%to find maximum row sum
a = X(1,1);
for i = 2:n
    if X(i,1)>a
        a = X(i,1);
    end
end
end