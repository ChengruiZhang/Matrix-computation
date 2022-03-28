function [x3_1, x3_2, x4_1, x4_2, k1, k2] = HW3_4
format rat
%% Q1
A = sym([1,2,3;2,3,5;3,4,7;4,5,11], 'f'); b = sym([1,1.5,3,6]','f'); delta_b = sym([0.1, 0, 0, 0]','f');
[q, r] = qr(A, 0);
x1 = inv(r)*q'*b;
err = norm(A * x1 - b,2)^2;
disp(err)

%% Q3
disp('************  Q3  ***********')
b = [1,1.5,3,6]'; delta_b = [0.1, 0, 0, 0]';
A3 = [1,2,2;2,2,2;3,3,3;1,1,0];  A_p3 = (A3'*A3)\(A3');
x3_1 = A_p3 * b; x3_2 = A_p3 * (b+delta_b);
k1 = norm(A3, 2)*norm(A_p3, 2);
fprintf('The results of Q3 for b and (b + delta_b) are:\n'); disp([x3_1, x3_2]);

%% Q4
disp('************  Q4  ***********')
A4 = [1,1,1;1,2,4;1,3,9;1,4,16]; 
A_p4 = (A4'*A4)\(A4');
x4_1 = A_p4 * b; x4_2 = A_p4 * (b+delta_b);
k2 = norm(A4, 2)*norm(A_p4, 2);
fprintf('The results of Q4 for b and (b + delta_b) are:\n'); disp([x4_1, x4_2]);

end

%% cholesky decomposition
function [A] = m_chol(A)
n = length(A);
for i = 1:n
    for k = 1:i-1
       A(i,i) = A(i,i) - A(k,i)*A(k,i); 
    end
    if A(i,i) <= 0, disp('A is not positive definite.'); exit; end
    A(i,i) = sqrt(A(i,i));
    for j = i+1:n
        for k = 1:i-1, A(i,j)= A(i,j) - A(k,i)*A(k,j); end
        A(i,j) = A(i,j)/A(i,i);
    end
end
A = triu(A); % get upper triangular part of matrix
A = A';
end


%% Gauss func
function x = Gaussian(A, b)
A2 = [A,b]; n = length(A); x = sym(zeros(n, 1),'f');
for j = 1 : n
    [~, max_loc] = max(A2(j:end,j)); % find the location of maximum value of column j
    max_loc = max_loc + j - 1;
    for i = j : n % only j:n rows
        if i == max_loc, continue; end
        A2(i,:) = A2(i,:) - A2(i,j)/A2(max_loc, j) * A2(max_loc, :);
    end
    A2(max_loc, :) = A2(max_loc, :) / A2(max_loc, j);
    tmp = A2(j, :); A2(j,:) = A2(max_loc,:); A2(max_loc, :) = tmp;
end
%% back substitution
b = A2(:,n+1); A = A2(:,1:n);
for i = n:-1:1
    if i == n, x(n, 1) = b(n); continue; end
    x(i,1) = b(i) - A(i,i+1:n)*x(i+1:n,1);
end
end
