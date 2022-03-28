% just use [x1, x2] = HW3_2([],[]) to get the results of gradient decent
% method and normal equation respectively
function [x1, x2, A, y] = HW3_2(A, y)
%% preparation
m = 1e3; n = 210;
if ~ischar(A)
    A = zeros(m, n); y = zeros(m, 1);
    fid = fopen('data.txt'); tline = fgetl(fid); i = 1;
    while ischar(tline)
        tmp = regexp(tline, '\s+', 'split'); tmp = str2num(char(tmp(1:n))); A(i, :) = tmp; tline = fgetl(fid);i = i + 1;
    end; fclose(fid);

    fid = fopen('label.txt'); tline = fgetl(fid); i = 1;
    while ischar(tline)
        y(i, 1) = str2double(tline); tline = fgetl(fid); i = i + 1;
    end; fclose(fid);
end
%% Q1
disp('************  Q1  ***********')
gamma = 1e-5; x1 = zeros(n, 1); err = 1e4; iteration = 1e2; i = 1; % tmp = Inf;
while err > gamma
    delta = 2 * (y' - (A * x1)') * A; % or delta = 2*A'*(A*x-y); (mn + m) + m + mn
    x1 = x1 + (gamma * delta)'; % n
    err = (norm(y-A*x1, 2))^2; i = i + 1; % disp(err); 
%     if err>=tmp, break; end
%     tmp = err;
    if i>iteration, break; end
end
% fprintf('The gradient norm and '); 
fprintf('The gradient norm of normal equation method is: %e\n', norm(2 * (y' - (A * x1)') * A,2));
fprintf('The loss f(x) of gradient descent method is: %e\n', err);

%% Q2 
disp('************  Q2  ***********')
G = m_chol(A'*A); % lower tri-angular n^2*m + 1/3*n^3
z = zeros(n,1); b = A'*y; z(1, 1) = b(1)/G(1,1);
for i = 2: n % forward n^2
    z(i,1) = (b(i, 1) - G(i, 1:i-1)*z(1:i-1,1))/G(i, i);
end
x2 = Gaussian(G', z); % backward n^2
err2 = (norm(y-A*x2, 2))^2;
fprintf('The gradient norm of normal equation method is: %e\n', norm(2 * (y' - (A * x2)') * A,2));
fprintf('The loss f(x) of normal equation method is: %e\n', err2);
% [q,r] = qr(A, 0);
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
A2 = [A,b]; n = length(A); x = zeros(n, 1);
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
