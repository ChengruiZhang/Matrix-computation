
%% main func
function [y,x] = HW2_6
y = zeros(4, floor(1000/50)-1);
for i = 100 : 50 : 1000
    disp(i)
    [y(:, floor(i/50)-1),x] = main(i);
end
plot(100:50:1000, y(1,:), 'LineWidth', 1)
hold on, grid on
plot(100:50:1000, y(2,:), 'LineWidth', 1), plot(100:50:1000, y(3,:), 'LineWidth', 1),plot(100:50:1000, y(4,:), 'LineWidth', 1)
xlabel('Size of matrix A'), ylabel('time/s')
legend('Inverse', 'Cramer', 'Gaussian', 'LU\_factorization')
end


function [time, x] = main(n)
time = zeros(4,1); x = zeros(n,4);
A = randi([0,1],n,n)/n; b = randn(n,1); 
A = A'*A + eye(n,n);
tic; x(:,1) = inverse(A, b); time(1) = toc;
tic; x(:,2) = Cramer(A, b); time(2) = toc;
tic; x(:,3) = Gaussian(A, b); time(3) = toc;
tic; x(:,4) = LU_self(A, b); time(4) = toc;
end


%% inv func
function x = inverse(A, b)
n = length(A); A2 = [A, eye(n, n)];
for j = 1 : n
    [~, max_loc] = max(A2(j:end,j)); % find the location of maximum value of column j
    max_loc = max_loc + j - 1;
    for i = 1 : n
        if i == max_loc, continue; end
        A2(i,:) = A2(i,:) - A2(i,j)/A2(max_loc, j) * A2(max_loc, :);
    end
    A2(max_loc, :) = A2(max_loc, :) / A2(max_loc, j);
    tmp = A2(j, :); A2(j,:) = A2(max_loc,:); A2(max_loc, :) = tmp;
end
x = A2(:, n+1:end) * b;
end


%% cramer func
function x = Cramer(A, b)
n = length(A); x = zeros(n,1);
for i = 1 : n
    tmp = A; tmp(:,i) = b;
    x(i) = det(tmp)/det(A);
end
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
    tmp = A2(j, :); A2(j,:) = A2(max_loc, :); A2(max_loc, :) = tmp;
end
%% back substitution
b = A2(:,n+1); A = A2(:,1:n);
for i = n:-1:1
    if i == n, x(n, 1) = b(n); continue; end
    x(i,1) = b(i) - A(i,i+1:n)*x(i+1:n,1);
end
end


%% LU func
function x = LU_self(A, b)
n = length(A); x = zeros(n, 1); A2 = A; P = zeros(n,n); P1 = (1:n)'; L = eye(n,n);
for j = 1 : n
    [~, max_loc] = max(A2(j:end,j)); % find the location of maximum value of column j
    max_loc = max_loc + j - 1;
    tmp1 = A2(j, :); A2(j, :) = A2(max_loc, :); A2(max_loc, :) = tmp1;
    tmp2 = P1(j, :); P1(j, :) = P1(max_loc, :); P1(max_loc, :) = tmp2;
    for i = j + 1 : n
        a = A2(i,j)/A2(j,j); A2(i, j:end) =  A2(i, j:end) - a * A2(j, j:end); A2(i, j) = a;
    end
end
for i = 2:n
    L(i,1:i-1) = A2(i, 1:i-1); A2(i, 1:i-1) = 0;
end
for i = 1 : n
    P(i, P1(i)) = 1;
end

%% forward substitution
U = A2; y = zeros(n,1); b = P * b;
for i = 1: n
    if i == 1, y(1, 1) = b(1); continue; end
    y(i,1) = b(i) - L(i, 1:i-1)*y(1:i-1,1);
end

%% backward substitution
x = Gaussian(U, y);

end
