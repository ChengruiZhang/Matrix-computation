function [err1, err2, err3] = HW5_3(A1, b1, x1, A2, b2, x2)
%% initialization
close all
% A1 = load('.//data_problem3//A1.txt'); b1 = load('.//data_problem3//b1.txt'); x1 = load('.//data_problem3//x1.txt');
% A2 = load('.//data_problem3//A2.txt'); b2 = load('.//data_problem3//b2.txt'); x2 = load('.//data_problem3//x2.txt');
N = 500; eps = 1e-8; x0_1 = rand(10, 1); x0_2 = rand(1000, 1);

%% Q1 
% Jacobi
n = length(x1); [x_J1, ~, num_J1] = iteration(A1, b1, x0_1, N, eps, 0, 1); 
n = length(x2); [x_J2, ~, num_J2] = iteration(A2, b2, x0_2, N, eps, 0, 1); 
% Gauss
n = length(x1); [x_G1, ~, num_G1] = iteration(A1, b1, x0_1, N, eps, 0, 2); 
n = length(x2); [x_G2, ~, num_G2] = iteration(A2, b2, x0_2, N, eps, 0, 2); 

%% Q2
% SOR
% w = 0.05:0.05:1.95; Q2(A1, b1, x1, A2, b2, x2, N, eps, w);

%% Q3
% Jacobi
err1 = zeros(2, max(num_J1, num_J2)); 
for i = 1 : num_J1, err1(1, i) = norm(x_J1(:, i)-x1, 2); end
for i = 1 : num_J2, err1(2, i) = norm(x_J2(:, i)-x2, 2); end
figure, hold on, grid on
plot(0:max(num_J1, num_J2)-1, err1(1, :)); plot(0:max(num_J1, num_J2)-1, err1(2, :))
legend('Data-10', 'Data-1000'); xlabel('iteration\_num'); ylabel('ERROR'); title('Jacobi')
% Gauss
err2 = zeros(2, max(num_G1, num_G2)); 
for i = 1 : num_G1, err2(1, i) = norm(x_G1(:, i)-x1, 2); end
for i = 1 : num_G2, err2(2, i) = norm(x_G2(:, i)-x2, 2); end
figure, hold on, grid on
plot(0:max(num_G1, num_G2)-1, err2(1, :)); plot(0:max(num_G1, num_G2)-1, err2(2, :))
legend('Data-10', 'Data-1000'); xlabel('iteration\_num'); ylabel('ERROR'); title('Gauss')
% SOR, let w = 1
n = length(x1); [x_S1, ~, num_S1] = iteration(A1, b1, x0_1, N, eps, 0.5, 3); 
n = length(x2); [x_S2, ~, num_S2] = iteration(A2, b2, x0_2, N, eps, 0.5, 3); 

err3 = zeros(2, max(num_S1, num_S2)); 
for i = 1 : num_S1, err3(1, i) = norm(x_S1(:, i)-x1, 2); end
for i = 1 : num_S2, err3(2, i) = norm(x_S2(:, i)-x2, 2); end
figure, hold on, grid on
plot(0:max(num_S1, num_S2)-1, err3(1, :)); plot(0:max(num_S1, num_S2)-1, err3(2, :))
legend('Data-10', 'Data-1000'); xlabel('iteration\_num'); ylabel('ERROR'); title('SOR of Q3: w = 0.5')

end

%% Q2
function Q2(A1, b1, x1, A2, b2, x2, N, eps, w)
iter = zeros(2, length(w)); k = 1;
for i = w
    disp(k)
    n = length(x1); [~, ~, iter(1, k)] = iteration(A1, b1, rand(n, 1), N, eps, i, 3); 
    n = length(x2); [~, ~, iter(2, k)] = iteration(A2, b2, rand(n, 1), N, eps, i, 3); 
    k = k + 1;
end
plot(w, iter(1, :)), hold on, grid on
plot(w, iter(2, :))
legend('Data-10', 'Data-1000'); xlabel('w'); ylabel('iteration\_num'); title('SOR of Q2');
end


%% main
function [X, err, k] = iteration(A, b, x0, N, eps, w, t)
x = x0; err = zeros(1, N); X = zeros(length(x0), N+1); X(:, 1) = x0;
for k = 1 : N
    tmp = x;
    switch(t)
        case 1, x = Jacobi(x, A, b);
        case 2, x = Gauss(x, A, b);
        case 3, x = SOR(x, A, b, w);
    end
    err(k) = norm(x - tmp, 2); X(:, k+1) = x;
    if err(k) < eps, break; end
end
if k < N, err(k+1:end) = []; X(:, k+1:end) = []; end
end

function x = Jacobi(x0, A, b)
% x = diag(diag(A)) \ ((-triu(A, 1) - tril(A, -1))*x0 + b);
x = x0;
for i = 1 : length(x0)
    x(i) = b(i);
    for j = 1 : length(x0)
        if i == j, continue; end
        x(i) = x(i) - A(i, j) * x0(j);
    end
    x(i) = x(i) / A(i, i);
end
end

function x = Gauss(x0, A, b)
x = x0;
for i = 1 : length(x0)
    x(i) = b(i);
    try
        for j = 1 : i - 1
            x(i) = x(i) - A(i, j) * x(j);
        end
    catch
    end
    for j = i+1 : length(x0)
        x(i) = x(i) - A(i, j) * x0(j);
    end
    x(i) = x(i) / A(i, i);
end
end

function x = SOR(x0, A, b, w)
tmp = Gauss(x0, A, b);
x = (1 - w) * x0 + w * tmp;
end