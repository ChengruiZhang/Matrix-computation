% input matrix A, miu (non-eigenvalue of A), maximum iteration num and err; 
% and output the (lambda1, q1), (lambda2, q2) by power iteration and 
% inverse iteration respectively.
function [v, d, lambda1, q, lambda2, q2] = HW4_5% (A, miu, max_iter, err)
close all
alpha = 2; beta = 8; % A = sym([0, alpha; beta, 0], 'f');
A = [0, alpha; beta, 0]; max_iter = 200; err = 1e-10; miu = 1;
n = length(A);

%% Q1
disp('**********Q1**********')
[v, d] = eig(A); 
fprintf('The eigenvalue of A is:\n'); disp(d); 
fprintf('The corresponding eigenvector of A is:\n'); disp(v);

%% Q2
disp('**********Q2**********')
q = randn(n, 1) + randn(n, 1)*1i;
z = A * q; q = z / norm(z, 2); 
lambda1 = conj(q)' * A * q;
tmp = lambda1; iteration1 = 1; a0 = real(tmp); b0 = imag(tmp); PI = zeros(max_iter, 1); PI(1) = abs(tmp);

while 1
    z = A * q; q = z / norm(z, 2);
    lambda1 = conj(q)' * A * q; 
    a1 = real(lambda1); b1 = imag(lambda1); 
    PI(iteration1+1) = abs(lambda1);
    if norm(tmp - lambda1) <= err
        fprintf('The iteration number of Power Iteration is: %d\n', iteration1); 
        break; 
    end
    if iteration1 >= max_iter
        fprintf('The iteration number of Power Iteration exceeds the maximum: %d\n', max_iter); 
        break; 
    end
    tmp = lambda1; iteration1 = iteration1 + 1; a0 = real(tmp); b0 = imag(tmp); 
end
plot(1:iteration1+1, PI(1:iteration1+1))
xlabel('iteration'); ylabel('norm(\lambda^k - \lambda^{k-1})'); title('Power Iteration');
grid on

%% Q3
disp('**********Q3**********')
q2 = randn(n, 1) + randn(n, 1)*1i; 
I = eye(n,n); 
z = (A - miu * I) \ q2; q2 = z / norm(z, 2);
lambda2 = q2' * A * q2;
tmp2 = lambda2; iteration2 = 1; 
a0 = real(tmp2); b0 = imag(tmp2); 
II = zeros(max_iter, 1); II(1) = abs(tmp2);

while 1
    z = (A - miu * I) \ q2; q2 = z / norm(z, 2);
    lambda2 = q2' * A * q2; 
    a1 = real(lambda2); b1 = imag(lambda2); 
    II(iteration2+1) = abs(lambda2);
    if norm(tmp2 - lambda2) <= err
        fprintf('The iteration number of Inverse Iteration is: %d\n', iteration2); 
        break; 
    end
    if iteration2 >= max_iter
        fprintf('The iteration number of Inverse Iteration exceeds the maximum: %d\n', max_iter); 
        break; 
    end
    tmp2 = lambda2; 
    iteration2 = iteration2 + 1; 
    a0 = real(tmp2); b0 = imag(tmp2);
end

figure, plot(1:iteration2+1, II(1:iteration2+1))
xlabel('iteration'); ylabel('norm(\lambda^k - \lambda^{k-1})'); title('Inverse Iteration');
grid on

end