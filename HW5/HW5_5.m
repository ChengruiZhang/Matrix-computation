function Cr = HW5_5
close all
load('.\\data_problem4\\data1\\data1.mat');
[D, N] = size(data1); tau = 150; K = min(D, N);
[~, S, ~] = svd(data1); singular_value = diag(S);

%% Q1
s = 0;
for d = 1 : K - 1
    s = s + singular_value(K-d+1)^2;
    if s > tau, break, end
end
d_star = K - d + 1;
plot(1:K, singular_value.^2), xlabel('index'), ylabel('squared singular values of A'), grid on
figure, plot(d_star+1:K, singular_value(d_star + 1:end).^2), xlabel('index'), ylabel('squared singular values of A when tau = 150'), grid on

%% Q2_1
s_all = sum(singular_value.^2); fd = zeros(1, length(singular_value));
for d = 1 : K - 1
    tmp = singular_value(d + 1 : end);
    fd(d) = sum(tmp.^2) / s_all;
end
figure, stairs(1:length(fd), fd), xlabel('index'), ylabel('variance explained ratio'); grid on

%% Q2_2
tau2 = [0.1, 0.05, 0.02, 0.005]; s_all = sum(singular_value.^2); Cr = zeros(1, length(tau2));
for i = 1 : length(tau2)
    for d = K - 1 : -1 : 0
        tmp = singular_value(d + 1 : end);
        s2 = sum(tmp.^2) / s_all;
        if s2 > tau2(i), break, end
    end
    d_star = d + 1;
    Cr(i) = d_star * (D + N + 1) / (D * N);
end

end