function HW5_4
close all
%% Q1
d = 10; load('.\\data_problem4\\data1\\data1.mat');
Q1(data1, d);

%% Q2
load('.\\data_problem4\\all_data.mat'); % data_train, data_test, Y_label_test, Y_label_train
[~, ~, ~] = Q2(data_train, 10, 1);
[~, ~, ~, ~] = Q4_2_3(data_test, 10, 1);
[~, ~, ~, ~] = Q4_2_3(data_test, 2, 1);

%% Q3
d = 2 : 8; err = zeros(1, length(d)); Num = size(data_test, 2);
for i = 1 : length(d)
    [~, ~, ~, Proj_train] = Q4_2_3(data_train, d(i), 0);
    [~, ~, ~, Proj_test] = Q4_2_3(data_test, d(i), 0);
    for j = 1 : Num
        pos = K_NN(Proj_train, Proj_test(:, j));
        if Y_label_train(pos) ~= Y_label_test(j), err(i) = err(i) + 1; end
    end
end
err = err / Num * 100; disp(err)
end

function Q1(data1, d)
meanface = mean(data1, 2); subplot(1,5,1), imshow(reshape(meanface, [48,42])); title('meanface')
X = data1 - meanface; [U, S, ~] = svd(X);
sig_value = diag(S); [~, pos] = sort(sig_value); % m from small to large
eigenface = U(:, pos(end:-1:end-d+1)); eigenface = mapminmax(eigenface', 0, 1)';
subplot(1,5,2), imshow(reshape(eigenface(:,1), [48,42])), title('1-eigenface')
subplot(1,5,3), imshow(reshape(eigenface(:,2), [48,42])), title('2-eigenface')
subplot(1,5,4), imshow(reshape(eigenface(:,3), [48,42])), title('3-eigenface')
subplot(1,5,5), imshow(reshape(eigenface(:,10), [48,42])), title('10-eigenface')
end

function [meanface, U, Y] = Q2(X, d, isplot) % data1: train; data2: test
meanface = mean(X, 2); [U, S, ~] = svd(X - meanface); Ud = U(:, 1:d); m = diag(S);
% sig_value = diag(S); [m, pos] = sort(sig_value); % m from small to large
% eigenface = U(:, pos(end:-1:end-d+1)); 
eigenface = U(:, 1:d); eigenface = mapminmax(eigenface', 0, 1)';
% Q4_2_1
if isplot
    figure, subplot(1,5,1), imshow(reshape(meanface, [48,42])); title('meanface')
    subplot(1,5,2), imshow(reshape(eigenface(:,1), [48,42])), title('1-eigenface')
    subplot(1,5,3), imshow(reshape(eigenface(:,2), [48,42])), title('2-eigenface')
    subplot(1,5,4), imshow(reshape(eigenface(:,3), [48,42])), title('3-eigenface')
    subplot(1,5,5), imshow(reshape(eigenface(:,10), [48,42])), title('10-eigenface')
end
% Q4_2_2
if isplot
    figure, plot(1:length(m), m); xlabel('index'); ylabel('singular value'); grid on 
end
Ud = U(:, 1:d); Y = Ud' * (X - meanface);
end

function [meanface, U, Y, Proj] = Q4_2_3(X, d, isplot)
% Q4_2_3
meanface = mean(X, 2); [U, ~, ~] = svd(X - meanface);
% sig_value_test = diag(S_test); [m, pos] = sort(sig_value_test); Ud_test = U_test(:, pos(end:-1:end-d+1));
Ud_test = U(:, 1:d); Y = Ud_test' * (X - meanface); 
Proj = meanface + Ud_test * Y; % Proj = mapminmax(Proj', 0, 1)';
if isplot
    figure, subplot(1,4,1), imshow(reshape(Proj(:,2), [48,42])), title('1-individual')
    subplot(1,4,2), imshow(reshape(Proj(:,6), [48,42])), title('2-individual')
    subplot(1,4,3), imshow(reshape(Proj(:,14), [48,42])), title('3-individual')
    subplot(1,4,4), imshow(reshape(Proj(:,20), [48,42])), title('4-individual')
end
end

% data: D*N, x: D*1; return the position of closest data
function pos = K_NN(data, x)
err = inf; pos = 1;
[~, N] = size(data);
for i = 1 : N
    tmp = norm(x - data(:, i));
    if tmp < err
        err = tmp; pos = i;
    end
end
end