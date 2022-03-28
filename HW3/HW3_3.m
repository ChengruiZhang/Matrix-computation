function [q3_4_1, q3_4_2] = HW3_3(eps)
format rat
%% Q1 classical. This is for question 3.1, if u focus on 3.4, u can ignore this part.
n = 4; a = sym(zeros(n, n),'f');
a(:, 1) = sym([1,2,3,4]','f'); a(:, 2) = sym([2,3,4,5]','f'); 
a(:, 4) = sym([3,4,5,6]','f'); a(:, 3) = sym([3,5,7,11]','f');

q3_1 = sym(zeros(n,n),'f'); qt = sym(zeros(n,n),'f'); 
qt(:,1) = a(:, 1); q3_1(:,1) = qt(:,1)/norm(qt(:,1), 2);
for i = 2:3
    for j = 1 : i-1
        qt(:, i) = qt(:, i) + (q3_1(:, j)' * a(:, i)) * q3_1(:, j);
    end
    q3_1(:, i) = a(:, i) - qt(:, i);
    q3_1(:, i) = q3_1(:, i)/norm(q3_1(:, i), 2);
end

%% Q3 classical. This is for question 3.4
% eps = 1e-9;
disp('************  Q3.4-CGS  ***********')
n = 3; a = [1, eps, eps; 1, eps, 0; 1, 0, eps]';
q3_4_1 = zeros(n,n); qt = zeros(n,n); 
qt(:,1) = a(:, 1); q3_4_1(:,1) = qt(:,1)/norm(qt(:,1), 2);
for i = 2:n
    for j = 1 : i-1
        qt(:, i) = qt(:, i) + (q3_4_1(:, j)' * a(:, i)) * q3_4_1(:, j);
    end
    q3_4_1(:, i) = a(:, i) - qt(:, i);
    q3_4_1(:, i) = q3_4_1(:, i)/norm(q3_4_1(:, i), 2);
end
fprintf('The ||Q^TQ-I|| of CGS method at eps = %e is: %e\n', [eps, norm(q3_4_1'*q3_4_1 - eye(3), 'fro')]);
fprintf('The Q1 is:\n'); disp(q3_4_1);
%% Q3 modified
disp('************  Q3.4-MGS  ***********')
v = zeros(n,n); q3_4_2 = zeros(n, n);
for j = 1 : n
    v(:, j)=a(:, j);
end
for j = 1 : n
    q3_4_2(:, j) = v(:, j)/norm(v(:, j), 2); 
    for  k = j + 1 : n
        v(:, k)=v(:, k)-(q3_4_2(:, j)'*v(:, k))*q3_4_2(:, j); 
    end
end
fprintf('The ||Q^TQ-I|| of MGS method at eps = %e is: %e\n', [eps, norm(q3_4_2'*q3_4_2 - eye(3), 'fro')]);
fprintf('The Q2 is:\n'); disp(q3_4_2);
end