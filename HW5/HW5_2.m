function [sol1, sol2] = HW5_2
close all
%% Q1
A = [0, 1, 0, 0;0 0 2 0; 0 0 0 3; 1/6000 0 0 0];
sol1 = Q1(A); % sol1: cell, {lambda, {U, S, V}}

%% Q2
B = [3, 2; 1, 4; 0.5, 0.5];
sol2 = Q2(B); % sol2: cell {U, S, V}

%% Q3 
Q3(30)

end

%% Q1
function sol1 = Q1(A)
lambda = eig(A); [U, S, V] = svd(A); sol1 = {lambda, {U,S,V}};
end


%% Q2
function sol2 = Q2(A)
ATA = A' * A; [V, D] = eig(ATA);
[m, n] = size(A); S = zeros(m, n);
for i = 1 : min(m, n), S(i, i) = sqrt(D(i, i)); end
[Q, R] = qr(A*V);
[G1, B] = eliminate(R); % QR = QGB, B is diagonal
G2 = mksame(S, B); % S = GB
U = Q * G1 * G2; sol2 = {U, S, V};
end

% find G to eliminate the upper-triangular matrix A to diagonal matrix G*A
function [G, A] = eliminate(A)
[m, n] = size(A); G = eye(m, m);
for i = 2 : min(m, n)
    tmp = eye(m, m);
    for j = 1 : i - 1
        tmp(j, i) = A(j, i)/A(i, i);
        A(j, i) = 0;
    end
    G = G*tmp;
end
end

% find G to make the diagonal matrix A and B be the same, A = G * B
function G = mksame(A, B)
[m, n] = size(A); G = eye(m, m);
for i = 1 : min(m, n)
    if A(i, i) ~= B(i, i), G(i,i) = A(i, i)/B(i, i);end
end
end


%% Q3
function Q3(N)
s1 = zeros(1, N); s2 = zeros(1, N);
for m = 1 : N
    for i = 1 : m, A = triu(ones(m, m), 1) + eye(m, m)*0.1; end % create A
    s1(m) = min(svd(A));
    sol2 = Q2(A); tmp = diag(sol2{2}); 
    s2(m) = abs(min(tmp));
end
semilogy(1:30, s1), hold on, grid on
semilogy(1:30, s2)
legend('SVD', 'Algorithm 1')
xlabel('m'); ylabel('The smallest singular value'); title('Comparison of Algorithm 1 and SVD');
end