function D = Delta_Matrices(A, B)
% Calculates the matrices required to reduce a n x n system of N x N 
% eigenvalue problems to n single eigenvalue problems.
%
% Input:
%
% A: n N x N matrices, shape [N N n]
% B: n^2 N x N matrices, shape [N N n n]
%
% Output:
%
% D: n+1 N^n x N^n matrices, D_0, D_1, ... D_n, shape [N^n, N^n, n+1]
%
% -------------------------------------------------------------------------
% Note: system of equations:
%
% (A_1 - K_1 B_11 - K_2 B_12 - ... - K_n B_1n) x_1 = 0
% (A_2 - K_1 B_21 - K_2 B_22 - ... - K_n B_2n) x_2 = 0
%  ...       ...      ...             ....       ....
% (A_n - K_1 B_n1 - K_2 B_n2 - ... - K_n B_nn) x_n = 0
%
% is converted to:
%
% (D_0 - K_i D_i) x = 0
%
% where D_0 is the determinant of the B_ij system and D_i is the i^th
% Cramer determinant obtained by replacing the i^th column of the B_ij
% system with the A_i matrices. Here x is the Kronecker product of the x_i.
%
% -------------------------------------------------------------------------
% Note: this method does not scale well for large systems as calculating
% the n+1 determinants (Cramer's rule) required to solve the system
% requires summing the n! terms in each determinant. For matrices of
% scalars this can be more easily done by matrix factorisation, however
% here we're dealing with matrices of matrices. A faster method could
% likely be done by splitting the system into n x n matrices of scalars and
% using factorisation on these.

s = size(A);

if length(s) > 2; n = s(3); else; n = 1; end    % size of system
N = s(1);                                       % size of matrices

p = perms(1:n);                                 % calcuate all permutations

D = zeros(N^n, N^n, n+1);                       % assign empty Delta matrices
M = zeros(N, N, n);

C = @(i, j, k) (j==k)*A(:, :, i) + (j~=k)*B(:, :, i, j);

for i = 1:factorial(n)

    eps = Levi_Civita(p(i, :));                 % sign of permutation

    for j = 1:n
        M(:, :, j) = B(:, :, j, p(i, j));       % calculate all matrices in permutation
    end

    D(:, :, 1) = D(:, :, 1) + eps * multi_kron(M);

    for k = 1:n

            for j = 1:n
                M(:, :, j) = C(j, p(i, j), k);  % calculate all matrices in permutation
            end

        D(:, :, k+1) = D(:, :, k+1) + eps * multi_kron(M);

    end

end

end