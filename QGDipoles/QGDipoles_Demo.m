% Calculate some dipole solutions as examples for solving the Two-Parameter
% Inhomogeneous (Quadratic) Eigenvalue Problem.

% Note: the K values here are the square roots of the eigenvalues,
% therefore we calculate the eigenvalues, λ, and determine K as sqrt(λ).

addpath('..\..\N-Layer_QG\functions\') % path to QGDipoles (MATLAB version)
addpath('..\functions\')


% Let's build a 2-layer example:

N = 2; M = 6;
mu = [0 1];
lambda = 1./[1 1];

A = zeros(N*M); B = zeros(N*M);
Dn = @(n) kron(eye(M), diag(1:N == n));
c = [ones(N,1); zeros((M-1)*N, 1)]/4;
dn = @(n) kron((-1).^(0:M-1)', (1:N == n)');

for j = 0:M-1
    for k = 0:M-1
        A(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) A_func(x,lambda,mu),j,k);
        B(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) B_func(x,lambda,mu),j,k);
    end
end

B1 = Dn(1)*B; B2 = Dn(2)*B;
c1 = Dn(1)*c; c2 = Dn(2)*c; c0 = mu(1)*c1+mu(2)*c2;
d1 = dn(1); d2 = dn(2);

% Let's solve using root finding method:

Bcell = cell(2,1); Bcell{1} = B1; Bcell{2} = B2;
[K, a_j] = EVP_optim(A, Bcell, [c0 c1 c2], ...
                     [d1 d2], [4 5].^2, zeros(2*M,1));
K = sqrt(K); a_j = reshape(a_j,[N M])';

% Let's now try using 2-parameter EVP method. Let's convert the system to a
% 2-parameter system of 2 equations and then find the two K values using 
% the Delta matrix method. We know from the previous method and insights 
% from the corresponding physical problem that the solutions we want lie 
% between 3 and 4, are real, and are not infinite. We'll get rid of all
% other eigenvalues.

I = eye(N*M); O = zeros(N*M);
A0_1 = (d1'*(A\c0))*A;
A0_2 = (d2'*(A\c0))*A;
A1_1 = (d1'*(A\c1))*A - (d1'*(A\c0))*B1 + (c0*d1')*(A\B1);
A1_2 = (d2'*(A\c1))*A - (d2'*(A\c0))*B1 + (c0*d2')*(A\B1);
A2_1 = (d1'*(A\c2))*A - (d1'*(A\c0))*B2 + (c0*d1')*(A\B2);
A2_2 = (d2'*(A\c2))*A - (d2'*(A\c0))*B2 + (c0*d2')*(A\B2);
A11_1 = - (d1'*(A\c1))*B1 + (c1*d1')*(A\B1);
A11_2 = - (d2'*(A\c1))*B1 + (c1*d2')*(A\B1);
A22_1 = - (d1'*(A\c2))*B2 + (c2*d1')*(A\B2);
A22_2 = - (d2'*(A\c2))*B2 + (c2*d2')*(A\B2);
A12_1 = - (d1'*(A\c2))*B1 + (c2*d1')*(A\B1);
A12_2 = - (d2'*(A\c2))*B1 + (c2*d2')*(A\B1);
A21_1 = - (d1'*(A\c1))*B2 + (c1*d1')*(A\B2);
A21_2 = - (d2'*(A\c1))*B2 + (c1*d2')*(A\B2);

A1 = [A0_1 O O; O I O; O O I];
A2 = [A0_2 O O; O I O; O O I];
B11 = [A1_1 A11_1 A12_1; -I O I; O O I];
B21 = [A1_2 A11_2 A12_2; -I O I; O O I];
B12 = [A2_1 A21_1 A22_1; O -I O; -I -I O];
B22 = [A2_2 A21_2 A22_2; O -I O; -I -I O];

As = zeros(3*N*M, 3*N*M, 2);
Bs = zeros(3*N*M, 3*N*M, 2, 2);
As(:, :, 1) = A1; As(:, :, 2) = A2;
Bs(:, :, 1, 1) = B11; Bs(:, :, 1, 2) = B12;
Bs(:, :, 2, 1) = B21; Bs(:, :, 2, 2) = B22;

D = Delta_Matrices(As, Bs);

K1s = eig(D(:, :, 2), -D(:, :, 1));
K2s = eig(D(:, :, 3), -D(:, :, 1));

K1 = sqrt(K1s((~isinf(K1s) & (K1s < 16) & (K1s > 9) & (imag(K1s) == 0))));
K2 = sqrt(K2s((~isinf(K2s) & (K2s < 16) & (K2s > 9) & (imag(K2s) == 0))));

% A complication here is that there's several K_1 and K_2 values close to
% what we expect. This is because we're getting all K_1 values that are
% close to 4, regardless of the corresponding K_2 value. The ones we want
% are where both K_1 and K_2 are close to 4, which we expect to be a unique
% pair. However we can't identify this pair without pairing up all the
% eigenvalues first which requires more work.

% The values for K that we found using root finding (EVP_optim.m) are
% K ≈ (3.8002, 3.9499). We note that these values are among the K1 and K2
% that we found by solving the Delta matrix system, however without
% doing any more calculations, we're unable to pair them up.

% Let's now get the K2 corresponding to our K1:

i = find(abs(K1 - 3.8002) < 1e-4); K1 = K1(i);

K22 = eig(A1 + K1^2*B11, -B12);
K22 = sqrt(K22((~isinf(K22) & (K22 < 16) & (K22 > 9) & (imag(K22) == 0))));

j = find(abs(K22 - 3.9499) < 1e-4); K2 = K22(j);

% We can now get the eigenvector from our original system and verify it is
% perpendicular to the d_i.

% The conclusion is that the root finding approach is much better here. For
% larger systems (larger N or M) the linear algebra method is no longer
% feasible whereas the root finiding method scales well, even for systems
% of size O(100-1000).
