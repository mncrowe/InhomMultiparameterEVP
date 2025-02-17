% Plots a few Dipole solutions. Let's use 1-layer solutions here for
% simplicity. These plots relate to the N-Layer QG vortex problem which
% motivated this study into weird eigenvalue problems.

addpath('..\..\N-Layer_QG\functions\') % path to QGDipoles (MATLAB version)
addpath('..\..\N-Layer_QG\plotting\')  % path to QGDipoles (MATLAB version)
addpath('..\functions\')

% Create terms in eigenvalue problem:

M = 8;
A = zeros(M);
B = zeros(M);
c = [1/4; zeros(M-1, 1)];
d = (-1).^(0:M-1)';

for j = 0:M-1
    for k = 0:M-1
        A(j+1, k+1) = JJ_int(@(x) A_func(x, 1, 1), j, k);
        B(j+1, k+1) = JJ_int(@(x) B_func(x, 1, 1), j, k);
    end
end

% Solve eignevalue problem, (A-K^2*B)*a = (mu+K^2)*c s.t. d'*a = 0:

[K2, a2] = EVP(A, B, [c c], d, 1);

% Let's first find the K value closest to 4:

[~, i] = min(abs(sqrt(K2) - 4));
K = sqrt(K2(i)); a = a2(:, i);

% We can now calculate the associated vortex solution and plot it:

Nx = 512; Ny = 512;
Lx = 10; Ly = 10;
[x, y, k, l, grid] = create_domain([Nx Ny], [Lx Ly]);

[~, q] = CalcPsi(a, 1, 1, 1, 1, grid);
Plot_2D(q, x, y, [-1.5 1.5], 1)

% Let's now plot a solution for a different K, say K close to 7:

[~, i] = min(abs(sqrt(K2) - 7));
K = sqrt(K2(i)); a = a2(:, i);

[~, q] = CalcPsi(a, 1, 1, 1, 1, grid);
Plot_2D(q, x, y, [-1.5 1.5], 1)

% And the next K value (K close to 10):

[~, i] = min(abs(sqrt(K2) - 10));
K = sqrt(K2(i)); a = a2(:, i);

[~, q] = CalcPsi(a, 1, 1, 1, 1, grid);
Plot_2D(q, x, y, [-1.5 1.5], 1)

% Finally let's look at K close to 13.5:

[~, i] = min(abs(sqrt(K2) - 13.5));
K = sqrt(K2(i)); a = a2(:, i);

[~, q] = CalcPsi(a, 1, 1, 1, 1, grid);
Plot_2D(q, x, y, [-1.5 1.5], 1)

% We can see that the radial mode is increasing each time. Hence the value
% of K is related to the internal mode number (or wavenumber) of the
% vortex. With multiple layers, we can have a different modal structure in
% each layer.
