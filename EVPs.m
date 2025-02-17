% solve some weird EVPs

addpath("functions\")

% Generalised problem: (A - KB)a = 0, typically solved by premultiplying by
% A⁻¹ or B⁻¹.

% See 'functions/Gen_EVP.m'.

rng(2)

A = magic(3); B = randi(4, 3, 3);

[a, E] = eig(A, B);
K = sum(E);

disp('Generalised Problem:'); disp(' ')
disp('K = '); disp(K'); disp('a = '); disp(a)

% Generalised Inhomogeneous problem: (A - KB)a = c, d'a = 0, can
% premultiply by B⁻¹ and mutiply result by cd' to get [cd'B⁻¹A]a =
% (d'B⁻¹c)c hence c = [cd'B⁻¹A]a/(d'B⁻¹c) and we have Generalised problem:
% (A* - KB)a = 0 where A* = A - [cd'B⁻¹A]/(d'B⁻¹c). This allows us to find
% K using the generalised problem. However we need to invert original
% problem for a as length is now arbitrary. Can also premultiply by A⁻¹ if
% B is singular, or less well conditioned.

% Note: A*(A⁻¹c) = 0 so A* has zero determinant and a zero eigenvector, so
% there's at most N-1 solutions which work. The Nth 'solution' has an a
% which is not perpendicular to d.

% See 'functions/Inhom_EVP.m'.

rng(2)

A = randi(4, 3, 3); B = randi(4, 3, 3); c = [1; 2; 0]; d = [1; 0; 1];

As = A - (c*d')*(B\A)/(d'*(B\c));
K = eig(As, B)';
a = zeros(length(A));

for i = 1:length(K)
    a(:, i) = (A - K(i)*B)\c;
end

disp('Generalised Inhomogeneous Problem:'); disp(' ')
disp('K = '); disp(K'); disp('a = '); disp(a)

% Generalised Inhomogeneous Quadratic Problem: (A - KB)a = c₀ + Kc₁, d'a =
% 0. Again, we premultiply by B⁻¹ (or A⁻¹) and repeat the same steps. This
% time we get: (A - [(c₀ + Kc₁)d'B⁻¹A]/(d'B⁻¹(c₀ + Kc₁)) - KB)a = 0, which
% has a K in the denominator. We should multiply through to get:

% ([d'B⁻¹(c₀ + Kc₁)]A - [(c₀ + Kc₁)d'B⁻¹A] - [d'B⁻¹(c₀ + Kc₁)]KB)a = 0,

% which is the quadratic problem:

% (A₀ + KA₁ + K²A₂)a = 0,

% where:    A₀ = [d'B⁻¹c₀]A - [c₀d']B⁻¹A, 
%           A₁ = [d'B⁻¹c₁]A - [c₁d']B⁻¹A - [d'B⁻¹c₀]B,
%           A₂ = -[d'B⁻¹c₁]B.

% This quadratic problem can be written as:

% [A₀ O][ a] - K[-A₁ -A₂][ a] = [0]
% [O  I][Ka]    [ I   O ][Ka]   [0],

% which is a 2N x 2N Generalised EVP with 2N eigenvalues. The addition of a
% K in the RHS means that the Kernel is different for each K value so we 
% can often find N eigenvalues corresponding to N independent vectors. The
% rest will contain a 0 (as A₀ singular since A₀A⁻¹c₀ = 0) and a load of 
% values where:

% d'B⁻¹(c₀ + Kc₁) = 0 => K = -(d'B⁻¹c₀)/(d'B⁻¹c₁)

% so the denominator is 0 and we get division by 0 in the derivation.
% If c₀ and c₁ are proportional, we get the previous case with only N-1
% solutions. The reason we can get N (I think) is that the extra K provides
% a degree of freedom on the right which allows us to avoid the issue
% above, possibly via a div by zero in the denominator at the same time?
% e.g. as if c were zero in case 1 where we now we can get N solutions?

% See 'functions/EVP_K_imhom.m'.

rng(2)

A = randi(4, 3, 3); B = randi(4, 3, 3); O = zeros(3); I = eye(3);
c0 = [1; 2; 0]; c1 = [0; 1; -2]; d = [1; 0; 1];

A0 = (d'*(B\c0))*A - (c0*d')*(B\A);
A1 = (d'*(B\c1))*A - (c1*d')*(B\A) - (d'*(B\c0))*B;
A2 = -(d'*(B\c1))*B;

K = eig([A0 O; O I],[-A1 -A2; I O]);

a = zeros(length(A), length(A0));

for i = 1:length(K)
    a(:, i) = (A - K(i)*B)\(c0 + K(i)*c1);
end

p = abs(d'*a) < 1e-10;

disp('Generalised Inhomogeneous Quadratic Problem:'); disp(' ')
disp('K = '); disp(K); disp('a = '); disp(a)
disp('Satisfy system?'); disp(p)

% Two-Parameter (Generalised) Inhomogeneous Problem:
% (A - K₁B₁ - K₂B₂)a = c, d₁'a = 0, d₂'a = 0
% We now need two conditions as we have two extra degrees of freedom via
% the two eigenvalues. From here we can do the same steps as the
% 1-parameter case, premultiply by cdᵢ'A⁻¹ to get:

% -[cdᵢ']A⁻¹(K₁B₁ + K₂B₂) = (dᵢ'A⁻¹c)c, for i = 1,2,

% hence:

% c = [cdᵢ']A⁻¹(K₁B₁ + K₂B₂)/(dᵢ'A⁻¹c), for i = 1,2,

% so:

% (A - K₁B₁ᵢ - K₂B₂ᵢ)a = 0, for i = 1,2,

% where:

% B₁ᵢ = B₁ - [cdᵢ']A⁻¹B₁/(dᵢ'A⁻¹c), for i = 1,2,
% B₂ᵢ = B₂ - [cdᵢ']A⁻¹B₂/(dᵢ'A⁻¹c), for i = 1,2.

% This is a two-parameter eigenvalue system analogous to simultaneous
% equations but for matrices instead of numbers:

% (A - K₁B₁₁ - K₂B₂₁)a = 0,
% (A - K₁B₁₂ - K₂B₂₂)a = 0.

% Can convert this to a one parameter system using a tensor product
% (Kronecker product) to get:

% (D₁ - K₁D₀)x = 0,
% (D₂ - K₂D₀)x = 0,

% where:

% D₀ =  B₁₁ ⊗ B₂₂ - B₂₁ ⊗ B₁₂,
% D₁ =    A ⊗ B₂₂ - B₂₁ ⊗ A  ,
% D₂ =  B₁₁ ⊗ A   -   A ⊗ B₁₂.

% This can be visualised as intersecting curves in 2D space. We vary K₁ and
% calculate N values of K₂ = K₂(K₁), giving N curves in 2D space.
% Intersections correspond to solutions of the system. This is somewhat 
% complicated by complex branches, whenever two lines merge or one line
% splits that's due to complex eigenvalues with the same real part. We're
% only plotting for real K₁ and the real part of K₂. However, here all
% intersections happen on real branches so all solutions are real.

n = 4;

A = eye(n);
B11 = (diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B12 = (diag(n:-1:1,0) + 2*diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B21 = (diag(n:-1:1,0) - diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B22 = (-diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 + 2*diag(n-1:-1:1,-1)/2);

K1 = linspace(0, 2, 501);

K21 = zeros(length(A),length(K1));
K22 = zeros(length(A),length(K1));

for i = 1:length(K1)
    K21(:, i) = sort(eig(A - K1(i)*B11, B21),'ComparisonMethod','real');
    K22(:, i) = sort(eig(A - K1(i)*B12, B22),'ComparisonMethod','real');
end

warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart')
plot(K1, K21, 'r', K1, K22, 'b', 'LineWidth', 2);
xlim([0 2]); ylim([-1 1]);
grid; xlabel('K_1'); ylabel('K_2')

D0 = kron(B11, B22) - kron(B21, B12);
D1 = kron(A, B22) - kron(B21, A);
D2 = kron(B11, A) - kron(A, B12);

K1 = eig(D1, D0);
K2 = eig(D2, D0);

disp('Two-Parameter Inhomogeneous Problem:'); disp(' ')
disp('K_1 = '); disp(K1); disp('K_2 = '); disp(K2)
disp('But which K_1 goes with which K_2? (See plot).')

% We can also use the general method, coded in DeltaMatrices.m, for
% building the one parameter system.

As = zeros(n, n, 2);
As(:, :, 1) = A; As(:, :, 2) = A;

Bs = zeros(n, n, 2, 2);
Bs(:, :, 1, 1) = B11; Bs(:, :, 2, 1) = B12;
Bs(:, :, 1, 2) = B21; Bs(:, :, 2, 2) = B22;

D = Delta_Matrices(As, Bs);

K1s = eig(D(:, :, 2), D(:, :, 1));
K2s = eig(D(:, :, 3), D(:, :, 1));

disp(' '); disp('And using the general method we get the same:'); disp(' ')
disp('K_1 = '); disp(K1s); disp('K_2 = '); disp(K2s)

% Two-Parameter Inhomogeneous Quadratic Problem:
% (A - K₁B₁ - K₂B₂)a = c₀ + K₁c₁ + K₂c₂, d₁'a = 0, d₂'a = 0
% We now need two conditions as we have two extra degrees of freedom via
% the two eigenvalues. From here we can do the same steps as before, lets
% premultiply by (c₀ + K₁c₁ + K₂c₂)dᵢ'A⁻¹ and rearrange to get an
% expression for (c₀ + K₁c₁ + K₂c₂) which we can substitute into our
% equation and rearrange to get a two-parameter homogeneous quadratic
% problem:

% (A₀ᵢ + K₁A₁ᵢ + K₂A₂ᵢ + K₁²A₁₁ᵢ + K₁K₂[A₁₂ᵢ + A₂₁ᵢ] + K₂²A₂₂ᵢ)a = 0

% where:

% A₀ᵢ  = (dᵢ'A⁻¹c₀)A
% A₁ᵢ  = (dᵢ'A⁻¹c₁)A - (dᵢ'A⁻¹c₀)B₁ + [c₀dᵢ']A⁻¹B₁
% A₂ᵢ  = (dᵢ'A⁻¹c₂)A - (dᵢ'A⁻¹c₀)B₂ + [c₀dᵢ']A⁻¹B₂
% A₁₁ᵢ = - (dᵢ'A⁻¹c₁)B₁ + [c₁dᵢ']A⁻¹B₁
% A₂₂ᵢ = - (dᵢ'A⁻¹c₂)B₂ + [c₂dᵢ']A⁻¹B₂
% A₁₂ᵢ = - (dᵢ'A⁻¹c₂)B₁ + [c₂dᵢ']A⁻¹B₁
% A₂₁ᵢ = - (dᵢ'A⁻¹c₁)B₂ + [c₁dᵢ']A⁻¹B₂

% which may be rewritten as:

% ( [A₀ᵢ O O]      [A₁ᵢ A₁₁ᵢ A₁₂ᵢ]      [A₂ᵢ A₂₁ᵢ A₂₂ᵢ] ) [  a]
% ( [ O  I O] + K₁ [-I   O    I  ] + K₂ [ O   -I   O  ] ) [K₁a] = 0
% ( [ O  O I]      [ O   O    I  ]      [-I   -I   O  ] ) [K₂a]

% which is a system of 2 3N x 3N eigenvalue problems for two parameters. We
% can solve this by using the usual Kronecker approach. The resulting 
% system is of size (3N)^2 x (2N)^2. This is just about doable. For higher 
% numbers of parameters, things become difficult as matrices start getting
% really big. Compared with the usual Canonical form, we've added the I to
% (2, 3) and (3, 3) in the second matrix and the -I to (2, 2) and (2, 3) in
% the third matrix. This corresponds to K₁K₂a - K₂K₁a = 0 and removes zero
% rows ensuring that the matrix is non-singular (if the first rows are well
% behaved) or at least better conditioned.

% Three-Parameter Problem: Let's skip the conversion step from the 
% inhomogeneous problem to the system of quadratic problems to the system
% of larger linear problems and just look at the system of linear problems.
% In general this is:

% (A₁ - K₁B₁₁ - K₂B₁₂ - K₃B₁₃)a₁ = 0,
% (A₂ - K₁B₂₁ - K₂B₂₂ - K₃B₂₃)a₂ = 0,
% (A₃ - K₁B₃₁ - K₂B₃₂ - K₃B₃₃)a₃ = 0,

% which can be converted to three single parameter systems:

% (Δ₁ - K₁Δ₀)x = 0,
% (Δ₂ - K₂Δ₀)x = 0,
% (Δ₃ - K₃Δ₀)x = 0,

% where x = a₁ ⊗ a₂ ⊗ a₃ and:

%      | B₁₁ B₁₂ B₁₃ |
% Δ₀ = | B₂₁ B₂₂ B₂₃ |
%      | B₃₁ B₃₂ B₃₃ |⊗

% and Δᵢ is obtained by replacing column i in Δ₀ with (A₁, A₂, A₃)'. This
% quantity is calculated as per usual determinants except each
% multiplication is replaced by the Kronecker product ⊗. In the
% two-parameter case this reduces to the result above. Note that the i,j
% indices are flipped here relative to the example above. In the case of
% 1x1 matrices this reduces to solving a matrix problem, Mx = y, and the Δ
% matrices represent determinants and Gaussian elimination (Cramer's rule),
% i.e. inversion of the system via cofactors; equivalent to finding M⁻¹.

% See 'functions/DeltaMatrices.m'.

n = 3;

A1 = eye(n);
A2 = 2*eye(n);
A3 = -eye(n);

B11 = (diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B12 = (diag(n:-1:1,0) + 2*diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B13 = (diag(n:-1:1,0) + diag(n-1:-1:1,1)/4 + diag(n-1:-1:1,-1)/4);
B21 = (diag(n:-1:1,0) - diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B22 = (-diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 + 2*diag(n-1:-1:1,-1)/2);
B23 = (diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 - 2*diag(n-1:-1:1,-1)/2);
B31 = (2*diag(n:-1:1,0) + diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B32 = (diag(n:-1:1,0)/2 + diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);
B33 = (diag(n:-1:1,0)/4 + diag(n-1:-1:1,1)/2 + diag(n-1:-1:1,-1)/2);

As = zeros(n, n, 3);
As(:, :, 1) = A1; As(:, :, 2) = A2; As(:, :, 3) = A3;

Bs = zeros(n, n, 3, 3);
Bs(:, :, 1, 1) = B11; Bs(:, :, 2, 1) = B21; Bs(:, :, 3, 1) = B31;
Bs(:, :, 1, 2) = B12; Bs(:, :, 2, 2) = B22; Bs(:, :, 3, 2) = B32;
Bs(:, :, 1, 3) = B13; Bs(:, :, 2, 3) = B23; Bs(:, :, 3, 3) = B33;

D = Delta_Matrices(As, Bs);

K1s = eig(D(:, :, 2), D(:, :, 1));
K2s = eig(D(:, :, 3), D(:, :, 1));
K3s = eig(D(:, :, 4), D(:, :, 1));

disp('Three-Parameter Homogeneous Problem:'); disp(' ')
disp('K_1 = '); disp(K1s);
disp('K_2 = '); disp(K2s);
disp('K_3 = '); disp(K3s);
disp(['We have loads of eigenvalues ... ' ...
    'but which K_1 goes with which K_2 and which K_3?'])
disp('¯\_(ツ)_/¯')

% How do we determine which (K₁, K₂, K₃) values go together? Look at
% eigenvectors and match based on common eigenvalues. Or given the
% eigenvector a₁ we can form a simple linear algebra problem:

% (A₁a₁) - K₁(B₁₁a₁) - K₂(B₁₂a₁) - K₃(B₁₃a₁) = 0,

% The bracketed terms are vectors so we can take the first 2 (linearly
% independent) rows of this system and solve as a 2 x 2 eigenvalue problem
% for K₂ and K₃. This doesn't always work well if there's small errors in
% the eigenvalues or eigenvectors. We also need to find a₁ which requires
% unwrapping the kronecker product.

% n-Parameter Inhomogeneous Quadratic Problem:
%(A - Σ KⱼBⱼ)a = c₀ + Σ Kⱼcⱼ, dᵢ'a = 0, where i, j ∈ {1, 2, ..., n}
%
% General procedure:
%
% 1. Apply n conditions of the form dᵢ'a = 0 to get a system of n quadratic
%    homogeneous problems with n²+n+1 terms in each. This is an n-parameter
%    quadratic homogeneous eigenvalue problem.
%
% 2. Put each quadratic problem into 'Canonical Form', i.e. convert to a
%    linear problem of size (n+1)N x (n+1)N. We now have an n-parameter
%    (linear) homogeneous eigenvalue problem.
%
% 3. Convert this system to n 1-parameter problems using a tensor product
%    over the n spaces and forming the Δᵢ matrices as outlined above. The 
%    resulting system is of size [(n+1)N]ⁿ x [(n+1)N]ⁿ. Each of these
%    1-parameter problems may be solved independently to calculate the Kⱼ
%    values.
%
% 4. However, we don't know how to group our Kⱼ values as solutions of the
%    form K = (K₁, K₂, ...). Try to match based on eigenvectors. Or, solve
%    recursively: find K₁ then solve the n-1 parameter system for K₂, and
%    then solve the n-2 system etc.

% However, this isn't really a very useful way of going about this. If we
% have a 5 layer vortex model with 10 coefficients in each layer (this work
% was motivated by an n-layer QG vortex solution method) then we have n = 5
% and N = 50. Therefore, our final problem is of size 2.43e12 x 2.43e12,
% which is unfeasibly large for this solution method.

% Is there a better method?

% If you care about nice linear algebra ... maybe, there's some recent work
% on subspace methods but these don't appear to have (yet) yielded fast,
% scalable numerical methods.

% If you don't care about linear algebra ... yes, you just need a
% reasonable initial guess and can use gradient based root finding.

% Let's start by applying the conditions somehow. Geometrically, they say
% that the vector a must lie in a space perpendicular to all the dᵢ. We
% assume that the set of n vectors dᵢ span Rⁿ, which they must otherwise 
% we don't have enough conditions to find unique solutions. This space can
% be extended to Rᴺ by adding N-n new vectors eⱼ. The set {dᵢ, eⱼ} now
% forms a basis of Rᴺ and the eⱼ can be calculated using the Gram-Schmidt
% process. We can now write:

% a = Σ αⱼeⱼ, where j ∈ {1, 2, ..., N-n}.

% Adding the n Kᵢ values gives the N values we want to determine:

% x = (Kᵢ, αⱼ), for i ∈ {1, 2, .., n}, j ∈ {1, 2, ..., N-n}.

% We have N equations (as original system has N rows) so we can determine
% these N values as roots of the function:

% f(x) = (A - Σ KᵢBᵢ)[Σ αⱼeⱼ] - c₀ - Σ Kᵢcᵢ.

% We want to use a fast gradient based optimisation method (in 
% 'functions/EVP_optim.m' we use the MATLAB 'fsolve' function) so
% calculating the Jacobian of f will slow things down. Instead we want to
% calculate J analytically and pass it to the solver. We can do this
% easily:

% ∂f/∂Kᵢ =  -Bᵢ[Σ αⱼeⱼ] - cᵢ
% ∂f/∂αⱼ =  (A - Σ KᵢBᵢ)eⱼ

% Hence J = [∂f/∂Kᵢ ∂f/∂αⱼ] is a (square) matrix-valued function of x.

% As initial guess for x can be given by guessing the Kᵢ and projecting the
% guess for 'a' back onto the basis to find the corresponding αⱼ. We then
% hope the solver finds the right root. Experience suggests this works
% quite well and a reasonable guess for the Kᵢ is generally sufficient even
% when guessing 'a = 0'. Of course, the physical problem this relates to
% has real eigenvalues. This method would get a bit more complicated in the
% general complex case but would probably still work if you split into real
% and imaginary parts and work with the new (double-size) system.

% How accurate are these solutions? We can control the error in the
% root-finding method to increase error as required.

% See 'functions/EVP_optim.m'.

% See also:

% https://github.com/mncrowe/QGDipoles.jl
% https://github.com/mncrowe/QGDipoles.m

% where this method is used to solve the original n-layer vortex problem
% that motivated these eigenvalue based ramblings.

% Let's see an example:

N = 10; n = 3;

A = (diag(N:-1:1,0) + diag(N-1:-1:1,1)/2 + diag(N-1:-1:1,-1)/2);

B = cell(n, 1);
for i = 1:n
    B{i} = -(diag((N:-1:1)./(1:N),0) + diag((N-1:-1:1)./(1:N-1),1)/i ...
                                     + diag((N-1:-1:1)./(2:N),-1)/i);
end

c = 2*((1:N)'<=(1:n+1))./((1:N)'+(1:n+1));
d = -(-1).^(1:N)'./((1:N)'.^(1:n));

K0 = 4*ones(n, 1);
[K, a] = EVP_optim(A, B, c, d, K0);

M = A;
C = c(:, 1);
for i = 1:n
    M = M - K(i)*B{i};
    C = C + K(i)*c(:, i+1);
end

diff = M*a - C;

disp(' '); disp('n-Parameter Inhomogeneous Problem:'); disp(' ')
disp('K = '); disp(K);
disp('How good is this solution? Difference between LHS and RHS:')
disp(' '); disp(diff)
