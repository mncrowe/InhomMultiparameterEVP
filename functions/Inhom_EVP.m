function [K,x] = Inhom_EVP(A,B,c,d,N)
% Solve the inhomogeneous EVP: (A-KB)x = c, the additional condition of
% x'Nx = d (matrix N, scalar d) or d'x = N (vector d, scalar N) is imposed.
%
% Inputs:
%
% - A,B: square matrices
% - c: vector
% - d,N: either matrix N and scalar d such that x'Nx = d or vector d and
%        scalar N such that d'x = N

cond_tol = 1e-8;
M = length(A);
n_eigs = min(M,10);

if length(N) == M && length(d) == 1     % Case 1: x'Nx = d

    if rcond(N) < cond_tol; error('Method does not currently support singular N.'); end
    if d == 0
        error('Method does not support d = 0 for this problem.')
    else
        N = N/d;
    end
    K = polyeig(A*(N\A')-c*c',-(A*(N\B')+A'*(N\B)),B*(N\B'));

else

if length(N) == 1 && length(d) == M     % Case 2: x'd = N

    if N == 0

        % Case of x'd = 0:
        if rcond(B) > rcond(A)  % invert B if better conditioned
            K = Gen_EVP(-B,A-(c*d')*(B\A)/(d'*(B\c)),n_eigs);
        else                    % invert A otherwise
            K = Gen_EVP(-B+(c*d')*(A\B)/(d'*(A\c)),A,n_eigs);
        end

    else

        % Case of x'd = N ~= 0:
        K = Gen_EVP(-B,A-(c*d')/N,n_eigs);

    end    

else

    error('Size of d and N are inconsistent with either condition format.')

end

end

% Order eigenvalues:
K = sort(K,1,"descend");

% Calculate eigenvectors for all cases:
x = zeros(M,length(K));
for ik = 1:length(K)
    x(:,ik) = (A-K(ik)*B)\c;
end

end