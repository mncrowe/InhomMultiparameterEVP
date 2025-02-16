function K = multi_kron(M)
% Calculates M_1 ⊗ M_2 ⊗ ... ⊗ M_m for given M matrices
%
% M: matrices, array of size [N N m]

s = size(M);

if length(s) > 2

    K = M(:, :, 1);

    for i = 2:s(3)
        K = kron(K, M(:, :, i));
    end

else

    K = M;
    
end

end