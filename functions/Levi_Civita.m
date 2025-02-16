function eps = Levi_Civita(p)
% Calculates the sign of the given permutation
%
% p: permutation of the numbers 1, 2, ..., n

I = speye(length(p));
eps = det(I(:,p));

end