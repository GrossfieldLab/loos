% Translate matrix to center (assuming col-wise format & homogenous)

function [A] = center(A)


m = mean(A')';
m(4) = 0;

B = diag(m) * ones(size(A));
A = A - B;
