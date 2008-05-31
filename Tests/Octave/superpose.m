% Compute the 4x4 transformation matrix to superimpose X onto Y
% Assumes X & Y are col-wise coords and homogenous...

function [M] = superpose(X, Y)

%%% First, center the two...


A = X(1:3, :);
B = Y(1:3, :);

ac = mean(A')';
bc = mean(B')';

A = A - diag(ac)*ones(size(A));
B = B - diag(bc)*ones(size(B));


R = A*B';   % Correlation
[U,S,V] = svd(R);

% Flip if necessary...
if det(R) < 0
	U(:, 3) = -U(:,3);
endif

MM = U * V';

M = eye(4);
M(1:3, 1:3) = MM';

% Add translation in...
M(4,1:3) = bc(1:3);



