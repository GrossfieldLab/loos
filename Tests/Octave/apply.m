function [Y] = apply(X, M)

[b a] = size(X);
for i=1:a,
	Y(:,i) = M * X(:,i);
end
