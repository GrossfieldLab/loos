function [d] = rmsd(A, B)

X = A(1:3, :);
Y = B(1:3, :);

Z = X-Y;
[b a] = size(Z);

for i=1:a,
	z(i) = norm(Z(:, i));
end

z = z .* z;
d = sqrt(mean(z));
