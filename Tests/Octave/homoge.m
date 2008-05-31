% Pull matrix A into homegenous coords (& transpose if necessary)

function [B] = homoge(A)

[b a] = size(A);
if b > a
	A=A';
	[a,b] = swap(a,b);
endif

B = ones(b+1, a);
B(1:b, :) = A;

	