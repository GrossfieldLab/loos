function [M] = translate(v)

M=eye(4);
M(1:3, 4) = v';
