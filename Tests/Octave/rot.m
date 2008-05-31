% Rotate about an arbitrary axis.  Assumes LHS with pos Z away from origin

function [M] = rot(v, angle)

M=eye(4);
theta = pi * angle / 180.0;
c = cos(theta);
s = sin(theta);

M(1,1) = v(1) * v(1) * (1.0 - c) + c;
M(1,2) = v(1) * v(2) * (1.0 - c) - v(3) * s;
M(1,3) = v(1) * v(3) * (1.0 - c) + v(2) * s;

M(2,1) = v(1) * v(2) * (1.0 - c) + v(3) * s;
M(2,2) = v(2) * v(2) * (1.0 - c) + c;
M(2,3) = v(2) * v(3) * (1.0 - c) - v(1) * s;

M(3,1) = v(1) * v(3) * (1.0 - c) - v(2) * s;
M(3,2) = v(2) * v(3) * (1.0 - c) + v(1) * s;
M(3,3) = v(3) * v(3) * (1.0 - c) + c;
