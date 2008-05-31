% Rotate about a specified axis and angle

function [M] = rotate(axis, angle)

if axis == 'x'
	M = rot([1 0 0], angle);
elseif axis == 'y'
	M = rot([0 1 0], angle);
else
	M = rot([0 0 1], angle);
endif
