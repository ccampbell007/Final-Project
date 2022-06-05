function [circular_nodes] = circularcoordinates(R,base,N);

% This code develops the x-y coordinates based on the base of the circle,
% the radius, and the number of nodes

% Base Circle
theta = 2*pi/(N-1);              % Angle between points

% Establish x-y coordinates
n = 0;                       % Counting variable
for angle = -pi/2:theta:3*pi/2
    n = n+1;                 % Increment
x(n) = R*cos(angle);        % x coordinate
y(n) = R*sin(angle);        % y coordinate
end

% Shift Coordinates to have bottom at the base

for i = 1:N
    circular_nodes(2*i-1,:) = x(i);
    circular_nodes(2*i,:) = y(i) + R + base;
end

end