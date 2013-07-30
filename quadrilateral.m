function [X, Y] = quadrilateral(a1, a2, a3, a4, point1, point2)
%   function [X, Y] = quadrilateral(a1, a2, a3, a4)
% finds the quadrilateral formed by the intersections of the four lines
% with linear coefficient a1, a2, a3 and a4. 
% Lines with coefficient a1, a2 start from point 1. 
% Lines with coefficents a3, a4 start from point 2
% INPUTS: point1 = [0 0], 1x2 vector. 

x1 = point1(1); y1 = point1(2); 
b1 = y1 - a1*x1; 
b2 = y1 - a2*x1; 

x2 = point2(1); y2 = point2(2); 
b3 = y2 - a3*x2; 
b4 = y2 - a3*x2;

X = []; Y = [];
[x0, y0] = lineIntersect(a1, b1, a3, b3);
X = [X; x0]; Y = [Y; y0]; 
[x0, y0] = lineIntersect(a1, b1, a4, b4);
X = [X; x0]; Y = [Y; y0]; 
[x0, y0] = lineIntersect(a2, b2, a4, b4);
X = [X; x0]; Y = [Y; y0]; 
[x0, y0] = lineIntersect(a2, b2, a3, b3);
X = [X; x0]; Y = [Y; y0]; 

end
