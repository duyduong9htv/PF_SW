function h = plotline(a, x, x0, y0, varargin)
% function h = plotline(a,x, x0, y0)
% plots the line equation y = a*(x - x0) + y0 
% with linear coefficient a, going through (x0, y0). 
% Eg: x = linspace(-3, 3, 100); 
%figure; plotline(1, x, 2, 3, '--r'); 



y = a*(x - x0) + y0; 
plot(x, y, varargin{:}); 


end
