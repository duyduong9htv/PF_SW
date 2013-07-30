function prob = calP_y_given_x(theta_y, theta_x)
% function prob = calP_y_given_x(theta_y, theta_x)
% calculates the probability that a bearing measurement theta_y is
% recorded, given that theta_x is the true bearing. 
% In this simplistic model, the beamwidth, centered at theta_x is assumed
% to be 1 degree to both sides. Therefore: p(theta_y) (within beam) = 1/2,
% p(theta_y) (outside beam) = 0. 

% if abs(theta_y - theta_x) <= 1. 
%     prob = 1/2; 
% else prob = 0 ; 
% end


prob = normpdf(theta_y, theta_x, 1.5); 
    
end
