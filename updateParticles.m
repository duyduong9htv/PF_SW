function p = updateParticles(data_in, std_dev)
% function p = updateParticles(data_in, std_dev)
% updates the particle based on the motion model. In this case the motion
% model is a random 2-D Brownian drift with constant velocity. std_dev =
% velocity *delta_t (depending on how far apart in time the 2 measurements
% are) 
% INPUTS: data_in (N x D), N = number of particles, D = dimention (1, 2,
% 3...) depending on the number of coordinates to be tracked. 


X = data_in; %for easy notation 
N = size(X, 1); %number of particles; 
D = size(X, 2); %number of dimensions 

p = X + std_dev*randn(N, D); 

end
