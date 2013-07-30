%% coding out a particle filter to track the moving sperm whale based on
%% bearing measurement


load whaleAMay14loc_results.mat

initRange = 1200; 

r = [w.rcvLocs(1, 1)+initRange*sind(w.trueBearings(1)); ...
    w.rcvLocs(1, 2) + initRange*cosd(w.trueBearings(1))]'; %initial range 

r_all = []; 

%generate N particles, first assuming that they're Gaussian distributed
%around the initial guess 

sigma = 20; 
N = 1000; 
v = 20; %let velocity drift be about 2 m/s
r_particles = repmat(r, N, 1); 
r_all = r;  

r_est = r; 

for t = 2:length(w.trueBearings)
    y_new = w.trueBearings(t); 
    p_y_given_x = []; 
    delta_t = w.whaleTime(t) - w.whaleTime(t - 1); 
    std_dev = v * delta_t; 
    
    if std_dev == 0
        std_dev = 10; 
    end
    
    %prediction 
    r_particles = updateParticles(r_particles, std_dev); 
    
    %update/correction 
    p_y_given_x = []; 

    theta_x = calBearings(r_particles, w.rcvLocs(t, :)); 
    p_update = zeros(size(theta_x)); 
    p_predict = p_update; 

    for k = 1:N
        p_y_given_x = [p_y_given_x; calP_y_given_x(y_new, theta_x(k))]; 
        p_predict(k) = mvnpdf(r_particles(k, :)', r_est', std_dev*eye(2)); 
        p_update(k) = p_predict(k)*p_y_given_x(k); 
    end
    
    
   
    %resample the particles: 
    inds = randsampleDD(N, p_update); 
%     inds = randsample(1:N, N, true,p_update/sum(p_update));
    r_particles = r_particles(inds, :); 
    
    %estimate location based on distribution
    
    r_est = [0 0]; 
    for k = 1:N
        r_est = r_est + p_update(k)*r_particles(k, :)/sum(p_update); 
    end

    r_all = [r_all; r_est]; 
    disp(t); 
    
end
%% check 

figure; hold on; 
plot2dd(w.rcvLocs, '--k'); 
plot(w.xSATinst, w.ySATinst, '*'); 
plot2dd(r_all, 'r--+'); 
plot2dd(r_all(1, :), 'go'); 
%% 

plot2dd(r_particles, 'k.')

%% 


t = 2; 
y_new = w.trueBearings(2); 
p_y_given_x = []; 

theta_x = calBearings(r_particles, w.rcvLocs(2, :)); 
p_update = zeros(size(theta_x)); 
p_predict = p_update; 

for k = 1:N
    p_y_given_x = [p_y_given_x; calP_y_given_x(y_new, theta_x(k))]; 
    p_predict(k) = mvnpdf(r_particles(k, :)', r', sigma*eye(2)); 
    p_update(k) = p_predict(k)*p_y_given_x(k); 
end

r_est = [0 0]; 
for k = 1:N
r_est = r_est + p_update(k)*r_particles(k, :)/sum(p_update); 
end




%% r = [w.xSATinst(20); w.ySATinst(1)]; 
sigma = 20; 
N = 1000; 

for k = 1:N
    r_particles(:, k) = r + sigma*randn(2, 1); 
end

figure; plot2dd(w.rcvLocs, '--k'); hold on; 
plot2dd(r_particles', 'k.'); 
axis equal

theta = w.trueBearings(2); 

bearing_err = 1.5; %degree (so 3 degrees beamwidth) 

theta_new = w.trueBearings(2); 

thetas = calBearings(r_particles', w.rcvLocs(2, :)); 

for k = 1:length(thetas)
    
end





