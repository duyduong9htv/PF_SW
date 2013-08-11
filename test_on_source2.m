
load('/Users/dtran/Research/RLS_bearings_only/sourceLocalizationCodes/track_array_GPS_bearing_info_track570_3_new.mat')
rcv = track.receiverLocation; 
rcvTime = track.timeVector; 
measuredBearings = track.bearingEstimatesTrueNorth; 
arrayHeading = track.arrayHeading; 
measurementTime = rcvTime; 

w = TMA(rcv, rcvTime, measuredBearings,measurementTime, arrayHeading);

initRange = 12e3; 

r = [w.rcvLocs(1, 1)+initRange*sind(w.trueBearings(1)); ...
    w.rcvLocs(1, 2) + initRange*cosd(w.trueBearings(1))]'; %initial range 

r_all = []; 


%% plot geometry of problem 

figure; hold on; plot2dd(w.rcvLocs, '--k'); axis equal; 
plot2dd(track.sourceLocation, 'b*'); %source true locations 

plot2dd(r, 'r*'); %initial estimate 


%%

%generate N particles, first assuming that they're Gaussian distributed
%around the initial guess 

sigma = 20; 
N = 100; 
v = 0.5; %let velocity drift be about 2 m/s
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

figure; plot2dd(track.sourceLocation, '-r'); hold on; 
plot2dd(w.rcvLocs, '--k'); axis equal
plot2dd(r_all, '-*'); 