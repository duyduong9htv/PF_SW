%% Jul 30 
load('/Users/dtran/Research/RLS_bearings_only/sourceLocalizationCodes/track_array_GPS_bearing_info_track540_1_new.mat')

load('/Users/dtran/Research/RLS_bearings_only/sourceLocalizationCodes/track_array_GPS_bearing_info_track570_3_new.mat')
rcv = track.receiverLocation; 
rcvTime = track.timeVector; 
measuredBearings = track.bearingEstimatesTrueNorth; 
arrayHeading = track.arrayHeading; 
measurementTime = rcvTime; 

w = TMA(rcv, rcvTime, measuredBearings,measurementTime, arrayHeading);


%% initialize
%parameters 
k = 1; 
N = 12; %N^2 = number of particles 
obs = []; %observation 
ocs = []; 
locs = []; 
step = 2; 
halfbw= 1.5;%half beam-width 
 
disp(k); 
alpha = 90 - w.trueBearings(k); 

a11 = tand(alpha - halfbw); 
a12 = tand(alpha + halfbw); 


alpha = 90 - w.trueBearings(k+step); 
a21 = tand(alpha - halfbw); 
a22 = tand(alpha + halfbw); 

for a1 = linspace(a11, a12, N)
    b1 = w.rcvLocs(k, 2) - a1*w.rcvLocs(k, 1); 
    for a2 = linspace(a21, a22, N); 
        b2 = w.rcvLocs(k+step, 2) - a2*w.rcvLocs(k+step, 1); 
        [x0, y0] = lineIntersect(a1, b1, a2, b2);
        obs = [obs; x0 y0]; 
%              figure(6); plot(x0, y0, 'ko', 'linewidth', 1); 

    end   
end     

% plot2dd((obs), 'bo', 'linewidth', 1); 

particles = obs; %set the initial observations to be the particles 

% 
% 
% x0 = w.rcvLocs(1, 1) + 2e3; 
% y0 = w.rcvLocs(1, 2); 
% particles = repmat([x0 y0], N^2, 1);
locs = [locs; mean(particles)]; 
%% Recursive particle filtering 



for k =  2:130 %150 %2:850% 100 %850 %100 %850% 551:600 %880 %:1:878
%     N = 15

    %% prediction 
    dt = mean([w.whaleTime(k) w.whaleTime(k+step)]) ...
            - mean([w.whaleTime(k - 1), w.whaleTime(k+step -1 )]); 
    std_dev = dt*2; 
    particles = updateParticles(particles, std_dev); 
    
    %% obtaining observations 
   
    obs = []; %observation 
%     disp(k); 
    
    % first beam steered from first rcv location 
    alpha = 90 - w.trueBearings(k); 
    
    a11 = tand(alpha - halfbw); 
    a12 = tand(alpha + halfbw); 
    
    %second beam steered from second rcv location
    alpha = 90 - w.trueBearings(k+step); 
    a21 = tand(alpha - halfbw); 
    a22 = tand(alpha + halfbw); 
    
    for a1 = linspace(a11, a12, N)
        b1 = w.rcvLocs(k, 2) - a1*w.rcvLocs(k, 1); 
        for a2 = linspace(a21, a22, N); 
            b2 = w.rcvLocs(k+step, 2) - a2*w.rcvLocs(k+step, 1); 
            [x0, y0] = lineIntersect(a1, b1, a2, b2);
            if ddist([x0 y0], w.rcvLocs(k, :)) < 30e3 %constrained within some distance 
            obs = [obs; x0 y0]; 
            end 
%              figure(6); plot(x0, y0, 'ko', 'linewidth', 1); 
             
        end    
    end 

    inds1 = []; 
    for k1 = 1:size(obs, 1)
        if ddist(obs(k1, :), locs(end, :)) < 10e3
            inds1 = [inds1; k1]; 
        end
    end
    
    obs1 = obs(inds1, :); 
    disp(size(obs, 1));
    
    if (size(obs1, 1) < 0.1*N^2)
        particles = obs; 
        locs = [locs; mean(obs)]; 
    else 
        obs = obs1; 
        
        
    
        %% correction of predicted probability 
        p_update = []; 
        bearings = []; 
        tests = []; 
        for ii = 1:size(obs, 1)
            theta_x = calBearings(particles(ii, :), w.rcvLocs(k, :)); 

            %first beam steered towards particle ii 
            alpha = 90 - theta_x; 
            a11 = tand(alpha - halfbw); 
            a12 = tand(alpha + halfbw); 

            %second beam steered toward particle ii 

            theta_x = calBearings(particles(ii, :), w.rcvLocs(k+1, :)); 
            alpha = 90 - theta_x; 
            a21 = tand(alpha - halfbw); 
            a22 = tand(alpha + halfbw); 


            %find quadrilateral region defined by intersections of the 2 beams

            [X, Y] = quadrilateral(a11, a12, a21, a22, ...
                                    w.rcvLocs(k, :), w.rcvLocs(k+1, :));

            % find which observations are valid given the predicted random
            % drift model 
             test = inpolygon(obs(ii, 1), obs(ii, 2), X, Y);
             tests = [tests; test]; 

        %     theta_y = calBearings(obs(ii, :), w.rcvLocs(k, :)); 
        %     p = calP_y_given_x(theta_y, theta_x); 
        %     p_update = [p_update; p];     
        end

        disp(sum(tests)); pause(0.1); 

        %% find valid observations according to model 
        inds = tests(:) > 0;
        temp = particles(inds, :);
        inds = find(temp(:, 1) < 20e3 + w.rcvLocs(1, 1)); %one check point, find only whales on the right
        temp = temp(inds, :);
        inds = find(temp(:, 1) > w.rcvLocs(1, 1) - 20e3); 
        temp = temp(inds, :); 

    %     inds1 = []; 
    %     for k1 = 1:size(temp, 1)
    %         if ddist(temp(k1, :), locs(end, :)) < 2e3
    %             inds1 = [inds1; k]; 
    %         end
    %         
    %     end
    %     
    %     temp = temp(inds1, :); 
    %     


        M = size(temp, 1); 


        if M > 10
            resampled_indices = randsample(M, N^2, true); 
            particles = temp(resampled_indices, :); 
            locs = [locs; mean(temp)];
        else 

            d = []; 
            for m = 1:length(obs)
                d = [d; ddist(obs(m, :), locs(end, :))]; 
            end

            inds = find(d(:) < 1000); 

            %constrain the target not to move too much from previous location
            if ~isempty(inds)
                obs = obs(inds, :); 
            else
                obs = particles; 
            end 


    %         inds = find(obs(:, 1) > w.rcvLocs(1, 1)); 
    %         obs = obs(inds, :); 
            if size(obs, 1) == 1
               locs = [locs; obs];
            else 
                locs = [locs; mean(obs)]; 
            end 
            particles = obs; 
            M = size(particles, 1); 
            resampled_indices = randsample(M, N^2, true);
            particles = particles(resampled_indices, :); 
        end 

        bearings = [bearings; theta_x]; 
    end 
    
end

%% check 

figure; plot2dd(w.rcvLocs, '--k'); hold on; 
plot2dd(locs, 'ks');  plot2dd(locs, '--r'); axis equal; 
plot2dd(track.sourceLocation, '-b', 'linewidth', 3); 

