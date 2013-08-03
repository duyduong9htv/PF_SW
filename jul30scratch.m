%% Jul 30 

load whaleAMay14loc_results.mat


%% initialize
%parameters 
k = 1; 
N = 30; %N^2 = number of particles 
obs = []; %observation 
ocs = []; 
locs = []; 

halfbw= 1;%half beam-width 
 
disp(k); 
alpha = 90 - w.trueBearings(k); 

a11 = tand(alpha - halfbw); 
a12 = tand(alpha + halfbw); 


alpha = 90 - w.trueBearings(k+30); 
a21 = tand(alpha - halfbw); 
a22 = tand(alpha + halfbw); 

for a1 = linspace(a11, a12, N)
    b1 = w.rcvLocs(k, 2) - a1*w.rcvLocs(k, 1); 
    for a2 = linspace(a21, a22, N); 
        b2 = w.rcvLocs(k+30, 2) - a2*w.rcvLocs(k+30, 1); 
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



for k =  2:850 %2:850% 100 %850 %100 %850% 551:600 %880 %:1:878
%     N = 15
    dt = mean([w.whaleTime(k) w.whaleTime(k+30)]) ...
            - mean([w.whaleTime(k - 1), w.whaleTime(k+29)]); 
    std_dev = dt*2; 
    particles = updateParticles(particles, std_dev); 
   
    obs = []; %observation 
    disp(k); 
    
    % first beam steered from first rcv location 
    alpha = 90 - w.trueBearings(k); 
    
    a11 = tand(alpha - halfbw); 
    a12 = tand(alpha + halfbw); 
    
    %second beam steered from second rcv location
    alpha = 90 - w.trueBearings(k+30); 
    a21 = tand(alpha - halfbw); 
    a22 = tand(alpha + halfbw); 
    
    for a1 = linspace(a11, a12, N)
        b1 = w.rcvLocs(k, 2) - a1*w.rcvLocs(k, 1); 
        for a2 = linspace(a21, a22, N); 
            b2 = w.rcvLocs(k+30, 2) - a2*w.rcvLocs(k+30, 1); 
            [x0, y0] = lineIntersect(a1, b1, a2, b2);
            obs = [obs; x0 y0]; 
%              figure(6); plot(x0, y0, 'ko', 'linewidth', 1); 
             
        end   
    end     

    
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

        theta_x = calBearings(particles(ii, :), w.rcvLocs(k+30, :)); 
        alpha = 90 - theta_x; 
        a21 = tand(alpha - halfbw); 
        a22 = tand(alpha + halfbw); 
        
        
        %find quadrilateral region defined by intersections of the 2 beams

        [X, Y] = quadrilateral(a11, a12, a21, a22, ...
                                w.rcvLocs(k, :), w.rcvLocs(k+30, :));
        
        % find which observations are valid given the predicted random
        % drift model 
         test = inpolygon(obs(ii, 1), obs(ii, 2), X, Y);
         tests = [tests; test]; 
         
    %     theta_y = calBearings(obs(ii, :), w.rcvLocs(k, :)); 
    %     p = calP_y_given_x(theta_y, theta_x); 
    %     p_update = [p_update; p];     
    end

    %% find valid observations according to model 
    inds = tests(:) > 0;
    temp = particles(inds, :);
    inds = find(temp(:, 1) > w.rcvLocs(1, 1)); %one check point, find only whales on the right
    temp = temp(inds, :);
     
    M = size(temp, 1); 
    if M > 50
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
        locs = [locs; mean(obs)]; 
        particles = obs; 
        M = size(particles, 1); 
        resampled_indices = randsample(M, N^2, true);
        particles = particles(resampled_indices, :); 
    end 
    
    bearings = [bearings; theta_x]; 
    
end

%% check 

figure; plot2dd(w.rcvLocs, '--k'); hold on; 
plot2dd(locs, 'ks');  plot2dd(locs, '--r'); axis equal 


%% 

figure; hold on; 
for k = 1:850
    d = ddist(locs(k, :), w.rcvLocs(k, :)); 
    plot(w.whaleTime(k), d, '*'); 
end


for k = 1:850
    plot2dd(locs(k, :), 'g*'); pause(0.1); 
end


%% choose points 

figure; plot(mod(w.whaleTime(1:850)/60, 60), locs(:, 1), 'ko'); hold on; 

[xPolygon,yPolygon]=getline(gcf,'closed')
inx = find(inpolygon(w.whaleTime(1:850), locs(:, 1), xPolygon, yPolygon))
hold on; 
plot(w.whaleTime(inx), locs(inx, 1), 'ro'); 



figure; plot(w.whaleTime(1:850), locs(:, 2), 'ko'); hold on; 

[xPolygon,yPolygon]=getline(gcf,'closed')
iny = find(inpolygon(w.whaleTime(1:850), locs(:, 2), xPolygon, yPolygon))
hold on; 
plot(w.whaleTime(iny), locs(iny, 2), 'ro'); 



inds = intersect(inx, iny)


%% plot 3-D result 
% clear
load sw_smooth_x_y
figure;
% depths_all = movingAverage(t, depths_all, 30); 
x_origin = w.rcvLocs(1, 1); 
y_origin = w.rcvLocs(1, 2); 
plot3(w.rcvLocs(:, 1) - x_origin,...
        w.rcvLocs(:, 2) - y_origin, -65*ones(size(w.rcvLocs(:, 1))), '-k', 'linewidth', 2); 

plot3(w.rcvLocs(:, 1) - x_origin,...
        w.rcvLocs(:, 2) - y_origin, -160*ones(size(w.rcvLocs(:, 1))), '--b', 'linewidth', 2); 

plot3([0 0], [0 0], [-65 0], '--k')
plot3(w.rcvLocs(end, 1)*[1 1] - x_origin, w.rcvLocs(end, 2)*[1 1] - y_origin, [-65 0], '--k')   
    
hold on; 
plot3(x1 - x_origin, y1-y_origin, -depths_all, '--ro');

plot3(x1 - x_origin, y1 - y_origin, -160*ones(size(x1)), 'bo');
plot3(x1 - x_origin, y1 - y_origin, -160*ones(size(x1)), '--k', 'linewidth', 2);
grid on; 
box on; 

% for k = [1 30 31 150 310 311 340 341 450 540]
%     plot3([x1(k) x1(k)] - x_origin, [y1(k), y1(k)]-y_origin, [-depths_all(k) -160], 'k-')
% end

setFont(16); setFigureAuto; 
xlabel('Eastings (m)'); ylabel('Northings (m)'); 
zlabel('Depth (m)'); 
ztl = get(gca, 'ztick'); 
set(gca, 'zticklabel', -ztl); 

