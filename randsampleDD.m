function inds = randsampleDD(N, w)
% function inds = randsampleDD(N, w)
% a hack of the randsample.m file
% resamples the numbers from 1 to N following a weight distribution vector
% w

inds = [];
w = w/sum(w); 
w1 = round(w*N); 

for k = 1:length(w1)
    if w1(k) ~= 0 
        inds = [inds; repmat(k, w1(k), 1)]; 
    end
end


[w2, i] = sort(w1);
k = i(end-5 : end); 

while length(inds) < N
    draw = randsample(k, 1); 
    inds = [inds; draw]; 
end




end
