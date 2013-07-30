alpha = 90 - w.trueBearings(1); 

a11 = tand(alpha - 1); 
a12 = tand(alpha + 1); 

figure(1); hold on; 

x = linspace(w.rcvLocs(1, 1), w.rcvLocs(1, 1) + 10e3, 100); 
plotline(a11, x, w.rcvLocs(1, 1), w.rcvLocs(1, 2), '-r'); 
plotline(a12, x, w.rcvLocs(1, 1), w.rcvLocs(1, 2), '-r'); 


alpha = 90 - w.trueBearings(30); 

a21 = tand(alpha - 1); 
a22 = tand(alpha + 1); 

figure(1); 

x0 = w.rcvLocs(30, 1); 
y0 = w.rcvLocs(30, 2); 

x = linspace(w.rcvLocs(1, 1), w.rcvLocs(1, 1) + 10e3, 100); 
plotline(a21, x, x0, y0, '-k'); 
plotline(a22, x, x0, y0, '-k'); 


[x1, y1] = lineIntersect(a12, w.rcvLocs(