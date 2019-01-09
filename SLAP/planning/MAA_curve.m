function Y = MAA_curve (d1,d2,v1,v2,T, res)
%A curve Y with Minimum Absolute Acceleration joining an initial position
%(d1) and velocity (v1) with a final position (d2) and velocity (v2), in a
%fixed time T
%res is the resolution of the returned curve, 
    %i.e. we return res+1 points spanning [0 T]

D = d2-d1;

if abs(v1-v2)<1e-10
    %spend half your time accelerating and half decelerating
    a = 4*((D./T) - v1)./T;
    x1 = linspace(0, T/2,ceil(100*res/2+1));
    y1 = d1+ v1*x1 + 0.5*a*x1.^2;
    
    y2 = d2 - v2*x1 -0.5*a*x1.^2;
    y2 = fliplr(y2);
else
    c1 = (v1-v2);
    c2 = 2*(v2*T - D);
    c3 = T*D - 0.5*v1*(T.^2) -0.5*v2*(T.^2);
    
    t = (-c2 + sqrt(c2.^2 - 4*c1*c3))/(2*c1);

if t<0 || t>T  %the other quadratic solution is to be used; we first accelerate then decelerate
    t = (-c2 - sqrt(c2.^2 - 4*c1*c3))/(2*c1);
    a = -(v2-v1)/(2*t - T);
    
    t = round(100*res*(t/T))* T./(100*res); %discretify

    x1 = linspace(0, t, ceil((t/T)*100*res+1));
    y1 = d1+ v1*x1 - 0.5*a*x1.^2;
    
    
    x2 = linspace(0, T-t, ceil((1-t/T)*100*res+1));
    y2 = d2 - v2*x2 + 0.5*a*x2.^2;
    y2 = fliplr(y2);
else  %the correct quadratic solution was found first, we decelerate, then accelerate
    a = (v2-v1)/(2*t - T);
    
    t = round(100*res*(t/T))* T./(100*res); %discretify
    
    x1 = linspace(0, t, ceil((t/T)*100*res+1));
    y1 = d1+ v1*x1 + 0.5*a*x1.^2;
    
    
    x2 = linspace(0, T-t, ceil((1-t/T)*100*res+1));
    y2 = d2 - v2*x2 - 0.5*a*x2.^2;
    y2 = fliplr(y2);
end
end

Y = [y1(1:end-1) y2(1:end)]; %[y1(1:end-1) (y1(end)+y2(1))/2 y2(2:end)]
Y = Y(1:100:end);