    function [Vx,Vy] = T2V(T,Td,Xoff,Yoff,angle)
        Vx = T*cos(angle) + Td*sin(angle) + Xoff;
        Vy = -T*sin(angle) + Td*cos(angle) + Yoff;
    end