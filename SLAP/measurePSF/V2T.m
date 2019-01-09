    function [T,Td] = V2T(Vx,Vy,Xoff,Yoff,angle)
        T = (Vx - Xoff)*cos(angle) - (Vy-Yoff)*sin(angle);
        Td = (Vx - Xoff)*sin(angle) + (Vy-Yoff)*cos(angle);
    end