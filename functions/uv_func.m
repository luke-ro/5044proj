function [u,v] = uv_func(r,lm,R,u0,v0)
    u = f*dot((lm - r), R(:,1))/dot((lm - r), R(:,3)) + u0;
    v = f*dot((lm - r), R(:,2))/dot((lm - r), R(:,3)) + v0;
end

