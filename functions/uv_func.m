function [u,v] = uv_func(r,lm,R,u0,v0)
    f = 2080.7959; % [pixels]
    n = size(lm,2);
    if(n==1)
        u = f*dot((lm - r), R(:,1))/dot((lm - r), R(:,3)) + u0;
        v = f*dot((lm - r), R(:,2))/dot((lm - r), R(:,3)) + v0;
        return
    end
    
    ic = repmat(R(:,1),1,n);
    jc = repmat(R(:,2),1,n);
    kc = repmat(R(:,3),1,n);

    % do the u,v calc for a vector of many landmarks by specifying axis
    % to do the dot product along
    u = f*dot((lm - r), ic, 1)./dot((lm - r), kc, 1) + u0;
    v = f*dot((lm - r), jc, 1)./dot((lm - r), kc, 1) + v0;
    
end

