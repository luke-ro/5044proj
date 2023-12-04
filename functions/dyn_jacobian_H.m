function [J] = dyn_jacobian_H(r, l, R)
    r1 = r(1);
    r2 = r(2);
    r3 = r(3);

    l1 = l(1);
    l2 = l(2);
    l3 = l(3);

    i = R(:,1);
    j = R(:,2);
    k = R(:,3);

    i1 = i(1);
    i2 = i(2);
    i3 = i(3);

    j1 = j(1);
    j2 = j(2);
    j3 = j(3);

    k1 = k(1);
    k2 = k(2);
    k3 = k(3);

    J11 = ((l2 - r2)*(i2*k1-i1*k2) + (l3 - r3)*(i3*k1 - i1*k3))/(dot(l-r,k))^2;
    J21 = ((l2 - r2)*(j2*k1-j1*k2) + (l3 - r3)*(j3*k1 - j1*k3))/(dot(l-r,j))^2;
    J12 = ((l1 - r1)*(i1*k2-i2*k1) + (l3 - r3)*(i3*k2 - i2*k3))/(dot(l-r,k))^2;
    J22 = ((l1 - r1)*(j1*k2-j2*k1) + (l3 - r3)*(j3*k2 - j2*k3))/(dot(l-r,j))^2;
    J13 = ((l2 - r2)*(i2*k3-i1*k3) + (l1 - r1)*(i1*k3 - i3*k1))/(dot(l-r,k))^2;
    J23 = ((l2 - r2)*(j2*k3-j1*k3) + (l1 - r1)*(j1*k3 - j3*k1))/(dot(l-r,j))^2;
    
    J = [J11 J12 J13 0 0 0;
         J21 J22 J23 0 0 0];



end