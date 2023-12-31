function [J] = dyn_jacobian_H(x, l, R,const)
    r1 = x(1);
    r2 = x(2);
    r3 = x(3);

    r = x(1:3);

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

    J11 = const.f*(-dot((l-r),k)*i1+dot(l-r,i)*k1)/(dot(l-r,k))^2;
    J12 = const.f*(-dot((l-r),k)*i2+dot(l-r,i)*k2)/(dot(l-r,k))^2;
    J13 = const.f*(-dot((l-r),k)*i3+dot(l-r,i)*k3)/(dot(l-r,k))^2;
    
    J21 = const.f*(-dot((l-r),k)*j1+dot(l-r,j)*k1)/(dot(l-r,k))^2;
    J22 = const.f*(-dot((l-r),k)*j2+dot(l-r,j)*k2)/(dot(l-r,k))^2;
    J23 = const.f*(-dot((l-r),k)*j3+dot(l-r,j)*k3)/(dot(l-r,k))^2;
    
    J = [J11 J12 J13 0 0 0;
         J21 J22 J23 0 0 0];
end