function [u_phys0,v_phys0,p_phys0] = initial_value(xx,yy)

m1 = size(xx); m2 = size(yy);
n1 = m1(1,1); n2 = m2(1,2);
rho = 30;
u_phys0 = tanh(rho*(yy./(2*pi)-0.25)).*(yy<=pi)...
   + tanh(rho*(0.75-yy./(2*pi))).*(yy>pi);
v_phys0 = 0.05.*sin(xx);
p_phys0 = zeros(n1, n2);

end