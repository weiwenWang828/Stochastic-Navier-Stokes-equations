function [u_phys,v_phys] = solve_u_2pi(u,v,p,delta_t,dx,dy)
phat = fftn(p);
p_x_hat = dx .* phat;
p_y_hat = dy .* phat;
p_x = real(ifftn(p_x_hat));
p_y = real(ifftn(p_y_hat));


u_phys = u - 2 * pi * delta_t * p_x;
v_phys = v - 2 * pi * delta_t * p_y;


end