function p_phys = solve_possion(u,v,dx,dy,delta_t,coe,N)

uhat = fftn(u); vhat = fftn(v);
u_x_hat = dx .* uhat; v_y_hat = dy .* vhat;
u_x = real(ifftn(u_x_hat)); v_y = real(ifftn(v_y_hat));

p_right = u_x + v_y;

phat = fftn(p_right) ./ coe;
phat([1,N/2+1],[1,N/2+1]) = 0;

p_phys = real(ifftn(phat));

end