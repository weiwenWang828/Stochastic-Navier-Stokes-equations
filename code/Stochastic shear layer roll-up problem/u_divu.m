function [udivu1, udivu2] = u_divu(u,v,dx,dy)
uhat = fftn(u); vhat = fftn(v);
u_x_hat = dx .* uhat; u_y_hat = dy .* uhat;
v_x_hat = dx .* vhat; v_y_hat = dy .* vhat;
u_x = real(ifftn(u_x_hat)); u_y = real(ifftn(u_y_hat));
v_x = real(ifftn(v_x_hat)); v_y = real(ifftn(v_y_hat));

udivu1 = u .* u_x + v .* u_y;
udivu2 = u .* v_x + v .* v_y;

end