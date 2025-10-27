function [p_x, p_y] = partial(p,dx,dy)
phat = fftn(p);
p_x_hat = dx .* phat;
p_y_hat = dy .* phat;
p_x = real(ifftn(p_x_hat));
p_y = real(ifftn(p_y_hat));

end