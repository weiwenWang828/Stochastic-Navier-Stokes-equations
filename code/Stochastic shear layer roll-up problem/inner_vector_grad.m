function [inner] = inner_vector_grad(u1, v1, u2, v2, dx, dy, params)
u1hat = fftn(u1); v1hat = fftn(v1);
u2hat = fftn(u2); v2hat = fftn(v2);

u1_x_hat = dx .* u1hat; u1_y_hat = dy .* u1hat;
v1_x_hat = dx .* v1hat; v1_y_hat = dy .* v1hat;

u2_x_hat = dx .* u2hat; u2_y_hat = dy .* u2hat;
v2_x_hat = dx .* v2hat; v2_y_hat = dy .* v2hat;

u1_x = real(ifftn(u1_x_hat)); u1_y = real(ifftn(u1_y_hat));
v1_x = real(ifftn(v1_x_hat)); v1_y = real(ifftn(v1_y_hat));

u2_x = real(ifftn(u2_x_hat)); u2_y = real(ifftn(u2_y_hat));
v2_x = real(ifftn(v2_x_hat)); v2_y = real(ifftn(v2_y_hat));


inner = sum(sum(u1_x .* u2_x)) * params.hx * params.hy + sum(sum(u1_y .* u2_y)) * params.hx * params.hy +...
    sum(sum(v1_x .* v2_x)) * params.hx * params.hy + sum(sum(v1_y .* v2_y)) * params.hx * params.hy;


end