function [u_phys0,v_phys0,p_phys0] = initial_value(xx,yy)

u_phys0 = - 128 * ((xx + 1) / 2).^2 .* ((xx - 1) / 2).^2 .* yy .* ((yy - 1) / 2) .* ((yy + 1) / 2);
v_phys0 = 128 * ((yy + 1) / 2).^2 .* ((yy - 1) / 2).^2 .* xx .* ((xx - 1) / 2) .* ((xx + 1) / 2);
p_phys0 = xx;

end