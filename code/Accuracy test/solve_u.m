function [u_phys,v_phys] = solve_u(u,v,p,delta_t)

p_sptr_L = D2Lenphy_specShen(p, 1);
p_x_sptr = LDTL(p_sptr_L); p_x_phys = D2Lenphy_specShen(p_x_sptr, 2);
p_y_sptr = TLLDY(p_sptr_L); p_y_phys = D2Lenphy_specShen(p_y_sptr, 2);
u_phys = u - 2 * delta_t * p_x_phys;
v_phys = v - 2 * delta_t * p_y_phys;

end