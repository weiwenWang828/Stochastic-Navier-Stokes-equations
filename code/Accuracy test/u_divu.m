function [udivu1, udivu2] = u_divu(u, v)
u_sptr_L = D2Lenphy_specShen(u, 1);
v_sptr_L = D2Lenphy_specShen(v, 1);
u_x_sptr = LDTL(u_sptr_L); u_x_phys = D2Lenphy_specShen(u_x_sptr, 2);
u_y_sptr = TLLDY(u_sptr_L); u_y_phys = D2Lenphy_specShen(u_y_sptr, 2);
v_x_sptr = LDTL(v_sptr_L); v_x_phys = D2Lenphy_specShen(v_x_sptr, 2);
v_y_sptr = TLLDY(v_sptr_L); v_y_phys = D2Lenphy_specShen(v_y_sptr, 2);
udivu1 = u .* u_x_phys + v .* u_y_phys;
udivu2 = u .* v_x_phys + v .* v_y_phys;

end