function p_phys = solve_possion(u,v,Inn_Lh,inv_left,coe,V,delta_t)

u_sptr_L = D2Lenphy_specShen(u, 1);
v_sptr_L = D2Lenphy_specShen(v, 1);
u_x_sptr = LDTL(u_sptr_L); u_x_phys = D2Lenphy_specShen(u_x_sptr, 2);
v_y_sptr = TLLDY(v_sptr_L); v_y_phys = D2Lenphy_specShen(v_y_sptr, 2);
p_right = u_x_phys + v_y_phys;

p_right_sptr_L = D2Lenphy_specShen(p_right, 1);
p_right_phys = Inn_Lh * p_right_sptr_L * Inn_Lh'; 

m1 = size(p_right_phys); 
n1 = m1(1,1); n2 = m1(1,2);

trans_p_sptr = (inv_left * p_right_phys * inv_left') ./ coe;
trans_p_sptr(n1, n2) = 0;
p_sptr = V * trans_p_sptr * V'; % p new sptr in g_k(x)g_l(y); (n-2)x(n-2);
p_sptr_L = phiDTL(p_sptr); % p new sptr in L_k(x)L_l(y); n x n;
p_phys = D2Lenphy_specShen(p_sptr_L, 2); % p new phys;


end