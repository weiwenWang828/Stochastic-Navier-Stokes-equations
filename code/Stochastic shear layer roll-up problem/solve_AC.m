function u_phys = solve_AC(U_right,Inn_Lh,inv_left,coe,V)
U_right_sptr = D2Lenphy_specShen(U_right, 1); % the spectral coefficient of U_right in L_k(x)L_l(y);
U_right_phys = Inn_Lh * U_right_sptr * Inn_Lh'; 
trans_u_sptr = (inv_left * U_right_phys * inv_left') ./ coe;
u_sptr = V * trans_u_sptr * V'; % u new sptr in h_k(x)h_l(y); (n-2)x(n-2);
u_sptr_L = phiTL(u_sptr); % u new sptr in L_k(x)L_l(y); n x n;
u_phys = D2Lenphy_specShen(u_sptr_L, 2); % u new phys;


end