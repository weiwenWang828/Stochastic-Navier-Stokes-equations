function [u_numerical_ref,v_numerical_ref,p_numerical_ref,xi_numerical_ref,eta_numerical_ref,Q_wiener1,Q_wiener2,p_numerical_dt] = ...
    NS_reference(xx,yy,u_phys0,v_phys0,p_phys0,delta_t,T,N,E,V,E_Neu,V_Neu,Inn_Lh,inv_left,Inn_Lh_Neu,inv_left_Neu,params)
KK = T / delta_t;
xi_ref = zeros(ceil(KK) + 1); tt_ref = xi_ref;
eta_ref = xi_ref;    
Q_wiener1 = zeros(N, N, KK);
Q_wiener2 = Q_wiener1;
p_numerical_dt = zeros(N, N, KK);

u_phys_n_1 = u_phys0;    
v_phys_n_1 = v_phys0;
p_phys_n_1 = p_phys0;
xi_n_1 = 1; 
eta_n_1 = 1; 
xi_ref(1) = xi_n_1;
eta_ref(1) = eta_n_1;
tt_ref(1) = 0;

%% main loop

for ii = 1 : KK
    TT = ii * delta_t;

    coe1 = 1 + params.mu * 4 * delta_t * (E * ones(N - 2) + ones(N - 2) * E);
    coe2 = -2 * delta_t * (E_Neu * ones(N - 2) + ones(N - 2) * E_Neu);


    [g1, g2] = g(u_phys_n_1, v_phys_n_1);

    
    [p_x_phys_n_1, p_y_phys_n_1] = partial(p_phys_n_1);
    [udivu1, udivu2] = u_divu(u_phys_n_1, v_phys_n_1);
    [wiener1,wiener2] = wiener(xx,yy,delta_t,N,params);
    Q_wiener1(:, :, ii) = wiener1;
    Q_wiener2(:, :, ii) = wiener2;
    
    U1_right = u_phys_n_1 - 2 * delta_t * p_x_phys_n_1; 
    V1_right = v_phys_n_1 - 2 * delta_t * p_y_phys_n_1;

    tilde_u1_phys = solve_AC(U1_right,Inn_Lh,inv_left,coe1,V);
    tilde_v1_phys = solve_AC(V1_right,Inn_Lh,inv_left,coe1,V);

    U2_right = - 2 * delta_t * udivu1;
    V2_right = - 2 * delta_t * udivu2;

    tilde_u2_phys = solve_AC(U2_right,Inn_Lh,inv_left,coe1,V);
    tilde_v2_phys = solve_AC(V2_right,Inn_Lh,inv_left,coe1,V);

    U3_right = g1 .* wiener1;
    V3_right = g2 .* wiener2;

    tilde_u3_phys = solve_AC(U3_right,Inn_Lh,inv_left,coe1,V);
    tilde_v3_phys = solve_AC(V3_right,Inn_Lh,inv_left,coe1,V);



    p1_phys = solve_possion(tilde_u1_phys,tilde_v1_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,delta_t); 
    p2_phys = solve_possion(tilde_u2_phys,tilde_v2_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,delta_t); 
    p3_phys = solve_possion(tilde_u3_phys,tilde_v3_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,delta_t); 

    [u1_phys,v1_phys] = solve_u(tilde_u1_phys,tilde_v1_phys,p1_phys,delta_t);
    [u2_phys,v2_phys] = solve_u(tilde_u2_phys,tilde_v2_phys,p2_phys,delta_t);
    [u3_phys,v3_phys] = solve_u(tilde_u3_phys,tilde_v3_phys,p3_phys,delta_t);


    
    inner1 = 1 + params.mu * 4 * delta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u2_phys,tilde_v2_phys,N) +...
        inner_product(tilde_u2_phys,tilde_u2_phys,N) + inner_product(tilde_v2_phys,tilde_v2_phys,N);

    inner2 = params.mu * 4 * delta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u3_phys,tilde_v3_phys,N) +...
        inner_product(tilde_u2_phys,tilde_u3_phys,N) + inner_product(tilde_v2_phys,tilde_v3_phys,N);

    inner3 = inner2;

    inner4 = 1 + params.mu * 4 * delta_t * inner_vector_grad(tilde_u3_phys,tilde_v3_phys,tilde_u3_phys,tilde_v3_phys,N) +...
        inner_product(tilde_u3_phys,tilde_u3_phys,N) + inner_product(tilde_v3_phys,tilde_v3_phys,N);

    A = [inner1 inner2; inner3 inner4];

    b1 = xi_n_1 + 2 * delta_t * (inner_product(udivu1,tilde_u1_phys,N) + inner_product(udivu2,tilde_v1_phys,N));

    b2 = eta_n_1 - inner_product(U3_right,tilde_u1_phys-u_phys_n_1,N) -...
        inner_product(V3_right,tilde_v1_phys-v_phys_n_1,N) +...
        inner_product(U3_right,U3_right,N) + inner_product(V3_right,V3_right,N);

    b = [b1; b2];

    X = A \ b;

    xi_n = X(1);
    eta_n = X(2);

    u_phys_n = u1_phys + xi_n * u2_phys + eta_n * u3_phys;
    v_phys_n = v1_phys + xi_n * v2_phys + eta_n * v3_phys;
    p_phys_n = p1_phys + p_phys_n_1 + xi_n * p2_phys + eta_n * p3_phys;



    u_phys_n_1 = u_phys_n; v_phys_n_1 = v_phys_n; p_phys_n_1 = p_phys_n; xi_n_1 = xi_n; eta_n_1 = eta_n; % at t = t_n
    

    xi_ref(ii + 1) = xi_n_1;
    eta_ref(ii + 1) = eta_n_1;
    tt_ref(ii + 1) = TT;
    p_numerical_dt(:, :, ii) = p_phys_n;



end

u_numerical_ref = u_phys_n_1;
v_numerical_ref = v_phys_n_1;
p_numerical_ref = p_phys_n_1;
xi_numerical_ref = xi_n_1;
eta_numerical_ref = eta_n_1;


end






