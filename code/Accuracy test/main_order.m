clear;clc
format long


T = 0.2;
k0 = 1/12800; 
n = 1;
delta_t = k0;
KK = T / delta_t;

n_dt = 5;
delta_t_m = 1/200; 
k_m = T / delta_t_m;

N_p = 300; 

rng('shuffle');


params.mu = 1; 
params.J = 4;
params.epsilon = 0.0001;


N = 40; 
[x,w1]= legslb(N);
[y,w2]= legslb(N);
[xx,yy] = meshgrid(x,y);
M = basisM(N - 1); 
S = basisS(N - 1); 
[V, E] = eig(M \ S); 
inv_left  = inv(M * V); 
Inn_Lh = NNP2PL(N - 1); 

M_Neu = basisMN(N - 1); 
S_Neu = basisMNDN(N - 1); 
[V_Neu, E_Neu] = eig(M_Neu \ S_Neu);  
inv_left_Neu  = inv(M_Neu * V_Neu); 
Inn_Lh_Neu = NP2NPNL(N); 
Inn_Lh_Neu = Inn_Lh_Neu'; 

error_velocity_L2 = zeros(n_dt, N_p);
error_pressure_L2 = error_velocity_L2;
xi_mean = error_velocity_L2;
eta_mean = error_velocity_L2;
error_pressure_time = error_velocity_L2;
xi_numerical_ref_mean = zeros(N_p, 1);
eta_numerical_ref_mean = xi_numerical_ref_mean;


for n_p = 1 : N_p

    
    [u_phys0,v_phys0,p_phys0] = initial_value(xx,yy); % initial value;
    time_step_ref = zeros(n , 1);
    xi_ref = zeros(ceil(2^(n-1) * KK + 1), n); tt_ref = xi_ref;
    eta_ref = xi_ref;
    
    time_step = zeros(n_dt , 1);
    xi = zeros(2^(n_dt-1) * k_m + 1, n_dt); tt = xi;
    eta = xi;
    
    u_numerical = zeros(N, N, n_dt);
    v_numerical = u_numerical;
    p_numerical = u_numerical;
    xi_numerical = zeros(n_dt, 1);
    eta_numerical = xi_numerical;
    
        
    [u_numerical_ref,v_numerical_ref,p_numerical_ref,xi_numerical_ref,eta_numerical_ref,Q_wiener1,Q_wiener2,p_numerical_dt_ref] = ...
        NS_reference(xx,yy,u_phys0,v_phys0,p_phys0,delta_t,T,N,E,V,E_Neu,V_Neu,Inn_Lh,inv_left,Inn_Lh_Neu,inv_left_Neu,params);

    sdelta_t = delta_t_m;
    kk = k_m;
    xi_numerical_ref_mean(n_p, 1) = xi_numerical_ref;
    eta_numerical_ref_mean(n_p, 1) = eta_numerical_ref;


    
    for it = 1 : n_dt
        
        time_step(it , 1) = sdelta_t;
        TT = 0;
        u_phys_n_1 = u_phys0;    
        v_phys_n_1 = v_phys0;
        p_phys_n_1 = p_phys0;
        xi_n_1 = 1; 
        eta_n_1 = 1; 
        xi(1, it) = xi_n_1;
        eta(1, it) = eta_n_1;
        tt(1, it) = 0;
    
        
        num_groups = kk;         
        group_size = sdelta_t / delta_t;         
        sum_Q_wiener1 = zeros(N, N, num_groups); 
        sum_Q_wiener2 = sum_Q_wiener1;

        p_numerical_dt = zeros(N, N, num_groups);
        
        
        for k = 1 : num_groups            
            start_idx = (k - 1) * group_size + 1;  
            end_idx = k * group_size;            
            
            
            current_group1 = Q_wiener1(:, :, start_idx:end_idx);
            current_group2 = Q_wiener2(:, :, start_idx:end_idx);
            
            
            sum_Q_wiener1(:, :, k) = sum(current_group1, 3);
            sum_Q_wiener2(:, :, k) = sum(current_group2, 3);
        end
    

    
        for ii = 1 : kk
            TT = ii * sdelta_t;
    
            coe1 = 1 + 4 * params.mu * sdelta_t * (E * ones(N - 2) + ones(N - 2) * E);
            coe2 = -2 * sdelta_t * (E_Neu * ones(N - 2) + ones(N - 2) * E_Neu);
            [g1, g2] = g(u_phys_n_1, v_phys_n_1);
    
           
            
            [p_x_phys_n_1, p_y_phys_n_1] = partial(p_phys_n_1);
            [udivu1, udivu2] = u_divu(u_phys_n_1, v_phys_n_1);
 
            wiener1 = sum_Q_wiener1(:, :, ii);
            wiener2 = sum_Q_wiener2(:, :, ii);
            
            U1_right = u_phys_n_1 - 2 * sdelta_t * p_x_phys_n_1; 
            V1_right = v_phys_n_1 - 2 * sdelta_t * p_y_phys_n_1;
    
            tilde_u1_phys = solve_AC(U1_right,Inn_Lh,inv_left,coe1,V);
            tilde_v1_phys = solve_AC(V1_right,Inn_Lh,inv_left,coe1,V);
    
            U2_right = - 2 * sdelta_t * udivu1;
            V2_right = - 2 * sdelta_t * udivu2;
    
            tilde_u2_phys = solve_AC(U2_right,Inn_Lh,inv_left,coe1,V);
            tilde_v2_phys = solve_AC(V2_right,Inn_Lh,inv_left,coe1,V);
    
            U3_right = g1 .* wiener1;
            V3_right = g2 .* wiener2;
    
            tilde_u3_phys = solve_AC(U3_right,Inn_Lh,inv_left,coe1,V);
            tilde_v3_phys = solve_AC(V3_right,Inn_Lh,inv_left,coe1,V);
    
    

            p1_phys = solve_possion(tilde_u1_phys,tilde_v1_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,sdelta_t); % p1_{n+1}-p_{n};
            p2_phys = solve_possion(tilde_u2_phys,tilde_v2_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,sdelta_t); % p2_{n+1};
            p3_phys = solve_possion(tilde_u3_phys,tilde_v3_phys,Inn_Lh_Neu,inv_left_Neu,coe2,V_Neu,sdelta_t); % p3_{n+1};
    
            [u1_phys,v1_phys] = solve_u(tilde_u1_phys,tilde_v1_phys,p1_phys,sdelta_t);
            [u2_phys,v2_phys] = solve_u(tilde_u2_phys,tilde_v2_phys,p2_phys,sdelta_t);
            [u3_phys,v3_phys] = solve_u(tilde_u3_phys,tilde_v3_phys,p3_phys,sdelta_t);
            
            inner1 = 1 + params.mu * 4 * sdelta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u2_phys,tilde_v2_phys,N) +...
                inner_product(tilde_u2_phys,tilde_u2_phys,N) + inner_product(tilde_v2_phys,tilde_v2_phys,N);
    
            inner2 = params.mu * 4 * sdelta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u3_phys,tilde_v3_phys,N) +...
                inner_product(tilde_u2_phys,tilde_u3_phys,N) + inner_product(tilde_v2_phys,tilde_v3_phys,N);
    
            inner3 = inner2;
    
            inner4 = 1 + params.mu * 4 * sdelta_t * inner_vector_grad(tilde_u3_phys,tilde_v3_phys,tilde_u3_phys,tilde_v3_phys,N) +...
                inner_product(tilde_u3_phys,tilde_u3_phys,N) + inner_product(tilde_v3_phys,tilde_v3_phys,N);
     
    
            A = [inner1 inner2; inner3 inner4];
    
            b1 = xi_n_1 + 2 * sdelta_t * (inner_product(udivu1,tilde_u1_phys,N) + inner_product(udivu2,tilde_v1_phys,N));
    
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
            
        
            xi(ii + 1, it) = xi_n_1;
            eta(ii + 1, it) = eta_n_1;
            tt(ii + 1, it) = TT;
            p_numerical_dt(:, :, ii) = p_phys_n;
    
    
        end
   
        error_pressure_time(it, n_p) = sum(sum((k0 * sum(p_numerical_dt_ref, 3) - sdelta_t * sum(p_numerical_dt, 3)).^2 .* (w1 * (w2'))));
        
    
        kk = 2 * kk;
        sdelta_t = sdelta_t / 2;
    
    
        u_numerical(:,:,it) = u_phys_n_1;
        v_numerical(:,:,it) = v_phys_n_1;
        p_numerical(:,:,it) = p_phys_n_1;
        xi_numerical(it) = xi_n_1;
        eta_numerical(it) = eta_n_1;

        error1 = sum(sum((u_phys_n_1 - u_numerical_ref).^2 .* (w1 * (w2')))) +...
            sum(sum((v_phys_n_1 - v_numerical_ref).^2 .* (w1 * (w2'))));
        error2 = sum(sum((p_phys_n_1 - p_numerical_ref).^2 .* (w1 * (w2'))));
        error1_relative = error1 / (sum(sum((u_numerical_ref).^2 .* (w1 * (w2')))) +...
            sum(sum((v_numerical_ref).^2 .* (w1 * (w2')))));
        error2_relative = error2 / sum(sum((p_numerical_ref).^2 .* (w1 * (w2'))));
    

        fprintf('n_p=%d, it=%d, delta t=%d, xi=%.4f, eta=%.4f,\n error_u=%d, error_u_relative=%d, error_p=%d, error_p_relative=%d\n', n_p, it, sdelta_t, xi_n_1, eta_n_1, error1, error1_relative, error2, error2_relative)


    
    end
    
    for i = 1 : n_dt
        error_velocity_L2(i, n_p) = sum(sum((u_numerical(:,:,i) - u_numerical_ref).^2 .* (w1 * (w2')))) +...
            sum(sum((v_numerical(:,:,i) - v_numerical_ref).^2 .* (w1 * (w2'))));
        error_pressure_L2(i, n_p) =  sum(sum((p_numerical(:,:,i) - p_numerical_ref).^2 .* (w1 * (w2'))));
        xi_mean(i, n_p) = xi_numerical(i);
        eta_mean(i, n_p) = eta_numerical(i);
    
    end
 

end

Error_velocity_L2 = mean(error_velocity_L2, 2);
Error_sqrt_velocity_L2 = sqrt(Error_velocity_L2);
error_pressure_L2 = mean(error_pressure_time, 2);
Error_pressure_L2 = sqrt(error_pressure_L2);


save('Error_velocity_L2.mat', 'Error_velocity_L2')
save('Error_sqrt_velocity_L2.mat', 'Error_sqrt_velocity_L2')
save('Error_pressure_L2.mat', 'Error_pressure_L2')
save('xi_ref.mat', 'xi_numerical_ref_mean')
save('eta_ref.mat', 'eta_numerical_ref_mean')
save('xi.mat', 'xi_mean')
save('eta.mat', 'eta_mean')

