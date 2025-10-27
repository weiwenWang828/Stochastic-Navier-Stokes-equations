clear;clc
format long


T = 1.2;
delta_t = 1e-4;
T_m = T / delta_t;
n_dt = 1;
KK = ceil(T_m);
N_p = 300;


params.mu = 1e-4; 
params.J = 4;
params.epsilon = 0.001;



Nx = 128; Ny = Nx;
params.hx = 2*pi/Nx;
params.hy = 2*pi/Ny;
x = (0:Nx-1)*params.hx;
y = (0:Ny-1)*params.hy;
[xx,yy] = meshgrid(x,y);
kx = [0:Nx/2-1 0 -Nx/2+1:-1]; 
ky = [0:Ny/2-1 0 -Ny/2+1:-1]; 
[k_x,k_y] = meshgrid(kx,ky);
k2 = -(k_x.^2 + k_y.^2); 
dx = 1i * k_x; 
dy = 1i * k_y;

vor_numerical4 = zeros(Nx, Ny, N_p);
vor_numerical8 = vor_numerical4;
vor_numerical12 = vor_numerical4;
xi_numerical4 = zeros(1, N_p);
xi_numerical8 = xi_numerical4;
xi_numerical12 = xi_numerical4;
eta_numerical4 = zeros(1, N_p);
eta_numerical8 = eta_numerical4;
eta_numerical12 = eta_numerical4;



for n_p = 1 : N_p

    
    [u_phys0,v_phys0,p_phys0] = initial_value(xx,yy); 

    vor_hat = dx .* fftn(v_phys0) - dy .* fftn(u_phys0);
    vor0 = real(ifftn(vor_hat));
    
    time_step = zeros(n_dt , 1);
    xi = zeros(2^(n_dt-1) * KK + 1, n_dt); tt = xi;
    eta = xi;


    for it = 1 : n_dt
        
        time_step(it , 1) = delta_t;
        TT = 0;
        u_phys_n_1 = u_phys0;    
        v_phys_n_1 = v_phys0;
        p_phys_n_1 = p_phys0;
        xi_n_1 = 1; 
        eta_n_1 = 1; 
        xi(1, it) = xi_n_1;
        eta(1, it) = eta_n_1;
        tt(1, it) = 0;
    

    
        for ii = 1 : KK
            TT = ii * delta_t;
    
            coe1 = 1 - params.mu * 4 * pi^2 * delta_t * k2;
            coe2 = 2 * pi * delta_t * k2;
            [g1, g2] = g(u_phys_n_1, v_phys_n_1, params);

          
            [p_x_phys_n_1, p_y_phys_n_1] = partial(p_phys_n_1, dx, dy);
            [udivu1, udivu2] = u_divu(u_phys_n_1, v_phys_n_1, dx, dy);
            [wiener1,wiener2] = wiener(xx,yy,delta_t,Nx,params);

    
            
            U1_right = u_phys_n_1 - 2 * pi * delta_t * p_x_phys_n_1; 
            V1_right = v_phys_n_1 - 2 * pi * delta_t * p_y_phys_n_1;
    
            tilde_u1_phys = real(ifftn(fftn(U1_right)./coe1));
            tilde_v1_phys = real(ifftn(fftn(V1_right)./coe1));
    
            U2_right = - 2 * pi * delta_t * udivu1;
            V2_right = - 2 * pi * delta_t * udivu2;
    
            tilde_u2_phys = real(ifftn(fftn(U2_right)./coe1));
            tilde_v2_phys = real(ifftn(fftn(V2_right)./coe1));
    
            U3_right = g1 .* wiener1;
            V3_right = g2 .* wiener2;
    
            tilde_u3_phys = real(ifftn(fftn(U3_right)./coe1));
            tilde_v3_phys = real(ifftn(fftn(V3_right)./coe1));
    

            p1_phys = solve_possion(tilde_u1_phys,tilde_v1_phys,dx,dy,delta_t,coe2,Nx); 
            p2_phys = solve_possion(tilde_u2_phys,tilde_v2_phys,dx,dy,delta_t,coe2,Nx); 
            p3_phys = solve_possion(tilde_u3_phys,tilde_v3_phys,dx,dy,delta_t,coe2,Nx); 
        
            [u1_phys,v1_phys] = solve_u_2pi(tilde_u1_phys,tilde_v1_phys,p1_phys,delta_t,dx,dy);
            [u2_phys,v2_phys] = solve_u_2pi(tilde_u2_phys,tilde_v2_phys,p2_phys,delta_t,dx,dy);
            [u3_phys,v3_phys] = solve_u_2pi(tilde_u3_phys,tilde_v3_phys,p3_phys,delta_t,dx,dy);
    
    
     
            inner1 = 1 + params.mu * 4 * pi^2 * delta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u2_phys,tilde_v2_phys,dx,dy,params) +...
                sum(sum(tilde_u2_phys.^2)) * params.hx * params.hy + sum(sum(tilde_v2_phys.^2)) * params.hx * params.hy;
        
            inner2 = params.mu * 4 * pi^2 * delta_t * inner_vector_grad(tilde_u2_phys,tilde_v2_phys,tilde_u3_phys,tilde_v3_phys,dx,dy,params) +...
                sum(sum(tilde_u2_phys .* tilde_u3_phys)) * params.hx * params.hy + sum(sum(tilde_v2_phys .* tilde_v3_phys)) * params.hx * params.hy;
        
            inner3 = inner2;
        
            inner4 = 1 + params.mu * 4 * pi^2 * delta_t * inner_vector_grad(tilde_u3_phys,tilde_v3_phys,tilde_u3_phys,tilde_v3_phys,dx,dy,params) +...
                sum(sum(tilde_u3_phys.^2)) * params.hx * params.hy + sum(sum(tilde_v3_phys.^2)) * params.hx * params.hy;
        
            A = [inner1 inner2; inner3 inner4];
        
            b1 = xi_n_1 + 2 * pi * delta_t * (sum(sum(udivu1 .* tilde_u1_phys)) * params.hx * params.hy +...
                sum(sum(udivu2 .* tilde_v1_phys)) * params.hx * params.hy);
        
            b2 = eta_n_1 - sum(sum(U3_right .* (tilde_u1_phys - u_phys_n_1))) * params.hx * params.hy -...
                sum(sum(V3_right .* (tilde_v1_phys - v_phys_n_1))) * params.hx * params.hy +...
                sum(sum(U3_right.^2)) * params.hx * params.hy + sum(sum(V3_right.^2)) * params.hx * params.hy;

    
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

            vor_hat = dx .* fftn(v_phys_n_1) - dy .* fftn(u_phys_n_1);
            vor = real(ifftn(vor_hat));

    
            fprintf('n_p=%d, it=%d, i=%d, time=%d, delta t=%d, xi=%.4f, eta=%.4f\n', n_p, it, ii, TT, delta_t, xi_n_1, eta_n_1)


        end

    end

end

save('vor04.mat', 'vor_numerical4')
save('xi04.mat', 'xi_numerical4')
save('eta04.mat', 'eta_numerical4')

save('vor08.mat', 'vor_numerical8')
save('xi08.mat', 'xi_numerical8')
save('eta08.mat', 'eta_numerical8')

save('vor12.mat', 'vor_numerical12')
save('xi12.mat', 'xi_numerical12')
save('eta12.mat', 'eta_numerical12')
