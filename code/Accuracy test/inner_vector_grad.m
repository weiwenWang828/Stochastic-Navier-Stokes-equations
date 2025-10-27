function [inner] = inner_vector_grad(u1, v1, u2, v2, N)

B = N1basisL(N);
u1_sptr_L = D2Lenphy_specShen(u1, 1); u2_sptr_L = D2Lenphy_specShen(u2, 1);
v1_sptr_L = D2Lenphy_specShen(v1, 1); v2_sptr_L = D2Lenphy_specShen(v2, 1);
u1_x_sptr = LDTL(u1_sptr_L); u1_y_sptr = TLLDY(u1_sptr_L); 
u2_x_sptr = LDTL(u2_sptr_L); u2_y_sptr = TLLDY(u2_sptr_L); 
v1_x_sptr = LDTL(v1_sptr_L); v1_y_sptr = TLLDY(v1_sptr_L); 
v2_x_sptr = LDTL(v2_sptr_L); v2_y_sptr = TLLDY(v2_sptr_L); 
inner = B*((u1_x_sptr.*(u2_x_sptr))*B') + B*((u1_y_sptr.*(u2_y_sptr))*B') +...
    B*((v1_x_sptr.*(v2_x_sptr))*B') + B*((v1_y_sptr.*(v2_y_sptr))*B');

end


function [B]=N1basisL(N)
%  [M]=basisLD(N) 
%  basis: the inner product of L_{k}.

B=zeros(1,N);
for i=1:N
      B(1,i)=2/(2*(i-1)+1);
end
end