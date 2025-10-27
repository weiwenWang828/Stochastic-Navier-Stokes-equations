function [inner] = inner_product(f,g,N)
% input: value of f, g (value at (x_i, y_j))
% output: inner product
f_sptr_L = D2Lenphy_specShen(f, 1); % sptr in L_k(x)L_l(y); n x n;
g_sptr_L = D2Lenphy_specShen(g, 1);
inner = N1basisL(N)*((f_sptr_L.*(g_sptr_L))*N1basisL(N)');
end

function [B]=N1basisL(N)
%  [M]=basisLD(N) 
%  basis: the inner product of L_{k}.

B=zeros(1,N);
for i=1:N
      B(1,i)=2/(2*(i-1)+1);
end
end