function [M]=basisMN(N)
%  [M]=basisMN(N) returns the values between the basis of inner product(2d)
%  basis: homogeneous Neumann boundary condition.L_{k}-{k(k+1)/((k+2)(k+3))}L_{k+2}.
B=zeros(N-1);
for i=1:N-3
    B(i,i)=2/(2*(i-1)+1)+(((i-1)*i)^2/((i+1)*(i+2))^2)*2/(2*(i+1)+1);
    B(i,i+2)=-(((i-1)*i)/((i+1)*(i+2)))*2/(2*(i+1)+1);
end
for i=(N-2):(N-1)
    B(i,i)=2/(2*(i-1)+1)+(((i-1)*i)^2/((i+1)*(i+2))^2)*2/(2*(i+1)+1);
end
M=(B+B')-diag(diag(B));
end