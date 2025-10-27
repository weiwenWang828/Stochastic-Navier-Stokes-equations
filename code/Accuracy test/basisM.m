function [M]=basisM(N)
%  [M]=basisM(N) returns the values between the basis of inner product(2d)
%  basis: homogeneous Dirichlet boundary condition. L_{k}-L_{k+2}.
B=zeros(N-1);
for i=1:N-3
    B(i,i)=2/(2*(i-1)+1)+2/(2*(i+1)+1);
    B(i,i+2)=-2/(2*(i+1)+1);
end
B(N-2,N-2)=2/(2*(N-3)+1)+2/(2*(N-1)+1);B(N-1,N-1)=2/(2*(N-2)+1)+2/(2*(N)+1);
M=(B+B')-diag(diag(B));
end