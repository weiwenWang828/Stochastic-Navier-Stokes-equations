function [B]=phiTL(U)
%  [B]=phiTL(U) change the basis L_{k}-L_{k+2} to basis L.
size2=size(U);n=size2(1,1);
A=zeros(n+2,n);
for i=1:n
A(i,i)=1;
A(i+2,i)=-1;
end
B=A*U*A';
   