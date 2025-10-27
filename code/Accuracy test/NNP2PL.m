function [B]=NNP2PL(N)
%  [B]=NNP2PL(N) returns the inner product of   (L_{k}- L_{k+2}) and L_{k}
B=zeros(N-1,N+1);
for i=1:N-1
    B(i,i)=2/(2*(i-1)+1);B(i,i+2)=-2/(2*(i+1)+1);
end
end