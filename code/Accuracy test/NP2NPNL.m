function [B]=NP2NPNL(n)
%  [B]=NP2NPNL(N) returns the inner product of  L_{k} and L_{k}-{k(k+1)/((k+2)(k+3))}L_{k+2} 
N=n-1;
B=zeros(N+1,N-1);
B(1,1)=2;
B(2,2)=2/3;
for i=3:N-1
    B(i,i)=2/(2*(i-1)+1);B(i,i-2)=-((i-3)*(i-2))/((i-1)*(i))*2/(2*(i-1)+1);
end
B(N,N-2)=-((N-3)*(N-2))/((N-1)*(N))*2/(2*(N-1)+1);
B(N+1,N-1)=-((N-2)*(N-1))/((N)*(N+1))*2/(2*(N)+1);
end