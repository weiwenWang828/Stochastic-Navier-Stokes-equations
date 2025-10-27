function [M]=basisMNDN(N)
%  [M]=basisMND(N) returns the values between the basis of inner product(2d)
%  basis: homogeneous Neumann boundary condition.The frist derivate L_{k}-{k(k+1)/((k+2)(k+3))}L_{k+2}.
B=zeros(N-1);
for i=2:N-1
        i1=i-1; b=0;  n=i1-1;   
       while n>=0
        b=2*(n)+1+b;
        n=n-2;
        end
    for j=i:N-1 
        j1=j-1;
       %if (mod(i+j,2)==0)&&(j>i)
       % B(i,j)=(1-(i1*(i1+1))/((i1+2)*(i1+3)))*(1-(j1*(j1+1))/((j1+2)*(j1+3)))*2*b-((i1*(i1+1))/((i1+2)*(i1+3)))*(1-(j1*(j1+1))/((j1+2)*(j1+3)))*2*(2*(i1+1)+1);  %奇偶一致。
       if (i==j)
         B(i,j)=(1-(i1*(i1+1))/((i1+2)*(i1+3)))*(1-(j1*(j1+1))/((j1+2)*(j1+3)))*2*b+((i1*(i1+1))/((i1+2)*(i1+3)))^2*2*(2*(i1+1)+1);  %奇偶一致。
        end   
    end
end
M=(B+B')-diag(diag(B));
end