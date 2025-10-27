function [A]=TLLDY(U)
%  [A]=TLLDY(U) change the basis L' to basis L.
size2=size(U);n=size2(1,1);
A=zeros(n,n);
A(:,n-1)=(2*(n-2)+1).*U(:,n);
A(:,n-2)=(2*(n-3)+1).*U(:,n-1);
  for i=(n-3):-2:1
    A(:,i)=(2*(i-1)+1)*(A(:,i+2)/(2*(i+1)+1)+U(:,i+1));   
  end
  for i=(n-4):-2:1
    A(:,i)=(2*(i-1)+1)*(A(:,i+2)/(2*(i+1)+1)+U(:,i+1));   
  end 
  