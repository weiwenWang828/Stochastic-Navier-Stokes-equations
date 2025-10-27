%--SAV for SDE-----
% dX=-nabla V(x)dt+\sigma(x) dB(t)---

function [x1T,x2T,xiT,xi1T,etaT]=SAV(N)
%-----One SAV for the nonlinear term-----

sigma=@(x) sqrt(2);
T=20;
nablaV=@(x) x.^3-x; 
tau=T/N;
x0=rand(1)*10;
W=randn(N,1);

X=zeros(N+1,1);  
X1=X;X1(1)=x0; 
xi=X; xi(1)=1;

X2=zeros(N+1,1);
X2(1)=x0;
xi1=xi; eta=xi;
for k=1:N
 %---for 1 SAV---
 deltaW=W(k)*sqrt(tau);
 xi(k+1)=(xi(k)+sigma(X1(k))*nablaV(X1(k))*deltaW*tau+tau^2*nablaV(X2(k))^2)/(1+tau^2*nablaV(X1(k))^2);
 %xi(k+1)=(xi(k)+p1)/(1+p1+tau^2*nablaV(X1(k))^2);
 X1(k+1)=X1(k)-xi(k+1)*nablaV(X1(k))*tau+sigma(X1(k))*deltaW;

 %--for 2 SAVs----
 % mean-reverting parameter
 %p2=abs(-sigma(X2(k))*tau*deltaW*nablaV(X2(k)));

  A=[1+tau^2*nablaV(X2(k))^2  -sigma(X2(k))*tau*deltaW*nablaV(X2(k));
      -sigma(X2(k))*tau*nablaV(X2(k))*deltaW  1+(sigma(X2(k)))^2*(deltaW)^2];
  b=[xi1(k)+tau^2*nablaV(X2(k))^2 eta(k)+sigma(X2(k))^2*(deltaW)^2]';
  %u=A\b;
  u(1)=(A(2,2)*b(1)-A(1,2)*b(2))/(A(1,1)*A(2,2)-A(1,2)^2);
  u(2)=(A(1,1)*b(2)-A(1,2)*b(1))/(A(1,1)*A(2,2)-A(1,2)^2);
 xi1(k+1)=u(1); eta(k+1)=u(2); 
 %xi1(k+1)=(xi1(k)+sigma(X2(k))*nablaV(X2(k))*deltaW*tau+p2)/(1+p2+tau^2*nablaV(X2(k))^2);
 %eta(k+1)=(eta(k)+sigma(X2(k))*nablaV(X2(k))*deltaW*tau+p2)/(1+p2+(sigma(X2(k)))^2*(deltaW)^2);
 X2(k+1)=X2(k)-xi1(k+1)*nablaV(X2(k))*tau+eta(k+1)*sigma(X2(k))*deltaW;
end

x1T=X1(N+1);
x2T=X2(N+1);
xiT=xi(N+1);
xi1T=xi1(end);
etaT=eta(end);
end