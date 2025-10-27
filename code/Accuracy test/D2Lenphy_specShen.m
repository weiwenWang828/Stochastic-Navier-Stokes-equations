function v=D2Lenphy_specShen(F,var)
%2 Dim  if we compare some examples that ledisctran will not work good
N=size(F,1);[X,W]=legslb(N);%%%% THIS error make others diffcult
%%%% returns n Legendre-Gauss-Lobatto points with x(1)=-1, x(n)=1 [X,W]=legslb(N);
%%%%  if Legendre-Gauss points  then [X,W]=legs(N);
%FDLT(physics-spectral) var=1 or BDLT(spectral-physics) var=2 ; 

%CLA=Neu 表示Neumann齐次边界条件，CLA=Dir 表示Dirichlet齐次边界条件，  
  if var==2
FV=ledisctran(N,X,W,F,1);
v=(ledisctran(N,X,W,FV',1));
  elseif var==1
FV=ledisctran(N,X,W,F,0);
v=ledisctran(N,X,W,(FV'),0);
  end
end