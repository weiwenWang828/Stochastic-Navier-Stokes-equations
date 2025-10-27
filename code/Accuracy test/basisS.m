function [M]=basisS(N)
%  [M]=basisS(N) returns the values between the basis of inner product(2d)
%  basis: homogeneous Dirichlet boundary condition. The frist derivat of L_{k}-L_{k+2}.
j = 0 : (N-2);
M = diag(4 * j + 6);
end
