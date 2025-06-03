function [M,K,C] = matrices_1D(xnode,xipg,wpg,N,Nxi)
% [M,K,C] = matrices(xnode,xipg,wpg,N,Nxi)
%
% This function computes mass matrix M, diffusion matrix K and 
% convection matriz C, which are used to obtain the transient solution
%
% Input:
%   xnode : nodes' corrdinates
%   xipg, wpg : numerical quadrature
%   N, Nxi : shape functions and its derivatives on Gauss' points
%
 
% Total number of nodes and elements
npoin = size(xnode,2); 
nelem = npoin-1; 

% Number of Gauss points
ngaus = size(wpg,1);

M = spalloc(npoin, npoin,3*npoin);
K = spalloc(npoin, npoin,3*npoin);
C = spalloc(npoin, npoin,3*npoin);

% Loop on elements
for i=1:nelem
    a = xnode(i);
    b = xnode(i+1);
    h = b-a;
    weigth = wpg*h/2;
    isp = [i i+1]; % global number of the nodes of the element
    % Loop on Gauss points
    for ig = 1:ngaus
        N_ig = N(ig,:);
        Nx_ig = Nxi(ig,:)*2/h;
        w_ig = weigth(ig);
        x = (a+b)/2+h/2*xipg(ig); 
        % Assembly local contribution
        M(isp,isp) = M(isp,isp) + w_ig*N_ig'*N_ig;
        K(isp,isp) = K(isp,isp) + w_ig*Nx_ig'*Nx_ig;
%         C(isp,isp) = C(isp,isp) + w_ig*N_ig'*Nx_ig;
        C(isp,isp) = C(isp,isp) + w_ig*2*(1/x)*N_ig'*Nx_ig;
    end
 end
