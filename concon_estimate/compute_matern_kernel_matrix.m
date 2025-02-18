function KM = compute_matern_kernel_matrix(nu, kappa, Lambda, U)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
%   compute graph matern kernel matrix for one hemisphere

% INPUT:
%   nu: differentiability parameter, nu>0
%   kappa: kernel bandwidth
%   Lambda: maxk x 1 vector of graph Laplacian eigenvalues, 
%           maxk is the truncation number
%   U: nnds x maxk matrix, columns are eigenvectors of graph Laplacian, 
%      nnds is the number of nodes on the hemisphere

% OUTPUT:
%   KM: nnds x nnds matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of nodes
nnds = size(U,1);

% compute kernel eigenvalues
rho = ( 2*nu/kappa^2 + Lambda ).^(-nu-1); % maxk x 1
rho(rho < 1e-30) = 0;

% compute kernel matrix
KM = ( U.* repmat(rho',[nnds,1]) ) * U'; 
