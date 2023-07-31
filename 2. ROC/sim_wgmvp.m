function [ wg ] = sim_wgmvp( S,T )
% This function simulates  Nx1 vector of GMVP portfolio weights
% see Okhrin and Schmidt 2006 for details
% Input: S = covariance matrix      
%        T = estimation window length
% Output wg: Nx1 vector of GMVP weights

% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018

GMVP=@(x) (inv(x)*ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));

N = size(S,1);
l = ones(N,1);
wg0 = GMVP(S);
mu_w = wg0(1:end-1);
R = inv(S)-(inv(S)*l*l'*inv(S))/(l'*inv(S)*l);
var_w = (1/(T-N-1))*R(1:end-1,1:end-1)/(l'*inv(S)*l);

wx = mu_w +sqrt((T-N+1-2)/(T-N+1)).*chol(var_w)*mvtrnd(eye(N-1),T-N+1)';
wg = [wx;1-sum(wx)];

end
