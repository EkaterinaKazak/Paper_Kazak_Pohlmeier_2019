function [ r ] = rsimul( mu,s,H,w1,w2 )
% This function simulates out-of-sample portfolio returns
% Input: mu = mean vector
%        s = covariance matrix
%        H = out-of-sample horizon
%        w1 = vector of portfolio weights for the 1st strategy
%        w2 = vector of portfolio weights for the 2nd strategy
% Output r: Hx2 vector of portfolio returns
% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018
S = (s+s)./2;
mu_sim = [mu*w1,mu*w2];
S_s = [w1'*S*w1,w1'*S*w2;...
    w2'*S*w1,w2'*S*w2];
S_sim  = (S_s+S_s.')./2;
r = mvnrnd(mu_sim,S_sim,H);

end