function [ X ] = findxr( r,ratio, g,d00,lam )
% This function searches for an adjustment multiplicator for the return
% vector which generates the out-of-sample CE difference for 
% 1. d00 H0 (out-of-sample CE difference implied from the data)
% 2. d00 + 1% annualized difference
% 3. d00 + 5% annualized difference

% Input: r = matrix of returns
%        ratio = N/T
%        g = risk aversion coefficient
%        d00 = theoretical out-of-sample CE difference
%        lam = shrinkage parameter for Ridge Covariance Matrix estimator
% Output X: 3x1 vector of adjustment coefficients

% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018


GMVP=@(x) (inv(x)*ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));

N = size(r,2);
l = ones(N,1);
wn = (1/N)*l;
T = round(N.*(1./ratio));

%c =0.001:0.001:1.5; lc = length(c);
c =-2:0.001:2; lc = length(c);

d0 = zeros(lc,1);

parfor i = 1:lc
    
    mu0 = mean(r);
    rm = r-repmat(mu0,size(r,1),1);
    
    rm(:,3) = rm(:,3).*c(i);
    rm = rm +repmat(mu0,size(r,1),1);
    
    mu = mean(rm); S = cov(rm);
    Sr =  lam*eye(N) + (1-lam)*S;
    ce1 = GMVP(Sr)'*mu'-g*(GMVP(Sr)'*S*GMVP(Sr));
    ce2 = wn'*mu'-g*(wn'*S*wn);
    
    RR = inv(S)-(inv(S)*l*l'*inv(S))/(l'*inv(S)*l);
    R = RR./(l'*inv(S)*l);
    vw =  (1/(T-N-1))*R;
    
    %theoretical differences
    d0(i) = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;
    
end


AT = [d00, ((1+d00)^(12)+0.01)^(1/12)-1, ((1+d00)^(12)+0.05)^(1/12)-1];
fx = zeros(3,1);
for j = 1:3
    
    at = AT(j);
    [~, fx(j)] = min(abs(d0-at));
    
end
X = c(fx);
end

