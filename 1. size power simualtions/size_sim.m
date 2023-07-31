% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018

% This code generates null rejection probablities over a Monte-Carlo study
% for testing CE differences

close all; clear all; clc; warning off;

spmd
    warning('off');
end

% functions
GMVP=@(x) (inv(x)*ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));
CE = @(x,gamma) (mean(x)-gamma*var(x));
dif_ce = @(x,gamma) (CE(x(:,1),gamma)-CE(x(:,2),gamma));

load ('30.mat'); % Kenneth French data on monthly returns of 30 industry portfolios
r = data(:,2:end);

MC = 50000; % number of Monte-Carlo repetitions
B = 5000; % number of bootstrap iterations
alpha = 0.05; % nominal level for testing
HH =[100,500,1000]; % out-of-sample horizon
ratio = 0.25; % N/T ratio, to make things comparable

g = 0.5; % gamma/2 - risk aversion parameter for CE

% Parameters from the data set for an emoirical Monte Carlo
S = cov(r); % covariance matrix
mu = mean(r); % mean vector
N = size(r,2); % number of assets
l = ones(N,1); % iota
wn = (1/N)*l; % equally weighted portfolio weights
T = round(N.*(1./ratio)); % in-sample estimation window length

% Variance of GMVP weigths
RR = inv(S)-(S\l*l'/S)/(l'/S*l);
R = RR./(l'/S*l);
vw =  (1/(T-N-1))*R;
% Difference under H0
ce1 = GMVP(S)'*mu'-g*(GMVP(S)'*S*GMVP(S));
ce2 = wn'*mu'-g*(wn'*S*wn); 
d00 = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2; %H0 out-of-sample CE difference

% see function findx.m
x = findx( r,ratio, g,d00 );lx = length(x);
%space holder
d0 = zeros(lx,1); dhat = zeros(3,length(HH)); bias = dhat;
dm2 = dhat; dm1 = dhat;
bs2 = dhat; bs1 = dhat;
ba2 = dhat; ba1 = dhat;

for i = 1:3 % under 1. H0 2. H0+1% 3. H0+5%
   
    % size:
    mu0 = mean(r);
    rm = r-repmat(mu0,size(r,1),1);
    rm(:,3) = rm(:,3).*x(i);
    rm = rm +repmat(mu0,size(r,1),1);
    mu = mean(rm); S = cov(rm);
    
    ce1 = GMVP(S)'*mu'-g*(GMVP(S)'*S*GMVP(S));
    ce2 = wn'*mu'-g*(wn'*S*wn);
    
    RR = inv(S)-(inv(S)*l*l'*inv(S))/(l'*inv(S)*l);
    R = RR./(l'*inv(S)*l);
    vw =  (1/(T-N-1))*R;
    
    %theoretical differences
    d0(i) = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;
    
    for hh = 1:length(HH) % out-of-sample evaluation window
        
        H = HH(hh);
        dce = zeros(MC,1);
        t = dce; tb = dce;
        q2 = zeros(MC,2);   q1 = zeros(MC,1);
        qa2 = q2;  qa1 = q1;
        
        parfor m = 1:MC
            rr = mvnrnd(mu,S,T);
            wg = GMVP(cov(rr));
            rsim = rsimul(mu,S,H,wg,wn); % out-of-sample portfolio returns
            % Estimated CE difference
            dce(m) = dif_ce(rsim,g);
            sdce = deltalw_ce(g,rsim);
            t(m) = (dce(m)-(d00))/sdce;
            
            J = ceil(rand(H,B)*H); %bootstrapped index
            tb = zeros(B,1);bdif = tb;
            for b = 1:B
                br = rsim(J(:,b),:);
                bdif(b) = dif_ce(br,g);
                tb(b) = (bdif(b)-dce(m))/deltalw_ce(g,br);
            end
            % Bootstrapped quantiles
            q2(m,:) = quantile((bdif-dce(m)), [1-alpha/2,alpha/2]); 
            q1(m) = quantile((bdif-dce(m)), 1-alpha);
            qa2(m,:) = quantile(tb, [1-alpha/2,alpha/2]); 
            qa1(m) = quantile(tb, 1-alpha);
        end
        
        dhat(i,hh) = mean(dce);
        bias(i,hh) = dhat(i,hh)-d0(i);
        % Delta method
        dm2(i,hh) = mean(abs(t)>norminv(1-alpha/2,0,1));
        dm1(i,hh) = mean(t>norminv(1-alpha,0,1));
        
        % Power curve: bootstrap-t interval
        bs2(i,hh) = mean(((t)>qa2(:,1))|(t)<qa2(:,2));
        bs1(i,hh) = mean(t>qa1);
        
        % Power curve: percentile bootstrap 
        ba2(i,hh) = mean(((dce-d00)>q2(:,1))|(dce-d00)<q2(:,2));
        ba1(i,hh) = mean((dce-d00)>q1);
        
    end
    
end


