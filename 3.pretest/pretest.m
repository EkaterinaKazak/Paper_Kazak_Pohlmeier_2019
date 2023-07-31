% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018

% This code generates the out-of-sample CEs over a 1000 randomly drawn
% portfolios of size N for 5 portfolio strategies:
% 1. GMVP
% 2. 1/N
% 3. Pretest 
% 4. Pretest + alpha shrinking
% 5. Pretest + alpha smoothing

clc; clear;

GMVP=@(x) (inv(x)*ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));
CE = @(x,gamma) (mean(x)-gamma*var(x));
dif_ce = @(x,gamma) (CE(x(:,1),gamma)-CE(x(:,2),gamma));

load('100.mat');
data = data(:,2:end);

% parameters
MC = 1000; % number of randomly drawn portfolios
g = 0.5; % gamma/2 risk aversion
T = 60; % in-sample estimation window length
W = size(data,1);
% grid of alphas for pretesting
alpha = 0.01:0.01:0.99; la = length(alpha);

NN = [5,30,50];
TAB = zeros(5,length(NN)); VV = zeros(3,length(NN));
for nn = 1:length(NN) % portfolio sizes
    
    N = NN(nn); % number of assets
    l = ones(N,1); %iota
    wn = (1/N)*l; % equally weighted portfolio weights
    %space holder
    CEmat = zeros(5,MC); p1 = zeros(MC,1); p2 = zeros(MC,1); p3 = zeros(MC,1);
    parfor m =1:MC
        
        idx = randperm(100,N); %randomly drawn assets
        r = data(:,idx);
        S = cov(r); mu = mean(r);
        
        % in-sample CE testing
        a_star =  zeros(W-T-1,1); P1 =  zeros(W-T-1,1); P2 = zeros(W-T-1,1);  P3 = zeros(W-T-1,1);
        rg = zeros(W-T-1,1); rn = zeros(W-T-1,1); ra = zeros(W-T-1,1);
        rl = zeros(W-T-1,1); rl2 = zeros(W-T-1,1);
        for t = 2:W-T
            sample = r(t:t+T-1,:);
            
            %in-sample estimation
            wg = GMVP(cov(sample));
            CEgin = mean(sample)*wg-g*wg'*cov(sample)*wg;
            CEnin = mean(sample)*wn-g*wn'*cov(sample)*wn;
            rp_in = [sample*wg, sample*wn]; % in-sample portfolio returns
            
            din = CEgin-CEnin; % in-sample CE difference
            tin = din/deltalw_ce(g,rp_in); % test statistic
            
            CEain = zeros(la,1);
            for j = 1:la
                a = alpha(j);
                I  = tin > norminv(1-a,0,1);
                wa = I*wg + (1-I)*wn;
                % CEin-sample for a grid of nominal level alpha
                CEain(j) = mean(sample)*wa-g*wa'*cov(sample)*wa;
            end
            [~,idxx] = max(CEain);
            a_star(t) = alpha(idxx); % alpha^*
            I_star = tin>norminv(1-a_star(t),0,1);
            wa = I_star*wg + (1-I_star)*wn;
            P1(t) = I_star;
            
            %shrinking alpha to nominal 5% level
            lam = 0.5;
            al = a_star(t)*(1-lam) + lam*0.05;
            I_l = tin>norminv(1-al,0,1);
            wa_l = I_l*wg + (1-I_l)*wn;
            rl(t) =  r(t+T,:)*wa_l;
            P2(t) = I_l;
            % smoothing alpha over time
            al2 = a_star(t)*(1-lam) + lam*a_star(t-1);
            I_l2 = tin>norminv(1-al2,0,1);
            wa_l2 = I_l2*wg + (1-I_l2)*wn;
            rl2(t) =  r(t+T,:)*wa_l2;
            P3(t) = I_l2;
            
            % out-of-sample returns
            rg(t) = r(t+T,:)*wg; %GMVP
            rn(t) = r(t+T,:)*wn; %1/N
            ra(t) = r(t+T,:)*wa; % Pretest
            
        end
        p1(m) = mean(P1); p2(m) = mean(P2); p3(m) = mean(P2);
        % out-of-sample CE
        CEmat(:,m) = [CE(rg,g); CE(rn,g); CE(ra,g); CE(rl,g); CE(rl2,g)];
        
    end
    % Average CE and selection probabilty
    TAB(:,nn) = mean(CEmat,2); VV(:,nn) = [mean(p1);mean(p2);mean(p3)];
    
end
%results
res.TAB = TAB;
res.VV = VV;
res.g = g;
res.T = T;







