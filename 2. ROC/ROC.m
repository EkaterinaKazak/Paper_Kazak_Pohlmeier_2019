% The notation follows Kazak and Pohlmeier (2018)
% This version 09/2018

% This code generates ROC curves 

close all; clear all; clc; warning off;
spmd
    warning('off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             N/T ratio                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% functions
GMVP=@(x) (inv(x)*ones(size(x,1),1))/(ones(1,size(x,1))*inv(x)*ones(size(x,1),1));
CE = @(x,gamma) (mean(x)-gamma*var(x));
dif_ce = @(x,gamma) (CE(x(:,1),gamma)-CE(x(:,2),gamma));

load ('30.mat');
r = data(:,2:end);

MC = 50000;

H = 1000;
g = 0.5; % gamma/ 2
ratio = 0.01; %N/T ratio, to make things comparable

S = cov(r);
N = size(r,2);
l = ones(N,1);
wn = (1/N)*l;
T = round(N.*(1./ratio));

RR = inv(S)-(S\l*l'/S)/(l'/S*l);
R = RR./(l'/S*l);
vw =  (1/(T-N-1))*R;
mu = mean(r);
ce1 = GMVP(S)'*mu'-g*(GMVP(S)'*S*GMVP(S));
ce2 = wn'*mu'-g*(wn'*S*wn);
d00 = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;

x = findx( r,ratio, g,d00); lx = length(x);

A = 0:0.01:1; la = length(A);
d0 = zeros(lx,1); dhat = zeros(3,1);
dm2 = zeros(3,la); dm1 = dm2;

for i = 1:2
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
    
    dce = zeros(MC,1);
    t = dce; 
    tic
    parfor m = 1:MC
        
        wg = sim_wgmvp( S,T );
        rsim = rsimul(mu,S,H,wg,wn);
        
        dce(m) = dif_ce(rsim,g);
        sdce(m) = deltalw_ce(g,rsim);
        t(m) = (dce(m)-(d00))/sdce(m); 
    end
    toc
    dhat(i) = mean(sdce);
    % Delta method
    for ia = 1:la
        alpha = A(ia);
        dm2(i,ia) = mean(abs(t)>norminv(1-alpha/2,0,1));
        dm1(i,ia) = mean(t>norminv(1-alpha,0,1));
    end
    
    
end


s_nt01 = dm1(1,:); p_nt01 =  dm1(2,:); 
J_nt01 = max(p_nt01-s_nt01);
s_2 = dm2(1,:); p_2 =  dm2(2,:); 
J_2 = max(p_2-s_2);


ratio = 0.1; %N/T ratio, to make things comparable

S = cov(r);
N = size(r,2);
l = ones(N,1);
wn = (1/N)*l;
T = round(N.*(1./ratio));

RR = inv(S)-(S\l*l'/S)/(l'/S*l);
R = RR./(l'/S*l);
vw =  (1/(T-N-1))*R;
mu = mean(r);
ce1 = GMVP(S)'*mu'-g*(GMVP(S)'*S*GMVP(S));
ce2 = wn'*mu'-g*(wn'*S*wn);
d00 = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;

x = findx( r,ratio, g,d00); lx = length(x);

A = 0:0.01:1; la = length(A);
d0 = zeros(lx,1); dhat = zeros(3,1);
dm2 = zeros(3,la); dm1 = dm2;

for i = 1:2
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
    
    dce = zeros(MC,1);
    t = dce; 
    tic
    parfor m = 1:MC
        
        wg = sim_wgmvp( S,T );
        rsim = rsimul(mu,S,H,wg,wn);
        
        dce(m) = dif_ce(rsim,g);
        sdce(m) = deltalw_ce(g,rsim);
        t(m) = (dce(m)-(d00))/sdce(m); 
    end
    toc
    dhat(i) = mean(sdce);
    % Delta method
    for ia = 1:la
        alpha = A(ia);
        dm2(i,ia) = mean(abs(t)>norminv(1-alpha/2,0,1));
        dm1(i,ia) = mean(t>norminv(1-alpha,0,1));
    end
    
    
end

s_nt1 = dm1(1,:); p_nt1 =  dm1(2,:); 
J_nt1 = max(p_nt1-s_nt1);


figure(1)
plot(s_nt01,p_nt01,'r','linewidth',7)
hold on
plot(s_nt1,p_nt1,'k','linewidth',7)
hold on
plot(s_nt01,s_nt01,'k:','linewidth',3)
%legend({'NT = 0.01 (J = ',num2str(J_nt01),')N/T = 0.1'},'location','southeast')
legend(strcat('N/T = 0.01  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('N/T = 0.1    J = ',num2str(num2str(round(J_nt1,3)))),'location','southeast')
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Estimation Noise')
set(gca,'FontSize',50)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Side                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(s_nt01,p_nt01,'r','linewidth',7)
hold on
plot(s_2,p_2,'k','linewidth',7)
hold on
plot(s_nt01,s_nt01,'k:','linewidth',3)
legend(strcat('One-sided  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('Two-sided  J = ',num2str(num2str(round(J_2,3)))),'location','southeast')

%legend('One-sided','Two-sided','location','southeast')
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Side of the test')
set(gca,'FontSize',50)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Horizon                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 500;
g = 0.5; % gamma/ 2
ratio = 0.01; %N/T ratio, to make things comparable

S = cov(r);
N = size(r,2);
l = ones(N,1);
wn = (1/N)*l;
T = round(N.*(1./ratio));

RR = inv(S)-(S\l*l'/S)/(l'/S*l);
R = RR./(l'/S*l);
vw =  (1/(T-N-1))*R;
mu = mean(r);
ce1 = GMVP(S)'*mu'-g*(GMVP(S)'*S*GMVP(S));
ce2 = wn'*mu'-g*(wn'*S*wn);
d00 = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;

x = findx( r,ratio, g,d00); lx = length(x);

A = 0:0.01:1; la = length(A);
d0 = zeros(lx,1); dhat = zeros(3,1);
dm2 = zeros(3,la); dm1 = dm2;

for i = 1:2
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
    
    dce = zeros(MC,1);
    t = dce; 
    tic
    parfor m = 1:MC
        
        wg = sim_wgmvp( S,T );
        rsim = rsimul(mu,S,H,wg,wn);
        
        dce(m) = dif_ce(rsim,g);
        sdce(m) = deltalw_ce(g,rsim);
        t(m) = (dce(m)-(d00))/sdce(m); 
    end
    toc
    dhat(i) = mean(sdce);
    % Delta method
    for ia = 1:la
        alpha = A(ia);
        dm2(i,ia) = mean(abs(t)>norminv(1-alpha/2,0,1));
        dm1(i,ia) = mean(t>norminv(1-alpha,0,1));
    end
    
    
end


s_h = dm1(1,:); p_h =  dm1(2,:); 
J_h = max(p_h-s_h);

figure(3)
plot(s_nt01,p_nt01,'r','linewidth',7)
hold on
plot(s_h,p_h,'k','linewidth',7)
hold on
plot(s_nt01,s_nt01,'k:','linewidth',3)
%legend('H = 1000','H = 500','location','southeast')
legend(strcat('H = 1000  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('H = 500    J = ',num2str(num2str(round(J_h,3)))),'location','southeast')

xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Out-of-sample horizon')
set(gca,'FontSize',50)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          benchmark                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rati = 0.01;
H = 1000;

%[~,lam] = covMarket(r,-1);
lam = 0.007;%findlam(r,ratio,g);%0.0078;
S = cov(r);
Sr = lam*eye(N) + (1-lam)*S;
N = size(r,2);
l = ones(N,1);
wn = (1/N)*l;
T = round(N.*(1./ratio));

RR = inv(S)-(S\l*l'/S)/(l'/S*l);
R = RR./(l'/S*l);
vw =  (1/(T-N-1))*R;
mu = mean(r);
ce1 = GMVP(Sr)'*mu'-g*(GMVP(Sr)'*S*GMVP(Sr));
ce2 = wn'*mu'-g*(wn'*S*wn);
d00 = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;

x = findxr( r,ratio, g,d00,lam); lx = length(x);
x(1) = 1;
A = 0:0.01:1; la = length(A);
d0 = zeros(lx,1); dhat = zeros(3,1);
dm2 = zeros(3,la); dm1 = dm2; lxx = dhat;

for i = 1:2
    % size:
    mu0 = mean(r);
    rm = r-repmat(mu0,size(r,1),1);
    rm(:,3) = rm(:,3).*x(i);
    rm = rm +repmat(mu0,size(r,1),1);
    mu = mean(rm); S = cov(rm);%eye(N)*lam + (1-lam)*cov(rm);% c
    
    Sr = lam*eye(N) + (1-lam)*S;
    ce1 = GMVP(Sr)'*mu'-g*(GMVP(Sr)'*S*GMVP(Sr));
    ce2 = wn'*mu'-g*(wn'*S*wn);
    
    RR = inv(S)-(inv(S)*l*l'*inv(S))/(l'*inv(S)*l);
    R = RR./(l'*inv(S)*l);
    vw =  (1/(T-N-1))*R;
    
    %theoretical differences
    d0(i) = ce1-g*trace(S*vw)-g*mu*vw*mu'-ce2;
    
    dce = zeros(MC,1);
    t = dce; 
    tic
    %Sx = lam*eye(N) + (1-lam)*S;
    parfor m = 1:MC
        
        wg = sim_wgmvp( Sr,T );
        rsim = rsimul(mu,S,H,wg,wn);
        
        dce(m) = dif_ce(rsim,g);
        sdce(m) = deltalw_ce(g,rsim);
        t(m) = (dce(m)-(d00))/sdce(m); 
    end
    toc
    dhat(i) = mean(sdce);
    % Delta method
    for ia = 1:la
        alpha = A(ia);
        dm2(i,ia) = mean(abs(t)>norminv(1-alpha/2,0,1));
        dm1(i,ia) = mean(t>norminv(1-alpha,0,1));
    end
    
    
end
rs2 = dm2(1,:); rp2 = dm2(2,:); 
rs1 = dm1(1,:); rp1 = dm1(2,:); 
J_r = max(rp1-rs1);

figure(4)
plot(s_nt01,p_nt01,'r','linewidth',7)
hold on
plot(rs1,rp1,'k','linewidth',7)
hold on
plot(s_nt01,s_nt01,'k:','linewidth',3)
%legend('GMVP vs 1/N','Ridge vs 1/N','location','southeast')
legend(strcat('GMVP vs 1/N  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('Ridge vs 1/N   J = ',num2str(num2str(round(J_r,3)))),'location','southeast')

xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Benchmark')
set(gca,'FontSize',50)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,2,1)
plot(s_nt01,p_nt01,'r','linewidth',5)
hold on
plot(s_nt1,p_nt1,'k:','linewidth',5)
hold on
plot(s_nt01,s_nt01,'k--','linewidth',3)
%legend({'NT = 0.01 (J = ',num2str(J_nt01),')N/T = 0.1'},'location','southeast')
legend(strcat('N/T = 0.01  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('N/T = 0.1    J = ',num2str(num2str(round(J_nt1,3)))),'location','southeast')
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Estimation Noise')
set(gca,'FontSize',24)

subplot(2,2,2)
plot(s_nt01,p_nt01,'r','linewidth',5)
hold on
plot(s_2,p_2,'k:','linewidth',5)
hold on
plot(s_nt01,s_nt01,'k--','linewidth',3)
legend(strcat('One-sided  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('Two-sided  J = ',num2str(num2str(round(J_2,3)))),'location','southeast')

%legend('One-sided','Two-sided','location','southeast')
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Side of the test')
set(gca,'FontSize',24)


subplot(2,2,3)
plot(s_nt01,p_nt01,'r','linewidth',5)
hold on
plot(s_h,p_h,'k:','linewidth',5)
hold on
plot(s_nt01,s_nt01,'k--','linewidth',3)
%legend('H = 1000','H = 500','location','southeast')
legend(strcat('H = 1000  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('H = 500    J = ',num2str(num2str(round(J_h,3)))),'location','southeast')

xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Out-of-sample horizon')
set(gca,'FontSize',24)

subplot(2,2,4)
plot(s_nt01,p_nt01,'r','linewidth',5)
hold on
plot(rs1,rp1,'k:','linewidth',5)
hold on
plot(s_nt01,s_nt01,'k--','linewidth',3)
%legend('GMVP vs 1/N','Ridge vs 1/N','location','southeast')
legend(strcat('GMVP vs 1/N  J = ', num2str(num2str(round(J_nt01,3)))),...
    strcat('Ridge vs 1/N   J = ',num2str(num2str(round(J_r,3)))),'location','southeast')

xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$1-\hat{\beta}$','interpreter','latex')
title('Benchmark')
set(gca,'FontSize',24)
