function cmp = CM_stat_BCa_bootstrap_correlation(data1,data2,data1c,data2c)

% perform a permutation test on estimated Ys
% n_cond = 4;
% n_hem = 1;
% data1 = data_RA(:,n_cond,n_hem,1);
% data2 = data_RA(:,n_cond,n_hem,2);
% data1c = data_RAc(:,n_cond,n_hem,1);
% data2c = data_RAc(:,n_cond,n_hem,2);

if nargin == 4
    Nsim = 1000;
end
Nsub = size(data1,1);


r = corr(data1,data2,'type','Spearman');
rc = corr(data1c,data2c,'type','Spearman');


dr_dist = zeros(1,Nsim);
for n_sim = 1:Nsim
    selec = ceil(rand(1,Nsub)*Nsub);
    selec(selec == 0) = 1;
    
    rr = corr(data1(selec),data2(selec),'type','Spearman');
    rrc = corr(data1c(selec),data2c(selec),'type','Spearman');
    dr_dist(n_sim) = rr-rrc;
end
dr_dist = sort(dr_dist);
p = sum(dr_dist<0)/(Nsim-sum(dr_dist == 0));
p = min(p,1-p)*2;



% BCa bootstrap (Effron 14.3, page 184)
% Jacknife estimation
dr_JN = zeros(1,Nsub);
for n = 1:Nsub
    keep = [1:n-1 n+1:Nsub];
    r_JN = corr(data1(keep),data2(keep),'type','Spearman');
    rc_JN = corr(data1c(keep),data2c(keep),'type','Spearman');
    dr_JN(n) = r_JN-rc_JN;
end


p = 0.5;
dp = 0.25;
while p/dp < 1000
    alpha = p/2;

    % bias correction
    dr = r-rc;
    z0 = norminv(mean(dr_dist < dr));

    % acceleration factor
    anum = sum((mean(dr_JN)-dr_JN).^3);
    adenom = 6*sum((mean(dr_JN)-dr_JN).^2)^(3/2);
    a = anum/adenom;

    % estimate alphas
    alpha1 = normcdf(z0 + (z0 + norminv(alpha))/(1-a*(z0 + norminv(alpha))));
    alpha2 = normcdf(z0 + (z0 + norminv(1-alpha))/(1-a*(z0 + norminv(1-alpha))));
    
    % BCa confidence interval
    CI = prctile(dr_dist,[alpha1 alpha2]*100);
    
    p = p - prod(sign(CI))*dp;
    dp = dp/2;
end


alpha = 0.025;

% bias correction
dr = r-rc;
z0 = norminv(mean(dr_dist < dr));

% acceleration factor
anum = sum((mean(dr_JN)-dr_JN).^3);
adenom = 6*sum((mean(dr_JN)-dr_JN).^2)^(3/2);
a = anum/adenom;

% estimate alphas
alpha1 = normcdf(z0 + (z0 + norminv(alpha))/(1-a*(z0 + norminv(alpha))));
alpha2 = normcdf(z0 + (z0 + norminv(1-alpha))/(1-a*(z0 + norminv(1-alpha))));

% BCa confidence interval
CI = prctile(dr_dist,[alpha1 alpha2]*100);




cmp = [];
cmp.r = r;
cmp.rc = rc;
cmp.Nsim = Nsim;
cmp.dr_dist = dr_dist;
cmp.p = p;
cmp.CI95 = CI;


% % cmp = [];
% % mse_dec1 = ((stat1.yVy-stat1.yVy_val).^2)/mean((stat1.yVy).^2,1);
% % mse_dec2 = ((stat2.yVy-stat2.yVy_val).^2)/mean((stat2.yVy).^2,1);
% % % cmp.p_sr = signrank(sqrt(mse_dec1)-sqrt(mse_dec2));
% % cmp.p_sr2 = signrank((mse_dec1)-(mse_dec2));
% % % [tmp, cmp.p_t] = ttest(sqrt(mse_dec1)-sqrt(mse_dec2));
% % [tmp, cmp.p_t2] = ttest((mse_dec1)-(mse_dec2));
% % 
% % cmp.p = cmp.p_t2;
% % % cmp.p = cmp.p_sr2;





% 
% load p_cmp_all p p_sr p_sr2 mse1 mse2 p_t2 p_t
% % p_cmp_all = []; p = []; p_sr = []; p_sr2 = []; mse1 = []; mse2 = []; p_t = []; p_t2 = [];
% 
% p(end+1) = cmp.p;
% p_sr(end+1) = cmp.p_sr;
% p_sr2(end+1) = cmp.p_sr2;
% p_t(end+1) = cmp.p_t;
% p_t2(end+1) = cmp.p_t2;
% mse1(end+1) = mean(mse_dec1);
% mse2(end+1) = mean(mse_dec2);
% 
% save p_cmp_all p p_sr p_sr2 mse1 mse2 p_t2 p_t






% 
% 
% subplot(2,2,1:2); hist(dmse_dist,20);
% subplot(2,2,3); plot(sqrt(mse_dec1)-sqrt(mse_dec2),'o')
% subplot(2,2,4); plot((mse_dec1)-(mse_dec2),'o')
% subplot(2,2,3); plot(sqrt(mse_dec1),sqrt(mse_dec2),'o')
% subplot(2,2,4); plot((mse_dec1),(mse_dec2),'o')
% mean(sqrt(mse_dec1)-sqrt(mse_dec2))
% std(sqrt(mse_dec1)-sqrt(mse_dec2))/sqrt(Nsub)
% 
% cmp
% 
% pause






% 
% 
% 
% for n = 1:1000
%     tmp1 = randn(Nsub,1);
%     tmp2 = randn(Nsub,1);
%     p_sr(n) = signrank(tmp1-tmp2);
%     
%     
% 
%     dmse_dist = zeros(1,Nsim);
%     for n_sim = 1:Nsim
%         selec = ceil(rand(1,Nsub)*Nsub);
%         selec(selec == 0) = 1;
% 
%         dmse_dist(n_sim) = mean(tmp1(selec))-mean(tmp2(selec));
%     end
%     p(n) = sum(dmse_dist<0)/Nsim;
%     p(n) = min(p(n),1-p(n))*2;
% 
% end





% % % 
% % % if nargin == 2
% % %     Nsim = 1000;
% % % end
% % % Nsub = size(stat1.Y,2);
% % % 
% % % Y1 = stat1.Vy'*stat1.Y; Y1 = Y1/std(Y1);
% % % Y2 = stat2.Vy'*stat2.Y; Y2 = Y2/std(Y2);
% % % X1 = stat1.Vx'*stat1.X; X1 = X1/std(X1);
% % % X2 = stat2.Vx'*stat2.X; X2 = X2/std(X2);
% % % 
% % % dr_dist = zeros(1,Nsim);
% % % for n_sim = 0:Nsim    
% % %     r1 = Y1*X1'/(Nsub-1);
% % %     r2 = Y2*X2'/(Nsub-1);
% % %     
% % %     if n_sim
% % %         dr_dist(n_sim) = r1 - r2;
% % %     else
% % %         dr = r1 - r2;
% % %     end
% % %     
% % %     perm = find(randn(1,Nsub)>0);
% % %     [Y1(:,perm),Y2(:,perm)] = deal(Y2(:,perm),Y1(:,perm));
% % % end
% % % p = mean(abs(dr)<abs(dr_dist));
% % % 
% % % cmp = [];
% % % cmp.Nsim = Nsim;
% % % cmp.dr = dr;
% % % cmp.dr_dist = dr_dist;
% % % cmp.p = p;
