% generate the data
X = log10(subcoh(:,:)');
if 0 % data normalized by the CVC without noise
    X(1:8,:) = X(1:8,:)./repmat(X(9,:),8,1);
    X(10:17,:) = X(10:17,:)./repmat(X(18,:),8,1);
end
if 0 % data normalized by the CVC without noise
    X = mean(X([9 18],:));
end
% Y = randn(5,size(X,2));
% Y(3,:) = Y(3,:) + [0 16 28 42 12 35 0 0 0 0 0 0 0 0 0 0 0 0]*X;

left_handers = [16 20 21 28 36];
Y = scores_37_sub;

X(:,left_handers) = [];
Y(:,left_handers) = [];



% put the data into a single matrix
Nx = size(X,1); x = 1:Nx;
Ny = size(Y,1); y = Nx+1:Nx+Ny;
L = size(X,2);
data = cat(1,X,Y);

% standardize the data
Mdata = mean(data,2);
data = data-repmat(Mdata,size(data(1,:)));
SDdata = std(data,[],2);
data = diag(1./SDdata)*data;

% Define test and train subjects
Ltest = floor(L/10);
Ltest = 10; test = 1:Ltest;
Ltrain = L-Ltest; train = Ltest+1:Ltest+Ltrain;
Nrun = 50;
[tmp,tt] = sort(rand(Nrun,L),2); % tt stands for test and train

Nsim = 100;

% initialize the variables and start the main loop
alphax = 2.^(5:-1:-10);
alphay = 2.^(5:-1:-10);
CCAeval = zeros(1,Nsim+1);
CCAevalruns = zeros(Nrun,Nsim+1);
for perm_stat = 0:100
    disp(['Processing perm_stat = ' num2str(perm_stat) '/' num2str(Nsim)])
    if perm_stat
        order = randperm(L);
        data(y,:) = data(y,order);
    end
    
    allCCA = zeros(Nrun,length(alphax));
    allVx = zeros(Nx,Nrun,length(alphax));
    allVy = zeros(Ny,Nrun,length(alphax));
    for n_run = 1:Nrun
        % training of the CCA model
        int = tt(n_run,train);
        Ctrain = data(:,int)*data(:,int)'/length(int);
        int = tt(n_run,test);
        Ctest = data(:,int)*data(:,int)'/length(int);

        for n_alpha = 1:length(alphax)
            % Estime the regularized cov matrices
            C = Ctrain;
            Cxxreg = C(x,x)+alphax(n_alpha)*eye(Nx);
            Cyyreg = C(y,y)+alphay(n_alpha)*eye(Ny);
            
            % Estimate CCA and coef vector for x
            [Vx,D] = eig(Cxxreg\C(x,y)*(Cyyreg\C(x,y)'));
            [M,n_keep] = max(abs(diag(D)));
            Vx = Vx(:,n_keep);
            
            % Estimate CCA and coef vector for y
            [Vy,D] = eig(Cyyreg\C(x,y)'*(Cxxreg\C(x,y)));
            [M,n_keep] = max(abs(diag(D)));
            Vy = Vy(:,n_keep);
            
            % estimate whether V has to be flipped (positive or negative
            % correlation for the training data)
            this_corr = (Vx'*C(x,y)*Vy)/sqrt(Vx'*C(x,x)*Vx)/sqrt(Vy'*C(y,y)*Vy); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flip = sign(this_corr);
            Vx = Vx*flip;
            allVx(:,n_run,n_alpha) = Vx;
            allVy(:,n_run,n_alpha) = Vy;

            % testing of the CCA model to optimize the reg parameter
            % (alpha)
            C = Ctest;
            allCCA(n_run,n_alpha) = (Vx'*C(x,y)*Vy)/sqrt(Vx'*C(x,x)*Vx)/sqrt(Vy'*C(y,y)*Vy);
        end
    end
    CCA = mean(allCCA,1);
    [M,n_alpha] = max(CCA);
    
    % save the results (save more for the non-permuted data)
    CCAeval(perm_stat+1) = mean(allCCA(:,n_alpha));
    CCAevalruns(:,perm_stat+1) = allCCA(:,n_alpha);
    if ~perm_stat
        [Vx,S,tmp] = svd(allVx(:,:,n_alpha),'econ');
        [M,n_Vx] = max(diag(S));
        Vx = Vx(:,n_Vx);
        [Vy,S,tmp] = svd(allVy(:,:,n_alpha),'econ');
        [M,n_Vy] = max(diag(S));
        Vy = Vy(:,n_Vy);
    
        Vxeval = Vx;
        Vyeval = Vy;
        alphaxeval = alphax(n_alpha);
        alphayeval = alphay(n_alpha);
        CCAcurves = CCA;
    end
end











%%
% generate the data
X = subcoh(:,:)';
if 0 % data normalized by the CVC without noise
    X(1:8,:) = X(1:8,:)./repmat(X(9,:),8,1);
    X(10:17,:) = X(10:17,:)./repmat(X(18,:),8,1);
end
if 0 % data normalized by the CVC without noise
    X = mean(X([9 18],:));
end
if 0
    X = mean(X(1:8,:)./repmat(X(9,:),8,1));
    X(X>2) = 2;
end
if 0
    X = mean(X(10:17,:)./repmat(X(18,:),8,1));
end
if 0 % diff inf vs. non-inf, left hem
    X = (mean(X(5:8,:))-mean(X(1:4,:)))./(mean(X(5:8,:))+mean(X(1:4,:)));
end
if 0 % diff inf vs. non-inf, right hem
    X = (mean(X(14:17,:))-mean(X(10:13,:)))./(mean(X(14:17,:))+mean(X(10:13,:)));
end
if 0 % diff inf vs. non-inf, both hem
    X = (mean(X([5:8 14:17],:))-mean(X([1:4 10:13],:)))./(mean(X([5:8 14:17],:))+mean(X([1:4 10:13],:)));
end
if 0 % benefit of the video, right hem
    X = (mean(X([1:2:8],:))-mean(X([2:2:8],:)))./(mean(X([1:2:8],:))+mean(X([2:2:8],:)));
end
if 0 % benefit of the video, left hem
    X = (mean(X([10:2:17],:))-mean(X([11:2:17],:)))./(mean(X([10:2:17],:))+mean(X([11:2:17],:)));
end
if 1 % benefit of the video, both hem
    X = (mean(X([1:2:8 10:2:17],:))-mean(X([2:2:8 11:2:17],:)))./(mean(X([1:2:8 10:2:17],:))+mean(X([2:2:8 11:2:17],:)));
end

Y = scores_37_sub;
plot(X,Y,'.')
[r,p] = corr(X',Y')


if 1 % data normalized by the CVC without noise
    X = mean(X([9 18],:));
end



%
keep = 1:37;
keep([16 20 21 28 36]) = [];
X = (mean(subcohatt(:,:),2)' - mean(subcohglobal(:,:),2)');
X = (mean(subcohatt(:,:),2)' - mean(subcohglobal(:,:),2)')./(mean(subcohatt(:,:),2)' + mean(subcohglobal(:,:),2)');
Y = scores_37_sub;
[r,p] = corr(X(keep)',Y(keep)')




