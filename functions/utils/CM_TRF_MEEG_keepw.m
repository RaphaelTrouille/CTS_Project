function  [CM pred data_ds ref_sigs_ds ref_sigs_ds_filt tds] = CM_TRF_MEEG_keepw(CM,bad,data)

%  CM = CM_coh_MEG_ref_one_pass(CM)
%  CM = CM_coh_MEG_ref_one_pass(CM,bad)
%   Compute coherence with all kind of reference signals
%
%   input : 
%
%    these fields are needed in CM
%     infile :    name of the data fiff file
%
%    these fields are optionals
%     quantum :   number of time steps in the epochs (default : 1024)
%     Find :      index of frequencies to include (default : 1:quantum/2)
%     win :       type of window used (any window known by matlab, ex: 
%                 'bartlett', 'blackman', 'boxcar', 'chebwin',
%                 'hamming', 'hann', 'kaiser', 'triang') (default : boxcar)
%     r :         degree of overlap of epochs (default 5)
%     maxave :    upper limit on the number of averages
%     bads :      cell array containing the bad channels' index.
%                 ex : bads = {'MEG0111' 'MEG2213'}. Those channels will 
%                 not be considered in the artefact rejection. (default :
%                 no bad channels).
%     ref :       contains reference signals and filter parameters
%       .sig :    cell array with reference channels name. If more than one
%                 name is entered, take the Euclidian norma after filtering
%       .filt :   structure containing filter parameters
%         .par      : type of filter, e.g. {'high'  'low'}
%         .f_vect   : normalized cut-off frequencies, e.g. [1 195]/Fs*2
%         .Ws_vect  : normalized width frequencies, e.g. [0.5 5]/Fs*2
%       .rectif : 1 to rectify ref signal and 0 to do not rectify
%       .norm :   1 to normalize the ref signal and 0 to do not normalize 
%     CSD :       cross spectral density matrix, assagn whatever value if 
%                 its computation is required
%     surrog :    add to obtain surrogate-data-based stats. Note that using
%                 the default will lead to long calculation time and strict
%                 statistics. It might be better to use a priori frequency
%                 and sensor selection.
%       .Nsim :   number of simulations (default : 1000)
%       .sens :   sensor selection used for stat only (default : all MEG)
%       .Find :   frequency indices used for stat only (default : same as
%                 CM.Find
%       .alpha :  alpha value used in stat (default : 0.95 (p < 0.05))
%
%   output :
%
%     CM :        structure containing the result from the corticomuscular
%                 coherence
%
%   exemple :
% 
%     CM = [];
%     CM.infile = '/home/neuromag/Sarah/meg/MEG_0499/MEG_0499_CVC_lu_mc_lcfc.fif';
%     raw = fiff_setup_read_raw(CM.infile);
%     CM.win = 'boxcar';
%     CM.quantum = 2048;
%     Fs = 1000;
%     CM.tdeb = double(raw.first_samp);
%     CM.tfin = double(raw.last_samp);
%     CM.Find = 1:210;
% 
%     % ref: voice of the subject filtered around F0
%     Nsig = 1; % ref signal : sound file
%     CM.ref(Nsig).chan = {'MISC001'};
%     CM.ref(Nsig).filt.par = {'high' 'low'};
%     CM.ref(Nsig).filt.f_vect =  [50 330]/Fs*2;
%     CM.ref(Nsig).filt.Ws_vect = [5 5 ]/Fs*2;
%     CM.ref(Nsig).rectif = 1;
%     CM.ref(Nsig).norm = 1;  
%     % ref: voice of the subject filtered around F0
%     Nsig = 2; % ref signal : sound file
%     CM.ref(Nsig).chan = {'ECG061'};
%     CM.ref(Nsig).norm = 1;
%     CM.CSD = [];
% 
%     % surrogate-data-based stat
%     CM.surrog.Nsim = 1000;
%     CM.surrog.Find = 2:5;
%     CM.surrog.sens = 34:54;
% 
%     CM = CM_coh_MEG_ref_one_pass(CM);
% 
%     % Use fieldtrip to produce topoplots
%     cfg = [];
%     cfg.xparam = 'time';
%     % cfg.zlim = [0 0.16]
%     cfg.layout = 'neuromag306mag.lay';
%     ave = [];
%     ave.avg = reshape(repmat(CM.cohgrad(:,3,1)',3,1),306,1);
%     ave.fsample = 1000;
%     ave.time = 0;
%     ave.label = raw.info.ch_names(fiff_pick_types(raw.info,1,0))';
%     ave.dimord = 'chan_time';
%     figure; ft_topoplotER(cfg,ave);
%       
%     % Use make_evoked
%     make_evoked(CM.fxy(1:306,:,1)/max(max(CM.fxy(1:306,:,1)))*1e-10,[CM.infile(1:end-4) '_fxy.fif'],raw,-round(CM.quantum/2)+1,CM.nave,CM.ref.chan(1))
% 




% Parameters initialization
if ~isfield(CM,'Nfold')
    CM.Nfold = 10;
end

% pick MEEG signals
if ~isfield(CM,'picksMEEG')
    picksMEEG = 1:length(CM.label);
else
    picksMEEG = CM.picksMEEG;
end
Nmeeg = length(picksMEEG);
Nref = length(CM.ref);
from = 1;
to = size(data,2);

if isfield(CM,'tdeb')
    bad_suppl = ones(1,to-from+1);
    for k = 1:length(CM.tdeb)
        bad_suppl(CM.tdeb(k)-from+1:CM.tfin(k)-from+1) = 0;
    end
    bad = bad | bad_suppl;
end
bad = [1 bad 1];
bad = diff(bad);
CM.tdeb = find(bad == -1)+from-1;
CM.tfin = find(bad == 1)+from-2;


% ref signals
picksref = [];
F_h_ref = zeros(length(CM.ref),to-from+1);
Nref = length(CM.ref);
ref_sigs = zeros(size(cat(1,CM.ref.chan),1),to-from+1);
ref_inds = [];
for ref = 1:Nref
    ref_inds = [ref_inds ref*ones(1,size(CM.ref(ref).chan,1))];
    these_inds = find(ref_inds == ref);
    for n = 1:size(CM.ref(ref).chan,1)
        ref_sigs(these_inds(n),:) = conv(CM.ref(ref).chan(n,:),ones(1,CM.ds)/CM.ds,'same');
    end
    
    if 0 % high pass filter the reference signal
        F_h_y = cosine_filter(size(ref_sigs,2),CM.filt.par(1),CM.filt.f_vect(1),CM.filt.Ws_vect(1));
        Fref_sigs = fft(ref_sigs,[],2);
        Fref_sigs = Fref_sigs.*repmat(F_h_y,size(ref_sigs(:,1)));
        ref_sigs = real(ifft(Fref_sigs,[],2));
    end
end

% clip and downsample
Fs = CM.Fs;
data_ds = [];
ref_sigs_ds = [];
ep_ind = [];
F_h_y = cosine_filter(size(data,2),CM.filt.par,CM.filt.f_vect,CM.filt.Ws_vect);
Fdata = fft(data,[],2);
Fdata = Fdata.*repmat(F_h_y,size(data(:,1)));
data = real(ifft(Fdata,[],2));
for n = 1:length(CM.tdeb)
    data_ds = cat(2,data_ds,data(:,CM.tdeb(n)-from+1:CM.ds:CM.tfin(n)-from+1));
    ref_sigs_ds = cat(2,ref_sigs_ds,ref_sigs(:,CM.tdeb(n)-from+1:CM.ds:CM.tfin(n)-from+1));
    ep_ind = cat(2,ep_ind,n*ones(1,size(data_ds,2)-length(ep_ind)));
end
SDdata = std(data_ds,[],2);
data_ds = data_ds./repmat(SDdata,size(data_ds(1,:)));
SDref_sigs = std(ref_sigs_ds,[],2);
ref_sigs_ds./repmat(SDref_sigs,size(ref_sigs_ds(1,:)));
Fs = Fs/CM.ds;

% also compute band-pass filtered sound envelope
F_h_y = cosine_filter(length(ref_sigs_ds),CM.filt.par(2),CM.filt.f_vect(2)*CM.ds,CM.filt.Ws_vect(2)*CM.ds);
Fref_sigs_ds = fft(ref_sigs_ds,[],2);
Fref_sigs_ds = Fref_sigs_ds.*repmat(F_h_y,size(ref_sigs_ds(:,1)));
ref_sigs_ds_filt = real(ifft(Fref_sigs_ds,[],2));


% [U,S,V] = svd(data_ds,'econ');
% CM.incl = find(diag(S)>max(diag(S))/20);
% data_ds_sav = data_ds;
% data_ds = V(:,CM.incl)';
% Nmeeg = length(CM.incl);


CM.SDdata = SDdata;
CM.SDref_sigs = SDref_sigs;
CM.Npts = size(Fref_sigs_ds,2);



pred = zeros([size(data_ds') Nref]);
if any(CM.map == 1)
    % simple Nfold-fold cross-validation
    Nfold = CM.Nfold;
    L = round(linspace(1,size(ref_sigs_ds,2)+1,Nfold+1));
    CM.r = zeros(Nmeeg,length(CM.lambda),Nfold,Nref);
    for n_val = 1:Nfold
        % indicate discontinuities arrising due to splitting
        this_ep_ind = ep_ind;
        this_ep_ind(L(n_val+1):end) = this_ep_ind(L(n_val+1):end)+1;
        
        % define test and val periods
        keep_test = [1:L(n_val)-1 L(n_val+1):L(end)-1];
        keep_val = L(n_val):L(n_val+1)-1;

        for n_l = 1:length(CM.lambda)
            tic;
            disp(['Running n_val=' num2str(n_val) '/' num2str(Nfold) ' n_l=' num2str(n_l) '/' num2str(length(CM.lambda))])

            % Fit and the model to training data and validate it on val data
            for n_ref = 1:Nref
                these_inds = find(ref_inds==n_ref);
                [w,t,C] = mTRFtrain_ep(ref_sigs_ds(these_inds,keep_test)',data_ds(:,keep_test)',Fs,1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,CM.lambda(n_l),this_ep_ind(keep_test));
                [pred(keep_val,:,n_ref),CM.r(:,n_l,n_val,n_ref)] = mTRFpredict_ep(ref_sigs_ds(these_inds,keep_val)',data_ds(:,keep_val)',w,Fs,1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,C,this_ep_ind(keep_val));
            end
        end
    end
    
    % find lambda for which max r across sensors is maximal
    M = sort(mean(CM.r,3),1,'descend');
    M = mean(M(1:round(Nmeeg/20),:,:,:));
    [CM.rmax,n_l] = max(M,[],2);
    
    CM.n_l = squeeze(n_l);
    CM.rmax = squeeze(CM.rmax);

    % evaluate the TRFs
    CM.w = [];
    for n_ref = 1:Nref
        these_inds = find(ref_inds==n_ref);
        [CM.w{n_ref},CM.t,CM.C{n_ref}] = mTRFtrain_ep(ref_sigs_ds(these_inds,:)',data_ds',Fs,1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,CM.lambda(CM.n_l(n_ref)),ep_ind);
    end
    
    pred = permute(pred,[2 1 3]);


    % Fit the model to training + val data and test on nest data at
    % optimal lambda
    CM.rval = zeros(Nmeeg,Nref,Nfold);
    if Nmeeg == 204
        CM.rval = zeros(Nmeeg/2,Nref,Nfold);
    end


    for n_val = 1:Nfold
        keep_val = L(n_val):L(n_val+1)-1;
        keep_test = 1:L(end)-1;
        keep_test(keep_val) = [];
        for n_ref = 1:Nref
            for n_chan = 1:Nmeeg
                if Nmeeg == 204
                    if mod(n_chan,2)
                %                         EE = sum((data_ds(n_chan+[0 1],keep_val)'-pred(keep_val,n_chan+[0 1])).^2);
                %                         PP = sum((data_ds(n_chan+[0 1],keep_val)').^2);
                %                         CM.mseval((n_chan+1)/2,n_ref,n_val) = EE*CM.SDdata(n_chan+[0 1]).^2/(PP*CM.SDdata(n_chan+[0 1]).^2);

                        % evaluate co
                        cx = sin(pi/20:pi/20:pi);
                        cy = cos(pi/20:pi/20:pi);
                        dat = data_ds(n_chan+[0 1],keep_test);
                        pr = pred(n_chan+[0 1],keep_test,n_ref);
                        r = -2;
                        for k = 1:length(cx)
                            this_r = corr(([cx(k) cy(k)]*dat)',([cx(k) cy(k)]*pr)');
                            if this_r > r
                                r = this_r;
                                k_sav = k;
                            end
                        end

                        dat = data_ds(n_chan+[0 1],keep_val);
                        pr = pred(n_chan+[0 1],keep_val,n_ref);
                        k = k_sav;
                        CM.rval((n_chan+1)/2,n_ref,n_val) = corr(([cx(k) cy(k)]*dat)',([cx(k) cy(k)]*pr)');
                    end
                else
                    CM.rval(n_chan,n_ref,n_val) = corr(data_ds(n_chan,keep_val)',pred(n_chan,keep_val)');
                end
            end
        end
    end
end


CM.w = [];
if any(CM.map == -1)
%     % to boost computation speed
%     [U,S,V] = svd(data_ds,'econ');
%     CM.incl = find(diag(S)>max(diag(S))/1000);
%     data_ds = S(CM.incl,CM.incl)*V(:,CM.incl)';
%     Nmeeg = length(CM.incl);

    % nested cross-validation
    Nfold = CM.Nfold;
    L = round(linspace(1,size(ref_sigs_ds,2)+1,Nfold+1));
    CM.r_nest = zeros(Nref,Nfold);
    for n_nest = 1:Nfold
        disp(['Running n_nest=' num2str(n_nest) '/' num2str(Nfold)])
        
        if 0
            r = cell(1,Nfold);
            for n_val = 1:Nfold
                r{n_val} = parfor_val_loop_rev(ep_ind,L,n_val,CM.lambda,CM.t_pre,CM.t_buff,CM.t_post,ref_sigs_ds,data_ds,ref_sigs_ds_filt,Fs,Nfold)
            end
            r = cat(3,r{:});
        else
            r = zeros(Nref,length(CM.lambda),Nfold);
            for n_val = [1:n_nest-1 n_nest+1:Nfold]
                % indicate discontinuities arrising due to splitting
                this_ep_ind = ep_ind;
                this_ep_ind(L(n_nest+1):end) = this_ep_ind(L(n_nest+1):end)+1;
                this_ep_ind(L(n_val+1):end) = this_ep_ind(L(n_val+1):end)+1;

                for n_l = 1:length(CM.lambda)
                    tic;
                    disp(['Running n_nest=' num2str(n_nest) '/' num2str(Nfold) ' n_val=' num2str(n_val) '/' num2str(Nfold) ' n_l=' num2str(n_l) '/' num2str(length(CM.lambda))])

                    % Fit the model to training data
                    n1 = min(n_nest,n_val);
                    n2 = max(n_nest,n_val);
                    keep_test = [1:L(n1)-1 L(n1+1):L(n2)-1 L(n2+1):L(end)-1];
                    [w,t,C] = mTRFtrain_ep(ref_sigs_ds(:,keep_test)',data_ds(:,keep_test)',Fs,-1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,CM.lambda(n_l),this_ep_ind(keep_test));

                    % Test the model on validation data
                    keep_val = L(n_val):L(n_val+1)-1;
                    [pred,r(:,n_l,n_val)] = mTRFpredict_ep(ref_sigs_ds_filt(:,keep_val)',data_ds(:,keep_val)',w,Fs,-1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,C,this_ep_ind(keep_val));
                    if ismac
                        pause(toc)
                    end
                end
            end
        end

        % find optimal reg param
        r(:,:,n_nest) = [];
        r_mean = mean(r,3);
        [r_max,n_l_opt] = max(r_mean,[],2);
        CM.n_l_opt(n_nest,:) = n_l_opt;
        CM.r_val(:,:,:,n_nest) = r;

        % indicate discontinuities arrising due to splitting
        this_ep_ind = ep_ind;
        this_ep_ind(L(n_nest+1):end) = this_ep_ind(L(n_nest+1):end)+1;

        % Fit the model to training + val data and test on nest data at
        % optimal lambda
        keep_test = [1:L(n_nest)-1 L(n_nest+1):L(end)-1];
        keep_nest = L(n_nest):L(n_nest+1)-1;
        for n_ref = 1:Nref
            % fit
            [w,t,C] = mTRFtrain_ep(ref_sigs_ds(n_ref,keep_test)',data_ds(:,keep_test)',Fs,-1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,CM.lambda(n_l_opt(n_ref)),this_ep_ind(keep_test));
            CM.w(:,:,n_ref,n_nest) = w;
            
            % test
            [pred,CM.r_nest(n_ref,n_nest)] = mTRFpredict_ep(ref_sigs_ds_filt(n_ref,keep_nest)',data_ds(:,keep_nest)',w,Fs,-1,CM.t_pre-CM.t_buff,CM.t_post+CM.t_buff,C,this_ep_ind(keep_nest));
        end
    end
end




function r = parfor_val_loop_fwd(ep_ind,L,n_val,lambda,t_pre,t_buff,t_post,Nref,ref_inds,ref_sigs_ds,data_ds,Fs,Nmeeg,Nfold)

% indicate discontinuities arrising due to splitting
this_ep_ind = ep_ind;
this_ep_ind(L(n_val+1):end) = this_ep_ind(L(n_val+1):end)+1;

% define test and val periods
keep_test = [1:L(n_val)-1 L(n_val+1):L(end)-1];
keep_val = L(n_val):L(n_val+1)-1;

r = zeros(Nmeeg,length(lambda),1,Nref);
for n_l = 1:length(lambda)
    tic;
    disp(['Running n_val=' num2str(n_val) '/' num2str(Nfold) ' n_l=' num2str(n_l) '/' num2str(length(lambda))])

    % Fit and the model to training data and validate it on val data
    for n_ref = 1:Nref
        these_inds = find(ref_inds==n_ref);
        [w,t,C] = mTRFtrain_ep(ref_sigs_ds(these_inds,keep_test)',data_ds(:,keep_test)',Fs,1,t_pre-t_buff,t_post+t_buff,lambda(n_l),this_ep_ind(keep_test));
        [pred,r(:,n_l,1,n_ref)] = mTRFpredict_ep(ref_sigs_ds(these_inds,keep_val)',data_ds(:,keep_val)',w,Fs,1,t_pre-t_buff,t_post+t_buff,C,this_ep_ind(keep_val));
    end
end



function r = parfor_val_loop_rev(ep_ind,L,n_val,lambda,t_pre,t_buff,t_post,ref_sigs_ds,data_ds,ref_sigs_ds_filt,Fs,Nfold)

 % indicate discontinuities arrising due to splitting
this_ep_ind = ep_ind;
this_ep_ind(L(n_val+1):end) = this_ep_ind(L(n_val+1):end)+1;

for n_l = 1:length(lambda)
    tic;
    disp(['Running n_val=' num2str(n_val) '/' num2str(Nfold) ' n_l=' num2str(n_l) '/' num2str(length(lambda))])

    % Fit the model to training data
    keep_test = [1:L(n_val)-1 L(n_val+1):L(end)-1];
    [w,t,C] = mTRFtrain_ep(ref_sigs_ds(:,keep_test)',data_ds(:,keep_test)',Fs,-1,t_pre-t_buff,t_post+t_buff,lambda(n_l),this_ep_ind(keep_test));

    % Test the model on validation data
    keep_val = L(n_val):L(n_val+1)-1;
    [pred,r(:,n_l)] = mTRFpredict_ep(ref_sigs_ds_filt(:,keep_val)',data_ds(:,keep_val)',w,Fs,-1,t_pre-t_buff,t_post+t_buff,C,this_ep_ind(keep_val));
    if ismac
        pause(toc)
    end
end
