%% 0. ENVIRONMENT & CONFIGURATION
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
if ~exist('setup_environment', 'file')
    addpath(fullfile(scriptPath, '..'));
end

[meg_dir, subjects, deriv_dir, snd_dir, vid_dir] = setup_environment();

subjects_to_discard_from_analysis = {'Meg6666'}; % load all subjects to discard due to possible autism or low intelligence
% only do this script for the subjects with ICA and CORRECTED!! 
% these are all subject where ICA still needs to be done: 
subjects_not_to_include = {'Meg3925','Meg3949' ... % subj to remove? different CTS at differnt days 
   }; 

%% Coherence analysis
% subjects(find(e ==0)) = [];

% >>>>>>>>>> No matfile <<<<<<<<<<
%load(matfile,'CTS_coh_all','group','nave');
%n_rm = round(mean(nave(:,:,1:47),3)-330)';


for perm_stat = 0:1
    for n_sub = 1:length(subjects)
        tic
        sub_name = subjects(n_sub).name;
        % 1. Filter subjects to discard
        [skip, msg] = check_exclusion_criteria(sub_name, subjects_not_to_include, subjects_to_discard_from_analysis);
        if skip, fprintf('Skipping %s\n', sub_name, msg); continue; end

        % 2. Prepare matfile name in which data will be saved
        mat_path = fullfile(deriv_dir, 'coherence', sub_name);
        if ~exist(mat_path, 'dir'), mkdir(mat_path); end

        suffix = iif(perm_stat, '_main_coh_language_perm_stat', '_main_coh_language');
        save_file = fullfile(mat_path, [sub_name suffix '.mat']);

        % Retrieve subject's and sound folders and payed set, order, and vids
        subfold_no_date = fullfile(meg_dir,sub_name);
        date = dir(fullfile(subfold_no_date,'2*'));
        subfold =  fullfile(meg_dir,subjects(n_sub).name,date.name);

        CMall = [];
        
        % check if subject exists in data participants, if not - throw
        % warning and continue
        try 
            [set, order, vids] = get_set_order_vids(subfold);
            set_order = ['set' num2str(set) '_order' num2str(order)];
        catch
            warning('subject %s was not found in the datas participant file!', subjects(n_sub).name)
            continue
        end
        
         if 1 % ~exist(matfile)
         
             % try
             for n_vid = 1:length(vids) % loop over the different trials: video1, video2 and audio only
                 % retrieve file names (meg, eeg & sound)
                subj_files = get_subject_files(subfold, snd_dir, sub_name, set, order, n_vid);
        %         infos = get_subject_files(meg_dir, sub_name, n_vid);
        % 
        % 
        %             % ignore this meg file 
        %             if ~exist(megfile, 'file')
        %                 warning('no fif found for subject %s trial %s \n \t-- no file %s found', subjects(n_sub).name, trial, megfile)
        %                 continue
        %             end
        %         end
        %       
        %       ls(megfile)
                
                % retrive the timing of the conditions
                [t, vid_en_in_SiN, t_dis] = get_vid_timings(order,vids(n_vid));
                
                % in the case of audio only, we want an extra condition!! 
                if n_vid == 3
                    vid_en_in_SiN(1,:) = 2; 
                end
                
                
                % initialize the CM structures
                surrogmeg = [];
                surrogeeg = [];
                CMmeg = [];
                CMeeg = [];

              
                

                % >>>>>>>>>> if use_MEG ????? <<<<<<<<<<<
                for use_MEG = 1
                    % read wav sound signal

                    Yglb = load_WAV_audio(subj_files.snd_global);

                    if t(end) ~= length(Yglb.signal)/Yglb.Fs
                        error('sound file length and theoretical timings not compatible');
                    end
                    
                    % Load MISC stuff (signal, Fs, sensors(if MEG))
                    MISCorig = load_MISC(subj_files.meg_file, perm_stat);
                    
                    % realign the sound files
                    [dec, tds] = realign_sound_file(subj_files.sync_file, MISCorig, Yglb);
                    
                    % verification
                    [t1, t2, L, gof(n_sub, n_vid)] = alignement_verification(tds, dec, MISCorig, Yglb, perm_stat, false);
                    
                    % init CM struct
                    CM = init_CM(MISCorig.Fs, 2, 'boxcar');
                    % CM.tap = 3;
                    if exist('min_nave','var')
                        CM.maxave = min_nave(n_sub);
                    end         
                    
                    % 1. Prepare data sources: Follow same writting if need
                    % to add signals or just comment to discard
                    ref_sources = struct();
                    ref_sources(1).name = 'global sound';       ref_sources(1).data = load_WAV_audio(subj_files.snd_global);
                    ref_sources(2).name = 'attended speech';    ref_sources(2).data = load_WAV_audio(subj_files.snd_att);
                    %ref_sources(3).name = 'noise';              ref_sources(3).data = load_WAV_audio(subj_files.snd_noise);
                    
                    % Add lip movement (specific case of video processing)
                    Vmouth_FS = lips_apperture(vid_dir, vids(n_vid), ref_sources(2).data);
                    ref_sources(end+1).name = 'Mouth Surface';  ref_sources(end).data = Vmouth_FS;
                    
                    % 2. One loop of processing (Enveloppe -> Downsample -> Permute -> adding to CM) 
                    cf = get_mel_freqs(150, 7000, 31);
                    for i = 1:length(ref_sources)
                        % Extract enveloppe only if it is audio signal
                        if isstruct(ref_sources(i).data) && isfield(ref_sources(i).data, 'signal')
                            Yp = envelope_extraction(ref_sources(i).data, cf);
                        else
                            Yp = ref_sources(i).data;       % Already low frequency signal (enveloppe/mouvement)
                        end
                        % Allignement and Downsampling
                        MISC = downsample_signal(Yp, tds, t1, t2, size(MISCorig.signal));
                        % Apply permutation if necessary
                        if perm_stat
                            MISC = permute_signal_segments(MISC, MISCorig.Fs);
                        end
                        CM = add_CM_ref(CM, ref_sources(i).name, MISC);
                    end

                                        
                    % clear all sound variables (not needed anymore)
                    clear Y F_h_x Yds t1 t2 MISC MISCorig.signal Yp
                    
                    % define the frequencies of interest
                    
                    %                     % surrogate-data-based stat
                    %                     if ~perm_stat & isfield(CM,'maxave')
                    %                         CM.surrog.Nsim = 1000;
                    %                         CM.surrog.Find = 1:41;            % 9 à 17 si analyse bande 4-8 (laisser 2 pour bande beta)
                    %                         % CM.surrog.sens = sensors.picksMEG;
                    %                     end
                    %                     CM.CSD = 1;
                    

                    % Filter out bad
                    [CM, bad] = prepare_CM_data(CM, dec, L, subj_files.meg_file);
                    
                    % remove time points during which the actor was singing
                    
                    for n_dis = 1:size(t_dis,1)
                        n_dis_Y = t_dis(n_dis,:)*FS;
                        % audiowrite('test_sing.wav',Y(n_dis_Y(1)-FS:n_dis_Y(2)+FS),FS)
                        [tmp,n_dis_Yds] = min(abs(bsxfun(@minus,tds,n_dis_Y')),[],2);
                        set_to_bad = dec+n_dis_Yds(1):dec+n_dis_Yds(2);
                        set_to_bad = min(max(set_to_bad,1),length(bad));
                        bad(set_to_bad) = 1;
                    end
                    
                    
                    for n_cond = 1:6

                        % get the four different conditions:
                        switch n_cond
                            case 1
                                % video without noise
                                n_t = find(~any(repmat(cond_vid_en_in_SiN(:,1),size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                            case 2
                                % video with noise
                                n_t = find(~any(repmat(cond_vid_en_in_SiN(:,3),size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                            case 3
                                % picture without noise
                                n_t = find(~any(repmat(cond_vid_en_in_SiN(:,5),size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                            case 4
                                % picture with noise
                                n_t = find(~any(repmat(cond_vid_en_in_SiN(:,7),size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                            case 5 
                                % blinded without noise
                                n_t = find(~any(repmat([2;0;0;0],size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                            case 6
                                % blinded with noise
                                n_t = find(~any(repmat([2;1;1;1],size(vid_en_in_SiN(1,:))) - vid_en_in_SiN));
                        end
                        
                        if isempty(n_t)
                            % if condition is not present, skip and go to
                            % next condition
                            continue;
                        end
                        
                        CM.tdeb = double(raw.first_samp)+max(dec+t(n_t)*MISCorig.Fs,0);
                        CM.tfin = double(raw.first_samp)+min(dec+t(n_t+1)*MISCorig.Fs,double(raw.last_samp-raw.first_samp));
                        
                        if use_MEG
                            CMmeg = CM;
                            badmeg = bad;
                        else
                            CMeeg = CM;
                            badeeg = bad;
                        end
                        
                        %
                        %                     if is_morgane % use the same bad for both conditions
                        %                         % realign MEG and EEG data to use the same bad
                        %                         [dec,tds] = CM_CVC_realign_MISC_son(CMmeg.ref(min(2,length(CMmeg.ref))).chan,CMeeg.ref(min(2,length(CMeeg.ref))).chan,Fs,Fs);
                        %                         L = min(length(CMmeg.ref(1).chan)-dec*(dec>0),length(CMeeg.ref(1).chan(tds))+dec*(dec<0));
                        %                         t1 = dec*(dec>0)+(1:L);
                        %                         t2 = -dec*(dec<0)+(1:L);
                        %                         if 0
                        %                             plot(CMmeg.ref(min(2,length(CMmeg.ref))).chan(t1)); hold on; plot(CMeeg.ref(min(2,length(CMmeg.ref))).chan(tds(t2)),'g')
                        %                         end
                        %
                        %                         % refine the realignment
                        %                         p10 = round(length(t1)/10);
                        %                         p90 = round(length(t1)*9/10);
                        %                         t10 = mean(t2(tds(p10+(-10*Fs:10*Fs))));
                        %                         t90 = mean(t2(tds(p90+(-10*Fs:10*Fs))));
                        %                         t0 = t10-(t90-t10)/8;
                        %                         t100 = t90+(t90-t10)/8;
                        %                         t2 = round(linspace(t0,t100,length(t1)));
                        %                         t2 = min(max(t2,1),length(CMeeg.ref(1).chan));
                        %                         if 0
                        %                             plot(CMmeg.ref(min(2,length(CMmeg.ref))).chan(t1)); hold on; plot(CMeeg.ref(min(2,length(CMmeg.ref))).chan(t2),'g')
                        %                         end
                        %
                        %                         % combine the bad vectors
                        %                         badmeg(t1) = badmeg(t1) |  badeeg(t2);
                        %                         badeeg(t2) = badmeg(t1);
                        %                     end
                        
                        if use_MEG
                            CMmeg = CM_coh_MEG_ref_one_pass(CMmeg,badmeg);
                            if (n_vid == 2 & n_cond < 5) | (n_vid == 3 & n_cond > 4) 
                                CMmeg.maxave = CMmeg.nave - n_rm(n_cond);
                                CMmeg = CM_coh_MEG_ref_one_pass(CMmeg,badmeg);
                                CMmeg = rmfield(CMmeg,'maxave');
                            end
                            CMmeg = rmfield(CMmeg,{'fxy','spectre','ref'});
                            if CMmeg.nave == 0
                                CMmeg.Fxx(:) = 0;
                                CMmeg.Fyy(:) = 0;
                                CMmeg.Fxy(:) = 0;
                                CMmeg.Fgrads(:) = 0;
                                CMmeg.CSD(:) = 0;
                            end
                            
                            if (n_vid == 1 && n_cond == 1) || (isempty(CMall) && n_cond == 1)
                                CMall = CMmeg;
                            else
                                if size(CMall, 2) >= n_cond
                                    % merge
                                    CMall(n_cond) = CM_combine_results(CMall(n_cond),CMmeg);
                                else
                                    % make new row 
                                    CMall(n_cond) = CMmeg;
                                end
                            end
                            
                        else
                            % for EEG -- not used
                            CMeeg(n_cond) = CM_coh_EEG_ref_one_pass(CMeeg(n_cond),data,badeeg);
                        end
                    end
                end
            end
            surrog = [];
            if ~isempty(CMall)
                if isfield(CMall(1,1),'surrog')
                    for s1 = 1:size(CMall,1)
                        for s2 = 1:size(CMall,2)
                            %  >>>>>>>>>> Explanations ? <<<<<<<<<<
                            surrog(s1,s2) = CMall.surrog;
                            CMall.surrog = [];
                        end
                    end
                end
                gof_align = gof(n_sub,:);
                save(matfile,'CMall','surrog','gof_align')
            end
            % make_evoked(CM.fxy(1:306,:,1)/max(max(CM.fxy(1:306,:,1)))*1e-10,[CM.infile(1:end-4)'_fxy.fif'],raw,-round(CM.quantum/2)+1,CM.nave)
            %         catch
            %             failure(n_sub) = 1;
            %         end
        end
        toc
    end
end

