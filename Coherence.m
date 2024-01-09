num = num_sub;
raw_data = sprintf('data_sub_%02d.mat',num); 
load(raw_data);
matrix_coh_FDR = zeros(15,15);
n_mat = 0;

for i = 1:size(ft_form.trial,2)/2;
    null_dis_coh = cell(1,15);
    
    cfg = [];
    cfg.dataset     = raw_data;
    cfg.trl = [(i-1)*2*512+1,i*2*512,2*512];
    
    cfg.bpfreq      = [1 60]; %interest frequency range for pre-processing
    cfg.reref       = 'yes';
    cfg.refchannel  = chans;
    prepro_data_rest = ft_preprocessing(cfg);
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'fourier';
    cfg.channel = chans;
    cfg.tapsmofrq = 3;
    cfg.foilim = [13 30]; %beta --> interest frequency range for connectivity estimation
    freq_rest = ft_freqanalysis(cfg, prepro_data_rest);
    
    cfg = [];
    cfg.method = 'coh';
    cfg.complex = 'abs';
    coh = ft_connectivityanalysis(cfg, freq_rest);
    
    coh_spc = squeeze(coh.cohspctrm(:,:,:));
    
    coh_mat = nan(chan,chan);
    for num1 = 1:chan
        for num2 = 1:chan
            if num1 ~= num2
                coh_mat(num1,num2) = mean(coh_spc(num1,num2,:));
                coh_mat(num2,num1) = mean(coh_spc(num2,num1,:));
            end
        end
    end
    
    fc_raw = prepro_data_rest.trial{1, 1};
    
    
    % Permutation test
    for j = 1:n; % permutation for each channel
        
        FP_coh = nan(nsamps,n,n);
        shuffled = [1:j-1 j+1:n];
        
        for s = 1:nsamps; % Making spurious connectivity values as much as the number of samples for the test
            
            for ii = 1:14 % chuffling time-block of each channel except the choosen channel
                shuffled_chan = selected_chan(shuffled(ii));
                shaped_raw = reshape(fc_raw(shuffled_chan,:),1,block_size,nblocks);
                time_locking(shuffled_chan,:) = reshape(shaped_raw(:,:,randperm(nblocks)),1,size(fc_raw,2));
            end;
            
            time_locking(selected_chan(j),:) = fc_raw(selected_chan(j),:);
            
            prepro_data_rest.trial{1,1} = time_locking;
            
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.taper = 'dpss';
            cfg.output = 'fourier';
            cfg.channel = chans;
            cfg.tapsmofrq = 3;
            cfg.foilim = [13 30]; %beta
            freq_rest = ft_freqanalysis(cfg, prepro_data_rest);
            
            cfg = [];
            cfg.method = 'coh';
            cfg.complex = 'abs';
            per_coh = ft_connectivityanalysis(cfg, freq_rest);
            
            per_coh_spc = squeeze(per_coh.cohspctrm(:,:,:));
            
            per_coh_mat = nan(chan,chan);
            
            for num1 = 1:chan % making a spurious connectivity matrix
                for num2 = 1:chan
                    if num1 ~= num2
                        per_coh_mat(num1,num2) = mean(per_coh_spc(num1,num2,:));
                        per_coh_mat(num2,num1) = mean(per_coh_spc(num2,num1,:));
                    end
                end
            end
            FP_coh(s,:,:) = per_coh_mat;
            
        end;
        null_dis_coh{1,j} = FP_coh; % making a null-distribution of the connectivity values
    end;

    
    matrix_coh = zeros(15,15);
    mhtc = 'FDR';
    for c = 1:size(coh_mat,1);
        for r = c:size(coh_mat,1);
            if c ~= r
                null_c = null_dis_coh{1,c};
                null_r = null_dis_coh{1,r};
                null_all = cat(1,null_c,null_r);
                pval_p = empirical_pval(coh_mat,null_all);
                sig_p  = significance(pval_p,alpha,mhtc);
                clear null_c null_r null_all
                
                if sig_p(r,c) == 1;
                    matrix_coh(r,c) = coh_mat(r,c);
                end;
                if sig_p(c,r) == 1;
                    matrix_coh(c,r) = coh_mat(c,r);
                end
                clear pval_p sig_p
            end;
        end;
    end;
    
    matrix_coh_FDR = matrix_coh_FDR + matrix_coh;
    clear null_dis_coh;
    
    n_mat = n_mat +1;
end;

coh_mat_beta = matrix_coh_FDR./n_mat;
f_name = sprintf('coh_%02d.mat',num);
save(['~\coh' f_name],'coh_mat_beta');
