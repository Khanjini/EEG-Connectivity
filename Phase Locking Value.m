chans = {'F3';'FC3';'C5';'C3';'CP3';'P3';'O1';'Cz';'F4';'FC4';'C4';'C6';'CP4';'P4';'O2'};
selected_chan = [5,10,14,13,18,21,27,48,40,45,50,51,55,58,64];

chan = size(selected_chan,2);
nsamps = 250;
m = 512*2;
n = size(selected_chan,2);
block_size = 16;
nblocks = m/block_size;
alpha = 0.05;

num = num_sub;
raw_data = sprintf('data_sub_%02d.mat',num); 
load(raw_data);
matrix_plv_FDR = zeros(15,15);
n_mat = 0;

for i = 1:size(ft_form.trial,2)/2;
    null_dis_plv = cell(1,15);
    
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
    
    plv_mat = nan(chan,chan);
        
    for r = 1:chan;
        for c = 1:chan;
            if r ~= c;
                plv_mat(r,c) = squeeze(mean(abs(mean(exp(1i*(angle(freq_rest.fourierspctrm(:,r,:)) - angle(freq_rest.fourierspctrm(:,c,:)))),1))));
            end;
        end;
    end;
    
    fc_raw = prepro_data_rest.trial{1, 1};  
    
    %Permutation test
    for j = 1:n; % permutation for each channel
        
        FP_plv = nan(nsamps,n,n);
        shuffled = [1:j-1 j+1:n];
        
        for s = 1:nsamps; % Making spurious connectivity values as much as the number of samples for the test
            
            for ii = 1:14 % chuffling time-block of each channel except the choosen channel
                shuffled_chan = selected_chan(shuffled(ii));
                shaped_raw = reshape(fc_raw(shuffled_chan,:),1,block_size,nblocks);
                time_locking(shuffled_chan,:) = reshape(shaped_raw(:,:,randperm(nblocks)),1,size(fc_raw,2));
            end;
            
            time_locking(selected_chan(j),:) = fc_raw(selected_chan(j),:);
            
            prepro_data_rest.trial{1,1} = time_locking;
                
            freq_rest_per = ft_freqanalysis(cfg, prepro_data_rest);
            
            per_plv_mat = nan(chan,chan);
            
            for r = 1:chan; % making a spurious connectivity matrix
                for c = 1:chan;
                    if r ~= c;
                        per_plv_mat(r,c) = squeeze(mean(abs(mean(exp(1i*(angle(freq_rest_per.fourierspctrm(:,r,:)) -angle(freq_rest_per.fourierspctrm(:,c,:)))),1))));
                    end;
                end;
            end;
            FP_plv(s,:,:) = per_plv_mat;
            
        end;
        null_dis_plv{1,j} = FP_plv; % making a null-distribution of the connectivity values
    end;
      
    matrix_plv = zeros(15,15);
    mhtc = 'FDR';
    for c = 1:size(coh_mat,1);
        for r = c:size(coh_mat,1);
            if c ~= r
                null_c = null_dis_plv{1,c};
                null_r = null_dis_plv{1,r};
                null_all = cat(1,null_c,null_r);
                pval_p = empirical_pval(coh_mat,null_all);
                sig_p  = significance(pval_p,alpha,mhtc);
                clear null_c null_r null_all
                
                if sig_p(r,c) == 1;
                    matrix_plv(r,c) = coh_mat(r,c);
                end;
                if sig_p(c,r) == 1;
                    matrix_plv(c,r) = coh_mat(c,r);
                end
                clear pval_p sig_p
            end;
        end;
    end;
    
    matrix_plv_FDR = matrix_plv_FDR + matrix_plv;
    clear null_dis_coh;
    
    n_mat = n_mat +1;
end;

plv_mat_beta = matrix_plv_FDR./n_mat;
f_name = sprintf('plv_%02d.mat',num);
save(['~\plv' f_name],'plv_mat_beta');
