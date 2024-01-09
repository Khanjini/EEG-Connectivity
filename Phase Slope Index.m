chans = {'F3';'FC3';'C5';'C3';'CP3';'P3';'O1';'Cz';'F4';'FC4';'C4';'C6';'CP4';'P4';'O2'};
selected_chan = [5,10,14,13,18,21,27,48,40,45,50,51,55,58,64];

chan = size(selected_chan,2);
nsamps = 250;
m = 512*2;
n = size(selected_chan,2);
block_size = 16;
nblocks = m/block_size;

num = num_sub;
raw_data = sprintf('data_sub_%02d.mat',num);
load(raw_data);
matrix_psi_FDR = zeros(15,15);
n_mat = 0;

for i = 1:size(ft_form.trial,2)/2;
    null_dis_psi = cell(1,15);
    
    cfg = [];
    cfg.dataset     = raw_data;
    cfg.trl = [(i-1)*2*512+1,i*2*512,2*512];
    
    cfg.bpfreq      = [1 60];
    cfg.reref       = 'yes';
    cfg.refchannel  = chans%
    prepro_data_rest = ft_preprocessing(cfg);
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'fourier';
    cfg.channel = chans;
    cfg.tapsmofrq = 3;
    cfg.foilim = [13 30]; %beta
    freq_rest = ft_freqanalysis(cfg, prepro_data_rest);
    
    cfg = [];
    cfg.method = 'psi';
    cfg.bandwidth = 5;
    psi = ft_connectivityanalysis(cfg, freq_rest);
    
    psi_spc = squeeze(psi.psispctrm(:,:,:));
    
    psi_mat = zeros(chan,chan);
    for num1 = 1:chan
        for num2 = 1:chan
            if num1 ~= num2
                psi_mat(num1,num2) = mean(psi_spc(num1,num2,:));
                psi_mat(num2,num1) = mean(psi_spc(num2,num1,:));
            end
        end
    end
    
    matrix_psi_FDR = matrix_psi_FDR + psi_mat;
    n_mat = n_mat +1;
end;

psi_mat_beta = matrix_psi_FDR./n_mat;
f_name = sprintf('psi_%02d.mat',num);
save(['~\psi' f_name],'psi_mat_beta');
