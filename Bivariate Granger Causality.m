selected_chan = [5,10,14,13,18,21,27,48,40,45,50,51,55,58,64];

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'BIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
tstat     = 'F';     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
nperms    = 500;    % number of permutations for permutation test
bsize     = 16;     % permutation test block size: empty for automatic (uses model order)
acmaxlags  = morder;
MA_size = 5; 
down_size = 2;

num = num_sub ;
f_name = sprintf('s%02d.mat',num);
load(f_name);

data = eeg.rest(selected_chan,:,:);

data = data.* 524288/2^24; %gain value
bf = ft_preproc_bandpassfilter(data, eeg.srate, [8 13], 3, 'but'); %alpha

%bf = reref(bf,[]); % Commom Average Reference

bf = bf(:,[512*3+1:512*63],:);
rest_data = reshape(bf,15,512*2,30); % Segementation for 2 seconds

ntrials = size(rest_data,3); 

% Moving Average for [-2 2]
X = rest_data;
MA_X = zeros(size(X,1),size(X,2)-4,size(X,3));
for i = 3:size(X,2)-2;
    for c = 1:size(X,1);
        MA_X(c,i-2,:) = mean(X(c,[i-2,i-1,i,i+1,i+2],:));
    end;
end;

% Down Sizing for catching the trend
for c = 1:size(X,1);
    D(c,:,:) = downsample(MA_X(c,:,:),down_size);
end;

X = D;

matrix = zeros(size(X,1),size(X,1));

BIC = tsdata_to_infocrit(X,momax,icregmode);
[~,bmo_BIC] = min(BIC);

morder = bmo_BIC;

nobs = size(X,2);

for rr = 1:14;
    for cc = rr+1:15;
        
        [F,A,SIG] = GCCA_tsdata_to_pwcgc(X([rr,cc],:,:),morder,regmode);
        pval = mvgc_pval(F,morder,nobs,ntrials,1,1,size(X,1)-1,tstat);
        sig  = significance(pval,alpha,mhtc);
        rho = var_specrad(A);
        
        if rho < 1;
            FP =permtest_tsdata_to_pwcgc(X([rr,cc],:,:),morder,bsize,nperms,regmode,acmaxlags);
            pval_p = empirical_pval(F,FP);
            sig_p  = significance(pval_p,alpha,mhtc);
        end;
        
        if sig(1,2) == 1 & sig_p(1,2) == 1;
            matrix(rr,cc) = F(1,2);
        end;
        
        if sig(2,1) == 1 & sig_p(2,1) == 1;
            matrix(cc,rr) = F(2,1);
        end
    end;
end;

BVGC_alpha = matrix;

f_name = sprintf('BVGC_%02d.mat',num);
save(['~\BVGC' f_name],'BVGC_alpha');
