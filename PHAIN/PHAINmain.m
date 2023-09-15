function [outsig, snr_procedure] = PHAINmain(insig, mask, param, paramsolver, oracle)

% insig ........... input gapped signal
% mask ............ logical vector indicating the missing samples

% param
%   .a
%   .M
%   .w
%   .offset
%   .type 

% paramsolver
%   .sigma .... parameter of prox_f*
%   .tau ...... parameter of prox_g
%   .alpha .... step size
%   .lambda ... threshold
%   .delta .... stop criterion
%   .epsilon .. stop criterion 2
%   .x0 ....... starting point of primal variable
%   .y0 ....... starting point of dual variable
%   .maxit .... maximal number of iterations
%   .tol ...... tolerance of the relative norm of x between two iterations

% Date: 09/10/2021

%% input signal shortening

a = param.a;
M = param.M;
w = param.w;

N = length(insig);
s = find(~mask,1,'first');
f = find(~mask,1,'last');
[q, L] = shortenForDGT(w, a, s, f, offset(s, f, a, param.offset));
if L < dgtlength(L,a,M)
    L = dgtlength(L,a,M);
end
if q+L-1 <= N && q >= 1
    signal = insig(q:q+L-1);
    mask = mask(q:q+L-1);
    oracle = oracle(q:q+L-1);
else
    q = max(q,1);
    padding = zeros(L-(N-q+1),1);
    signal = [ insig(q:end); padding ];
    mask = [ mask(q:end); true(L-(N-q+1),1) ];
    oracle = [ oracle(q:end); padding ];
end

%% 

[win, ~] = generalizedCosWin(w, 'hanning');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);
diff_win = numericalDiffWin(tight_win);
    
zeroPhaseFlag = true;
rotateFlag = true;

[sigIdx,sumIdx,sumArray,ifftArray,rotIdx] = precomputationForFDGT(length(signal),w,a,M);
ana = @(signal) FDGT(signal,tight_win,sigIdx,M,rotIdx,zeroPhaseFlag);
syn = @(spec) invFDGT(spec,tight_win,sumIdx,sumArray,ifftArray,rotIdx,zeroPhaseFlag)*w;
dana = @(signal) FDGT(signal,diff_win,sigIdx,M,rotIdx,zeroPhaseFlag);
IF= @(signal) calcInstFreq(ana(signal),dana(signal),M,w,rotateFlag);
pcana = @(spec,ifcum) instPhaseCorrection(spec,ifcum,a,M);
pcsyn = @(iPCspec,ifcum) invInstPhaseCorrection(iPCspec,ifcum,a,M);

T = @(x) [x(:,1:end-1) - x(:,2:end)];
T_ = @(y) [y(:,1) (y(:,2:end)-y(:,1:end-1)) -y(:,end)];

J = @(x,ifcum) T(pcana(ana(x),ifcum));
J_ = @(spec,ifcum) syn(pcsyn(T_(spec),ifcum));
    
%% weighting

wts = 1;

%% reconstruction

soft = @(x,gamma) sign(x) .* max(abs(x)-gamma, 0);

if strcmp(param.type,'P')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(signal);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    paramsolver.x0 = signal;

    param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
    param.K = @(sig) J(sig,inst_freq);
    param.K_adj = @(spec) J_(spec,inst_freq);
    [y_hat, snr_procedure] = PDS(param,paramsolver, oracle, mask);

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'Po')

    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(oracle);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    paramsolver.x0 = signal;
    
    param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
    param.K = @(sig) J(sig,inst_freq);
    param.K_adj = @(spec) J_(spec,inst_freq);
    [y_hat, snr_procedure] = PDS(param,paramsolver, oracle, mask);
    
    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'RP')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(signal);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    y_hat = zeros(param.dim,1); % initial solution
    paramsolver.x0 = signal;

    snr_procedure = NaN(paramsolver.iter, paramsolver.maxit);

    % the outer cycle
    for iteration = 1:paramsolver.maxit
        y_old = y_hat;
        param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
        param.K = @(sig) J(sig,inst_freq);
        param.K_adj = @(spec) J_(spec,inst_freq);
        [y_hat, snr_procedure(:,iteration)] = PDS(param,paramsolver, oracle, mask);
        if norm(y_old - y_hat) < paramsolver.delta
            break;
        end
%         rwts = 1./(abs(param.K(y_hat)) + paramsolver.epsilon);
%         rwts = abs(param.K(y_hat));
        denom = movmean(abs(ana(y_hat)),[0,1],2);
%         denom = denom./vecnorm(denom,2);
        rwts = 1./(denom + paramsolver.epsilon);
%         rwts = rwts./max(rwts,[],1);
        rwts = rwts(:,1:end-1);
    end

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'RPo')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(oracle);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    y_hat = zeros(param.dim,1); % initial solution
    paramsolver.x0 = signal;

    snr_procedure = NaN(paramsolver.iter, paramsolver.maxit);

    % the outer cycle
    for iteration = 1:paramsolver.maxit
        y_old = y_hat;
        param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
        param.K = @(sig) J(sig,inst_freq);
        param.K_adj = @(spec) J_(spec,inst_freq);
        [y_hat, snr_procedure(:,iteration)] = PDS(param,paramsolver, oracle, mask);
        if norm(y_old - y_hat) < paramsolver.delta
            break;
        end
%         rwts = 1./(abs(param.K(y_hat)) + paramsolver.epsilon);  
%         rwts = abs(param.K(y_hat));
        denom = movmean(abs(ana(y_hat)),[0,1],2);
%         denom = denom./vecnorm(denom,2);
        rwts = 1./(denom + paramsolver.epsilon);
%         rwts = rwts./max(rwts,[],1);
        rwts = rwts(:,1:end-1);
    end

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'UP')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(signal);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    y_hat = zeros(param.dim,1); % initial solution
    paramsolver.x0 = signal;

    snr_procedure = NaN(paramsolver.iter, paramsolver.maxit);

    % the outer cycle
    for iteration = 1:paramsolver.maxit
        y_old = y_hat;
        param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
        param.K = @(sig) J(sig,inst_freq);
        param.K_adj = @(spec) J_(spec,inst_freq);
        [y_hat, snr_procedure(:,iteration)] = PDS(param,paramsolver, oracle, mask);
        if norm(y_old - y_hat) < paramsolver.delta
            break;
        end
        inst_freq = IF(y_hat);
    end

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'URP')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(signal);
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    y_hat = zeros(param.dim,1); % initial solution
    paramsolver.x0 = signal;

    snr_procedure = NaN(paramsolver.iter,paramsolver.maxit);

    % the outer cycle
    for iteration = 1:paramsolver.maxit
        y_old = y_hat;
        param.prox_f = @(x) soft(x,rwts.*wts*lambda/sigma);
        param.K = @(sig) J(sig,inst_freq);
        param.K_adj = @(spec) J_(spec,inst_freq);
        [y_hat, snr_procedure(:,iteration)] = PDS(param,paramsolver, oracle, mask);
        if norm(y_old - y_hat) < paramsolver.delta
            break;
        end
%         rwts = 1./(abs(param.K(y_hat)) + paramsolver.epsilon); 
%         rwts = abs(param.K(y_hat));
        denom = movmean(abs(ana(y_hat)),[0,1],2);
%         denom = denom./vecnorm(denom,2);
        rwts = 1./(denom + paramsolver.epsilon);
%         rwts = rwts./max(rwts,[],1);
        rwts = rwts(:,1:end-1);
        inst_freq = IF(y_hat);
        param.K = @(sig) J(sig,inst_freq); 
    end

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'UPL2')
    
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    inst_freq = IF(signal);
    y_hat = zeros(param.dim,1); % initial solution
    paramsolver.x0 = signal;

    snr_procedure = NaN(paramsolver.iter, paramsolver.maxit);

    % the outer cycle
    for iteration = 1:paramsolver.maxit
        y_old = y_hat;
        param.K = @(sig) J(sig,inst_freq);
        param.K_adj = @(spec) J_(spec,inst_freq);
        [y_hat, snr_procedure(:,iteration)] = ProximalGradient(param, paramsolver, oracle, mask);
        if norm(y_old - y_hat) < paramsolver.delta
            break;
        end
        inst_freq = IF(y_hat);
        param.K = @(sig) J(sig,inst_freq);
    end

    % solution
    restored = real(y_hat(1:L));

elseif strcmp(param.type,'PCTV')
    
    rwts = 1;
    param.prox_g = @(x) proj_time(x,mask,signal.*mask);
    param.dim = L;
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    paramsolver.x0 = signal;

    param.prox_f = @(x) soft(x, rwts.*wts*lambda/sigma);
    param.K = @(sig) T(ana(sig));
    param.K_adj = @(spec) syn(T_(spec));
    [y_hat, snr_procedure] = PDS(param,paramsolver, oracle, mask);

    % solution
    restored = real(y_hat(1:L));
    
end

%% putting the restored signal together
outsig = insig;
outsig(q:q+L-1) = restored;
outsig = outsig(1:N);