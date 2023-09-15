%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   PHASE-AWARE AUDIO INPAINTING (PHAIN)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
close all
clear
clc

rng(0)

%% paths

addpath(genpath('dataset'))
addpath('utils')
addpath(genpath('phase_correction'));
addpath('PHAIN');
addpath('results');

%% settings

Sounds = dir('dataset/EBU_SQAM/*.wav');
NN = length(Sounds);
data = cell(NN,1);
info = audioinfo(Sounds(1).name);
Fs = info.SampleRate;
for nn = 1:NN
    data{nn} = audioread(Sounds(nn).name);
end
clear audio info
     
gaps = 5:5:50; % [ms]
M = length(gaps);

N = 8; % # of gaps

sig_counter = 0;

methodLabels = {'B_PHAIN', 'B_PHAIN_ora', 'R_PHAIN', 'R_PHAIN_ora', 'UR_PHAIN', 'U_PHAIN'};

inner = 100;
outer = 10;

for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN, M);
end

SNR  = NaN(NN, M, N, length(methodLabels));  % SNR per gap
TIME = NaN(NN, M, N, length(methodLabels));  % execution time per gap
SNR_procedure = cell(NN, M, N, length(methodLabels));  % iterative behavior

%% parameters

winLen = 2800;
shiftLen = winLen/4;
FFTnum = winLen;

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;
param.offset = off;

paramsolver.epsilon = 1e-3;
paramsolver.tol = 5e-4;

paramsolver.tau = 0.25;
paramsolver.sigma = 1;
paramsolver.alpha = 1;

off = 'half';


%% inpainting

for nn = 1:NN

    signal = data{nn};

    for m = 1:M

        fprintf('\nSignal: %d / %d', nn, NN)
        fprintf('\nGap Length: %d [ms]\n\n', gaps(m))

        gapLength = gaps(m); % [ms]
        h = round(Fs*gapLength/1000); % [samples]
        full.length = length(signal);
        full.mask = true(full.length, 1);

        notDegraded = 0.5; % [s]
        segment.length = round((length(signal) - 2*notDegraded*Fs) / N);

        for i = 1:length(methodLabels)
            solution.(methodLabels{i}){nn, m} = signal;
        end

        for n = 1:N
            
            fprintf('Gap Number: %d / %d\n', n, N)
            idx = round(notDegraded*Fs) + ((n - 1)*segment.length+1:n*segment.length);

            s = round((winLen + 1) + rand()*(segment.length - 2*winLen - h));
            f = s + h - 1;
            segment.mask = true(segment.length, 1);
            segment.mask(s:f) = false;

            segment.data = signal(idx);
            segment.max = max(abs(segment.data));
            segment.data = segment.data/segment.max;
            segment.gapped = segment.data.*segment.mask;
            segment.gappedStore = segment.gapped;
            full.mask(idx) = segment.mask;

            % shortening the signal
            [q, ~, ~, ~, ~, ~, ~, ~, U, V, L] =...
                min_sig_supp_2(winLen, shiftLen, FFTnum, s, f, segment.length, 1, offset(s, f, shiftLen, off));
            origL = L;
            enoughL = ceil(L/lcm(shiftLen, FFTnum))*lcm(shiftLen, FFTnum);
            if L < enoughL
                L = enoughL;
            end
            Q = q + L - 1;
            if Q <= segment.length
                ssegment.mask = segment.mask(q:Q);
                ssegment.gapped = segment.gapped(q:Q);
                ssegment.data = segment.data(q:Q);
            else
                ssegment.mask = segment.mask(q:end);
                ssegment.gapped = segment.gapped(q:end);
                ssegment.data = segment.data(q:end);
                ssegment.mask(end+1:L) = true;
                ssegment.gapped(end+1:L) = 0;
                ssegment.data(end+1:L) = 0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%% B-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%

            fprintf('B-PHAIN...\n')
            tic
            param.type = 'B';
            paramsolver.I = 1000;

            [segment.solution, SNR_procedure{nn, m, n, 1}] = ...
                PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);

            solution.B_PHAIN{nn, m}(idx) = segment.solution*segment.max;
            TIME(nn, m, n, 1) = toc;

        end

    end

end


%% segment processing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Phain (ora.)                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(15)
        tic
        fprintf('\nStarting Phain(oracle)...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for ll = 1:numLambda
            paramsolver.lambda = lambdaArray(ll);
    
            % algorithm
            [segment.solution, snr_procedure{sig_counter,gap_counter,i,15,ll}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
    
            % updating the global solution
            solution.PhainOra{sig_counter,gap_counter,ll}(idxs) = segment.solution*segment.max;
    
            TIME(sig_counter,gap_counter,i,15,ll) = toc;
        end
    end % Phain(ora.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  R-Phain                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(16)
        tic
        fprintf('\nStarting R-Phain...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'RP';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);

            for ll = 1:numLambda
    
                paramsolver.lambda = lambdaArray(ll);
    
                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,16,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.RPhain{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,16,ll,pp) = toc;
    
            end
        end

    end % R-Phain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              R-Phain(ora.)                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(17)
        tic
        fprintf('\nStarting R-Phain(oracle)...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'RPo';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);
    
            for ll = 1:numLambda
    
                paramsolver.lambda = lambdaArray(ll);
    
                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,17,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.RPhainOra{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,17,ll,pp) = toc;

            end
        end

    end % R-Phain(ora.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  U-Phain                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(18)
        tic
        fprintf('\nStarting U-Phain...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'UP';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);

            for ll = 1:numLambda
                paramsolver.lambda = lambdaArray(ll);
        
                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,18,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.UPhain{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,18,ll,pp) = toc;
            end
        end

    end % U-Phain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                 UR-Phain                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(19)
        tic
        fprintf('\nStarting UR-Phain...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'URP';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);

            for ll = 1:numLambda
                paramsolver.lambda = lambdaArray(ll);

                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,19,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.URPhain{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,19,ll,pp) = toc;

            end
        end

    end % UR-Phain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              p-shrinkage                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(20)
        tic
        fprintf('\nStarting p-shrinkage...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'off';
        param.reweighting = 'off';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;

        for ll = 1:numLambda
            paramsolver.lambda = lambdaArray(ll);
    
            paramsolver.a = pForPSH;  % p=1: soft-thresholding, p->-inf: hard-thresholding
    
            % algorithm
            [segment.solution, snr_procedure{sig_counter,gap_counter,i,20,ll}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
    
            % updating the global solution
            solution.psh{sig_counter,gap_counter,ll}(idxs) = segment.solution*segment.max;
    
            TIME(sig_counter,gap_counter,i,20,ll) = toc;
        end
    end % p-shrinkage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Smooth-Hard                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(21)
        tic
        fprintf('\nStarting Smooth-Hard...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'off';
        param.reweighting = 'off';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,21}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.SH{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,21) = toc;
    end % Smooth-Hard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  SCAD                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(22)
        tic
        fprintf('\nStarting SCAD...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'off';
        param.reweighting = 'off';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;  % >2

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,22}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.SCAD{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,22) = toc;
    end % SCAD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Logarithmic                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(23)
        tic
        fprintf('\nStarting Logarithmic...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'off';
        param.reweighting = 'off';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;  % <= 1/lambda

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,23}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Log{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,23) = toc;
    end % logarithmic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           re-p-shrinkage                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(24)
        tic
        fprintf('\nStarting re-p-shrinkage...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'off';
        param.reweighting = 'on';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;  % p=1: soft-thresholding, p->-inf: hard-thresholding

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,24}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.repsh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,24) = toc;
    end % re-p-shrinkage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           re-Smooth-Hard                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(25)
        tic
        fprintf('\nStarting re-Smooth-Hard...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'off';
        param.reweighting = 'on';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,25}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.reSH{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,25) = toc;
    end % re-Smooth-Hard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               re-SCAD                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(26)
        tic
        fprintf('\nStarting re-SCAD...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'off';
        param.reweighting = 'on';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;  % >2

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,26}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.reSCAD{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,26) = toc;
    end % reSCAD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            re-Logarithmic                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(27)
        tic
        fprintf('\nStarting re-Logarithmic...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'off';
        param.reweighting = 'on';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;  % <= 1/lambda

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,27}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.reLog{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,27) = toc;
    end % relogarithmic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Phain-p                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(28)
        tic
        fprintf('\nStarting Phain-p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,28}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Phain_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,28) = toc;
    end % Phain_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Phain (ora.) -p                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(29)
        tic
        fprintf('\nStarting Phain(oracle)-p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,29}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.PhainOra_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,29) = toc;
    end % Phain(ora.)_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  R-Phain                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(30)
        tic
        fprintf('\nStarting R-Phain_p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'RP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,30}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhain_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,30) = toc;
    end % R-Phain_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            R-Phain(ora.)_p                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(31)
        tic
        fprintf('\nStarting R-Phain(oracle)_p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'RPo';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,31}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhainOra_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,31) = toc;
    end % R-Phain(ora.)_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                U-Phain_p                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(32)
        tic
        fprintf('\nStarting U-Phain_p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'UP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,32}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.UPhain_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,32) = toc;
    end % U-Phain_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               UR-Phain_p                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(33)
        tic
        fprintf('\nStarting UR-Phain_p...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'psh';
        param.phainsw = 'on';
        param.phaintype = 'URP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = pForPSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,33}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.URPhain_p{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,33) = toc;
    end % UR-Phain_p


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Phain-sh                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(34)
        tic
        fprintf('\nStarting Phain-sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,34}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Phain_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,34) = toc;
    end % Phain_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Phain (ora.) -sh                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(35)
        tic
        fprintf('\nStarting Phain(oracle)-sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,35}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.PhainOra_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,35) = toc;
    end % Phain(ora.)_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               R-Phain_sh                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(36)
        tic
        fprintf('\nStarting R-Phain_sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'RP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,36}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhain_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,36) = toc;
    end % R-Phain_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            R-Phain(ora.)_sh                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(37)
        tic
        fprintf('\nStarting R-Phain(oracle)_sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'RPo';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,37}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhainOra_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,37) = toc;
    end % R-Phain(ora.)_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                U-Phain_sh                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(38)
        tic
        fprintf('\nStarting U-Phain_sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'UP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,38}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.UPhain_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,38) = toc;
    end % U-Phain_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               UR-Phain_sh                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(39)
        tic
        fprintf('\nStarting UR-Phain_sh...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'sh';
        param.phainsw = 'on';
        param.phaintype = 'URP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSH;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,39}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.URPhain_sh{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,39) = toc;
    end % UR-Phain_sh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Phain-scad                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(40)
        tic
        fprintf('\nStarting Phain-scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,40}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Phain_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,40) = toc;
    end % Phain_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Phain (ora.) -scad                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(41)
        tic
        fprintf('\nStarting Phain(oracle)-scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,41}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.PhainOra_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,41) = toc;
    end % Phain(ora.)_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               R-Phain_scad                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(42)
        tic
        fprintf('\nStarting R-Phain_scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'RP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,42}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhain_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,42) = toc;
    end % R-Phain_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            R-Phain(ora.)_scad                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(43)
        tic
        fprintf('\nStarting R-Phain(oracle)_scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'RPo';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,43}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhainOra_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,43) = toc;
    end % R-Phain(ora.)_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                U-Phain_scad                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(44)
        tic
        fprintf('\nStarting U-Phain_scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'UP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,44}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.UPhain_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,44) = toc;
    end % U-Phain_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               UR-Phain_scad                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(45)
        tic
        fprintf('\nStarting UR-Phain_scad...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'scad';
        param.phainsw = 'on';
        param.phaintype = 'URP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForSCAD;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,45}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.URPhain_scad{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,45) = toc;
    end % UR-Phain_scad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Phain-log                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(46)
        tic
        fprintf('\nStarting Phain-log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,46}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Phain_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,46) = toc;
    end % Phain_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Phain (ora.) -log                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(47)
        tic
        fprintf('\nStarting Phain(oracle)-log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,47}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.PhainOra_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,47) = toc;
    end % Phain(ora.)_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               R-Phain_log                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(48)
        tic
        fprintf('\nStarting R-Phain_log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'RP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,48}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhain_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,48) = toc;
    end % R-Phain_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            R-Phain(ora.)_log                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(49)
        tic
        fprintf('\nStarting R-Phain(oracle)_log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'RPo';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,49}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhainOra_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,49) = toc;
    end % R-Phain(ora.)_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                U-Phain_log                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(50)
        tic
        fprintf('\nStarting U-Phain_log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'UP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,50}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.UPhain_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,50) = toc;
    end % U-Phain_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               UR-Phain_log                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(51)
        tic
        fprintf('\nStarting UR-Phain_log...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'log';
        param.phainsw = 'on';
        param.phaintype = 'URP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = 1;

        paramsolver.a = aForLOG;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,51}] = nonconvexMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.URPhain_l{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,51) = toc;
    end % UR-Phain_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                    Cos                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(52)
        tic
        fprintf('\nStarting Cosine thresholding...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'off';
        param.reweighting = 'off';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa; 

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,52}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,52) = toc;
    end % Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                 re-Cos                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(53)
        tic
        fprintf('\nStarting re-Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'off';
        param.reweighting = 'on';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 5;
        paramsolver.tau = 0.2;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa; 

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,53}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.reCos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,53) = toc;
    end % re-Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Phain-Cos                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(54)
        tic
        fprintf('\nStarting Phain-Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,54}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.Phain_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,54) = toc;
    end % Phain_Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Phain (ora.) -Cos                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(55)
        tic
        fprintf('\nStarting Phain(oracle)-Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'Po';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,55}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.PhainOra_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,55) = toc;
    end % Phain(ora.)_Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               R-Phain_Cos                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(56)
        tic
        fprintf('\nStarting R-Phain_Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'RP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,56}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhain_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,56) = toc;
    end % R-Phain_Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            R-Phain(ora.)_Cos                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(57)
        tic
        fprintf('\nStarting R-Phain(oracle)_Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'RPo';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,57}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.RPhainOra_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,57) = toc;
    end % R-Phain(ora.)_Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                U-Phain_Cos                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(58)
        tic
        fprintf('\nStarting U-Phain_Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'UP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,58}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.UPhain_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,58) = toc;
    end % U-Phain_Cos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               UR-Phain_Cos                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(59)
        tic
        fprintf('\nStarting UR-Phain_Cos...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.phainsw = 'on';
        param.phaintype = 'URP';

        paramsolver.maxit = 10;
        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.iter = 100;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;
        paramsolver.lambda = lambdaForCos;

        paramsolver.kappa = kappa;

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,59}] = cosMain(segment.gapped, segment.mask, param, paramsolver, segment.data);

        % updating the global solution
        solution.URPhain_Cos{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;

        TIME(sig_counter,gap_counter,i,59) = toc;
    end % UR-Phain_Cos


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                A-SPAIN LEARNED                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    if turnon(60)
        tic
        fprintf('Starting A-SPAIN LEARNED...\n');

        ASPAINBasisOptparam.algorithm = 'aspain_LEARNED';
        ASPAINBasisOptparam.w = w;
        ASPAINBasisOptparam.a = a;
        ASPAINBasisOptparam.wtype = 'hann';
        ASPAINBasisOptparam.M = M;
        ASPAINBasisOptparam.Ls = L;
        ASPAINBasisOptparam.mask = ssegment.mask;
        ASPAINBasisOptparam.gwindow = gabwin(ASPAINBasisOptparam.wtype, ASPAINBasisOptparam.a, ASPAINBasisOptparam.M);
        ASPAINBasisOptparam.gwindow = normalize(ASPAINBasisOptparam.gwindow,'peak'); % peak-normalization of the analysis window
        ASPAINBasisOptparam.gdual = gabdual(ASPAINBasisOptparam.gwindow,ASPAINBasisOptparam.a,ASPAINBasisOptparam.M); 

        % solver settings
        ASPAINBasisOptparamsolver.s = 1; % increment of k
        ASPAINBasisOptparamsolver.r = 1; % every r-th iteration increment k by s   
        ASPAINBasisOptparamsolver.epsilon = 0.001; % stopping criterion of termination function
        ASPAINBasisOptparamsolver.maxit = ceil(floor(ASPAINBasisOptparam.w/2+1)*ASPAINBasisOptparamsolver.r/ASPAINBasisOptparamsolver.s)*2; % maximum number of iterations
        ASPAINBasisOptparamsolver.store_snr = 0; 
        ASPAINBasisOptparamsolver.store_obj = 0;
        
        %% basis optimization

        q1=q-ASPAINBasisOptparam.w*5;
        Q1=Q+ASPAINBasisOptparam.w*5;
        
        
        if Q1 > length(segment.gapped_store) 
                segment.gapped_store = [segment.gapped_store ; zeros(Q1-length(segment.gapped_store) ,1)];
        end
      

        % prepare learning neighborhood
        coeff_TF = dgtreal(segment.gapped_store((max(q1,1)):Q1),ASPAINBasisOptparam.gwindow,ASPAINBasisOptparam.a,ASPAINBasisOptparam.M,'timeinv');
        Mprime=size(coeff_TF,1);



        mitte = round((Q1-q1)/(2*ASPAINBasisOptparam.a)+min(q1,0)/ASPAINBasisOptparam.a);
        offs = round((Q-q)/(2*ASPAINBasisOptparam.a));
        learn_offs=(2*ASPAINBasisOptparam.w/ASPAINBasisOptparam.a); %determines size of learning neighborhood



        X_tr = [coeff_TF(2:(Mprime-1),(max((mitte-offs-learn_offs),1)):(max((mitte-offs),1))) coeff_TF(2:(Mprime-1),(mitte+offs):(mitte+offs+learn_offs))];

        [Basis_small,sparsity_init,sparsity_final] = basis_opt_new(X_tr,1,2^(-20));
        Basis=eye(Mprime);
        Basis(2:(Mprime-1),2:(Mprime-1))=Basis_small;

        
        % algorithm   
        [recovered_ASPAIN_LEARNED, snr_procedure{sig_counter,gap_counter,i,60}, ~] = a_spain_learned(ssegment.gapped, ASPAINBasisOptparam, ASPAINBasisOptparamsolver, ssegment.data,Basis);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_ASPAIN_LEARNED(1:origL);

        % updating the global solution
        solution.ASPAIN_LEARNED{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
    
        TIME(sig_counter,gap_counter,i,60) = toc;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                iPCTV_L2                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(61)
        tic
        fprintf('\nStarting iPCTV_L2...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'UPL2';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;

        paramsolver.alpha = 1/4;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);

            for ll = 1:numLambda
                paramsolver.lambda = lambdaArray(ll);
        
                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,61,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.iPCTV_L2{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,61,ll,pp) = toc;
            end
        end

    end % iPCTV_L2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                PCTV                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(62)
        tic
        fprintf('\nStarting PCTV...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'PCTV';

        paramsolver.epsilon = 1e-3;
        paramsolver.delta = 0.01;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for pp = 1:numIterProp

            paramsolver.iter = inner(pp);
            paramsolver.maxit = outer(pp);

            for ll = 1:numLambda
                paramsolver.lambda = lambdaArray(ll);
        
                % algorithm
                [segment.solution, snr_procedure{sig_counter,gap_counter,i,62,ll,pp}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
        
                % updating the global solution
                solution.PCTV{sig_counter,gap_counter,ll,pp}(idxs) = segment.solution*segment.max;
        
                TIME(sig_counter,gap_counter,i,62,ll,pp) = toc;
            end
        end

    end % PCTV

% audiowrite("CP_"+num2str(w)+"_"+num2str(a)+".wav", solution.CP{sig_counter,gap_counter,ll}, fs, BitsPerSample = 32);

%% figures
fulltime = (1:length(signal))/fs;

% full time
if strcmp(gradual.type,'analysis')
    gradualweighting = CPweighting;
else
    gradualweighting = DRweighting;
end
leg = {'original',...
       'syn. $$\ell_1$$ (\texttt{none})',...
       ['$$\ell_1$$'],...
       ['syn. $$\ell_1$$ (\texttt{', DRweighting, '})'],...
       ['ana. $$\ell_1$$ (\texttt{', CPweighting, '})'],...
       'syn. $$\ell_1$$ (\texttt{iterative})',...   
       '\texttt{reweighting-}$$\ell_1$$',... 
       [gradual.type(1:3), '. $$\ell_1$$ (\texttt{', gradualweighting '}) + gradual'],...
       [TDCparam.type(1:3), '. $$\ell_1$$ (\texttt{', TDCweighting '}) + tdc'],...
       'S-SPAIN H',...
       'S-SPAIN OMP',...
       'A-SPAIN',...
       'OMP',...
       'Janssen',...
       'Phain',...
       'Phain(oracle)',...
       'R-Phain',...
       'R-Phain(oracle)',...
       'U-Phain',...
       'UR-Phain',...
       'p-shrinkage',...
       'Smooth-Hard',...
       'SCAD',...
       'Logarithmic',...
       're-p-shrinkage',...
       're-Smooth-Hard',...
       're-SCAD',...
       're-Logarithmic',...
       'Phain-p',...
       'Phain(oracle)-p',...
       'R-Phain-p',...
       'R-Phain(oracle)-p',...
       'U-Phain-p',...
       'UR-Phain-p',...
       'Phain-sh',...
       'Phain(oracle)-sh',...
       'R-Phain-sh',...
       'R-Phain(oracle)-sh',...
       'U-Phain-sh',...
       'UR-Phain-sh',...
       'Phain-scad',...
       'Phain(oracle)-scad',...
       'R-Phain-scad',...
       'R-Phain(oracle)-scad',...
       'U-Phain-scad',...
       'UR-Phain-scad',...
       'Phain-log',...
       'Phain(oracle)-log',...
       'R-Phain-log',...
       'R-Phain(oracle)-log',...
       'U-Phain-log',...
       'UR-Phain-log',...
       'Cos',...
       'reCos',...
       'Phain-Cos',...
       'Phain(oracle)-Cos',...
       'R-Phain-Cos',...
       'R-Phain(oracle)-Cos',...
       'U-Phain-Cos',...
       'UR-Phain-Cos',...
       'A-SPAIN LEARNED',...
       'iPCTV_L2',...
       'PCTV'
       };

fields = fieldnames(solution);
   
% figure()
% hold on
% plot(fulltime,signal)
% for i = 1:sum(turnon)
%     plot(fulltime,solution.(fields{i}));
% end
% legend(leg([true turnon]),'interpreter','latex')
% xlabel('time [s]')
% title(sprintf('signal: %s, gap length: %d ms, full signal',signame,gap_length),'interpreter','none')

%%
% 
% % only the gaps
% figure()
% hold on
% plot(signal(~full.mask))
% for i = 1:sum(turnon)
%     plot(solution.(fields{i}){sig_counter,gap_counter}(~full.mask));
% end
% legend(leg([true turnon]),'interpreter','latex')
% title(sprintf('signal: %s, gap length: %d ms, only the gaps',signame,gap_length),'interpreter','none')

%% SNRs (and ODGs)
% fprintf('\nSNRs and ODGs computed from all the gaps at once:\n')

for pp = 1:numIterProp
    for ll = 1:numLambda
        for i = 1:sum(turnon)
            restoredsignal = solution.(methodLabels{i}){sig_counter,gap_counter,ll,pp};
            yves = signal(~full.mask);
            olhae = restoredsignal(~full.mask);
            gowon = zeros(n,1);
            for ii = 1:n
                gowon(ii) = snr_n(yves(1+h*(ii-1):h*ii),olhae(1+h*(ii-1):h*ii));
            end
            snr_separate(sig_counter,gap_counter,:,i,ll,pp) = gowon;
            snr_part(sig_counter,gap_counter,i,ll,pp) = mean(gowon);
            snr_full(sig_counter,gap_counter,i,ll,pp) = snr_n(yves,olhae);
        end
    end
end
    
for pp = 1:numIterProp
    for ll = 1:numLambda
        for i = 1:sum(turnon)
            restoredsignal = solution.(methodLabels{i}){sig_counter,gap_counter,ll,pp};
            [~, ~, ODG_full(sig_counter,gap_counter,i,ll,pp), ~] = audioqual(signal,restoredsignal,fs);
        end
    end
end


LINE_Notify("Audio Inpainting Done")

%% 

snr_vec = squeeze(median(snr_separate,[1,3]));
figure()
plot(snr_vec)
legend(leg([false turnon]),'interpreter','latex')
xlabel('gap length 5:5:100 [ms]')

odg_vec = squeeze(median(ODG_full,1));
figure()
plot(odg_vec)
legend(leg([false turnon]),'interpreter','latex')
xlabel('gap length 5:5:100 [ms]')

%% for lambdas

gaplen = 8;
heejin = squeeze(median(snr_separate(:,gaplen,:,:,:),3));
figure('Position', [393 79.4000 1.5032e+03 970.4000]), tiledlayout(2,5)
for iii = 1:10
    nexttile
    semilogx(lambdaArray,squeeze(heejin(iii,:,:)))
    if iii == 10
        legend(leg([false turnon]),'interpreter','latex')
    end
end

%% for iter props

% heejin = squeeze(median(snr_separate,[1,3]));
% figure('Position', [393 79.4000 1.5032e+03 970.4000]), tiledlayout(2,2)
% for jj = 1:4
%     nexttile
%     plot(squeeze(heejin(:,jj,:)))
%     title(methodLabels(jj))
%     legend(['inner: ', num2str(inner(1)), ' outer: ', num2str(outer(1))],['inner: ', num2str(inner(2)), ' outer: ', num2str(outer(2))],...
%         ['inner: ', num2str(inner(3)), ' outer: ', num2str(outer(3))], ['inner: ', num2str(inner(4)), ' outer: ', num2str(outer(4))])
% end

%% 
% 
% R = squeeze(snr_procedure(:,8,:,6));
% R_procedure = 0;
% for nn = 1:8
%     for n = 1:10
%         if n+nn == 2
%             R_procedure = R{n,nn};
%         else
%             R_procedure = [R_procedure,R{n,nn}];
%         end
%     end
% end
% 
% J = squeeze(snr_procedure(:,8,:,13));
% J_procedure = 0;
% for nn = 1:8
%     for n = 1:10
%         if n+nn == 2
%             J_procedure = J{n,nn};
%         else
%             J_procedure = [J_procedure,J{n,nn}];
%         end
%     end
% end
% 
% A = squeeze(snr_procedure(:,8,:,11));
% A_procedure = 0;
% for nn = 1:8
%     for n = 1:10
%         if n+nn == 2
%             A_procedure = A{n,nn};
%         else
%             A_procedure = [A_procedure,A{n,nn}];
%         end
%     end
% end
% 
% P = squeeze(snr_procedure(:,8,:,15));
% P_procedure = 0;
% for nn = 1:8
%     for n = 1:10
%         if n+nn == 2
%             P_procedure = P{n,nn};
%         else
%             P_procedure = [P_procedure,P{n,nn}];
%         end
%     end
% end
% 
% U = squeeze(snr_procedure(:,8,:,18));
% U_procedure = 0;
% for nn = 1:8
%     for n = 1:10
%         if n+nn == 2
%             U_procedure = U{n,nn};
%         else
%             U_procedure = [U_procedure,U{n,nn}];
%         end
%     end
% end
% 
% figure('Position',[158.6000 41.8000 1.0496e+03 1.0288e+03]), tiledlayout(1,2)
% nexttile
% plot(R_procedure(:,150:200))
% nexttile
% plot(U_procedure(:,150:200))
% 
% figure('Position',[158.6000 41.8000 1.0496e+03 1.0288e+03]), tiledlayout(1,2)
% nexttile
% plot(J_procedure(:,1:50))
% nexttile
% plot(A_procedure(:,1:50))

%%

D = datestr(now,'yyyy_mm_dd_HH');
cd Results
mkdir(D)
cd ../.

save(['Results/',D,'/result.mat'],'fields','methodLabels','leg','turnon','snr_separate','snr_full','snr_part','ODG_full','TIME',"snr_procedure");

%% 

for nn = [2, 5, 8]
    signame = sigs{nn};
    signal = eval(signame);
    for n = [5, 7, 10]
        audiowrite(['Results/', D, '/true_sig', num2str(nn), '_gap', num2str(n), '.wav'], signal, fs, BitsPerSample = 32);
        audiowrite(['Results/', D, '/ASPAIN_sig', num2str(nn), '_gap', num2str(n), '.wav'], solution.ASPAIN{nn,n}, fs, BitsPerSample = 32);
        audiowrite(['Results/', D, '/UPHAIN_sig', num2str(nn), '_gap', num2str(n), '.wav'], solution.UPhain{nn,n}, fs, BitsPerSample = 32);
        audiowrite(['Results/', D, '/UPHAIN_p_sig', num2str(nn), '_gap', num2str(n), '.wav'], solution.UPhain_p{nn,n}, fs, BitsPerSample = 32);
        audiowrite(['Results/', D, '/UPHAIN_L2_sig', num2str(nn), '_gap', num2str(n), '.wav'], solution.iPCTV_L2{nn,n}, fs, BitsPerSample = 32);
    end
end
