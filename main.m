%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      INPAINTING METHODS COMPARISON                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 09-Oct-2021~ revised for thesis
% 10-Feb-2022~ revised for journal
 
close all
clear
clc

rng(0)

% global conj_atoms;
% conj_atoms = true;

%% paths

addpath('reweighted l1 relaxation');
addpath('SPAIN');
addpath('Janssen');
addpath('OMP');
addpath('ltfat');
addpath('Toolbox');
addpath('PemoQ_1.4.0_free');
addpath('Phain');
addpath('Results');
addpath('Nonconvex');
addpath('basisopt');
addpath(genpath('cvx'));

ltfatstart

%% settings

load('EBU_SQAM.mat');

sigs = { 'a08_violin',...
         'a16_clarinet',...
         'a18_bassoon',...
         'a25_harp',...
         'a35_glockenspiel',...
         'a41_celesta',...
         'a42_accordion',...
         'a58_guitar_sarasate',...
         'a60_piano_schubert',...
         'a66_wind_ensemble_stravinsky'
         };
     
gaps = 5:5:50; % [ms]

% control, which methods to use in the default settings
turnon = logical([ 0,... %  1. DR
                   0,... %  2. CP
                   0,... %  3. wDR
                   0,... %  4. wCP
                   0,... %  5. reDR
                   0,... %  6. reCP
                   0,... %  7. gradual
                   0,... %  8. tdc
                   0,... %  9. SSPAIN H
                   0,... % 10. SSPAIN OMP
                   0,... % 11. ASPAIN
                   0,... % 12. OMP
                   0,... % 13. Janssen
                   ...
                   0,... % 14. Phain
                   0,... % 15. Phain(ora.)
                   0,... % 16. R-Phain
                   0,... % 17. R-Phain(ora.)
                   0,... % 18. U-Phain
                   0,... % 19. UR-Phain
                   ...
                   0,... % 20. p-shrinkage
                   0,... % 21. Smooth-Hard
                   0,... % 22. SCAD
                   0,... % 23. log-arithmic
                   ...
                   0,... % 24. re-p-shrinkage
                   0,... % 25. re-Smooth-Hard
                   0,... % 26. re-SCAD
                   0,... % 27. re-log-arithmic
                   ...
                   0,... % 28. Phain p-shrinkage
                   0,... % 29. Phain(ora.) p-shrinkage
                   0,... % 30. R-Phain p-shrinkage
                   0,... % 31. R-Phain(ora.) p-shrinkage
                   1,... % 32. U-Phain p-shrinkage
                   0,... % 33. UR-Phain p-shrinkage
                   ...
                   0,... % 34. Phain Smooth-Hard
                   0,... % 35. Phain(ora.) Smooth-Hard
                   0,... % 36. R-Phain Smooth-Hard
                   0,... % 37. R-Phain(ora.) Smooth-Hard
                   0,... % 38. U-Phain Smooth-Hard
                   0,... % 39. UR-Phain Smooth-Hard
                   ...
                   0,... % 40. Phain SCAD
                   0,... % 41. Phain(ora.) SCAD
                   0,... % 42. R-Phain SCAD
                   0,... % 43. R-Phain(ora.) SCAD
                   0,... % 44. U-Phain SCAD
                   0,... % 45. UR-Phain SCAD
                   ...
                   0,... % 46. Phain logarithmic
                   0,... % 47. Phain(ora.) logarithmic
                   0,... % 48. R-Phain logarithmic
                   0,... % 49. R-Phain(ora.) logarithmic
                   0,... % 50. U-Phain logarithmic
                   0,... % 51. UR-Phain logarithmic
                   ...
                   0,... % 52. Cosine thresholding
                   0,... % 53. re-Cosine thresholding
                   0,... % 54. Phain Cos
                   0,... % 55. Phain(ora.) Cos
                   0,... % 56. R-Phain Cos
                   0,... % 57. R-Phain(ora.) Cos
                   0,... % 58. U-Phain Cos
                   0,... % 59. UR-Phain Cos
                    ...
                   0,... % 60. ASPAIN_LEARNED
                   0,... % 61. iPCTV_L2
                   0,... % 62. PCTV
                   ]);

fprintf('\nAvailable signals:\n')
for i = 1:length(sigs)
    fprintf('  %2d.  %s\n',i,sigs{i})
end

fprintf('\nAvailable gap lengths:\n')
for i = 1:length(gaps)
    fprintf('  %2d.  %2d ms\n',i,gaps(i))
end

prompt = '\nInitialization with speech-shaped noise (0:No, 1:Yes): ';
init_speechnoise = input(prompt);

prompt = '\nLambda varies? (0:No, 1:Yes): ';
lambdaTurnOn = input(prompt);

prompt = '\nModify the propotion between the inner and outer iteration? (0:No, 1:Yes): ';
iterPropTurnOn = input(prompt);

prompt = '\nChoose signal number (1-10). Vector input accepted: ';
signums = input(prompt);

prompt = '\nChoose gap length number (1-10). Vector input accepted: ';
gapnums = input(prompt);

prompt = '\nChoose number of gaps per one signal (recommended max. 10): ';
n = input(prompt);

fprintf('\nControl the algorithms to be used:\n')
prompt = '  Use defaults? (0/1): ';
defaults = input(prompt);   

sig_counter = 0;

methodLabels = {'DR',...
                'CP',...
                'wDR',...
                'wCP',...
                'reDR',...
                'reCP',...
                'gradual',...
                'tdc',...
                'SSPAIN_H',...
                'SSPAIN_OMP',...
                'ASPAIN',...
                'OMP',...
                'Janssen',...
                'Phain',...
                'PhainOra',...
                'RPhain',...
                'RPhainOra',...
                'UPhain',...
                'URPhain',...
                'psh',...
                'SH',...
                'SCAD',...
                'Log',...
                'repsh',...
                'reSH',...
                'reSCAD',...
                'reLog',...
                'Phain_p',...
                'PhainOra_p',...
                'RPhain_p',...
                'RPhainOra_p',...
                'UPhain_p',...
                'URPhain_p',...
                'Phain_sh',...
                'PhainOra_sh',...
                'RPhain_sh',...
                'RPhainOra_sh',...
                'UPhain_sh',...
                'URPhain_sh',...
                'Phain_scad',...
                'PhainOra_scad',...
                'RPhain_scad',...
                'RPhainOra_scad',...
                'UPhain_scad',...
                'URPhain_scad',...
                'Phain_l',...
                'PhainOra_l',...
                'RPhain_l',...
                'RPhainOra_l',...
                'UPhain_l',...
                'URPhain_l',...
                'Cos',...
                'reCos',...
                'Phain_Cos',...
                'PhainOra_Cos',...
                'RPhain_Cos',...
                'RPhainOra_Cos',...
                'UPhain_Cos',...
                'URPhain_Cos',...
                'ASPAIN_LEARNED',...
                'iPCTV_L2',...
                'PCTV'
                };
methodLabels = methodLabels(turnon);

pForPSH  = 1.5;    % p = 1: soft-thresholding, p -> -inf: hard-thresholding
aForSH   = 0.01;   % 
% aForSCAD = 2.5;     % a > 2
aForSCAD = [2, 2.5];  % [2, integer>2] is default
aForLOG  = 0.9;   % a < 1/lambda
kappa    = 1;
lambdaForCos = 1;

if lambdaTurnOn
    lambdaArray = logspace(-2,2,10);
else
    lambdaArray = 1;
end
numLambda = length(lambdaArray);

if iterPropTurnOn
    inner = [100, 200, 100, 100, 50];
    outer = [ 20,  10,  10,   5, 10];
else
    inner = 100;
    outer = 10;
end
numIterProp = length(inner);

for i = 1:sum(turnon)
    solution.(methodLabels{i}) = cell(numel(signums),numel(gapnums),numLambda,numIterProp);
end

snr_separate = NaN(numel(signums),numel(gapnums),n,sum(turnon),numLambda,numIterProp);
snr_part = NaN(numel(signums),numel(gapnums),sum(turnon),numLambda,numIterProp);  % Average SNR of each gap 
snr_full = NaN(numel(signums),numel(gapnums),sum(turnon),numLambda,numIterProp);  % SNR of a whole signal
ODG_full = NaN(numel(signums),numel(gapnums),sum(turnon),numLambda,numIterProp);  % ODG of a whole signal (because of its length, one cannot compute ODG of each gap)
TIME = NaN(numel(signums),numel(gapnums),n,numel(turnon),numLambda,numIterProp);  % execution time per one gap
snr_procedure = cell(numel(signums),numel(gapnums),n,numel(turnon),numLambda,numIterProp);  % SNR values in each iteration

%% inpainting
for signum = signums
    
sig_counter = sig_counter + 1;
gap_counter = 0;    

for gapnum = gapnums

gap_counter = gap_counter + 1;    
    
fprintf('\nSignal number: %d (no. %d of the chosen %d)',signum,sig_counter,length(signums));
fprintf('\nGap length: %d ms (no. %d of the chosen %d)\n\n',gaps(gapnum),gap_counter,length(gapnums));
    
%% loading signal                                              
signame = sigs{signum};
signal = eval(signame);

fprintf('Signal: %s\n',signame);
fprintf('Sampling rate: %d Hz\n',fs);

%% transform setting
% window length approximately 64 ms + divisible by 4
% w = 2800;
% a = w/4;
w = 2048;
a = w/4;
M = w;

% frame for methods based on convex relaxation
F = frametight(frame('dgt',{'hann',w},a,M,'timeinv'));

% parameter of the transform for methods based on non-convex heuristics
DFTred = 4;

%% signal degradation
gap_length  = gaps(gapnum);
h           = round(fs*gap_length/1000); % gap length in samples 
full.length = length(signal);
full.mask   = true(full.length,1);

%% segment setting
not_degraded   = 0.5; % length of reliable part at the start / end of the signal (is seconds)
segment.length = round((length(signal)-2*not_degraded*fs)/n);

%% some global parameters 
off           = 'half';        % 'full' / 'half' / 'none'
gradual.type  = 'analysis';    % 'analysis' / 'synthesis'
TDCparam.type = 'analysis';    % 'analysis' / 'synthesis'

%% some not so global parameters
DRweighting   = 'norm';        % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
CPweighting   = 'energy';      % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
TDCweighting  = 'energy';      % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'

%% initializing solutions

for pp = 1:numIterProp
    for ll = 1:numLambda
        for i = 1:sum(turnon)
            solution.(methodLabels{i}){sig_counter,gap_counter,ll,pp} = signal;
        end
    end
end

%% segment processing
soft = @(x,gamma) sign(x) .* max(abs(x)-gamma, 0);
for i = 1:n
    fprintf('\nGap number %d of %d.\n',i,n);
    idxs = round(not_degraded*fs)+((i-1)*segment.length+1:i*segment.length);
    
    % degrading signal
    fprintf('Degrading signal.\n');      
    s = (w+1) + rand()*(segment.length-2*w-h); % gap start in samples
    s = round(s);
    f = s + h - 1; % gap end in samples
    segment.mask = true(segment.length,1); % indicator of the reliable samples
    segment.mask(s:f) = 0;
    segment.data = signal(idxs);
    segment.max = max(abs(segment.data));
    segment.data = segment.data/segment.max;
    segment.gapped = segment.data.*segment.mask; % lost samples set to 0 
    segment.gapped_store = segment.gapped; 
    full.mask(idxs) = segment.mask;

    if init_speechnoise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for initialization with speech-shaped
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% noise
        l_part = segment.data(s - h : s - 1);
        r_part = segment.data(f + 1 : f + h);
    
        % average_mag = (abs(fft(l_part)) + abs(fft(r_part)))/2;
        average_mag = abs(fft(l_part));
    
        if mod(h, 2) % odd
            endid = floor(h/2) + 1;
            rand_angle = -pi + 2*pi*rand(endid, 1);
            % rand_angle = angle(fft(l_part));
            speechNoise_fft = average_mag(1:endid) .* exp(1i * rand_angle(1:endid));
            speechNoise = ifft([0; speechNoise_fft(2:end); flip(conj(speechNoise_fft(2:end)))]);
        else % even
            endid = h/2 + 1;
            rand_angle = -pi + 2*pi*rand(endid, 1);
            % rand_angle = angle(fft(l_part));
            speechNoise_fft = average_mag(1:endid) .* exp(1i * rand_angle(1:endid));
            speechNoise = ifft([0; speechNoise_fft(2:end - 1); real(speechNoise_fft(end)); flip(conj(speechNoise_fft(2:end-1)))]);
        end

        segment.gapped(s:f) = speechNoise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
    end
    
    % shortening the signal
    [q,~,~,~,~,~,~,~,U,V,L] = min_sig_supp_2(w,a,M,s,f,segment.length,1,offset(s,f,a,off));
    origL = L;
    if L < framelength(F,L)
        L = framelength(F,L);
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
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Douglas-Rachford                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(1)
        fprintf('Starting Douglas-Rachford...\n')
        
        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.DR{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Chambolle-Pock                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(2)
        tic
        fprintf('Starting Chambolle-Pock...\n');
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = false;
        paramsolver.maxit = 1000;

        for ll = 1:numLambda
            paramsolver.lambda = lambdaArray(ll);
    
            % algorithm
            [segment.solution, snr_procedure{sig_counter,gap_counter,i,2,ll}] = reweighted(segment.gapped, segment.mask, param, paramsolver, [], [], segment.data);
    
            % updating the global solution
            solution.CP{sig_counter,gap_counter,ll}(idxs) = segment.solution(1:segment.length)*segment.max;
            TIME(sig_counter,gap_counter,i,2,ll) = toc;
           
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       weighted Douglas-Rachford                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(3)
        fprintf('Starting weighted Douglas-Rachford...\n')

        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = DRweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.wDR{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        weighted Chambolle-Pock                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(4)
        fprintf('Starting weighted Chambolle-Pock...\n');
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = CPweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.wCP{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                iteratively reweighted Douglas-Rachford                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(5)
        fprintf('Starting iteratively reweighted Douglas-Rachford...\n');
        
        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;

        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = true;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], RWparamsolver);

        % updating the global solution
        solution.reDR{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 iteratively reweighted Chambolle-Pock                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(6)
        tic
        fprintf('Starting iteratively reweighted Chambolle-Pock...\n');
        
        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = true;

        for ll = 1:numLambda
        paramsolver.lambda = lambdaArray(ll);

        % algorithm
        [segment.solution, snr_procedure{sig_counter,gap_counter,i,6,ll}] = reweighted(segment.gapped, segment.mask, param, paramsolver, [], RWparamsolver,segment.data);

        % updating the global solution
        solution.reCP{sig_counter,gap_counter,ll}(idxs) = segment.solution(1:segment.length)*segment.max;
        TIME(sig_counter,gap_counter,i,6,ll) = toc;
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                gradual                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(7) && strcmp(gradual.type,'analysis')
        tic
        fprintf('Starting gradual...\n');
        
        gapindexes = find(~ssegment.mask);
        gradual.s = gapindexes(1);
        gradual.f = gapindexes(end);
        gradual.r = ceil(h/4);
        gradual.mask = ssegment.mask;
        gradual.gapped = ssegment.gapped;
        gradual.atoms = syn_atoms(F);

        % not changing parameters
        CPparam.g = @(x) 0; % it is actually the indicator function... 
        CPparam.K = @(x) frana(F,x);
        CPparam.K_adj = @(x) frsyn(F,x);
        CPparam.dim = length(frsyn(F,frana(F,ssegment.gapped)));
        
        % solver settings
        CPparamsolver = [];
  
        % the main cycle
        while gradual.s <= gradual.f
            % computing the weights for l1 relaxation
            gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, CPweighting);

            % changing parameters
            CPparam.f = @(x) norm(gradual.wts.*x,1);
            CPparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);  
            CPparam.prox_g = @(x,gamma) proj_time(x,gradual.mask,gradual.gapped);

            % algorithm
            [ recovered_gradual, ~, ~ ] = ChambollePock(CPparam,CPparamsolver);
            recovered_gradual = real(recovered_gradual);

            % s, f update
            gradual.s = gradual.s + gradual.r;
            gradual.f = gradual.f - gradual.r;
            gradual.mask = true(L,1);
            gradual.mask(gradual.s:gradual.f) = false;
            gradual.gapped = gradual.mask.*recovered_gradual;
        end

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_gradual(1:origL);

        % updating the global solution
        solution.gradual{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;
        toc
    end  
    
    if turnon(7) && strcmp(gradual.type,'synthesis')
        tic
        fprintf('Starting gradual...\n');
        
        gapindexes = find(~ssegment.mask);
        gradual.s = gapindexes(1);
        gradual.f = gapindexes(end);
        gradual.r = ceil(h/4);
        gradual.mask = ssegment.mask;
        gradual.gapped = ssegment.gapped;
        gradual.atoms = syn_atoms(F);

        % not changing parameters        
        DRparam.g = @(x) 0; % it is actually the indicator function...
        c = frana(F,ssegment.gapped);
        DRparam.dim = length(c);
        
        % solver settings
        DRparamsolver = [];
  
        % the main cycle
        while gradual.s <= gradual.f
            % computing the weights for l1 relaxation
            gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, DRweighting);

            % changing parameters
            DRparam.f = @(x) norm(gradual.wts.*x,1);
            DRparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);
            c = frana(F,gradual.gapped);
            DRparam.prox_g = @(x,gamma) x - frana(F,gradual.mask.*frsyn(F,x)) + c;

            % algorithm
            [ x_hat, ~, ~ ] = DouglasRachford(DRparam,DRparamsolver);
            x_hat = DRparam.prox_g(x_hat);
            recovered_gradual = frsyn(F,x_hat);
            recovered_gradual = real(recovered_gradual);

            % s, f update
            gradual.s = gradual.s + gradual.r;
            gradual.f = gradual.f - gradual.r;
            gradual.mask = true(L,1);
            gradual.mask(gradual.s:gradual.f) = false;
            gradual.gapped = gradual.mask.*recovered_gradual;
        end

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_gradual(1:origL);

        % updating the global solution
        solution.gradual{sig_counter,gap_counter}(idxs) = segment.solution*segment.max;
        toc
    end      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 direct time-domain compensation (tdc)                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if turnon(8)      
        tic
        fprintf('Starting direct time-domain compensation (tdc)...\n');

        % (1) taking the basic reconstruction of the original gap
        
        % transform settings
        param.type = TDCparam.type;
        param.F = F;
        param.offset = off;   
        param.weighting = TDCweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.tdc(idxs) = segment.solution(1:segment.length)*segment.max;
        
        % (2) direct time domain compensation (working with the whole signal)
        
        % ensuring only the single gap
        TDCmask = true(length(solution.tdc),1);
        TDCmask(idxs) = segment.mask;
        
        % TDC parameters
        TDCparam.F = F;
        TDCparam.offset = off;
        TDCparam.weighting = TDCweighting;
        TDCparam.reweighting = false;
        TDCparamsolver.gaps = 4;
        TDCparamsolver.segs = 10;      
        TDCparamsolver.shift = w/2;
        TDCparamsolver.lens = h/4;
                
        % compensation
        solution.tdc{sig_counter,gap_counter} = tdc( solution.tdc, TDCmask, TDCparam, [], [], TDCparamsolver );    
        toc
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               S-SPAIN H                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(9)
        tic
        fprintf('Starting S-SPAIN H...\n');
        
        SPAINparam.algorithm = 'sspain';
        SPAINparam.w = w;
        SPAINparam.a = a;
        SPAINparam.wtype = 'hann';
        SPAINparam.M = M;
        SPAINparam.Ls = L;
        SPAINparam.mask = ssegment.mask;
        SPAINparam.F = frame('dft');
        SPAINparam.F.redundancy = DFTred;
        SPAINparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPAINparam.F.redundancy-1),1)]);
        SPAINparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPAINparam.F.redundancy);

        % solver settings
        SPAINparamsolver.s = 1; % increment of k
        SPAINparamsolver.r = 1; % every r-th iteration increment k by s   
        SPAINparamsolver.epsilon = 0.1; % stopping criterion of termination function
        SPAINparamsolver.maxit = ceil(floor(SPAINparam.w*SPAINparam.F.redundancy/2+1)*SPAINparamsolver.r/SPAINparamsolver.s); % maximum number of iterations
        SPAINparamsolver.store_snr = 0; 
        SPAINparamsolver.store_obj = 0;
        SPAINparamsolver.f_update = 'H';

        % algorithm    
        [recovered_SSPAIN, ~, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_SSPAIN(1:origL);

        % updating the global solution
        solution.SSPAIN_H{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;      
        toc
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              S-SPAIN OMP                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(10)
        tic
        fprintf('Starting S-SPAIN OMP...\n');
        
        SPAINparam.algorithm = 'sspain';
        SPAINparam.w = w;
        SPAINparam.a = a;
        SPAINparam.wtype = 'hann';
        SPAINparam.M = M;
        SPAINparam.Ls = L;
        SPAINparam.mask = ssegment.mask;
        SPAINparam.F = frame('dft');
        SPAINparam.F.redundancy = DFTred;
        SPAINparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPAINparam.F.redundancy-1),1)]);
        SPAINparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPAINparam.F.redundancy);

        % solver settings
        SPAINparamsolver.s = 1; % increment of k
        SPAINparamsolver.r = 1; % every r-th iteration increment k by s   
        SPAINparamsolver.epsilon = 0.1; % stopping criterion of termination function
        SPAINparamsolver.maxit = ceil(floor(SPAINparam.w*SPAINparam.F.redundancy/2+1)*SPAINparamsolver.r/SPAINparamsolver.s); % maximum number of iterations
        SPAINparamsolver.store_snr = 0; 
        SPAINparamsolver.store_obj = 0;
        SPAINparamsolver.f_update = 'OMP';

        fprintf('OMP chosen to compute the f-update. The algorithm may take up to several minutes.\n');

        % algorithm    
        [recovered_SSPAIN, ~, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_SSPAIN(1:origL);

        % updating the global solution
        solution.SSPAIN_OMP{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;      
        toc
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                A-SPAIN                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(11)
        tic
        fprintf('Starting A-SPAIN...\n');
        
        SPAINparam.algorithm = 'aspain';
        SPAINparam.w = w;
        SPAINparam.a = a;
        SPAINparam.wtype = 'hann';
        SPAINparam.M = M;
        SPAINparam.Ls = L;
        SPAINparam.mask = ssegment.mask;
        SPAINparam.F = frame('dft');
        SPAINparam.F.redundancy = DFTred;
        SPAINparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPAINparam.F.redundancy-1),1)]);
        SPAINparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPAINparam.F.redundancy);

        % solver settings
        SPAINparamsolver.s = 1; % increment of k
        SPAINparamsolver.r = 1; % every r-th iteration increment k by s   
        SPAINparamsolver.epsilon = 0.1; % stopping criterion of termination function
        SPAINparamsolver.maxit = ceil(floor(SPAINparam.w*SPAINparam.F.redundancy/2+1)*SPAINparamsolver.r/SPAINparamsolver.s); % maximum number of iterations
        SPAINparamsolver.store_snr = 1; 
        SPAINparamsolver.store_obj = 0;
           
        % algorithm    
        [recovered_ASPAIN, snr_procedure{sig_counter,gap_counter,i,11}, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_ASPAIN(1:origL);

        % updating the global solution
        solution.ASPAIN{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;  
        TIME(sig_counter,gap_counter,i,11) = toc;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  OMP                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(12)
        tic
        fprintf('Starting OMP...\n');
        
        OMPparam.w = w;
        OMPparam.a = a;
        OMPparam.wtype = 'hann';
        OMPparam.M = M;
        OMPparam.Ls = L;
        OMPparam.mask = ssegment.mask;
        OMPparam.F = frame('dft');
        OMPparam.F.redundancy = DFTred;
        OMPparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(OMPparam.F.redundancy-1),1)]);
        OMPparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/OMPparam.F.redundancy);

        % solver settings
        OMPparamsolver.crit = 'mse';    % termination function
        OMPparamsolver.epsilon = 0.1/w;  % stopping criterion of termination function
        OMPparamsolver.maxit = OMPparam.w*OMPparam.F.redundancy;
        OMPparamsolver.store_snr = 0; 
        OMPparamsolver.store_obj = 0;
        OMPparamsolver.sol = 'qr'; % algorithm to compute projection in OMP
        
        % algorithm    
        recovered_OMP = omp_segmentation(ssegment.gapped,OMPparam,OMPparamsolver);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_OMP(1:origL);

        % updating the global solution
        solution.OMP{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
        toc
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                Janssen                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(13)
        tic
        fprintf('Starting Janssen...\n');
        
        Janssenparam.w = w;
        Janssenparam.a = a;
        Janssenparam.wtype = 'hann';
        Janssenparam.Ls = L;
        Janssenparam.mask = ssegment.mask;

        % solver settings
        Janssenparamsolver.Nit = 50; % number of iterations
        % Janssenparamsolver.p = 2*a;    
    
        % algorithm    
        [recovered_Janssen, snr_procedure{sig_counter,gap_counter,i,13}] = janssen(ssegment.gapped,Janssenparam,Janssenparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_Janssen(1:origL);

        % updating the global solution
        solution.Janssen{sig_counter,gap_counter}(idxs) = segment.solution(1:segment.length)*segment.max;
        TIME(sig_counter,gap_counter,i,13) = toc;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                 Phain                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(14)
        tic
        fprintf('\nStarting Phain...\n');
        
        % parameters of the main cycle
        param.a = a;
        param.M = M;
        param.w = w;
        param.offset = off;
        param.type = 'P';

        paramsolver.epsilon = 1e-3;
        paramsolver.iter = 1000;
        paramsolver.tol = 5e-4;
        
        paramsolver.sigma = 1;
        paramsolver.tau = 0.25;
        paramsolver.alpha = 1;

        for ll = 1:numLambda
            paramsolver.lambda = lambdaArray(ll);
    
            % algorithm
            [segment.solution, snr_procedure{sig_counter,gap_counter,i,14,ll}] = propMain(segment.gapped, segment.mask, param, paramsolver, segment.data);
    
            % updating the global solution
            solution.Phain{sig_counter,gap_counter,ll}(idxs) = segment.solution*segment.max;
    
            TIME(sig_counter,gap_counter,i,14,ll) = toc;
        end

    end % Phain

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


end % gap num

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

end % gapnum
end % signum

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
