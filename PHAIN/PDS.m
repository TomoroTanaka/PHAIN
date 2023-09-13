function [x_hat, snr_procedure] = PDS(param,paramsolver, oracle, mask)

% inputs:
%      param
%         .K ........ linear function K(x)
%         .K_adj .... linear function K*(x)
%         .prox_f ... proximal operator of f(x), used as
%                     prox_f(arg, parameter)
%         .prox_g ... proximal operator of g(x), used as
%                     prox_g(arg, parameter)
%         .dim ...... length of x
%      paramsolver
%         .sigma .... parameter
%         .tau ...... parameter
%         .alpha .... step size
%         .x0 ....... starting point of primal variable
%         .y0 ....... starting point of dual variable
%         .maxit .... maximal number of iterations
%         .tol ...... tolerance of the relative norm of x between two
%                     iterations
%
% outputs:
%     x_hat ......... solution
%
% Date: 17/02/21

%% default settings
% if ~isfield(paramsolver,'sigma')
%     paramsolver.sigma = 5;
% end
% if ~isfield(paramsolver,'tau')
%     paramsolver.tau = 0.2;
% end
% if ~isfield(paramsolver,'alpha')
%     paramsolver.alpha = 1;
% end
% if ~isfield(paramsolver,'x0')
%     paramsolver.x0 = zeros(param.dim,1);
% end
if ~isfield(paramsolver,'y0')
    paramsolver.y0 = zeros(size(param.K(zeros(param.dim,1))));
end
% if ~isfield(paramsolver,'iter')
%     paramsolver.iter = 100;
% end
% if ~isfield(paramsolver,'tol')
%     paramsolver.tol = 5e-4;
% end

%% initialization
x_old = paramsolver.x0;
x = paramsolver.x0;
v = paramsolver.y0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;

snr_procedure = NaN(paramsolver.iter,1);

%% algorithm
i = 1;
rel_norm = inf;
% while i <= paramsolver.iter && (rel_norm > paramsolver.tol || i < 10) 
while i <= paramsolver.iter
    
    p = param.prox_g(x - tau*sigma*param.K_adj(v));
    r = v + param.K(2*p - x);
    q = r - param.prox_f(r);
    
    x = x + alpha*(p - x);
    v = v + alpha*(q - v);
    
    i = i + 1;
    rel_norm = norm(x - x_old)/norm(x_old);
    snr_procedure(i - 1) = snr_n(oracle(~mask),x(~mask));
    
    x_old = x;
end

%% output
x_hat = x;