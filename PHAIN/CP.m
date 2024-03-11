function [x, snr_procedure] = CP(param, paramsolver, oracle, mask)

%% initialization

x = paramsolver.x0;
u = paramsolver.u0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;

snr_procedure = NaN(paramsolver.I, 1);

%% iteration

for i = 1:paramsolver.I
    
    p = param.proj(x - tau*sigma*param.L_adj(u));
    v = u + param.L(2*p - x);
    q = v - param.prox(v);
    
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    
    snr_procedure(i) = calcSNR(oracle(~mask), x(~mask));

end