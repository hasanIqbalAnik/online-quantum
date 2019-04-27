% Regret minimization using Regularized follow the leader algorithm

% define regret: sum of all losses from prediction minus sum of all losses
% when the true state was the worst for us
% r_t = sum(loss(prediction)) - minimum psi from C_n (sum(loss(measurement * psi)))

% this loss function must be convex and L-Lipschitz

% derivative of absolute value function is -1 if x < 0 or +1 if x > 0,
% undefined if x = 0

% Cn is the set of all trace-1 positive semi-definite complex matrices of
% dimension 2^n


[Es, bs, rho, L] = generate_dt_regret(100);


[xtra, T] = size(Es); % number of rounds
n = 1;




eta = .2; % random eta
curpred = (2 ^ (- n)) * eye(2 ^ n);
pd = RFTL(Es, bs, T, eta, n, rho, curpred) % best prediction
rho % true state

calculate_regret(bs, Es, rho, pd, T)
sqrt(T * n)


eta = sqrt((log(2) * n) / (2 * T * (L^2)));
curpred = (2 ^ (- n)) * eye(2 ^ n);
pd = RFTL(Es, bs, T, eta, n, rho, curpred) % best prediction
rho % true state

calculate_regret(bs, Es, rho, pd, T)
2 * 1 * sqrt(2 * log(2))