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