eps = .3;

[Es, bs, rho, L] = generate_dt_regret(100, eps);

[xtra, T] = size(Es); % number of rounds
n = 1;

curpred = (2 ^ (- n)) * eye(2 ^ n);
num_mistake = 0;

for i=1:T
    v1 = trace(Es{i} *curpred);
    v2 = trace(Es{i} * rho);
    if abs(trace(Es{i} *curpred) - trace(Es{i} * rho)) > eps
        num_mistake = num_mistake + 1;
        curpred = RFTL(Es, bs, i, 1, n, rho, curpred);
    end
end

num_mistake
(1/eps^2)
curpred
rho



