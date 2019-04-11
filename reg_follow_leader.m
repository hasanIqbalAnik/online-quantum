% Regret minimization using Regularized follow the leader algorithm

% define regret: sum of all losses from prediction minus sum of all losses
% when the true state was the worst for us
% r_t = sum(loss(prediction)) - minimum psi from C_n (sum(loss(measurement * psi)))

% this loss function must be convex and L-Lipschitz

% derivative of absolute value function is -1 if x < 0 or +1 if x > 0,
% undefined if x = 0

% Cn is the set of all trace-1 positive semi-definite complex matrices of
% dimension 2^n

T = 30; % number of rounds
[Es, bs, rho] = generate_data(1000);

pd = RFTL(Es, bs, T, .4, 1, rho)
rho

calculate_regret(bs, Es, pd, rho, T)


% this regret should be bounded by O(L*sqrt(T * n)) for L-1 and L-2 losses
upper_bound = 2 * sqrt(2 * log(2) * 20)

function sv = sum_val(t, lpt, phi)
sv = 0;
for i = 1:t
    sv = sv + trace(lpt * phi);
end
end


function pred = RFTL(Es, bs, T, eta, n, rho)
w1 = (2 ^ (- n)) * eye(2 ^ n);

for i = 1:T
    
    lpt = 2 * (trace(Es{i} * w1) - bs(i)) * Es{i};
 
    cvx_begin quiet
    variable x(2, 2) semidefinite;
    minimize(eta * sum_val(i, lpt, x) - quantum_entr(x))
    subject to
    trace(x) == 1
    cvx_end
    
    w1 = x;
 
end
pred = w1;
end

function total_regret = calculate_regret(bs, Es, psi, omega, T)
% bs = {set of supposed measurement outcomes}, Es={set of measurement
% operators}, psi = true state, omega = predicted state, T = number of
% rounds

total_regret = 0;
for i = 1:T
    total_regret = total_regret + (trace(Es{i} * omega) - bs(i))^2 - (trace(Es{i} * psi) - bs(i))^2;
end
end

