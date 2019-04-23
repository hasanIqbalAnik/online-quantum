% Regret minimization using Regularized follow the leader algorithm

% define regret: sum of all losses from prediction minus sum of all losses
% when the true state was the worst for us
% r_t = sum(loss(prediction)) - minimum psi from C_n (sum(loss(measurement * psi)))

% this loss function must be convex and L-Lipschitz

% derivative of absolute value function is -1 if x < 0 or +1 if x > 0,
% undefined if x = 0

% Cn is the set of all trace-1 positive semi-definite complex matrices of
% dimension 2^n


[Es, bs, rho, L] = generate_data(500);
[xtra, T] = size(Es); % number of rounds
n = 1;


eta = sqrt((log(2) * n) / (2 * T * (L^2)));
pd = RFTL(Es, bs, T, eta, n, rho) % best prediction
rho % true state

calculate_regret(bs, Es, pd, rho, T)

% this regret should be bounded by O(L*sqrt(T * n)) for L-1 and L-2 losses
upper_bound = 2 * L * sqrt(2 * log(2) * T)



function sv = sum_val(t, meas_ops, pred_st, loss_q, phi)
sv = 0;
for i = 1:t
    lpt = 2 * (trace(meas_ops{i} * pred_st) - loss_q(i)) * meas_ops{i};    
    sv = sv + trace(lpt * phi);
end
end


function pred = RFTL(Es, bs, T, eta, n, rho)
w1 = (2 ^ (- n)) * eye(2 ^ n);
for i = 1:T
    
    cvx_begin quiet
    variable x(2, 2) semidefinite;
    minimize(eta * sum_val(i, Es, w1, bs, x) - quantum_entr(x))
    subject to
    trace(x) == 1
    cvx_end
    
    w1 = x;
 
end
pred = w1;
end

function cl = calculate_loss(Es, bs, x, T)
total = 0;
for i = 1:T
    total = total + (trace(Es{i} * x) - bs(i));
end
cl = total;
end

function total_regret = calculate_regret(bs, Es, psi, omega, T)
% bs = {set of supposed measurement outcomes}, Es={set of measurement
% operators}, psi = true state, omega = predicted state, T = number of
% rounds

cvx_begin quiet
variable x(2, 2) semidefinite;
minimize(calculate_loss(Es, bs, x, T))
subject to
trace(x) == 1
cvx_end

pred_loss = calculate_loss(Es, bs, omega, T);
min_loss = calculate_loss(Es, bs, x, T);

total_regret = pred_loss - min_loss;

end

