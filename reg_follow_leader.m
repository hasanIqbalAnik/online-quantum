% Regret minimization using Regularized follow the leader algorithm

% define regret: sum of all losses from prediction minus sum of all losses
% when the true state was the worst for us
% r_t = sum(loss(prediction)) - minimum psi from C_n (sum(loss(measurement * psi)))

% Cn is the set of all trace-1 positive semi-definite complex matrices of
% dimension 2^n


function total_regret = calculate_regret(bs, Es, psi, omega, T)
% bs = {set of supposed measurement outcomes}, Es={set of measurement
% operators}, psi = true state, omega = predicted state, T = number of
% rounds

    total_regret = 0;
    for i = 1:T
        total_regret = total_regret + (abs(trace(Es{i} * omega) - bs(i)) - abs(trace(Es{i} * psi) - bs(i)));
    end
end