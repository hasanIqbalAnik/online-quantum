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