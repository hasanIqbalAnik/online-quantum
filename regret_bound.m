% Regret minimization using Regularized follow the leader algorithm

% define regret: sum of all losses from prediction minus sum of all losses
% when the true state was the worst for us
% r_t = sum(loss(prediction)) - minimum psi from C_n (sum(loss(measurement * psi)))

% this loss function must be convex and L-Lipschitz

% derivative of absolute value function is -1 if x < 0 or +1 if x > 0,
% undefined if x = 0

% Cn is the set of all trace-1 positive semi-definite complex matrices of
% dimension 2^n

regrets = [];
bounds = [];

% for i = 1:10
%     [Es, bs, rho, L] = generate_dt_regret(100);
% 
%     [xtra, T] = size(Es); % number of rounds
%     n = 1;
% 
% 
%     eta = .2; % random eta
%     curpred = (2 ^ (- n)) * eye(2 ^ n); % current random prediction
%     pd = RFTL(Es, bs, T, eta, curpred) % best prediction
%     rho % true state
% 
%     regrets = [regrets calculate_regret(bs, Es, rho, pd, T)] % calculate regret
%     bounds = [bounds sqrt(T * n)]; % theoretical bound
%     
% end


for i = 1:10

    [Es, bs, rho, L] = generate_dt_regret(100);

    [xtra, T] = size(Es); % number of rounds
    n = 1;

eta = sqrt((log(2) * n) / (2 * T * (L^2))); % eta from theorem
curpred = (2 ^ (- n)) * eye(2 ^ n); % current prediction
pd = RFTL(Es, bs, T, eta, curpred); % best prediction
rho; % true state

regrets = [regrets calculate_regret(bs, Es, rho, pd, T)] % calculate regret
bounds = [bounds 2 * 1 * sqrt(2 * log(2) * T)]; % theoretical bound
    
end

plot(regrets)
hold on
plot(bounds)
xlabel('rounds')
ylabel('bounds')
