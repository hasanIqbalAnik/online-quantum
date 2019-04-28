mistakes = [];
upper_bounds = [];

for eps = .2:.01:.33
    
[Es, bs, rho, L] = generate_dt_regret(100, eps); % generate necessary data
[xtra, T] = size(Es); % number of rounds
n = 1; % number of qubits
curpred = (2 ^ (- n)) * eye(2 ^ n); % maximally mixed prediction
num_mistake = 0; % cound number of mistakes

for i=1:T
    if abs(trace(Es{i} *curpred) - trace(Es{i} * rho)) > eps
        num_mistake = num_mistake + 1;
        curpred = RFTL(Es, bs, i, 1, curpred);
    end
end


mistakes = [mistakes num_mistake];
upper_bounds = [upper_bounds 1/(eps^2)];

end

mistakes
upper_bounds
plot(mistakes);
hold on
plot(upper_bounds);

xlabel('eps values')
ylabel('bounds')
