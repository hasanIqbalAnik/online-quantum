tracedist = [];
upper_bounds = [];
X = [];
%k=3;
eps = .2
%for eps = .01:.01:.33
for k = 1:1:10
    [rho, Es] = generateData(k, eps);
    hyp = LearningWPostSelection(k, Es, rho);
    
    [xtra, T] = size(Es); % number of rounds
    
    traced = trace(hyp - rho);
    
    tracedist = [tracedist traced];
    upper_bounds = [upper_bounds 4*sqrt(k*eps)];
    %X = [X eps];
    X = [X k]
end

tracedist
upper_bounds

plot(X,tracedist);
hold on
plot(X,upper_bounds);

%xlabel('eps values')
xlabel('k values')
ylabel('trace distance')