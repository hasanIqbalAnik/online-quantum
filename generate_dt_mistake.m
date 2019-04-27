function [Es, bs, rho, L] = generate_data(s, rndmin, rndsd, eps)

% Generate test state
psi = sqrt(5 / 7) * [1; 0] + sqrt(2 / 7) * [0; 1];  
rho = psi * psi'; % density matrix

% Generate measurement operators (E_i, i = 1.......T)
% They should be Hermitian and eigenvalues would be in [0, 1]

Es = {}; % cell array to store matrices
ctr = 1;

    for i = 0:s
        m = rand(2) + 1j * rand(2); % random complex matrix
        E = (m + m') / 2; % random hermitian matrix

        es = eig(E);
        if all(es >= 0) && all(es <= 1)  % eigenvalues should be in [0, 1]
            Es{ctr} = E;
            eig(E);
            ctr = ctr + 1;
        end
    end
% finished generating measurement operators

% Realizable case:
% Generate b_t = Tr(E_t * rho), these values are in eps/3 >= bt >= -(eps/3)
% range

    
    bs = [];
    for j = 1:ctr - 1
        
        v1 = trace(Es{j} * rho);
        a = -((eps/3) - v1);
        b = (eps/3) - v1;
        bt = (b-a).*rand(1,1) + a; % random uniform sampling from that range
        bs = [bs bt];
    end
% finished generating errors b_t

% Generate Loss functions until l_T which is  l_t(z) = |z - b_t|
% This z value is our prediction from the algorithm and b_t is defined before.

% Non-realizable case:
% Generate b_t arbitrarily



L = 1;

end