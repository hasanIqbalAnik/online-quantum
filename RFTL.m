
function pred = RFTL(Es, bs, T, eta, curpred)
for i = 1:T
    
    cvx_begin quiet
    variable x(2, 2) semidefinite;
    minimize(eta * sum_val(i, Es, curpred, bs, x) - quantum_entr(x))
    subject to
    trace(x) == 1
    cvx_end
    
    curpred = x;
end
pred = curpred;
end
