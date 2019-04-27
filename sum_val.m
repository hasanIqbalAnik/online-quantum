function sv = sum_val(t, meas_ops, pred_st, loss_q, phi)
sv = 0;
for i = 1:t
    lpt = 2 * (trace(meas_ops{i} * pred_st) - loss_q(i)) * meas_ops{i};    
    sv = sv + trace(lpt * phi);
end
end
