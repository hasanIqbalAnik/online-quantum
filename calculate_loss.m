function cl = calculate_loss(Es, bs, x, T)
total = 0;
for i = 1:T
    total = total + (trace(Es{i} * x) - bs(i));
end
cl = total;
end