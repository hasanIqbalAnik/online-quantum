function [State, Measurements] = generateData(k, eps)
    psi = sqrt(5 / 7) * [1; 0] + sqrt(2 / 7) * [0; 1]; 

    u = psi;
    x=u(:).'/norm(u);
    yz=null(x).';
    xyz=[x;yz];  %The rows of this matrix are the axes
    
    if xyz(:, 1) == psi
        psi_n = xyz(:,2);
    else
        psi_n = xyz(:,1);
    end
    Measurements  = {1:k};
    State = psi*psi';
    for i = 1:k
        r = eps.*rand(1,1);
        meas = sqrt(1-r)*psi+sqrt(r)*psi_n;
        Measurements{i} = meas*meas';
    end
    
end
