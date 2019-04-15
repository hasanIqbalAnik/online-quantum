function HypothesisState = LearningWPostSelection(k, Es, rho)
    Op = eye(2);
    for idx = [1,k]
        M = PostSelection(rho, Es{idx});
        Op = M * Op;
    end
    HypothesisState = Op * rho;
end
