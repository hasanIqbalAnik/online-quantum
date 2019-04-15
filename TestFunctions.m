function passed = TestFunctions(k, eps, times)
    for i = 1:times
        [rho, Es] = generateData(k, eps);
        hyp = LearningWPostSelection(k, Es, rho);
        pass = trace(hyp - rho) <= 4*sqrt(k*eps);
        if pass == 0
            passed = 0;
            return 
        end
    end
    passed = 1;
end
