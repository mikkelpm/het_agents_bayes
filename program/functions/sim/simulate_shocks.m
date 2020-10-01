function shocks = simulate_shocks(M_, T, numrep)

    % Simulate normally distributed shocks

    cholSigma = chol(M_.Sigma_e);
    shocks = randn(numrep*T,M_.exo_nbr)*cholSigma;
    
    % Reshape if more the one repetition is requested
    if numrep>1
        shocks = permute(reshape(shocks, [numrep T M_.exo_nbr]), [2 3 1]);
    end

end