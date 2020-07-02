function Sigma_Z = cov_smpl(Sigma_X)

    % Compute var-cov matrix of Z=(X|E=0,X^2|E=0,X^3|E=0,X|E=1,X^2|E=1,X^3|E=1)',
    % where X|E=0 ~ N(0,Sigma_X(1)), X|E=1 ~ N(0,Sigma_X(2)) (approximately)
    
    % https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Higher_moments
    global nEpsilon nMeasure aggEmployment
    
    Sigma_Z = zeros(nMeasure*nEpsilon);
    
    % Diagonal
    Sigma_Z(1,1) = Sigma_X(1)/(1-aggEmployment); % Var(X|E=0)/P(E=0)
    Sigma_Z(2,2) = 2*Sigma_X(1)^2/(1-aggEmployment); % Var(X^2|E=0)/P(E=0)
    Sigma_Z(3,3) = 15*Sigma_X(1)^3/(1-aggEmployment); % Var(X^3|E=0)/P(E=0)
    Sigma_Z(nMeasure+1,nMeasure+1) = Sigma_X(2)/aggEmployment; % Var(X|E=1)/P(E=1)
    Sigma_Z(nMeasure+2,nMeasure+2) = 2*Sigma_X(2)^2/aggEmployment; % Var(X^2|E=1)/P(E=1)
    Sigma_Z(nMeasure+3,nMeasure+3) = 15*Sigma_X(2)^3/aggEmployment; % Var(X^3|E=1)/P(E=1)
    
    % Non-zero upper off-diagonal elements
    Sigma_Z(1,3) = 3*Sigma_X(1)^2/(1-aggEmployment); % Cov(X,X^3|E=0)/P(E=0)
    Sigma_Z(nMeasure+1,nMeasure+3) = 3*Sigma_X(2)^2/aggEmployment; % Cov(X,X^3|E=1)/P(E=1)
    
    % Copy upper off-diagonal part to lower part
    lower_ind = tril(ones(size(Sigma_Z)));
    Sigma_Z_t = Sigma_Z';
    Sigma_Z(lower_ind==1) = Sigma_Z_t(lower_ind==1);

end