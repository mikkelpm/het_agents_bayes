function Sigma_Z = cov_smpl(Sigma_X)

    % Compute var-cov matrix of Z=(X_1,X_2,X_1^2,X_1*X_2,X_2^2)',
    % where (X_1,X_2)' ~ N(0,Sigma_X)
    
    % https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Higher_moments
    
    Sigma_Z = zeros(5);
    
    % Diagonal
    Sigma_Z(1,1) = Sigma_X(1,1); % Var(X_1)
    Sigma_Z(2,2) = Sigma_X(2,2); % Var(X_2)
    Sigma_Z(3,3) = 2*Sigma_X(1,1)^2; % Var(X_1^2)
    Sigma_Z(4,4) = Sigma_X(1,1)*Sigma_X(2,2) + Sigma_X(1,2)^2; % Var(X_1*X_2)
    Sigma_Z(5,5) = 2*Sigma_X(2,2)^2; % Var(X_2^2)
    
    % Non-zero upper off-diagonal elements
    Sigma_Z(1,2) = Sigma_X(1,2); % Cov(X_1,X_2)
    Sigma_Z(3,4) = 2*Sigma_X(1,1)*Sigma_X(1,2); % Cov(X_1^2,X_1*X_2)
    Sigma_Z(3,5) = 2*Sigma_X(1,2)^2; % Cov(X_1^2,X_2^2)
    Sigma_Z(4,5) = 2*Sigma_X(2,2)*Sigma_X(1,2); % Cov(X_1*X_2,X_2^2)
    
    % Copy upper off-diagonal part to lower part
    lower_ind = tril(ones(size(Sigma_Z)));
    Sigma_Z_t = Sigma_Z';
    Sigma_Z(lower_ind==1) = Sigma_Z_t(lower_ind==1);

end