function Sigma_Z = cov_smpl(mom_X)

    % Compute var-cov matrix of Z=(X,X^2,X^3)', E[X] = 0
    % Fisher (1930, PLMS)
    
    % For moment-based methods up to the 3rd moments. 
    % NEED to change for moment-based methods with higher order moments.
    
    Sigma_Z = zeros(length(mom_X)/2);
    
    % Diagonal
    Sigma_Z(1,1) = mom_X(2); % Var(X)
    Sigma_Z(1,2) = mom_X(3); % Cov(X,X^2)
    Sigma_Z(1,3) = mom_X(4)-3*mom_X(2)^2; % Cov(X,X^3)
    Sigma_Z(2,2) = mom_X(4)-mom_X(2)^2; % Var(X^2)
    Sigma_Z(2,3) = mom_X(5)-4*mom_X(3)*mom_X(2); % Cov(X^2,X^3)
    Sigma_Z(3,3) = mom_X(6)-6*mom_X(4)*mom_X(2)-mom_X(3)^2+9*mom_X(2)^3; % Var(X^3)
    
    % Copy upper off-diagonal part to lower part
    lower_ind = tril(ones(size(Sigma_Z)));
    Sigma_Z_t = Sigma_Z';
    Sigma_Z(lower_ind==1) = Sigma_Z_t(lower_ind==1);

end