function Sigma_Z = cov_smpl(Sigma_X)

    % Compute var-cov matrix of Z=(X|E=0,X^2|E=0,X^3|E=0,X|E=1,X^2|E=1,X^3|E=1)',
    % where X|E=0 ~ N(0,Sigma_X(1)), X|E=1 ~ N(0,Sigma_X(2)) (approximately)
    
    % https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Higher_moments
    % OK for nMeasure = 2 or 3, NEED to change if nMeasure takes other values
    
    global nMeasure aggEmployment
    
    v_coef = [1 2 6]';
    Sigma_Z_aux = v_coef(1:nMeasure).*repmat(Sigma_X,nMeasure,1).^repmat((1:nMeasure)',1,length(Sigma_X))...
        ./[1-aggEmployment aggEmployment];
    Sigma_Z = diag(Sigma_Z_aux(:));

end