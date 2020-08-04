% Cov matrix of sample moments based on s.s.

% Moments E[lambda^j]
mom_lambda = exp((1:nMeasure*2)*mu_l + 0.5*(1:nMeasure*2).^2*(-2*mu_l)); 

% Higher-order moments of HH assets
for eepsilon = 0:1
    % prepare the parameters
    moment = nan(1,nMeasure);
    measureCoefficient = nan(1,nMeasure);
    for i_Measure = 1:nMeasure
        moment(i_Measure) = eval(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)]);
        measureCoefficient(i_Measure) = eval(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)]);
    end
    
    % unconstrained assets
    moment_aux = moment;
    moment_aux(1) = 0;
    g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux'));
    normalization = integral(g, aaBar, Inf);
    for i_mom = nMeasure+1:nMeasure*2
        mMoments(eepsilon+1,i_mom) = integral(@(a) g(a).*(a-moment(1)).^i_mom, aaBar, Inf)/normalization;
    end
end

% Higher-order moments of HH income
smpl_m(:,4) = r^4*mMoments(:,4)*mom_lambda(4)...
    +4*r^3*mMoments(:,3).*smpl_m(:,1)*(mom_lambda(4)-mom_lambda(3))...
    +6*r^2*mMoments(:,2).*smpl_m(:,1).^2*(mom_lambda(4)-2*mom_lambda(3)+mom_lambda(2))...
    +smpl_m(:,1).^4*(mom_lambda(4)-4*mom_lambda(3)+6*mom_lambda(2)-3); % Fourth central moment
smpl_m(:,5) = r^5*mMoments(:,5)*mom_lambda(5)...
    +5*r^4*mMoments(:,4).*smpl_m(:,1)*(mom_lambda(5)-mom_lambda(4))...
    +10*r^3*mMoments(:,3).*smpl_m(:,1).^2*(mom_lambda(5)-2*mom_lambda(4)+mom_lambda(3))...
    +10*r^2*mMoments(:,2).*smpl_m(:,1).^3*(mom_lambda(5)-3*mom_lambda(4)+3*mom_lambda(3)-mom_lambda(2))...
    +smpl_m(:,1).^5*(mom_lambda(5)-5*mom_lambda(4)+10*mom_lambda(3)-10*mom_lambda(2)+4); % Fifth central moment
smpl_m(:,6) = r^6*mMoments(:,6)*mom_lambda(6)...
    +6*r^5*mMoments(:,5).*smpl_m(:,1)*(mom_lambda(6)-mom_lambda(5))...
    +15*r^4*mMoments(:,4).*smpl_m(:,1).^2*(mom_lambda(6)-2*mom_lambda(5)+mom_lambda(4))...
    +20*r^3*mMoments(:,3).*smpl_m(:,1).^3*(mom_lambda(6)-3*mom_lambda(5)+3*mom_lambda(4)-mom_lambda(3))...
    +15*r^2*mMoments(:,2).*smpl_m(:,1).^4*(mom_lambda(6)-4*mom_lambda(5)+6*mom_lambda(4)-4*mom_lambda(3)+mom_lambda(2))...
    +smpl_m(:,1).^6*(mom_lambda(6)-6*mom_lambda(5)+15*mom_lambda(4)-20*mom_lambda(3)+15*mom_lambda(2)-5);  % Sixth central moment
M_.H(1,1) = ssigmaMeas^2;
M_.H(2:4,2:4) = cov_smpl(smpl_m(1,:))/N_micro/(1-aggEmployment);
M_.H(5:7,5:7) = cov_smpl(smpl_m(2,:))/N_micro/aggEmployment;