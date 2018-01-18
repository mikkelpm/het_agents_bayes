function [residual, g1, g2, g3] = example1_new_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(6, 1);
T14 = y(8)^(1+params(6));
T30 = y(5)*exp(y(9))/(exp(y(12))*y(11));
T45 = exp(y(7))*y(1)^params(3);
T46 = y(8)^(1-params(3));
T47 = T45*T46;
lhs =y(5)*params(5)*T14;
rhs =(1-params(3))*y(4);
residual(1)= lhs-rhs;
lhs =y(6);
rhs =params(1)*T30*(params(3)*exp(y(12))*y(10)+y(6)*(1-params(4)));
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T47;
residual(3)= lhs-rhs;
lhs =y(6);
rhs =exp(y(9))*(y(4)-y(5))+(1-params(4))*y(1);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =params(2)*y(2)+params(7)*y(3)+params(8)*x(it_, 1);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(2)*params(7)+params(2)*y(3)+params(9)*x(it_, 2);
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 14);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-(1-params(3)));
  g1(1,5)=params(5)*T14;
  g1(1,8)=y(5)*params(5)*getPowerDeriv(y(8),1+params(6),1);
  g1(2,10)=(-(params(1)*T30*params(3)*exp(y(12))));
  g1(2,5)=(-(params(1)*(params(3)*exp(y(12))*y(10)+y(6)*(1-params(4)))*exp(y(9))/(exp(y(12))*y(11))));
  g1(2,11)=(-(params(1)*(params(3)*exp(y(12))*y(10)+y(6)*(1-params(4)))*(-(y(5)*exp(y(9))*exp(y(12))))/(exp(y(12))*y(11)*exp(y(12))*y(11))));
  g1(2,6)=1-params(1)*T30*(1-params(4));
  g1(2,9)=(-(params(1)*T30*(params(3)*exp(y(12))*y(10)+y(6)*(1-params(4)))));
  g1(2,12)=(-(params(1)*((params(3)*exp(y(12))*y(10)+y(6)*(1-params(4)))*(-(y(5)*exp(y(9))*exp(y(12))*y(11)))/(exp(y(12))*y(11)*exp(y(12))*y(11))+T30*params(3)*exp(y(12))*y(10))));
  g1(3,4)=1;
  g1(3,1)=(-(T46*exp(y(7))*getPowerDeriv(y(1),params(3),1)));
  g1(3,7)=(-T47);
  g1(3,8)=(-(T45*getPowerDeriv(y(8),1-params(3),1)));
  g1(4,4)=(-exp(y(9)));
  g1(4,5)=exp(y(9));
  g1(4,1)=(-(1-params(4)));
  g1(4,6)=1;
  g1(4,9)=(-(exp(y(9))*(y(4)-y(5))));
  g1(5,2)=(-params(2));
  g1(5,7)=1;
  g1(5,3)=(-params(7));
  g1(5,13)=(-params(8));
  g1(6,2)=(-params(7));
  g1(6,3)=(-params(2));
  g1(6,9)=1;
  g1(6,14)=(-params(9));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,2744);
end
end
end
end
