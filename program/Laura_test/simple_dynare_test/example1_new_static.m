function [residual, g1, g2, g3] = example1_new_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 6, 1);

%
% Model equations
%

T14 = y(5)^(1+params(6));
T37 = exp(y(4))*y(3)^params(3);
T38 = y(5)^(1-params(3));
T39 = T37*T38;
lhs =y(2)*params(5)*T14;
rhs =(1-params(3))*y(1);
residual(1)= lhs-rhs;
lhs =y(3);
rhs =params(1)*(y(1)*params(3)*exp(y(6))+y(3)*(1-params(4)));
residual(2)= lhs-rhs;
lhs =y(1);
rhs =T39;
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(3)*(1-params(4))+exp(y(6))*(y(1)-y(2));
residual(4)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(2)+y(6)*params(7)+params(8)*x(1);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(4)*params(7)+y(6)*params(2)+params(9)*x(2);
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-(1-params(3)));
  g1(1,2)=params(5)*T14;
  g1(1,5)=y(2)*params(5)*getPowerDeriv(y(5),1+params(6),1);
  g1(2,1)=(-(params(1)*params(3)*exp(y(6))));
  g1(2,3)=1-params(1)*(1-params(4));
  g1(2,6)=(-(params(1)*y(1)*params(3)*exp(y(6))));
  g1(3,1)=1;
  g1(3,3)=(-(T38*exp(y(4))*getPowerDeriv(y(3),params(3),1)));
  g1(3,4)=(-T39);
  g1(3,5)=(-(T37*getPowerDeriv(y(5),1-params(3),1)));
  g1(4,1)=(-exp(y(6)));
  g1(4,2)=exp(y(6));
  g1(4,3)=1-(1-params(4));
  g1(4,6)=(-(exp(y(6))*(y(1)-y(2))));
  g1(5,4)=1-params(2);
  g1(5,6)=(-params(7));
  g1(6,4)=(-params(7));
  g1(6,6)=1-params(2);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,216);
end
end
end
end
