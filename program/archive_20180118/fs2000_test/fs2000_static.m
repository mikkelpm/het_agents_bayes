function [residual, g1, g2, g3] = fs2000_static(y, x, params)
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

residual = zeros( 15, 1);

%
% Model equations
%

T43 = (-(params(6)/(1-params(6))));
T55 = exp((params(3)+x(1))*(-params(1)));
T58 = y(7)^params(1);
T60 = y(9)^(-params(1));
T70 = y(9)^(1-params(1));
T71 = T58*T55*(1-params(1))*y(2)*params(2)*T70;
T114 = y(7)^(params(1)-1);
T120 = T70*params(1)*exp((-params(1))*(params(3)+log(y(4))))*T114+(1-params(7))*exp((-(params(3)+log(y(4)))));
T121 = y(2)*params(2)*T120;
T122 = T121/(y(1)*y(2)*y(3));
lhs =y(14);
rhs =exp(params(3)+x(1));
residual(1)= lhs-rhs;
lhs =log(y(1));
rhs =(1-params(5))*log(params(4))+log(y(1))*params(5)+x(2);
residual(2)= lhs-rhs;
residual(3) = (-y(2))/(y(1)*y(2)*y(3))+T122;
lhs =y(5);
rhs =y(10)/y(9);
residual(4)= lhs-rhs;
residual(5) = y(10)/y(9)+T43*y(2)*y(3)/(1-y(9));
lhs =y(6);
rhs =y(2)*(1-params(1))*T55*T58*T60/y(5);
residual(6)= lhs-rhs;
residual(7) = 1/(y(2)*y(3))-T71/(y(2)*y(3)*y(1)*y(10));
lhs =y(3)+y(7);
rhs =T70*T55*T58+y(7)*(1-params(7))*exp((-(params(3)+x(1))));
residual(8)= lhs-rhs;
lhs =y(2)*y(3);
rhs =y(1);
residual(9)= lhs-rhs;
lhs =y(1)-1+y(8);
rhs =y(10);
residual(10)= lhs-rhs;
lhs =y(4);
rhs =exp(x(1));
residual(11)= lhs-rhs;
lhs =y(13);
rhs =T55*T58*T70;
residual(12)= lhs-rhs;
lhs =y(11);
rhs =y(14)*y(13)/y(13);
residual(13)= lhs-rhs;
lhs =y(12);
rhs =y(1)/y(14);
residual(14)= lhs-rhs;
lhs =y(15);
rhs =T122;
residual(15)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(15, 15);

  %
  % Jacobian matrix
  %

T158 = (y(1)*y(2)*y(3)*params(2)*T120-T121*y(1)*y(3))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3));
T208 = y(2)*params(2)*(T70*T114*params(1)*exp((-params(1))*(params(3)+log(y(4))))*(-params(1))*1/y(4)+(1-params(7))*exp((-(params(3)+log(y(4)))))*(-(1/y(4))))/(y(1)*y(2)*y(3));
T218 = y(2)*params(2)*T70*params(1)*exp((-params(1))*(params(3)+log(y(4))))*getPowerDeriv(y(7),params(1)-1,1)/(y(1)*y(2)*y(3));
T219 = getPowerDeriv(y(7),params(1),1);
T236 = getPowerDeriv(y(9),1-params(1),1);
T239 = y(2)*params(2)*params(1)*exp((-params(1))*(params(3)+log(y(4))))*T114*T236/(y(1)*y(2)*y(3));
  g1(1,14)=1;
  g1(2,1)=1/y(1)-params(5)*1/y(1);
  g1(3,1)=(-((-y(2))*y(2)*y(3)))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3))+(-(y(2)*y(3)*T121))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3));
  g1(3,2)=((-(y(1)*y(2)*y(3)))-(-y(2))*y(1)*y(3))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3))+T158;
  g1(3,3)=(-((-y(2))*y(1)*y(2)))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3))+(-(T121*y(1)*y(2)))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3));
  g1(3,4)=T208;
  g1(3,7)=T218;
  g1(3,9)=T239;
  g1(4,5)=1;
  g1(4,9)=(-((-y(10))/(y(9)*y(9))));
  g1(4,10)=(-(1/y(9)));
  g1(5,2)=T43*y(3)/(1-y(9));
  g1(5,3)=T43*y(2)/(1-y(9));
  g1(5,9)=(-y(10))/(y(9)*y(9))+T43*y(2)*y(3)/((1-y(9))*(1-y(9)));
  g1(5,10)=1/y(9);
  g1(6,2)=(-(T60*T58*(1-params(1))*T55/y(5)));
  g1(6,5)=(-((-(y(2)*(1-params(1))*T55*T58*T60))/(y(5)*y(5))));
  g1(6,6)=1;
  g1(6,7)=(-(T60*y(2)*(1-params(1))*T55*T219/y(5)));
  g1(6,9)=(-(y(2)*(1-params(1))*T55*T58*getPowerDeriv(y(9),(-params(1)),1)/y(5)));
  g1(7,1)=(-((-(T71*y(2)*y(3)*y(10)))/(y(2)*y(3)*y(1)*y(10)*y(2)*y(3)*y(1)*y(10))));
  g1(7,2)=(-y(3))/(y(2)*y(3)*y(2)*y(3))-(y(2)*y(3)*y(1)*y(10)*T70*T58*T55*(1-params(1))*params(2)-T71*y(3)*y(1)*y(10))/(y(2)*y(3)*y(1)*y(10)*y(2)*y(3)*y(1)*y(10));
  g1(7,3)=(-y(2))/(y(2)*y(3)*y(2)*y(3))-(-(T71*y(2)*y(1)*y(10)))/(y(2)*y(3)*y(1)*y(10)*y(2)*y(3)*y(1)*y(10));
  g1(7,7)=(-(T70*T55*(1-params(1))*y(2)*params(2)*T219/(y(2)*y(3)*y(1)*y(10))));
  g1(7,9)=(-(T58*T55*(1-params(1))*y(2)*params(2)*T236/(y(2)*y(3)*y(1)*y(10))));
  g1(7,10)=(-((-(T71*y(2)*y(1)*y(3)))/(y(2)*y(3)*y(1)*y(10)*y(2)*y(3)*y(1)*y(10))));
  g1(8,3)=1;
  g1(8,7)=1-((1-params(7))*exp((-(params(3)+x(1))))+T70*T55*T219);
  g1(8,9)=(-(T55*T58*T236));
  g1(9,1)=(-1);
  g1(9,2)=y(3);
  g1(9,3)=y(2);
  g1(10,1)=1;
  g1(10,8)=1;
  g1(10,10)=(-1);
  g1(11,4)=1;
  g1(12,7)=(-(T55*T70*T219));
  g1(12,9)=(-(T55*T58*T236));
  g1(12,13)=1;
  g1(13,11)=1;
  g1(13,14)=(-1);
  g1(14,1)=(-(1/y(14)));
  g1(14,12)=1;
  g1(14,14)=(-((-y(1))/(y(14)*y(14))));
  g1(15,1)=(-((-(y(2)*y(3)*T121))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3))));
  g1(15,2)=(-T158);
  g1(15,3)=(-((-(T121*y(1)*y(2)))/(y(1)*y(2)*y(3)*y(1)*y(2)*y(3))));
  g1(15,4)=(-T208);
  g1(15,7)=(-T218);
  g1(15,9)=(-T239);
  g1(15,15)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],15,225);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],15,3375);
end
end
end
end
