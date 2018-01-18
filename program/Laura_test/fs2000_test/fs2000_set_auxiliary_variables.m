function y = fs2000_set_auxiliary_variables(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

y(15)=y(2)*params(2)*(y(9)^(1-params(1))*params(1)*exp((-params(1))*(params(3)+log(y(4))))*y(7)^(params(1)-1)+(1-params(7))*exp((-(params(3)+log(y(4))))))/(y(1)*y(2)*y(3));
