%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'example1_new';
M_.dynare_version = '4.5.3';
oo_.dynare_version = '4.5.3';
options_.dynare_version = '4.5.3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('example1_new.log');
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.exo_names_long = 'e';
M_.exo_names = char(M_.exo_names, 'u');
M_.exo_names_tex = char(M_.exo_names_tex, 'u');
M_.exo_names_long = char(M_.exo_names_long, 'u');
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'b');
M_.endo_names_tex = char(M_.endo_names_tex, 'b');
M_.endo_names_long = char(M_.endo_names_long, 'b');
M_.endo_partitions = struct();
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'std_e');
M_.param_names_tex = char(M_.param_names_tex, 'std\_e');
M_.param_names_long = char(M_.param_names_long, 'std_e');
M_.param_names = char(M_.param_names, 'std_u');
M_.param_names_tex = char(M_.param_names_tex, 'std\_u');
M_.param_names_long = char(M_.param_names_long, 'std_u');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 6;
M_.param_nbr = 9;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'y'};
options_.varobs(2)  = {'c'};
options_.varobs_id = [ 1 2  ];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('example1_new_static');
erase_compiled_function('example1_new_dynamic');
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 0 4 10;
 0 5 11;
 1 6 0;
 2 7 0;
 0 8 0;
 3 9 12;]';
M_.nstatic = 1;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.NNZDerivatives = [26; -1; -1];
M_.params( 3 ) = 0.36;
alpha = M_.params( 3 );
M_.params( 2 ) = 0.95;
rho = M_.params( 2 );
M_.params( 7 ) = 0.025;
tau = M_.params( 7 );
M_.params( 1 ) = 0.99;
beta = M_.params( 1 );
M_.params( 4 ) = 0.025;
delta = M_.params( 4 );
M_.params( 6 ) = 0;
psi = M_.params( 6 );
M_.params( 5 ) = 2.95;
theta = M_.params( 5 );
M_.params( 8 ) = 0.009;
std_e = M_.params( 8 );
M_.params( 9 ) = 0.009;
std_u = M_.params( 9 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 1.08068253095672;
oo_.steady_state( 2 ) = 0.80359242014163;
oo_.steady_state( 5 ) = 0.29175631001732;
oo_.steady_state( 3 ) = 11.08360443260358;
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 6 ) = 0;
oo_.exo_steady_state( 1 ) = 0;
oo_.exo_steady_state( 2 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
options_.irf = 0;
options_.order = 1;
var_list_ = char();
info = stoch_simul(var_list_);
save('example1_new_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('example1_new_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('example1_new_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('example1_new_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('example1_new_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('example1_new_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('example1_new_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
