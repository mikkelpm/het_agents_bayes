tic
dynare fs2000 noclearall;
t = toc;
disp(['Preprocess: ' num2str(t) 's'])

rhos = 0.8:0.05:0.95;
for i = 1:length(rhos)
    tic
    set_param_value('rho',rhos(i));
    dynare_estimation(var_list_);
    t = toc;
    disp(['rho = ' num2str(rhos(i)) ': ' num2str(t) 's'])
end