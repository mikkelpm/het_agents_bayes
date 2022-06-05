doc_fig_path = pwd;
addpath(doc_fig_path)
cd ../../program

v_model_name = {'hh','firm'};
n_model = length(v_model_name);

addpath(genpath('./functions'));
for i_model = 1:n_model
    model_name = v_model_name{i_model};
    addpath(genpath(['./' model_name '_model']));
    
    plot_depfun(['run_mcmc_' model_name],[model_name '_filemap'])
    
    rmpath(genpath(['./' model_name '_model']));
end
rmpath(genpath('./functions'));
rmpath(doc_fig_path)