% Macro Only model (likelihood type 2) doesn't depend on micro sample size (N),
% so we can simulate it once (with tag '_N1000') and copy its output file across different Ns

for i_rep = 1:n_rep
    
    v_file = cell(n_liktype,1);
    for i_type = plot_liktypes
        % Model file name and results
        the_filename = sprintf('%s%s%s%d%s%02d',model_name,subspec,...
            '_liktype',plot_liktypes(i_type),'_',plot_reps(i_rep));
        v_file{i_rep} = fullfile(results_folder, the_filename);
    end
    
    the_filename_N1000 = sprintf('%s%s%s%d%s%02d',model_name,'_N1000',...
        '_liktype',plot_liktypes(2),'_',plot_reps(i_rep));
    the_file_N1000 = fullfile(results_folder, the_filename_N1000);
    
    if isfile(strcat(v_file{1}, '.mat')) && ~isfile(strcat(v_file{2}, '.mat'))...
            && isfile(strcat(the_file_N1000, '.mat'))
        copyfile(strcat(the_file_N1000, '.mat'), strcat(v_file{2}, '.mat'));
    end
    
end