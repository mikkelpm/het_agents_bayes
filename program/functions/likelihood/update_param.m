function update_param(param, param_names)

    % Update global parameters

    for i=1:length(param)
        evalin('caller', sprintf('%s%s%f%s', param_names{i}, '=', param(i), ';'));
    end

end