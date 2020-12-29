function print_param(param, param_names, str)

    % Print current parameters
    
    param_names_comma = strjoin(param_names,',');
    fprintf(['%s%s%s%s' repmat('%6.4f ',1,length(param)) '%s\n'], str, ' [', param_names_comma, '] = [',param,']');

end