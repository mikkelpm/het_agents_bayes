function varargout = load_mat(filename)

    % Load from .mat file with suffix
    
    global mat_suff; % File suffix
    the_file = strcat(filename, mat_suff);
    
    varargout = cell(1,nargout);
    
    if nargout==0
        evalin('caller', ['load ' the_file]);
    else
        varargout{1} = load(the_file);
    end

end