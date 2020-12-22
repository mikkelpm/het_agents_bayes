function save_mat(filename, varargin)

    % Save specified variables to .mat file with suffix
    
    global mat_suff; % Suffix
    evalin('caller', ['save ' strcat(filename, mat_suff, '.mat') ' ' strjoin(varargin) ';']);

end