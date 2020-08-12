function str_cell = str_add_numbers(str, numbers)

    % Add list of numbers to end of string

    str_cell = cellfun(@(x) sprintf('%s%d', str, x), num2cell(numbers), 'UniformOutput', false)';
    
end