function data = read_data(filename)
    data = readtable(filename, 'Delimiter', 'tab' ,'VariableNamingRule','modify');
    if width(data) == 3
        data = renamevars(data, ["t__s_", "T__C_", "m__g_"], ["t", "T", "m"]);
    else 
        data = renamevars(data, ["t__s_", "T__C_"], ["T", "t"]);
    end
end
