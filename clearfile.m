function clearfile = clearfile(file)
fileID = fopen(file, 'w');
fclose(fileID)
end

