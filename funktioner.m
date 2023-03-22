fid = fopen('small_surface_1.txt','rt');
C = textscan(fid, '%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','	', 'HeaderLines',1);
fclose(fid);
P=cell2mat(C)
plot(C{1,1},C{1,2})