function write_data_file_total_biomass(T)

time = num2cell([0,T]);
time{1} = 'time';
p = num2cell(zeros(length(time),1)');
p{1} = 'p';
p{end} = 1;
c = [time;p]';

%   Write into file
fid = fopen('cellgrowth_total_biomass.csv', 'w') ;
fprintf(fid, '%s,', c{1,1:1}) ;
fprintf(fid, '%s\n', c{1,2}) ;
fclose(fid) ;
dlmwrite('cellgrowth_total_biomass.csv', c(2:end,:), '-append') ;

end
