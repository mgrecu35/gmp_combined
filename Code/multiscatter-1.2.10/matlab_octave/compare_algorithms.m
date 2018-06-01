% This matlab/octave script compares different versions of 
% the algorithm on the same profile

file = 'benchmark_profile50a_aniso.in';
%file = 'calipso_cirrus_aerosol_example.in';
%file = 'hogan2008fig8.in'
%file = 'benchmark_profile_irregular5_AD.in';
%file = 'cloudsat_radar_scenario3.in';
multiscatter = '../bin/multiscatter';

% Note that the fortran code requires all entries to be present,
% but some *.in files omit some entries
do_fortran = 1;

general_options = '';
options = {'-algorithms fast tdts',...
	   '-algorithms fast none',...
	   '-algorithms none tdts',...
	   '-algorithms fast tdts -optimize-wide-angle-gates',...
	   '-algorithms none tdts -optimize-wide-angle-gates'}
styles = {'k','r','b', 'ko', 'bo','gx'};

names = {'Full','Small-angle only','Wide-angle only','Optimized wide-angle','Optimized wide-angle only','Fortran version'};

for iopt = 1:length(options)
  command_line = [multiscatter ' ' general_options ' ' options{iopt} ...
		  ' ../examples/' file ' 2>/dev/null'];
  disp(['Executing ' command_line]);
  [status, data{iopt}] = unix(command_line);
  data{iopt} = str2num(data{iopt});
end

if do_fortran
  % Fortran version
  iopt = iopt+1;
  command_line = ['../bin/multiscatterf < ../examples/' file ' 2>/dev/null'];
  disp(['Executing ' command_line]);
  [status, data{iopt}] = unix(command_line);
  data{iopt} = str2num(data{iopt});
end


clf
for iopt = 1:length(data)
  semilogx(data{iopt}(:,5), data{iopt}(:,2), styles{iopt});
  hold on;
end
legend(names,2);
