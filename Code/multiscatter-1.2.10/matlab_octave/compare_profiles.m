% This matlab/octave script compares different profiles.

% Compare a profile on a regular grid with the same profile but
% on an irregular grid; this tests the part of the algorithm
% that interpolates irregular grids on to regular ones for the 
% wide-angle part
files_in = {'regular_profile.in', 'irregular_profile.in'};

% Compare the same profile but at different resolutions; this
% checks the numerical schemes 
%files_in = {'regular_profile.in', 'regular_profile2.in'};

% Perform the following algorithm combinations on each input profile
options_in = {'-algorithms single none', '-algorithms fast tdts'};

files = {files_in{1}, files_in{1}, files_in{2}, files_in{2}};
multiscatter = '../bin/multiscatter';

general_options = '';
options = {options_in{1}, options_in{2}, options_in{1}, options_in{2}};
names = {'Profile 1 single','Profile 1 full','Profile 2 single','Profile 2 full'};
styles = {'k','b','r--','m--'};

for ifile = 1:length(files)
	      command_line = [multiscatter ' ' general_options ' ' options{ifile} ...
			       ' ../examples/' files{ifile} ' 2>/dev/null'];
  disp(['Executing ' command_line]);
  [status, data{ifile}] = unix(command_line);
  data{ifile} = str2num(data{ifile});
end

clf
for ifile = 1:length(files)
  semilogx(data{ifile}(:,5), data{ifile}(:,2), styles{ifile});
  hold on;
end
legend(names);
