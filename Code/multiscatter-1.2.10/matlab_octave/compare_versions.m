% This matlab/octave script compares different versions of 
% the code on the same profile

is_radar = 0;
file = 'benchmark_profile50_aniso.in';
%file = 'benchmark_profile50.in';
%file = 'calipso_cirrus_aerosol_example.in';
%file = 'hogan2008fig8.in'
%file = 'benchmark_profile_irregular5_AD.in';
%file = 'cloudsat_radar_scenario3.in'; is_radar = 1;
file = 'cloudsat_radar_scenario3a.in'; is_radar = 1;
multiscatter_new = '../bin/multiscatter';
%multiscatter_old = '../../multiscatter-1.1.26/multiscatter';
%old_data_dir = '../../multiscatter-1.1.26';
multiscatter_old = '../../multiscatter-1.2.7/bin/multiscatter';
old_data_dir = '../examples';

if ~is_radar
  options_new = {'-algorithms fast tdts','-algorithms fast none'};
  options_old = {'-algorithms fast tdts','-algorithms fast none'};
%  options_old = {'-fast-qsa','-qsa-only -fast-qsa'};
  names = {'Full new','Full old', 'Small-angle new', 'Small-angle old',...
	   'Wide-angle new','Wide-angle old'};
else
%  options_new = {'-algorithms fast tdts'};
%  options_old = {'-no-forward-lobe'};
  options_new = {'-algorithms single tdts','-algorithms single none'};
%  options_old = {'-no-forward-lobe','-single-only'};
  options_old = {'-algorithms single tdts','-algorithms single none'};
  names = {'Full new','Full old','Single new','Single old',...
	   'Wide-angle new','Wide-angle old'};
end

styles = {'k','r','b', 'kx', 'ro','bs'};


for iopt = 1:length(options_new)
  command_line = [multiscatter_new ' ' options_new{iopt} ...
		  ' ../examples/' file ' 2>/dev/null'];
  disp(['Executing ' command_line]);
  [status, data_new{iopt}] = unix(command_line);
  data_new{iopt} = str2num(data_new{iopt});

  command_line = [multiscatter_old ' ' options_old{iopt} ...
		  ' ' old_data_dir '/' file ' 2>/dev/null'];
  disp(['Executing ' command_line]);
  [status, data_old{iopt}] = unix(command_line);
  data_old{iopt} = str2num(data_old{iopt});
end

clf
for iopt = 1:length(data_new)
  semilogx(data_new{iopt}(:,5), data_new{iopt}(:,2), styles{iopt});
  hold on;
  semilogx(data_old{iopt}(:,5), data_old{iopt}(:,2), styles{iopt+3});
end
semilogx(data_new{1}(:,5)-data_new{2}(:,5), data_new{1}(:,2),styles{3});
semilogx(data_old{1}(:,5)-data_old{2}(:,5), data_old{1}(:,2),styles{6});

legend(names,2);
