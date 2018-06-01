function data = multiscatter(options, wavelength, alt, rho_div, rho_fov, ...
                             range, ext, radius, S, ext_air, ssa, g, ...
                             ssa_air, droplet_fraction, ...
                             pristine_ice_fraction, ...
                             bscat_AD, bscat_air_AD);
% MULTISCATTER  Perform multiple scattering calculation for radar or lidar
%
%  This function calculates the apparent backscatter coefficient
%  accounting for multiple scattering, calling the "multiscatter" C
%  code.  It can be called in a number of different ways.
%
%  For lidar multiple scattering including both small-angle (SA) and
%  wide-angle (WA) scattering:
%
%    data = multiscatter(options, wavelength, alt, rho_div, rho_fov, ...
%                        range, ext, radius, S, ext_air, ssa, g, ...
%                        ssa_air, droplet_fraction, pristine_ice_fraction)
%
%  where the input variables are:
%    options: string containing space-separated settings for the code -
%      for the default options use the empty string ''
%    wavelength: instrument wavelength in m
%    alt: instrument altitude in m
%    rho_div: half-angle 1/e beam divergence in radians
%    rho_fov: half-angle receiver field-of-view in radians (or vector of FOVs)
%    range: vector of ranges at which the input data are located, in
%       m; the first range should be the closest to the instrument
%    ext: vector of extinction coefficients, in m-1
%    radius: vector of equivalent-area radii, in m
%    S: vector of backscatter-to-extinction ratios, in sr
%    ext_air: vector of molecular extinction coefficients, in m-1
%    ssa: vector of particle single-scatter albedos
%    g: vector of scattering asymmetry factors
%    ssa_air: vector of molecular single-scatter albedos
%    droplet_fraction: vector of fractions (between 0 and 1) of the particle
%       backscatter that is due to droplets - only affects small-angle
%       returns
%    pristine_ice_fraction: vector of fractions (between 0 and 1) of the
%       particle backscatter that is due to pristine ice - only affects
%       small-angle returns
%    bscat_AD: adjoint input for backscatter values
%    bscat_air_AD: adjoint input for the backscatter of air
%
%  For radar multiple scattering use the same form but with options
%  containing the string '-no-forward-lobe'
%
%  If droplet_fraction and pristine_ice_fraction are omitted then they are
%  assumed zero, i.e. the near-backscatter phase function is assumed
%  isotropic.  If "ssa_air" is also omitted it is assumed to be 1 if
%  wavelength < 1e-6 m and 0 if wavelength > 1e-6. If "ssa" and "g" are
%  omitted then the WA part of the calculation is not performed. If
%  "ext_air" is omitted then it is assumed to be zero.
%
%  The "options" string contains a space-separated list of arguments
%  passed directly to the executable.  Valid ones to use from matlab are:
%
% General options
%  -single-only    Single scattering only
%  -qsa-only       Don't include wide-angle multiple scattering
%  -wide-only      Only wide-angle multiple scattering
%  -hsrl           Output particulate and air backscatter separately
%  -gaussian-receiver
%                  Receiver is Gaussian rather than top-hat shaped
%  -adjoint        Output the adjoint as well
%  -numerical-jacobian
%                  Output the Jacobian instead
%  -jacobian       Approximate but faster Jacobian: only for small-angle
%  -ext-only       Only calculate the Jacobian with respect to extinction
%
% Options for quasi-small-angle (QSA) multiple scattering calculation
%  -explicit n     Use an explicit model with n orders of scattering
%  -fast-qsa       Use fast O(N) QSA model
%  -wide-angle-cutoff <theta>
%                  Forward scattering at angles greater than <theta> radians
%                  are deemed to escape, a crude way to deal with a problem
%                  associated with aerosols
%  -approx-exp     Appriximate the exp() function call for speed
%
% Options for wide-angle multiple-scattering
%  -no-forward-lobe
%                  Radar-like phase function behaviour: use single-scattering
%                  rather than QSA
%  -simple-2s-coeffts
%                  Use the simple upwind Euler formulation (can be unstable
%                  for high optical depth)
%  -ssa-scales-forward-lobe
%                  Single-scattering albedo less than unity reduces the
%                  forward scattering lobe as well as wide-angle scattering
%  -num-samples m  Output m samples, allowing sampling of signals appearing
%                  to originate below ground
%
%  The output is written to the structure "data" containing the member
%  "bscat" which is the apparent backscatter coefficient. It also
%  includes "range", "ext" and "radius", which are the same as the input
%  values. 
%
%  If the "-adjoint" option is specified then the adjoint variables
%  ext_AD, ssa_AD, g_AD and ext_bscat_ratio_AD will also be returned in
%  data based on the input variables bscat_AD and bscat_air_AD.
%
%  Note that this function may need to be edited to indicate where the
%  "multiscatter" executable file is located.

% Location of multiscatter executable
multiscatter_exec = '../bin/multiscatter';

% Check existence of executable
if ~exist(multiscatter_exec, 'file')
  error([multiscatter_exec ' not found']);
end

% If fewer than 9 arguments, show the help
if nargin < 9
  if nargin == 1
    if strcmp(options, 'executable')
      data = multiscatter_exec;
      return
    end
  end
  help multiscatter
  return
end

% Useful vectors
z = zeros(size(range));
o = ones(size(range));

% Allow for S and radius to be a single number
if length(S) == 1
  S = S.*o;
end
if length(radius) == 1
  radius = radius.*o;
end
% Interpret some of the options
nfov = length(rho_fov);
is_hsrl = (length(strfind(options, '-hsrl')) > 0);
is_adjoint = (length(strfind(options, '-adjoint')) > 0);
is_jacobian = (length(strfind(options, '-jacobian')) > 0);

% Open the input file and write the comments and the first line
fid = fopen('MSCATTER.DAT', 'w');
fprintf(fid, '# This file contains an input profile of atmospheric properties for\n');
fprintf(fid, '# the multiscatter code. The comand-line should be of the form:\n');
fprintf(fid, '#   ./multiscatter [options] [input_file] > [output_file]\n');
fprintf(fid, '# or\n');
fprintf(fid, '#   ./multiscatter [options] < [input_file] > [output_file]\n');
if ~isempty(options) 
  fprintf(fid, ['# Suitable options for this input file are:\n#   ' options '\n']);
else
  fprintf(fid, ['# No options are necessary for this particular file.\n']);
end
fprintf(fid, '# The file format consists of any number of comment lines starting\n');
fprintf(fid, '# with "#", followed by a line of 5 data values:\n');
fprintf(fid, '#   1. number of points in profile\n');
fprintf(fid, '#   2. instrument wavelength (m)\n');
fprintf(fid, '#   3. instrument altitude (m)\n');
fprintf(fid, '#   4. transmitter 1/e half-width (radians)\n');
fprintf(fid, '#  5+. receiver half-widths (radians)\n');
fprintf(fid, '#      (1/e half-width if the "-gaussian-receiver" option\n');
fprintf(fid, '#       is specified, top-hat half-width otherwise)\n'); 
fprintf(fid, '# The subsequent lines contain 4 or more elements:\n');
fprintf(fid, '#   1. range above ground, starting with nearest point to instrument (m)\n');
fprintf(fid, '#   2. extinction coefficient of cloud/aerosol only (m-1)\n');
fprintf(fid, '#   3. equivalent-area particle radius of cloud/aerosol (m)\n');
fprintf(fid, '#   4. extinction-to-backscatter ratio of cloud/aerosol (sterad)\n');
fprintf(fid, '#   5. extinction coefficient of air (m-1) (default 0)\n');
fprintf(fid, '#   6. single scattering albedo of cloud/aerosol\n');
fprintf(fid, '#   7. scattering asymmetry factor of cloud/aerosol\n');
fprintf(fid, '#   8. single scattering albedo of air (isotropic scattering)\n');
fprintf(fid, '#   9. fraction of cloud/aerosol backscatter due to droplets\n');
fprintf(fid, '#  10. fraction of cloud/aerosol backscatter due to pristine ice\n');
fprintf(fid, '# Note that elements 6-8 correspond to the wide-angle\n');
fprintf(fid, '# multiple-scattering calculation and if this is omitted then only\n');
fprintf(fid, '# the small-angle multiple-scattering calculation is performed.\n');
fprintf(fid, '# For more help on how to run the code, type:\n');
fprintf(fid, '#   ./multiscatter -help\n');
fprintf(fid, '%d %g %g %g', ...
        [length(range) wavelength alt rho_div]);
fprintf(fid, ' %g', rho_fov);
fprintf(fid, '\n');
% Set variables that may be missing
if nargin < 15
  pristine_ice_fraction = z;
  if nargin < 14
    droplet_fraction = z;
    if nargin < 13
      if wavelength > 1e-6
        ssa_air = z;
      else
        ssa_air = o;
      end
    end
  end
end

if nargin >= 12
  % Allow g, ssa, ssa_air and *_fraction to be single values
  if length(g) == 1
    g = g.*o;
  end
  if length(ssa) == 1
    ssa = ssa.*o;
  end
  if length(ssa_air) == 1
    ssa_air = ssa_air.*o;
  end
  if length(droplet_fraction) == 1
    droplet_fraction = droplet_fraction.*o;
  end
  if length(pristine_ice_fraction) == 1
    pristine_ice_fraction = pristine_ice_fraction.*o;
  end
  if nargin > 15
    % Adjoint
    if ~is_hsrl 
      bscat_air_AD = zeros(size(bscat_AD));
    end
    for ii = 1:length(range)
      fprintf(fid, '%g ', ...
              [range(ii) ext(ii) radius(ii) S(ii) ext_air(ii) ssa(ii) g(ii) ...
               ssa_air(ii) droplet_fraction(ii) pristine_ice_fraction(ii) ...
               bscat_AD(ii,:) bscat_air_AD(ii,:)]');
      fprintf(fid, '\n');
    end

  else
    % Wide-angle calculation
    fprintf(fid, '%g %g %g %g %g %g %g %g %g %g\n', ...
            [range(:) ext(:) radius(:) S(:) ext_air(:) ssa(:) g(:) ...
             ssa_air(:) droplet_fraction(:) pristine_ice_fraction(:)]');
  end
elseif nargin >= 10
  % Small-angle calculation only, with specified molecular extinction
  fprintf(fid, '%g %g %g %g %g\n', ...
          [range(:) ext(:) radius(:) S(:) ext_air(:)]');
else
  % Small-angle calculation only, with no molecular extinction
  fprintf(fid, '%g %g %g %g\n', ...
          [range(:) ext(:) radius(:) S(:)]');
end

fclose(fid);

% Run executable and read output
command_line = [multiscatter_exec ' ' options ' < MSCATTER.DAT 2> STDERR.TXT'];
disp(['! ' command_line]);

[status, output] = unix(command_line);
if status
  output
  unix('cat STDERR.TXT');
  error('An error occurred in running the executable');
end

try
  output = str2num(output);
catch
  output
  unix('cat STDERR.TXT');
  error('Problem interpretting output as a table of numbers');
end

if ~is_jacobian
  % Create a structure containing the output variables
  data.range = output(:,2);
  data.ext = output(:,3);
  data.radius = output(:,4);
  index = 4;
  data.bscat = output(:,index+[1:nfov]);
  index = 4+nfov;
  if is_hsrl
    data.bscat_air = output(:,index+[1:nfov]);
    index = index+nfov;
  end
  if is_adjoint
    data.ext_AD = output(:,index+1);
    data.ssa_AD = output(:,index+2);
    data.g_AD = output(:,index+3);
    data.ext_bscat_ratio_AD = output(:,index+4);
  end
else
  n = length(range);
  if is_hsrl
    data.d_bscat_d_ext = output(1:n,1:end/2);
    data.d_bscat_air_d_ext = output(1:n,end/2+1:end);
  else
    data.d_bscat_d_ext = output(1:n,1:end);
  end

  if size(output,1) > n
    if is_hsrl
      data.d_bscat_d_ssa = output(n+[1:n],1:end/2);
      data.d_bscat_air_d_ssa = output(n+[1:n],end/2+1:end);
      data.d_bscat_d_g = output(2*n+[1:n],1:end/2);
      data.d_bscat_air_d_g = output(2*n+[1:n],end/2+1:end);
      data.d_bscat_d_radius = output(3*n+[1:n],1:end/2);
      data.d_bscat_air_d_radius = output(3*n+[1:n],end/2+1:end);
      data.d_bscat_d_ext_bscat_ratio = output(4*n+1, 1:end/2);
    else
      data.d_bscat_d_ssa = output(n+[1:n],1:end);
      data.d_bscat_d_g = output(2*n+[1:n],1:end);
      data.d_bscat_d_radius = output(3*n+[1:n],1:end);
      data.d_bscat_d_ext_bscat_ratio = output(4*n+1,1:end);
    end
  end
end

