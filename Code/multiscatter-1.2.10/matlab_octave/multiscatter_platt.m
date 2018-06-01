function bscat = multiscatter_platt(drange, ext, S, ext_air, mu)
% MULTISCATTER_PLATT  Lidar multiple scattering using Platt's approximation
%   bscat = multiscatter_platt(range, ext, S, ext_air, mu)
% where
%   drange is the distance between each point (m)
%   ext is the extinction coefficient (m-1)
%   S is the extinction-to-backscatter ratio (sr)
%   ext_air is the extinction coefficient of the air (m-1)
%   mu is the "multiple scatter coefficient" between 1 (single
%     scattering) and 0.5 (all forward scattered photons remaining in
%     receiver field-of-view)
%
%   bscat is the apparent backscatter coefficient (m-1 sr-1)

  if nargin < 5
    help multiscatter_platt
    return;
  end
  
  % Calculate effective total extinction
  ext_eff = mu.*ext + ext_air;
  % Calculate cumulative optical depth
  od = cumsum([0 ext_eff(1:end-1).*drange]);
  % Calculate the effective transmission to each gate, accounting for
  % attenuation within the gate in question, and forward scattered
  % photons remaining within the receiver field-of-view
  trans = exp(-2.*od).*(1-exp(-2.*ext_eff.*drange)) ...
          ./(2.*ext_eff.*drange);
  % Correctly deal with gates containing no extinction
  index = find(~finite(trans));
  trans(index) = exp(-2.*od(index));
  % Calculate the apparent backscatter
  bscat = trans.*(ext./S + ext_air./(8.*pi./3));
  
