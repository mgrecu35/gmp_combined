% Script to show the use of the multiple scattering code for calculating
% multiple field-of-view lidar returns, and the corresponding adjoint.

% THOR
alt = 7000; % Instrument altitude (m)
wavelength = 532e-9; % (m)
fov = 0.5*0.001*[0.000 0.840;...
                 1.029 1.681;...
                 1.681 3.361;...
                 3.361 6.723;...
                 6.723 13.40;...
                 13.40 26.72;...
                 26.72 53.40;...
                 53.40 106.7];
rho_div = 1e-5; % Half-angle beam divergence (radians)
rho_fov = unique(fov(2:end)); % Half-angle field-of-view (radians)

footprints = rho_fov.*alt.*2; % Receiver footprint on the ground (m)
drange = 50; % Distance between range gates (m)

disp(['Receiver footprints on the ground: ' num2str(footprints) ' m']);

% Open figure with the appropriate size to be printed nicely
figure(1)
set(gcf,'units','inches',...
        'paperposition',[0.5 0.5 5.5 7],'position',[0.5 0.5 5.5 7],...
        'defaultaxesfontsize',13,'defaulttextfontsize',13)
clf


% Liquid cloud
range = 2500:-drange:0; % (m)
ext = zeros(size(range)); % Extinction coefficient (m-1)
radius = 10e-6.*ones(size(range)); % (m)
index = find(range >= 1000 & range < 2000); % Cloud 1-2 km
ext(index) = 0.03;
% asymmetry and single-scatter albedo for 10-micron droplets:
g = [0.862617];
ssa = [1.0];
droplet_frac = ones(size(ext));
pristine_ice_frac = zeros(size(ext));

S = 18.0.*ones(size(range)); % Extinction-backscatter ratio (sr)
ext_air = 1.6e-6.*exp(-range./8000).*8.*pi./3; % Molecular extinction (m-1)
ssa_air = ones(size(ext_air)); % Single-scatter albedo of air

% Adjoint
bscat_AD = zeros(length(ext), length(rho_fov));
bscat_AD(fix(end*3/4), 2) = 1.0;

% Options
options = '';
waoptions = [' ' options];
adoptions = ['-adjoint ' options];

optical_depth = sum(ext).*abs(median(diff(range)));
disp(['Optical depth: ' num2str(optical_depth)]);
  
sa_only = multiscatter([options ' -algorithms fast none'], ...
                       wavelength, alt, rho_div, rho_fov, ...
                       range, ext, radius, S, ext_air);

wa = multiscatter(waoptions, wavelength, alt, rho_div, rho_fov, ...
                  range, ext, radius, S, ext_air, ssa(1), g(1), ...
                  ssa_air, droplet_frac, pristine_ice_frac);

adjoint = multiscatter(adoptions, wavelength, alt, rho_div, rho_fov, ...
                  range, ext, radius, S, ext_air, ssa(1), g(1), ...
                  ssa_air, droplet_frac, pristine_ice_frac, bscat_AD);

% Single scattering
single_bscat = multiscatter_platt(drange, ext, S, ext_air, 1.0);
  
% Platt approximation with mu=0.5
platt_bscat = multiscatter_platt(drange, ext, S, ext_air, 0.5);
  
% Plot the result

semilogx(single_bscat,wa.range./1000,'k--',...
         'linewidth',3,'color',[1 1 1].*0.6);
hold on
semilogx(platt_bscat, range./1000,'k',...
         'color',[1 1 1].*0.6,'linewidth',3);
semilogx(wa.bscat, wa.range./1000,'linewidth',1);

axis([1.e-9 1.e-2 0 2.5]);
%ylim(fliplr(range([1 end]))./1000);
xlabel('Apparent backscatter (m^{-1} sr^{-1})');
ylabel('Height (km)');
set(gca,'xtick',10.^[-9:-2]);
legend('Single scattering','Platt \eta = 0.5',...
       '{\itr_e} = 5 \mum','{\itr_e} = 10 \mum',...
       4);
title(['Liquid cloud, optical depth = ' num2str(optical_depth)]);

jacobian = multiscatter([options ' -numerical-jacobian -ext-only'], wavelength, alt, rho_div, rho_fov, ...
                        range, ext, radius, S, ext_air, ssa(1), g(1), ...
                        ssa_air, droplet_frac, pristine_ice_frac);
