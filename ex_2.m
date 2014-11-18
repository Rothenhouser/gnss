% Constants
% Speed of light (m/s)
c = 299792458;
% Earth rotation rate (rad/s)
earth_rot_rate = 7.2921151467E-5;
% A priori receiver coordinates (m)
wank_xr = 4235956.688;
wank_yr = 834342.467;
wank_zr = 4681540.682;
% Recorded data
wank_c1 = importdata('WANK_C1');
wank_satt = importdata('WANK_SATT');
wank_xs_raw = importdata('WANK_SATX');
wank_ys_raw = importdata('WANK_SATY');
wank_zs_raw = importdata('WANK_SATZ');
epochs = importdata('Epochs.txt');
epochs = epochs(:,1);
%% 3 (a): Correct satellite positions for Earth rotation during light
%% travel time.
% Raw geometric distance:
rho_sr_raw = sqrt((wank_xs_raw - wank_xr).^2 ...
+ (wank_ys_raw - wank_yr).^2 + (wank_zs_raw - wank_zr).^2);
dOmega = rho_sr_raw * earth_rot_rate / c;
% Correct the coordinates.
wank_xs = wank_xs_raw .* cos(dOmega) + wank_ys_raw .* sin(dOmega);
wank_ys = - wank_xs_raw .* sin(dOmega) + wank_ys_raw .* cos(dOmega);
wank_zs = wank_zs_raw;
% Recompute the geometric distance.
rho_sr = sqrt((wank_xs - wank_xr).^2 ...
+ (wank_ys - wank_yr).^2 + (wank_zs - wank_zr).^2);
%% 3 (b) Compute the derivates of the observation equation.
% Todo: change this back to c? Need to adjust calculations below!
dPdt = 1;
cs =[1 1 1 1 1 1 1] * dPdt;
dPdx = -(wank_xs - wank_xr) ./ rho_sr;
dPdy = -(wank_ys - wank_yr) ./ rho_sr;
dPdz = -(wank_zs - wank_zr) ./ rho_sr;
% Correct pseudorange data with satellite clock correction.
c1_corrected = wank_c1 + (c * wank_satt); % changed sign before bracket
%% Solve the normal equations.
for i=1:length(epochs);
    % Setup the grand design matrix A for every epoch.
    A(:,1)=dPdx(i,:);
    A(:,2)=dPdy(i,:);
    A(:,3)=dPdz(i,:);
    A(:,4)=cs;
    % Solve normal equation.
    deltay(i,:) = c1_corrected(i,:) - rho_sr(i,:); %#ok<*SAGROW>
    atransp = transpose(A);
    N = atransp*A;
    % deltap are the deviations of computed from a priori coordinates. 
    deltap(:,i) = (N \ atransp) * transpose(deltay(i,:));
    % Compute the residuals.
    epsilon(:,i) = transpose(deltay(i,:)) - (A * deltap(:,i));
    m0(:,i) = sqrt(transpose(epsilon(:,i)) * epsilon(:,i) ./ 3);
    % Formal errors
    Q = inv(N);
    sigmax(:,i) = m0(:,i) * Q(1,1);
    sigmay(:,i) = m0(:,i) * Q(2,2);
    sigmaz(:,i) = m0(:,i) * Q(3,3);
    sigmat(:,i) = m0(:,i) * Q(4,4);
end
%% Calculate station coordinates for each epoch.
% Add or subtract the difference???
wank_gps_coords = [wank_xr + deltap(1,:); wank_yr + deltap(2,:); ... 
    wank_zr + deltap(3,:)]; 
%% For ex 2.2: Plot satellite positions to perhaps explain why y accuracy 
%% is better than x?
% Drift may be because troposphere not yet corrected and not constant??

%% B: 2.1 Single set of coordinates.
% Simply average coordinates?
wank_single_gps_coords = [mean(wank_gps_coords(1,:)); ...
    mean(wank_gps_coords(2,:)); mean(wank_gps_coords(3,:))];
% "Again give coordinate corrections and the formal errors."

%% Calculate coordinates with tropospheric correction.
% Calculate cos(z), the tropospheric zenith angle for each epoch.
% Need a single set of station coordinates
wank_xr_gps = wank_single_gps_coords(1,:);
wank_yr_gps = wank_single_gps_coords(2,:);
wank_zr_gps = wank_single_gps_coords(3,:);
% OR use the a priori coordinates
% wank_xr_gps = wank_xr
xyz = sqrt((wank_xr_gps).^2 + (wank_yr_gps).^2 + (wank_zr_gps).^2);
rho_sr = sqrt((wank_xs - wank_xr_gps).^2 ...
+ (wank_ys - wank_yr_gps).^2 + (wank_zs - wank_zr_gps).^2);
cos_z = (wank_xr_gps./xyz) .* (wank_xs - wank_xr_gps)./rho_sr ...
    + (wank_yr_gps./xyz) .* (wank_ys - wank_yr_gps)./rho_sr ...
    + (wank_zr_gps./xyz) .* (wank_zs - wank_zr_gps)./rho_sr;
% Tropospheric correction:
tropo_corr = 2.3 ./ cos_z;
% Apply correction.
c1_corrected = wank_c1 + (c * wank_satt) + tropo_corr; % Or subtract???
% Redo the coordinate calculation.
% Call least-squares again as function (?)


