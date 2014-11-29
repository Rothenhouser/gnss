clear all
close all

% Constants
% Speed of light (m/s)
c = 299792458;
% Time step between epochs
dt = 30.0;
% Number of satellites
ns = 7;
% Earth rotation rate (rad/s)
earth_rot_rate = 7.2921151467E-5;

epochs = importdata('data/Epochs.txt');
epochs = epochs(:,1);

%% Station specific data.
% A priori receiver coordinates (m)
wank_xr = 4235956.688;
wank_yr = 834342.467;
wank_zr = 4681540.682;
% Recorded data
wank_c1 = importdata('data/WANK_C1');
wank_satt = importdata('data/WANK_SATT');
wank_xs_raw = importdata('data/WANK_SATX');
wank_ys_raw = importdata('data/WANK_SATY');
wank_zs_raw = importdata('data/WANK_SATZ');

zugs_xr = 4246098.549;
zugs_yr = 824269.097;
zugs_zr = 4675790.018;
% Recorded data
zugs_c1 = importdata('data/ZUGS_C1');
zugs_satt = importdata('data/ZUGS_SATT');
zugs_xs_raw = importdata('data/ZUGS_SATX');
zugs_ys_raw = importdata('data/ZUGS_SATY');
zugs_zs_raw = importdata('data/ZUGS_SATZ');


stn = 'zugs';
if strcmp(stn, 'wank');
    stn_xr = wank_xr;
    stn_yr = wank_yr;
    stn_zr = wank_zr;
    stn_c1 = wank_c1;
    stn_satt = wank_satt;
    stn_xs_raw = wank_xs_raw;
    stn_ys_raw = wank_ys_raw;
    stn_zs_raw = wank_zs_raw;
elseif strcmp(stn, 'zugs');
    stn_xr = zugs_xr;
    stn_yr = zugs_yr;
    stn_zr = zugs_zr;
    stn_c1 = zugs_c1;
    stn_satt = zugs_satt;
    stn_xs_raw = zugs_xs_raw;
    stn_ys_raw = zugs_ys_raw;
    stn_zs_raw = zugs_zs_raw;        
end

%% 3 (a): Correct satellite positions for Earth rotation during light
%% travel time.
[stn_xs, stn_ys, stn_zs, rho_sr] = correctLightTravelTime(stn_xs_raw, ...
    stn_ys_raw, stn_zs_raw, stn_xr, stn_yr, stn_zr, earth_rot_rate, c);

%% 3 (b) Compute the derivates of the observation equation.
% Todo: change this back to c? Need to adjust calculations below!
dPdt = 1;
dPdx = -(stn_xs - stn_xr) ./ rho_sr;
dPdy = -(stn_ys - stn_yr) ./ rho_sr;
dPdz = -(stn_zs - stn_zr) ./ rho_sr;
% Correct pseudorange data with satellite clock correction.
c1_corrected = stn_c1 + (c * stn_satt); % changed sign before bracket

%% Calculate station coordinates for each epoch.
[deltap, epsilon, sigmax, sigmay, sigmaz, sigmat] = ...
    calcEpochCoordinates(dPdx, dPdy, dPdz, dPdt, c1_corrected, ...
    rho_sr, epochs);
stn_gps_coords = [stn_xr + deltap(1,:); stn_yr + deltap(2,:); ...
    stn_zr + deltap(3,:)];
%% For ex 2.2: Plot satellite positions to perhaps explain why y accuracy
%% is better than x?
% Drift may be because troposphere not yet corrected and not constant??

%% B: 2.1 Single set of coordinates.
[deltap, epsilon, sigma] = calcSingleCoordinates(dPdx, ...
    dPdy, dPdz, dPdt, c1_corrected, rho_sr, epochs);
stn_single_gps_coords = [stn_xr + deltap(1,:); stn_yr + deltap(2,:); ...
    stn_zr + deltap(3,:)];

% "Again give coordinate corrections and the formal errors."

%% Calculate coordinates with tropospheric correction.
c1_corr_tropo = correctTroposphere(stn_c1, stn_satt, ...
    stn_single_gps_coords, stn_xs, stn_ys, stn_zs, c);
% Redo the coordinate calculation.
% Call least-squares again as function (?)

%% Part B, 2.3
% All the corrections need to be applied to the pseudorange observations
% before the LS-alogrithm is executed???
% Compute satellite inertial velocities:
delta_epoch = 30.0;
for i = 1:length(epochs)-1;
    % uses 2point finite differences -> no derivative for last epoch!
    inert_vel_x(i,:) = (stn_xs(i+1,:) - stn_xs(i,:)) / dt;
    inert_vel_y(i,:) = (stn_ys(i+1,:) - stn_ys(i,:)) / dt;
    inert_vel_z(i,:) = (stn_zs(i+1,:) - stn_zs(i,:)) / dt;
    % Correct for earth's rotation
    inert_vel_x_corrected(i,:) = inert_vel_x(i,:) - earth_rot_rate .* stn_ys(i,:);
    inert_vel_y_corrected(i,:) = inert_vel_y(i,:) + earth_rot_rate .* stn_xs(i,:);
    inert_vel_z_corrected(i,:) = inert_vel_z(i,:);
end
% Discard last epoch from sattelite postions to have same matrix dimensions
stn_xs_short = stn_xs(1:end-1,:);
stn_ys_short = stn_ys(1:end-1,:);
stn_zs_short = stn_zs(1:end-1,:);
% Compute relativistic correction according to given formula
for i = 1:ns;
    delta_tk(:,i) = -2 * (stn_xs_short(:,i) .* inert_vel_x_corrected(:,i) ...
    + stn_ys_short(:,i) .* inert_vel_y_corrected(:,i) ...
    + stn_zs_short(:,i) .* inert_vel_z_corrected(:,i)) / c^2;
    % relativistic correction in meters
    delta_tk_m(:,i) = delta_tk(:,i) .* c;
end
