function [stn_xs, stn_ys, stn_zs, rho_sr] = correctLightTravelTime(stn_xs_raw, ...
    stn_ys_raw, stn_zs_raw, stn_xr, stn_yr, stn_zr, earth_rot_rate, c)
% Raw geometric distance:
rho_sr_raw = sqrt((stn_xs_raw - stn_xr).^2 ...
+ (stn_ys_raw - stn_yr).^2 + (stn_zs_raw - stn_zr).^2);
dOmega = rho_sr_raw * earth_rot_rate / c;
% Correct the coordinates.
stn_xs = stn_xs_raw .* cos(dOmega) + stn_ys_raw .* sin(dOmega);
stn_ys = - stn_xs_raw .* sin(dOmega) + stn_ys_raw .* cos(dOmega);
stn_zs = stn_zs_raw;
% Recompute the geometric distance.
rho_sr = sqrt((stn_xs - stn_xr).^2 ...
+ (stn_ys - stn_yr).^2 + (stn_zs - stn_zr).^2);
end