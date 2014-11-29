function [tropo_corr] = correctTroposphere(c1, satt, single_gps_coords, ...
                                          xs, ys, zs, c)
% Calculate cos(z), the tropospheric zenith angle for each epoch.
% Need a single set of station coordinates
xr_gps = single_gps_coords(1,:);
yr_gps = single_gps_coords(2,:);
zr_gps = single_gps_coords(3,:);
% OR use the a priori coordinates
% xr_gps = xr
xyz = sqrt((xr_gps).^2 + (yr_gps).^2 + (zr_gps).^2);
rho_sr = sqrt((xs - xr_gps).^2 ...
+ (ys - yr_gps).^2 + (zs - zr_gps).^2);
cos_z = (xr_gps./xyz) .* (xs - xr_gps)./rho_sr ...
    + (yr_gps./xyz) .* (ys - yr_gps)./rho_sr ...
    + (zr_gps./xyz) .* (zs - zr_gps)./rho_sr;
% Tropospheric correction:
tropo_corr = 2.3 ./ cos_z;
end