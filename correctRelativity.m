function [c1_corr_relativ] = correctRelativity(epochs, stn_xs, stn_ys, stn_zs, earth_rot_rate, c)
% Time step between epochs
dt = 30.0;
% Number of satellites
ns = length(stn_xs(1,:));
for i = 1:length(epochs)-1;
    % uses 2point finite differences -> no derivative for last epoch!
    inert_vel_x(i,:) = (stn_xs(i+1,:) - stn_xs(i,:)) / dt;
    inert_vel_y(i,:) = (stn_ys(i+1,:) - stn_ys(i,:)) / dt;
    inert_vel_z(i,:) = (stn_zs(i+1,:) - stn_zs(i,:)) / dt;
    % Correct for earth's rotation
    inert_vel_x_corr(i,:) = inert_vel_x(i,:) - earth_rot_rate .* stn_ys(i,:);
    inert_vel_y_corr(i,:) = inert_vel_y(i,:) + earth_rot_rate .* stn_xs(i,:);
    inert_vel_z_corr(i,:) = inert_vel_z(i,:);
end
% Discard last epoch from sattelite postions to have same matrix dimensions
stn_xs_short = stn_xs(1:end-1,:);
stn_ys_short = stn_ys(1:end-1,:);
stn_zs_short = stn_zs(1:end-1,:);
% Compute relativistic correction according to given formula
for i = 1:ns;
    delta_tk(:,i) = -2 * (stn_xs_short(:,i) .* inert_vel_x_corr(:,i) ...
	+ stn_ys_short(:,i) .* inert_vel_y_corr(:,i) ...
	+ stn_zs_short(:,i) .* inert_vel_z_corr(:,i)) / c^2;
	% relativistic correction in meters
    c1_corr_relativ(:,i) = delta_tk(:,i) .* c;
end
