function [c1_corr_relativ] = correctRelativity(epochs, stn_xs, stn_ys, stn_zs, earth_rot_rate, c)
% Time step between epochs
dt = 30.0;
% Number of satellites
ns = length(stn_xs(1,:));
inert_vel_x(1,:) = (stn_xs(2,:) - stn_xs(1,:)) / dt;
inert_vel_y(1,:) = (stn_ys(2,:) - stn_ys(1,:)) / dt;
inert_vel_z(1,:) = (stn_zs(2,:) - stn_zs(1,:)) / dt;
inert_vel_x(length(epochs),:) = (stn_xs(length(epochs),:) - stn_xs(length(epochs)-1,:)) / dt;
inert_vel_y(length(epochs),:) = (stn_ys(length(epochs),:) - stn_ys(length(epochs)-1,:)) / dt;
inert_vel_z(length(epochs),:) = (stn_zs(length(epochs),:) - stn_zs(length(epochs)-1,:)) / dt;

for i = 2:length(epochs) - 1;
    % uses 2point finite differences
    inert_vel_x(i,:) = (stn_xs(i+1,:) - stn_xs(i-1,:)) / (2*dt);
    inert_vel_y(i,:) = (stn_ys(i+1,:) - stn_ys(i-1,:)) / (2*dt);
    inert_vel_z(i,:) = (stn_zs(i+1,:) - stn_zs(i-1,:)) / (2*dt);   
end
for i = 1:length(epochs);
    % Correct for earth's rotation
    inert_vel_x_corr(i,:) = inert_vel_x(i,:) - earth_rot_rate .* stn_ys(i,:);
    inert_vel_y_corr(i,:) = inert_vel_y(i,:) + earth_rot_rate .* stn_xs(i,:);
    inert_vel_z_corr(i,:) = inert_vel_z(i,:);
end
% Compute relativistic correction according to given formula
for i = 1:ns;
    delta_tk(:,i) = -2 * (stn_xs(:,i) .* inert_vel_x_corr(:,i) ...
	+ stn_ys(:,i) .* inert_vel_y_corr(:,i) ...
	+ stn_zs(:,i) .* inert_vel_z_corr(:,i)) / c^2;
	% relativistic correction in meters
    c1_corr_relativ(:,i) = delta_tk(:,i) .* c;
end
