function [deltap, epsilon, sigmax, sigmay, sigmaz, sigmat] = ...
    calcEpochCoordinates(dPdx, dPdy, dPdz, dPdt, c1_corrected, ...
    rho_sr, epochs)
cs =[1 1 1 1 1 1 1] * dPdt;
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
    Q = N \ eye(size(N)); %changed, so singularity warning doesn't show
    sigmax(:,i) = m0(:,i) * Q(1,1);
    sigmay(:,i) = m0(:,i) * Q(2,2);
    sigmaz(:,i) = m0(:,i) * Q(3,3);
    sigmat(:,i) = m0(:,i) * Q(4,4);
end
end