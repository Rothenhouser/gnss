function [deltap, epsilon, sigma] = calcSingleCoordinates(dPdx, ...
    dPdy, dPdz, dPdt, c1_corrected, rho_sr, epochs)

ns = length(rho_sr(1,:));
% Setup the grand design matrix A -> estimate only one set of coordinates
A = zeros(length(epochs) * ns, length(epochs) + 3);
it = 1;
for i = 1:ns:length(epochs)*ns;
    for j = 0:ns-1;
        A(i + j,1) = dPdx(it,j + 1);
        A(i + j,2) = dPdy(it,j + 1);
        A(i + j,3) = dPdz(it,j + 1);
    end
    A(i:i + 6, it + 3) = dPdt;
    it = it + 1; 
end
% setup observation matrix
c1_corrected_resh = reshape(c1_corrected,[],1);
rho_sr_resh = reshape(rho_sr,[],1);
deltay = c1_corrected_resh - rho_sr_resh;

% Solve normal equation.
atransp = transpose(A);
N = atransp*A;
% deltap are the deviations of computed from a priori coordinates. 
deltap = (N \ atransp) * deltay;
% Compute the residuals.
epsilon = deltay - (A * deltap);
m0 = sqrt(transpose(epsilon) * epsilon ./ 3);


% Formal errors
Q = N \ eye(size(N));  %changed, so singularity warning doesn't show 
for i = 1:size(N);
    sigma(i) = m0 * sqrt(Q(i,i));
end;
sigma = transpose(sigma);

end