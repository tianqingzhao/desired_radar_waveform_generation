function [R, alpha1] = radar_cross_correlation(Pd, a, power, M)
    % radar beamforming, design radar waveform convariance matrix
    % paper: "Transmit beamforming for MIMO radar systems using signal cross-correlation"
    m = size(a, 2);
    cvx_begin quiet
    variable R(M, M) hermitian semidefinite
    variable alpha1 nonnegative
    expression u(m)
    for i = 1:m
        u(i) = alpha1*Pd(i) - a(:,i)'*R*a(:,i);
    end
    minimize norm(u, 2)
    subject to
%         trace(R) == power;  % sum power constrant
        diag(R) == (power/M) * ones(M, 1);  % constant of every antenna power is equal
    cvx_end
end