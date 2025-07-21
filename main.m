close all;
clear;
clc;

% parameter 
M = 15; % transmit antenna number
power = 1;  % power budget
L = 1024;   % frame


C = 181;
a = zeros(M, C);    % steering vector
angle_space_rad = linspace(-pi/2, pi/2, C);
for i = 1:M
    for j = 1:C
        a(i, j) = exp(1j * pi * (i-1) * sin(angle_space_rad(j)));
    end
end

angle_target = [-40, 0, 40];    % target angle
bandwidth = 10;
Pd = zeros(C, 1);   % desired beampattern
angle_space_deg = linspace(-90, 90, C);
for i = 1:length(angle_target)
    angle_target_i = angle_target(i);
    for c = 1:C
        if (abs(angle_target_i-angle_space_deg(c)) <= (bandwidth/2))
            Pd(c) = 1;
        end 
    end 
end

tic
[Rd, alpha1] = radar_cross_correlation(Pd, a, power, M);    % solve SDP model (E.q. 14)
toc
disp(['SDP runtime:', num2str(toc)]);

tic
s = signal_bpsk_ori(M);
s = sqrt(power / M) * s;
S = zeros(M^2, size(s,2)/2);
for i = 1:size(s,2)/2
    temp = s(:,i)*s(:,i)';
    S(:,i) = temp(:);
end
r = Rd(:);
p = fcls(real(r'), real(S));    % solve CLS model (E.q. 17), because for N=2, the sampled point is complex conjugate, so the dimension of p is 1/2*N^M
p_area = zeros(size(s,2)+1, 1);
p_area(1) = 0;
for i = 2:size(s,2)+1
    if (i <= (length(p)+1))
        p_area(i) = p_area(i-1) + p(i-1)/2;
    else
        p_area(i) = p_area(i-1) + p(size(s,2)+2-i)/2;
    end
end
Xd = zeros(M, L);  % reference radar signal Xd
for i = 1:L
    tempi = rand();
    indexi = find(p_area > tempi);
    Xd(:,i) = s(:,indexi(1)-1);
end
toc
disp(['NNLS runtime:', num2str(toc)]);
Rxrd = Xd * Xd' / L;

BP_d = abs(alpha1 * Pd);
BP_radar = zeros(C, 1);
BP_radar_synthesis = zeros(C, 1);
for c = 1:C
    BP_radar(c) = abs(a(:,c)' * Rd * a(:,c));
    BP_radar_synthesis(c) = abs(a(:,c)' * Rxrd * a(:,c));  
end

figure;
plot(angle_space_deg, BP_d, 'LineWidth', 2, 'color', 'k'); hold on;
plot(angle_space_deg, BP_radar, 'LineWidth', 2, 'color', 'b'); hold on;
plot(angle_space_deg, BP_radar_synthesis, 'LineWidth', 2, 'color', 'g'); hold on;
xlabel('\theta');
ylabel('Amplitude');
legend( 'desired beampattern', 'optimal covariance matrix R', 'reference radar signal Xd'); grid on;








