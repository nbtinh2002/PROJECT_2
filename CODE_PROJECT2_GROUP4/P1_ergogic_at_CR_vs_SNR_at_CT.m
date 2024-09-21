clear all; 
clc;
%Thiet lap thong so mo phong
    n = 1e6;        %So mau 
    R = 60;         %Sl phan tu phan xa RISc 
    pc_dB = 0:1:30; %SNR tai CT

% Chuyen SNR tu dB sang thang do tuyen tinh
    pc_linear = 10.^(pc_dB/10);

% Thiet lap mang luu tru Ergodic
    ergodic_rates_analysis = zeros(size(pc_linear));
    
%Phan tich Ergodic
    for i = 1:length(pc_linear)
        Ec  = 1 + (R * (16 - pi^2)) / 16 + (R * pi / 4)^2 + (R * pi * sqrt(pi)) / 4;
        Ec2 = (R * pi / 4)^4 + 6 * (R * pi / 4)^2 * (R * (16 - pi^2) / 16) + 3 * (R * (16 - pi^2) / 16)^2 ...
            + 2 * sqrt(pi) * ((R * pi / 4)^3 + 3 * (R * pi / 4) * (R * (16 - pi^2) / 16)) ...
            + 6 * (R * (16 - pi^2) / 16 + (R * pi / 4)^2) + 3 * R * pi * sqrt(pi) + 2;
        kc = (Ec^2) / (Ec2 - Ec^2);
        thetac = (Ec2 - Ec^2) / Ec;

        % Tinh toc do Ergodic tai CR
        ergodic_rates_analysis(i) = (1/log(2)) * (log(1 + kc * pc_linear(i) * thetac)...
                                  - (kc * (pc_linear(i) * thetac).^2) / (2 * (1 + kc * pc_linear(i) * thetac).^2));
    end
    % Ve do thi cho cac ket qua
    figure;
    plot(pc_dB, ergodic_rates_analysis, 'k-.','LineWidth', 1.5);
    xlabel('p_c(in dB))');
    ylabel('R^C_L (in bps/Hz)');
    title('Tỷ lệ Ergodic tại CR với SNR tại CT(p_C)');
    grid on;
    ax = gca;
    ax.XMinorGrid = 'on';
    
    pc_dB = [0 3 7 14 18 22 27 30];
    pc_linear = 10.^(pc_dB/10);
    ergodic_rates_simulated = zeros(size(pc_linear));

%Mo phong Ergodic
    for i = 1:length(pc_linear)
        ergodic_rate_sum = 0;
        for j = 1:n
        % Tao kenh Rayleigh theo phan phoi Gaussian
            h1 = (randn() + 1i * randn()) / sqrt(2);
            h2 = (randn(R, 1) + 1i * randn(R, 1)) / sqrt(2);
            h3 = (randn(R, 1) + 1i * randn(R, 1)) / sqrt(2);

            Theta_c = diag(exp(1i * rand(R,1) * 2 * pi));% Ma tran dich pha cua RISc
            gamma_L_C = pc_linear(i) * abs(h1 + h3' * Theta_c * h2).^2;% SINR tai CT
            ergodic_rate_sum = ergodic_rate_sum + log2(1 + gamma_L_C);% Ergodic tai CR
        end
        ergodic_rates_simulated(i) = ergodic_rate_sum / n;% Ergodic tai CR
    end
    hold on;
    plot(pc_dB, ergodic_rates_simulated, 'o-', 'LineWidth', 1);
    legend('Phân tích','Mô phỏng');