clear all;
clc;
%Thiet lap thong so mo phong
    n = 1e6;        %So mau 
    R = 60;         %Sl phan tu phan xa RISc
    S = 30;         %Sl phan tu phan xa RISd
    M = 15;         %Sl angten DR
    sig_sq = 1;     %Phuong sai nhieu
    phi = 0.5;      %Ty le phan bo cong suat D2D
    Pd = 10;        %Tong cong suat mang D2D
    pc_dB =0:1:30;                                     % SNR tai CT
    pd_dB = phi*Pd/(sig_sq^2);                         % SNR tai DT
    pj_dB = ((1 - phi) * Pd) / ((M - 1) * (sig_sq^2)); % SNR tai DR

%Chuyen SNR tu dB sang thang do tuyen tinh
    pc_linear = 10.^(pc_dB/10);
    pd_linear = 10^(pd_dB/10);
    pj_linear = 10^(pj_dB/10);
    
%Thiet lap mang luu tru Ergodic
    ergodic_rates_analysis = zeros(size(pc_linear));

%Phan tich Ergodic
    for i = 1:length(pc_linear)
        EU1 = pc_linear(i)/(R+1) + pd_linear/S + pj_linear/(M-1);
        VU1 = (pc_linear(i)/(R+1))^2 + (pd_linear/S)^2 + (pj_linear/(M-1))^2;
        EU2 = pd_linear/S + pj_linear/(M-1);
        VU2 = (pd_linear/S)^2 + (pj_linear/(M-1))^2;
        EU3 = pc_linear(i)/(R+1) + pj_linear/(M-1);
        VU3 = (pc_linear(i)/(R+1))^2 + (pj_linear/(M-1))^2;

        % Tinh toc do Ergodic tai E
        ergodic_rates_analysis(i) = (1/log(2)) * (log((1 +EU1)/(1+EU2)) ...
                                    -VU1/(2*(1+EU1)^2) + VU2/(2*(1+EU2)^2));
    end    
%Ve do thi cho cac ket qua    
    figure;
    plot(pc_dB, ergodic_rates_analysis, 'k-.','LineWidth', 1.5);
    xlabel('p_c(in dB))');
    ylabel('R^C_E (in bps/Hz)');
    title('Tỷ lệ Ergodic tại E với SNR tại CT(p_C)');
    grid on;
    ax = gca;
    ax.XMinorGrid = 'on';
    hold on;

pc_dB = [0 3 7 14 18 22 27 30];
pc_linear = 10.^(pc_dB/10);
ergodic_rates_simulated = zeros(size(pc_linear));

%Mo phong Ergodic
    for i = 1:length(pc_linear)
        ergodic_rate_sum = 0;
        for j = 1:n
        % Tao kenh Rayleigh theo phan phoi Gaussian
            h2 = (randn(R, 1) + 1i * randn(R, 1)) / sqrt(2);
            h4 = (randn() + 1i * randn()) / sqrt(2);
            h5 = (randn(R, 1) + 1i * randn(R, 1)) / sqrt(2);
            h6 = (randn(S, 1) + 1i * randn(S, 1)) / sqrt(2);
            h8 = (randn(S, 1) + 1i * randn(S, 1)) / sqrt(2);
            h9  = (randn(M-1,1) + 1i * randn(M-1,1)) / sqrt(2);
            h10 = (randn(M-1,1) + 1i * randn(M-1,1)) / sqrt(2);
            h11 = (randn(M-1,1) + 1i * randn(M-1,1)) / sqrt(2);
            
            Theta_c = diag(exp(1i * rand(R, 1) * 2 * pi));% Ma tran dich pha cua RISc
            Theta_d = diag(exp(1i * rand(S, 1) * 2 * pi));% Ma tran dich pha cua RISd
            Omega = ones(M - 1, M - 1);%Ma tran Beamforming
            gamma_E_C = (pc_linear(i) * (abs(h4 + h5' * Theta_c * h2)).^2)...
                        /((pd_linear * abs(h8' * Theta_d * h6).^2) + (pj_linear *(norm(h9' * Omega))^2) + 1);      
            ergodic_rate_sum = ergodic_rate_sum + log2(1 + gamma_E_C);% Ergodic tai E
        end
        ergodic_rates_simulated(i) = ergodic_rate_sum / n;% Ergodic tai CR
    end    
hold on;
    plot(pc_dB, ergodic_rates_simulated, 'o-', 'LineWidth', 1);
    legend('Phân tích','Mô phỏng');