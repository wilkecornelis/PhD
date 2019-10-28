close all;
clear all;

%% setup variables
lambda = 0.3;
k = 2*pi/lambda;
d = lambda/2;                   % Inter-element spacing in meter

B = 10e6;                   % Integration bandwidth in Hz
tau = 0.01;                % Integration time in s
Lsignal = B * tau;          % Number of time samples
kB = 1.38e-23;              % Boltzmann constant
Tsys = 35;                  % System temperature in kelvin

%% Import EEP from FEKO model

[labels1,a1,phi0] = readColData('data\phi0.txt',2,1,1); %Antenna 1 E-field at 1GHz
[labels1,a1,phi45] = readColData('data\phi45.txt',2,1,1); %Antenna 1 E-field at 1GHz
[labels1,a1,phi90] = readColData('data\phi90.txt',2,1,1); %Antenna 1 E-field at 1GHz


%% Tile signals

% Construct tile coordinate vector (this is the coordinate of the reference of each tile)
M = 1;          % Number of elements on a side of a tile
P = 576/(M^2);  % Number of tiles (or receive paths)

snr_range = linspace(0.01,10,1000)./1000;
for snr_ind = 30:40
    tic
    % matrix of source signals
    s_t = (randn(1, Lsignal) + 1i * randn(1, Lsignal));
    
    % construct vector of direction cosines of the sources, given source
    % declination measured from broadside
    s_decl_range = linspace(1,85,85);
    
    % true (complex) gain values
    g_true = ones(P, 1) + 0.1i * randn(P, 1);
    % ensure proper phase referencing
    g_true = g_true / (g_true(1) / abs(g_true(1)));
    % amplitude normalization to facilitate comparison
    g_true = g_true / mean(abs(g_true));
    
    SNR = snr_range(snr_ind);
    % Add gain noise and thermal noise to tile signals
    noise = sqrt(1 / SNR) * (randn(P, Lsignal) + 1i * randn(P, Lsignal));
    sigma_n = (1/Lsignal)*(noise*noise');
    n_true = diag(sigma_n);
    
    s_tile = g_true.*s_t + noise;
    sigma_s = ((1/Lsignal)*(s_t*s_t'))*ones(P,P);
    
    w_ref = ones(P,1);
    y = w_ref'*s_tile;
    
    %     r1_s = (1/Lsignal)*(s_tile*y');
    r2 = diag((1/Lsignal)*(s_tile*s_tile'));
    
    %% Start calibration
    for mc_ind = 1:100
        g_est = ones(P,1);
        sigma_cal = 2;
        g_prev = 0;

        Niter = 10;
        for iter = 1:Niter
            iter;
            % update reference beam weight vector
            w_ref = (ones(P, 1)./ conj(g_est));

            % update crosscorrelation vector
            y_t = w_ref'*s_tile;
            r1 = (1/Lsignal)*(s_tile*y_t');

            % Compile "y-vector"
            R = [r1;r2];

            % Construct submatrices of the A-matrix using the previous gain estimate
            a1 = (g_est'*w_ref)*eye(P)*sigma_cal;
            a2 = diag(w_ref);
            a3 = diag((conj(g_est)))*sigma_cal;
            a4 = eye(P);

            % compile A-mat
            A = [a1 a2; a3 a4];
            A_inv = inv(A);

            % Solve unknown gain and system noise vector ("x-vector")
            x_est = A\R;

            % update gain estimates
            g_est = (x_est(1:P));% + g_prev)/2;

            % ensure proper phase referencing
            g_est = g_est / (g_est(1) / abs(g_est(1)));
            % perform amplitude normalization to facilitate comparison
            g_est = g_est / mean(abs(g_est));

            g_prev = g_est;

            %             g_est_store(1,iter) = abs(g_est(2));
        end
        %     plot(g_est_store)
        %     hold on;

        g_mag(mc_ind,snr_ind) = mean(abs(abs(g_true)- abs(g_est)));
        g_phase(mc_ind,snr_ind) = mean(abs(angle(g_true)- angle(g_est)));
       
    end    
    %         n_est(:,snr_ind) = x_est(P+1:2*P);
    
    %         term1 = w_ref.*n_true;
    %         term2 = diag(g_tilde*sigma_s*g_tilde');
    %         term3 = g_tilde*sigma_s*g_tilde'*w_ref;
    %         term4 = sigma_n*w_ref;
    
    %     mean(abs(term1));
    %     mean(abs(term2));
    %     mean(abs(term3));
    %     mean(abs(term4));
    
    %         term1_store(mc_ind,snr_ind) = mean(term1);
    %         term2_store(mc_ind,snr_ind) = mean(term2);
    %         term3_store(mc_ind,snr_ind) = mean(term3);
    %         term4_store(mc_ind,snr_ind) = mean(term4);
    toc 
end


%
% save('SNR_data\term1.mat','term1_store');
% save('SNR_data\term2.mat','term2_store');
% save('SNR_data\term3.mat','term3_store');
% save('SNR_data\term4.mat','term4_store');

% plot(g_est_store)







