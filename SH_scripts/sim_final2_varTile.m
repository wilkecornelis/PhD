close all;
clear all;

%% setup variables
lambda = 0.3;
k = 2*pi/lambda;
d = lambda/2;                   % Inter-element spacing in meter

B = 10e6;                   % Integration bandwidth in Hz
tau = 0.02;                % Integration time in s
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

ind =1;
for m = 0:sqrt(P)-1
    for n = 0:sqrt(P)-1
        D_tile(ind,:) = M*d.*[m n];
        ind = ind +1;
    end
end

% plot(D_tile(:,1),D_tile(:,2),'*')
%
% for i = 1:length(D_tile)
%     text(D_tile(i,1),D_tile(i,2),int2str(i));
% end

% Construct element coordinate vector at tile level (MxM tile)
%
% for tile_ind = 1:P
%     ind =1;
%     for m = 0:M-1
%         for n = 0:M-1
%             D_elem(ind,:) = d.*[m n] + D_tile(tile_ind,:);
%             ind = ind +1;
%         end
%     end
%
% end

% matrix of source signals
s_t = (randn(9, Lsignal) + 1i * randn(9, Lsignal));

% construct vector of direction cosines of the sources, given source
% declination measured from broadside
s_decl_range = linspace(1,85,85);

% true (complex) gain values
g_true = ones(P, 1);% + 0.1i * randn(P, 1);
% ensure proper phase referencing
g_true = g_true / (g_true(1) / abs(g_true(1)));
% amplitude normalization to facilitate comparison
g_true = g_true / mean(abs(g_true));

for decl_ind = 1:1
    tic
    decl_ind
    s_decl = s_decl_range(decl_ind);
    s_decl = 89;
    
    phi_range = linspace(0,315,8);
    l_s = zeros(9,2);
    l_s(1,:) = [0 0];
    l_s(2:9,:) = [sin(s_decl*pi/180).*cos(phi_range.'*pi/180) sin(s_decl*pi/180).*sin(phi_range.'*pi/180)];
    
    % Value of EEP at source locations (9 sources)
    EEP = zeros(9,1);
%     EEP(1,1) = 1;
%     EEP(2,1) = phi0(s_decl+91);
%     EEP(3,1) = phi45(s_decl+91);
%     EEP(4,1) = phi90(s_decl+91);
%     EEP(5,1) = phi45(s_decl+91);
%     EEP(6,1) = phi0(91-s_decl);
%     EEP(7,1) = phi45(91-s_decl);
%     EEP(8,1) = phi90(91-s_decl);
%     EEP(9,1) = phi45(91-s_decl);
    
    EEP(1,1) = 1;
    EEP(2,1) = 1;
    EEP(3,1) = 1;
    EEP(4,1) = 1;
    EEP(5,1) = 1;
    EEP(6,1) = 1;
    EEP(7,1) = 1;
    EEP(8,1) = 1;
    EEP(9,1) = 1;
    
    % Calculate total signal for each tile
    s_tile = zeros(P,Lsignal);
    s_tile_i = zeros(P,Lsignal);
    if M>1
        for tile_ind = 1:P
            ind =1;
            s_elem = zeros(M^2,Lsignal);
            s_i = zeros(M^2,Lsignal);
            for m = 0:M-1
                for n = 0:M-1
                    % total antenna signal with calib source
                    EEP(1,1) = 1;
                    D_elem = d.*[m n] + D_tile(tile_ind,:);
                    w_s = (exp(1j*k*l_s*D_elem.'));
                    s_elem(ind,:) = sum(EEP.*w_s.*s_t);
                    
                    % total antenna signal without calib source for calculation of SIR
                    EEP(1,1) = 0;
                    D_elem = d.*[m n] + D_tile(tile_ind,:);
                    w_s = (exp(1j*k*l_s*D_elem.'));
                    s_i(ind,:) = sum(EEP.*w_s.*s_t);
                    
                    ind = ind +1;
                end
            end
            s_tile(tile_ind,:) = sum(s_elem);%;/(sqrt(M^2)); % Not exactly sure yet why need to divide by sqrt(M^2)
            s_tile_i(tile_ind,:) = sum(s_i);
        end
    end
    
    % The all-digital case
    if M==1
        for ind = 1:P
            % total antenna signal with calib source
            EEP(1,1) = 1;
            w_s = (exp(1j*k*l_s*D_tile(ind,:).'));
            s_tile(ind,:) = sum(EEP.*w_s.*s_t);
            
            % total antenna signal without calib source for calculation of SIR
            EEP(1,1) = 0;
            w_s = (exp(1j*k*l_s*D_tile(ind,:).'));
            a_i(ind,:) = (w_s);
            s_tile_i(ind,:) = sum(EEP.*w_s.*s_t);           
        end
    end
    
    % Add gain noise and thermal noise to tile signals
    noise = sqrt(1 / SNR) * (randn(N^2, Lsignal) + 1i * randn(N^2, Lsignal));
    
    % Check term1
%     a_i(:,1)=[];
%     
%     term1 = 0;
%     for i_ind = 1:8
%         term1 = term1 + a_i(2,i_ind)*(a_i(:,i_ind)'*ones(P,1));
%     end
%     abs(2*term1)
    
    s_tile = g_true.*s_tile + noise;
    
    sigma_i = (1/Lsignal)*(s_tile_i*s_tile_i');
    sigma_s = (1/Lsignal)*(s_tile*s_tile');
    
%     term2 = sigma_i*ones(P,1);
%     abs(term2(2))
    
    w_ref = ones(P,1);
    y = w_ref'*s_tile;
    
    %     r1_s = (1/Lsignal)*(s_tile*y');
    r2 = diag((1/Lsignal)*(s_tile*s_tile'));
    
    %% Start calibration
%     g_est = ones(P,1);
%     sigma_cal = 2;
%     g_prev = 0;
%     
%     Niter = 20;
%     for iter = 1:Niter
%         % update reference beam weight vector
%         w_ref = (ones(P, 1)./ conj(g_est));
%         
%         % update crosscorrelation vector
%         y_t = w_ref'*s_tile;
%         r1 = (1/Lsignal)*(s_tile*y_t');
%         
%         % Compile "y-vector"
%         R = [r1;r2];
%         
%         % Construct submatrices of the A-matrix using the previous gain estimate
%         a1 = (g_est'*w_ref)*eye(P)*sigma_cal;
%         a2 = diag(w_ref);
%         a3 = diag((conj(g_est)))*sigma_cal;
%         a4 = eye(P);
%         
%         % compile A-mat
%         A = [a1 a2; a3 a4];
%         A_inv = inv(A);
%         
%         % Solve unknown gain and system noise vector ("x-vector")
%         x_est = A\R;
%         
%         % update gain estimates
%         g_est = (x_est(1:P) + g_prev)/2;
%         
%         % ensure proper phase referencing
%         g_est = g_est / (g_est(1) / abs(g_est(1)));
%         % perform amplitude normalization to facilitate comparison
%         g_est = g_est / mean(abs(g_est));
%         
%         g_prev = g_est;
%         
%         %         g_mag(:,iter) = abs(abs(g_true)-abs(g_est));
%         %         g_phase(:,iter) = abs(angle(g_true)-angle(g_est));
%         g_est_store(:,iter) = g_est;
%         
%     end
%     Rxx(:,decl_ind) = r2;
%     Rxy(:,decl_ind) = r1;
%     Rxy_i(:,decl_ind) = sigma_i*ones(P,1);
%     g_mag(:,decl_ind) = abs(abs(g_true)- abs(g_est));
%     g_phase(:,decl_ind) = abs(angle(g_true)- angle(g_est));
%     toc
end
%
% save('data_varTile\g_mag_M4.mat','g_mag');
% save('data_varTile\g_phase_M4.mat','g_phase');
% save('data_varTile\Rxx_M4.mat','Rxx');
% save('data_varTile\Rxy_M4.mat','Rxy');
% save('data_varTile\Rxy_i_M4.mat','Rxy_i');










