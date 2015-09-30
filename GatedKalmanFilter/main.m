% 代码入口：Acc,Gyro,Mag,TimeLine,Angle
%%
% *************************************************************************
% 0)基本数据的设置
% *************************************************************************

clear ; close all; clc;
load data.mat
% Acc=Acc(1:1:1000);
% Gyro=Gyro(1:1:1000);
% Mag=Mag(1:1:1000);
% TimeLine=TimeLine(1:1:1000);
% 引入函数包
wag = wagLibrary2;
% 是否画图
showPlots = 'yes'; 
file_name='yaw';%还不清楚用处
j=1;
% 采样率
f=0;

%% 
% *************************************************************************
% 1)读入校准后的数据（Workspace中需要有Acc,Gyro,Mag,TimeLine,Angle数据）
% *************************************************************************

% -------------------------------------------------------------------------
% 1.3) Initialize Monte Carlo Simulation variables.
% -------------------------------------------------------------------------
% The for loop which iterates through all the signals in the folder is
% actually a Monte Carlo Simulation, since we will apply all the algorithms
% over all the signals, store the results and, at the end of the loop, we
% will average the results. Here we initizalize all the vectors where the
% optimal parameters of each algorithm as well as the resulting RMSE are
% stored. 

% Gated Kalman Filter algorithm: It has four parameters: 'alpha1',
% 'alpha2', 'beta1' and 'beta2':
opt_alphas1_GKF = 0;
opt_alphas2_GKF = 0;
opt_betas1_GKF = 0;
opt_betas2_GKF = 0;
rmse_GKF = 0;

% -------------------------------------------------------------------------
% 1.4) Monte Carlo Simulation.
% -------------------------------------------------------------------------

% ---------------------------------------------------------------------
% 1.3) Check the name of the Euler rotation angle.
% ---------------------------------------------------------------------
% Since the mechanical device to which the MIMU is attached only allows
% the rotation around a single Euler angle, each one of the files will
% contain random rotations around just one Euler angle. The name of the
% .CSV file contains the name of the Euler angle so we can extract it.
if strfind(file_name, 'pitch')
    selectedAngle = 'pitch';
elseif strfind(file_name, 'roll')
    selectedAngle = 'roll';
elseif strfind(file_name, 'yaw')
    selectedAngle = 'yaw';
end

% Build an iterator which starts by j = 1 and increases with the value 
% of the variable 'iter'. This is done to avoid using 'iter - 2' 
% everytime we want to store a value in a vector.
%     j = iter - 2;

% ---------------------------------------------------------------------
% 1.4) Calibrate raw data. 
% ---------------------------------------------------------------------
% Once data are loaded, we proceed to calibrate them, i.e. we transform
% the raw units [0-1023] into meaningful physical units (g, deg/s and 
% Gauss).      

% Acc是校准过的角速度
axC=Acc(:,1);
ayC=Acc(:,2);
azC=Acc(:,3);

% Truncate noise (hardcore style).
axC(abs(axC(:)) < 1e-6) = 0;      
ayC(abs(ayC(:)) < 1e-6) = 0;      
azC(abs(azC(:)) < 1e-6) = 0;  

% Calibrate MAGNETIC FIELD. 
hxC=Mag(:,1);
hyC=Mag(:,2);
hzC=Mag(:,3);                     

% Calibrate ANGULAR RATE. 
gxC=Gyro(:,1);
gyC=Gyro(:,2);
gzC=Gyro(:,3);

% Compensate remaining deviation (hardcore style). This is done since
% in our case the MIMU was static and already stable at the beginning
% of the data gathering process.
gxC = gxC - gxC(1);          
gyC = gyC - gyC(1);         
gzC = gzC - gzC(1);  

% Transform units from deg/s to rad/s.
gxC = gxC * pi/180;          
gyC = gyC * pi/180;         
gzC = -gzC * pi/180;         

% Build the acceleration, angular rate and magnetic field matrices.
Accelerometer = [axC ayC azC]; 
Gyroscope = [gxC gyC -gzC];    
Magnetometer = [hxC' hyC' hzC']; 

% Build the time vector.
time = TimeLine;
f=1/(time(2)-time(1));
    
%% 2) ATTITUDE ESTIMATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ---------------------------------------------------------------------
%  |
%  |_ Once all data are properly compensated and calibrated, we proceed
%     to estimate the orientation of the body which is being monitored. 
%     We will use different methods since the goal is to carry out a 
%     comparative study among them to determine their performance and 
%     accuracy under different motion conditions. 

%% --------------------------------------------------------------------
% 2.1) Projection of gravity and magnetic field. 
% ---------------------------------------------------------------------
% Compute pitch, roll and yaw. Pitch and roll are computed using the
% gravity projections and yaw is computed using the magnetic field
% projections. 
[pitchAcc, rollAcc] = wag.pitch_roll_decomp(axC, ayC, azC, 'rad'); 
[yawAccM, magXh, magYh] = wag.yaw_decomp(pitchAcc, rollAcc, hxC,...
    hyC, hzC);

% Compute RMSE with respect to the reference angle and plot signals. 
% The RMSE is computed from the 300th sample. This is done to allow 
% slower algorithms reach their convergence region and therefore allow
% for a fairer comparison between methods. 
if strcmpi(selectedAngle, 'pitch')
    acc_angle = pitchAcc;
elseif strcmpi(selectedAngle, 'roll')
    acc_angle = rollAcc;
elseif strcmp(selectedAngle, 'yaw')
    acc_angle = yawAccM;
end
    
% Plot estimated orientation angle. 
% if strcmpi(showPlots, 'yes')
%     wag.create_figures('other','Acc+Mag. Projections', time, ...
%         180 / pi * acc_angle, '--y'); 
% end
if strcmpi(showPlots, 'yes')
    figure(1);
    plot(time,acc_angle*180/pi,'b.')
    legend('Acc+Mag')
end
% a=[yawAccM,rollAcc,pitchAcc];
% plot(time,a*180/pi);

%% --------------------------------------------------------------------
% 2.7) Detection of motion intensity (low/high).
% ---------------------------------------------------------------------
% Now we will apply a series of motion intensity detectors to determine
% the degree of intensity of the motion. In this case we will make two
% distinctions: smooth and intense motion. The detectors build a binary
% marker which indicates whether the motion is smoot ('0') or intense
% ('1'). 

% 2.7.1) Initialize parameters of algorithm.

% MBGT (window size and decision threshold).
lwin_mbgt = 14;       threshold_mbgt = 1.4;

% MBCUSUM (window size and decision threshold).
lwin_mbcusum = 18;    threshold_mbcusum = 2.8e-6;   

% FSD (window size, decision threshold, overlapping and normalization
% factor). 
lwin_fsd = 20;  threshold_fsd = 6;  shift_fsd = 19; lambda = 30;

% LTSD (window size, decision threshold and overlapping).
lwin_ltsd = 14;       threshold_ltsd = 4;   shift_ltsd = 2;

% AMD (window size and decision threshold).
lwin_amd = 86;        threshold_amd = 0.0011;

% AMVD (window size and decision threshold).
lwin_amvd = 16;       threshold_amvd = 0.0173;

% ARED (window size and decision threshold).
lwin_ared = 8;        threshold_ared = 1;

% SHOD (window size and decision threshold).
lwin_shod = 19;       threshold_shod = 1;

% Define input signal.
input_signal = sqrt(axC .^ 2 + ayC .^ 2 + azC .^ 2);

% Computation of intensity markers.
[marker_mbgt, T_mbgt, marker_mbcusum, T_mbcusum, marker_fsd, ...
    T_fsd, T_fsd_expanded, marker_ltsd, T_ltsd, T_ltsd_expanded, ...
    marker_amd, T_amd, marker_amvd, T_amvd, marker_ared, T_ared, ...
    marker_shod, T_shod] = wag.build_int_markers(input_signal, axC, ...
    ayC, azC, gxC, gyC, gzC, threshold_mbgt, threshold_mbcusum, ...
    threshold_fsd, threshold_ltsd, lwin_mbgt, lwin_mbcusum, ...
    lwin_fsd, lwin_ltsd, shift_fsd, shift_ltsd, lambda, lwin_amd, ...
    threshold_amd, lwin_amvd, threshold_amvd, lwin_ared, ...
    threshold_ared, lwin_shod, threshold_shod);

% % Plot input signal and markers in order to decide which marker to
% % use. 
% if strcmpi(showPlots, 'yes')
%     detectors_figure = figure(2);
%     subplot(3, 3, 1)
%     plot(T_mbgt)
%     hold on
%     plot(threshold_mbgt * ones(1, length(T_mbgt)), 'r')
%     legend('Detector output (MBGT)', 'Detection threshold')
%     subplot(3, 3, 2)
%     plot(T_mbcusum)
%     hold on
%     plot(threshold_mbcusum * ones(1, length(T_mbcusum)), 'r')
%     legend('Detector output (MBCUSUM)', 'Detection threshold')
%     subplot(3, 3, 3)
%     plot(T_fsd_expanded)
%     hold on
%     plot(threshold_fsd * ones(1, length(T_fsd_expanded)), 'r')
%     legend('Detector output (FSD)', 'Detection threshold')
%     subplot(3, 3, 4)
%     plot(T_ltsd_expanded)
%     hold on
%     plot(threshold_ltsd * ones(1, length(T_ltsd_expanded)), 'r')
%     legend('Detector output (LTSD)', 'Detection threshold')
%     subplot(3, 3, 5)
%     plot(T_amd)
%     hold on
%     plot(threshold_amd * ones(1, length(T_amd)), 'r')
%     legend('Detector output (AMD)', 'Detection threshold')
%     subplot(3, 3, 6)
%     plot(T_amvd)
%     hold on
%     plot(threshold_amvd * ones(1, length(T_amvd)), 'r')
%     legend('Detector output (AMVD)', 'Detection threshold')
%     subplot(3, 3, 7)
%     plot(T_ared)
%     hold on
%     plot(threshold_ared * ones(1, length(T_ared)), 'r')
%     legend('Detector output (ARED)', 'Detection threshold')
%     subplot(3, 3, 8)
%     plot(T_shod)
%     hold on
%     plot(threshold_shod * ones(1, length(T_shod)), 'r')
%     legend('Detector output (SHOD)', 'Detection threshold')
% 
%     markers_figure = figure(3);
%     subplot(3, 3, 1)
%     plot(input_signal)
%     hold on
%     plot(marker_mbgt + 1,'r')
%     legend('Input signal','MBGT decision')
%     subplot(3, 3, 2)
%     plot(input_signal)
%     hold on
%     plot(marker_mbcusum + 1,'r')
%     legend('Input signal','MBCUSUM decision')
%     subplot(3, 3, 3)
%     plot(input_signal)
%     hold on
%     plot(marker_fsd + 1, 'r')
%     legend('Input signal','FSD decision')
%     subplot(3, 3, 4)
%     plot(input_signal)
%     hold on
%     plot(marker_ltsd + 1, 'r')
%     legend('Input signal','LTSD decision')
%     subplot(3, 3, 5)
%     plot(input_signal)
%     hold on
%     plot(marker_ared + 1, 'r')
%     legend('Input signal','ARED decision')
%     subplot(3, 3, 6)
%     plot(input_signal)
%     hold on
%     plot(marker_amvd + 1,'r')
%     legend('Input signal','AMVD decision')
%     subplot(3, 3, 7)
%     plot(input_signal)
%     hold on    
%     plot(marker_amd + 1, 'r')
%     legend('Input signal','AMD decision')
%     subplot(3, 3, 8)
%     plot(input_signal)
%     hold on       
%     plot(marker_shod + 1,'r')
%     legend('Input signal','SHOD decision')
% end

% Select motion intensity marker.
chosen_marker = marker_fsd;

%% --------------------------------------------------------------------
% 2.9) Gated Kalman Filter (DKF) I.
% ---------------------------------------------------------------------
% After building the motion intensity markers and selecting the most
% appropriate one, we will start to test the gating strategy on the
% same sensor fusion algorithms that were tested above. The goal is to
% check whether the gating strategy allows for an increment in the
% precision of the orientation estimates. 
% The first algorithm that we will check is, again, the standard Kalman
% Filter.

% 2.9.1) Definition of the variables of the optimization process.  
if strcmp(selectedAngle, 'pitch')
    gyro_GKF = gyC;
    obs_GKF = pitchAcc;
elseif strcmp(selectedAngle, 'roll')
    gyro_GKF = gxC;
    obs_GKF = rollAcc;
elseif strcmp(selectedAngle, 'yaw')
    gyro_GKF = gzC;
    obs_GKF = yawAccM;
end

% 2.9.2) Parameter optimization.

% Definition of the initial value of the parameters. 
p0_GKF = [100 100 0.1 0.1];


fprintf('--------------------GKF OPTIMIZATION----------------------\n')
fprintf('----------------------------------------------------------\n')

% Get the optimal parameters. 
opt_alpha1_GKF = 10;
opt_alpha2_GKF = 0;
opt_beta1_GKF = 1;
opt_beta2_GKF = 1;

% Save optimal parameters.
opt_alphas1_GKF(j) = opt_alpha1_GKF;
opt_alphas2_GKF(j) = opt_alpha2_GKF;
opt_betas1_GKF(j) = opt_beta1_GKF;
opt_betas2_GKF(j) = opt_beta2_GKF;

% 2.9.3) Computation of the orientation estimation using optimal 
%        parameters.
angle_GKF = wag.fusionGKF(gyro_GKF, obs_GKF, f, opt_alpha1_GKF, ...
    opt_alpha2_GKF, opt_beta1_GKF, opt_beta2_GKF, chosen_marker);

% % figure
% if strcmpi(showPlots, 'yes')
%     figure(1);
%     hold on
%     wag.create_figures('other', 'GKF', time, 180 / pi * ...
%         angle_GKF, '--k');
%     hold off
% end
if strcmpi(showPlots, 'yes')
    figure(1);
    hold on;
    plot(time,angle_GKF*180/pi,'r-')
    legend('GKF');
    hold off;
end

%% ------------------------------------------------------------------------
% 3) Compute mean RMSE and mean Parameters.
% -------------------------------------------------------------------------
% At the end of the Monte Carlo simulation we need to compute the mean and
% standard deviation of all the parameters and RMSEs.

mean_opt_alphas1_GKF = mean(opt_alphas1_GKF);
std_opt_alphas1_GKF = std(opt_alphas1_GKF);
mean_opt_alphas2_GKF = mean(opt_alphas2_GKF);
std_opt_alphas2_GKF = std(opt_alphas2_GKF);
mean_opt_betas1_GKF = mean(opt_betas1_GKF);
std_opt_betas1_GKF = std(opt_betas1_GKF);
mean_opt_betas2_GKF = mean(opt_betas2_GKF);
std_opt_betas2_GKF = std(opt_betas2_GKF);
mean_rmse_GKF = mean(rmse_GKF);
std_rmse_GKF = std(rmse_GKF);