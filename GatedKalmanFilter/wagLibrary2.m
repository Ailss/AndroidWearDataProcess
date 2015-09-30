function wag = wagLibrary

% FUNCTION WAGLIBRARY contains the references to all necesarry functions to 
% load, extract, calibrate, process and plot data gathered using Wagyromag.
% Additionally, it contains all the necessary functions for the study of
% gating in orientation estimation algorithms.

% Functions to estimate orientation using the acceleration and magnetic
% field decomposition.
wag.pitch_roll_decomp = @pitch_roll_decomp;
wag.yaw_decomp = @yaw_decomp;

% Gated Kalman filter functions.
wag.optimizeGKF = @optimizeGKF;
wag.fusionGKF = @fusionGKF;

% Auxiliary functions.
wag.create_figures = @ create_figures;
wag.build_int_markers = @build_int_markers;

end
% END OF WAGLIBRARY FUNCTION

function create_figures(order, method_name, x_axis, y_axis, format)

% FUNCTION CREATE_FIGURES Plots the orientation estimates. The function
% updates the legend automatically when a new signal is added to the plot.
% There is only a distinction between the first time the function is called
% and the subsequent executions. 
% 
% - INPUT PARAMETERS:
%   |_ 'order': Order in which the function is called. 
%               Possible values:
%                   - 'first': To indicate that the function is called for
%                   the first time.
%                   - 'other': Any other time but the first.
%   |_ 'method_name': Name of the orientation estimation algorithm which
%                     will be displayed in the legend.
%   |_ 'x_axis': Signal to be displayed in the X axis. Usually the time
%                signal.
%   |_ 'y_axis': Signal to be displayed in the Y axis. Usually the angle
%                signal.
%
% - OUTPUT PARAMETERS: This function has no output parameters.

% If the function is called for the first time, plot the signals normally.
if strcmp(order, 'first')
    plot(x_axis, y_axis, format, 'LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    legend(method_name);
end

% If the function has already been called at least once, plot the signals,
% get the names in the legend and update them. 
if strcmp(order, 'other')
    p2 = plot(x_axis, y_axis, format, 'LineWidth', 2);
    hold on
    ch = get(gca, 'children');
    ls = get(legend, 'String');
    if isempty(ls)
        legend(p2, method_name);
        xlabel('Time(s)')
        ylabel('Angle (deg)')
    else
        legend([p2 ch(1)], {char(ls), method_name});
        legend('-DynamicLegend');
        xlabel('Time(s)')
        ylabel('Angle (deg)')
    end 
end
end
% END OF CREATE_FIGURES FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function [pitch, roll] = pitch_roll_decomp(ax, ay, az, units)

% FUNCTION PITCH_ROLL_ACC Computes pitch and roll angles using the
% decomposition of the acceleration measured with a triaxial accelerometer. 
%
% - INPUT PARAMETERS:
%    |_'ax': vector containing the linear accelerations along the cartesian
%          X-axis.
%    |_'ay': vector containing the linear accelerations along the cartesian
%          Y-axis.
%    |_'az': vector containing the linear accelerations along the cartesian
%          Z-axis.
%    |_'units': string containing the units of angles to return. Set 'deg' 
%               or 'rad'.
%
% - OUTPUT PARAMETERS:
%    |_'pitch': vector containing the pitch calculated angle values given 
%               the 'ay' and 'az' components of linear accelerations. Pitch
%               is assumed to be the rotation angle about the Y-axis.
%    |_'roll':  vector containing the roll calculated angle values given 
%               the 'ax' and 'az' components of linear accelerations. Roll
%               is asumed to be the rotation angle about the X-axis.
%

% Check input parameters.
if (length(ax) ~= length(ay) || length(ay) ~= length(az) || ...
    length(az) ~= length(ax))
    error('Input vectors ax, ay and az must be the same length.');
end

if ~(strcmp(units, 'deg') || strcmp(units, 'rad'))
    error('Specified units are not correct. Set ''deg'' or ''rad''.');
end

% Compute pitch and roll. Pitch and roll values are obtained (in rad). Both
% pitch and roll values have different expressions depending on the 
% quadrant they are placed.
pitch = atan2(sqrt(ay .^ 2 + az .^ 2), -ax);
roll = atan2(az, ay);

% Roll and pitch values are transformed in degrees if so specified.
if strcmp(units, 'deg')
    roll = (180 / pi) .* roll;
    pitch = (180 / pi) .* pitch;
end

end
% END OF PITCH_ROLL_ACC FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function [yaw, magXh, magYh] = yaw_decomp(pitch, roll, magX, magY, magZ)

% FUNCTION YAW_DECOMP computes the yaw angle using pitch and roll angles 
% and the magnetic field gathered with a triaxial magnetometer. 
%
% - INPUT PARAMETERS: 
%    |_ 'pitch': vector containing the a priori known pitch angle values.
%    |_ 'roll': vector containing the a priori known roll angle values.
%    |_ 'magX': vector containing the a priori known magnetic field 
%               component along X-axis.
%    |_ 'magY': vector containing the a priori known magnetic field 
%               component along Y-axis.
%    |_ 'magZ': vector containing the a priori known magnetic field 
%               component along Z-axis.
%    |_ 'plotGraphics': string containing 'yes' or 'no' to specify if 
%                       results are to be plotted.
%
% - OUTPUT PARAMETERS:
%    |_ 'yaw': vector containing the calculated yaw angles, considering 
%              both local magnetic field and pitch and roll values.
%    |_ 'magXh': vector containing the projections of magX over the XY 
%                cartesian plane. These values are used to compute YAW 
%                output values.
%    |_ 'magYh': vector containing the projections of magY over the XY 
%                cartesian plane. These values are used to compute YAW 
%                output values.
%

% Check input parameters.
len = length(pitch);
if (length(roll) ~= len || length(magX) ~= len || length(magY) ~= ...
        len || length(magZ) ~= len)
    error('All input arguments must be the same length.');
end

% Compute Yaw Angle. Calculation of horizontal magnetic X coordinate 
% component (magXh) and horizontal magnetic Y coordinate component (magYh).
% This is done by derotating the measured magnetic field in axes X and Y by
% the pitch and roll angles of the body. That is, we find the projection of
% the magnetic field in axes X and Y in the XY plane. 
magXh = magX .* cos(pitch) + magY .* sin(roll) .* sin(pitch) - ...
        magZ .* cos(roll) .* sin(pitch);
magYh = magY .* cos(roll) + magZ .* sin(roll);
yaw = atan((-magYh) ./ magXh);

% Quadrant compensations. This part may need to be changed depending on the
% body reference frame of the MIMU. 
for i = 1 : length(yaw)
    if magXh(i) > 0 && magYh(i) == 0
        yaw(i) = pi / 2;
    end
    if magXh(i) < 0 && magYh(i) == 0
        yaw(i) = -pi / 2;
    end
    if magXh(i) < 0 && magYh(i) > 0
        yaw(i) = yaw(i) - pi;
    end
    if magXh(i) < 0 && magYh(i) < 0
        yaw(i) = yaw(i) + pi;
    end
    
    while yaw(i) > pi 
        yaw(i) = yaw(i) - 2 * pi; 
    end
    while yaw(i) < -pi 
        yaw(i) = yaw(i) + 2 * pi; 
    end
end

end
% END OF YAW_DECOMP FUNCTION
 
% % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [marker_mbgt, T_mbgt, marker_mbcusum, T_mbcusum, marker_fsd, ...
    T_fsd, T_fsd_expanded, marker_ltsd, T_ltsd, T_ltsd_expanded, ...
    marker_amd, T_amd, marker_amvd, T_amvd, marker_ared, T_ared, ...
    marker_shod, T_shod] = build_int_markers(input_signal, ax, ay, az, ...
    gx, gy, gz, threshold_mbgt, threshold_mbcusum, threshold_fsd, ...
    threshold_ltsd, lwin_mbgt, lwin_mbcusum, lwin_fsd, lwin_ltsd, ...
    shift_fsd, shift_ltsd, lambda, lwin_amd, threshold_amd, lwin_amvd, ...
    threshold_amvd, lwin_ared, threshold_ared, lwin_shod, threshold_shod)

% FUNCTION BUILD_INT_MARKERS computes a series of motion intensity markers
% applying the MBGT, MBCUSUM, FSD, LTSD, AMD, AMVD, ARED and SHOD
% algorithms. For more information about these algorithms check the
% following publication: 'Detection of (In)activity Periods in Human Body 
% Motion Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%
% - INPUT PARAMETERS:
%   |_ 'input_signal': Signal acting as the input for the MBGT, MBCUSUM,
%                      FSD and LTSD algorithms.
%   |_ 'ax': Calibrated acceleration measured along the X-axis.
%   |_ 'ay': Calibrated acceleration measured along the Y-axis.
%   |_ 'az': Calibrated acceleration measured along the Z-axis.
%   |_ 'gx': Calibrated angular rate measured along the X-axis.
%   |_ 'gy': Calibrated angular rate measured along the Y-axis.
%   |_ 'gz': Calibrated angular rate measured along the Z-axis.
%   |_ 'threshold_mbgt': Decision threshold of the MBGT algorithm.
%   |_ 'threshold_mbcusum': Decision threshold of the MBCUSUM algorithm.
%   |_ 'threshold_fsd': Decision threshold of the FSD algorithm.
%   |_ 'threshold_ltsd': Decision threshold of the LTSD algorithm.
%   |_ 'lwin_mbgt': Size of the sliding window of the MBGT algorithm.
%   |_ 'lwin_mbcusum': Size of the sliding window of the MBCUSUM algorithm.
%   |_ 'lwin_fsd': Size of the sliding window of the FSD algorithm.
%   |_ 'lwin_ltsd': Size of the sliding window of the LTSD algorithm.
%   |_ 'shif_fsd': Overlapping of the sliding window of the FSD algorithm.
%   |_ 'shif_ltsd': Overlapping of the sliding window of the LTSD 
%                   algorithm.
%   |_ 'lambda': Normalization factor of the MBCUSUM algorithm.
%   |_ 'lwin_amd': Size of the sliding window of the AMD algorithm.
%   |_ 'threshold_amd': Decision threshold of the AMD algorithm.
%   |_ 'lwin_amvd': Size of the sliding window of the AMVD algorithm.
%   |_ 'threshold_amvd': Decision threshold of the AMVD algorithm.
%   |_ 'lwin_ared': Size of the sliding window of the ARED algorithm.
%   |_ 'threshold_ared': Decision threshold of the ARED algorithm.
%   |_ 'lwin_shod': Size of the sliding window of the SHOD algorithm.
%   |_ 'threshold_shod': Decision threshold of the SHOD algorithm.

% - OUTPUT PARAMETERS:
%   |_ 'marker_mbgt': Motion intensity marker of the MBGT algorithm.
%   |_ 'T_mbgt': Decision signal of the MBGT algorithm to which the
%                threshold is compared.
%   |_ 'marker_mbcusum': Motion intensity marker of the MBCUSUM algorithm.
%   |_ 'T_mbcusum': Decision signal of the MBCUSUM algorithm to which the
%                   threshold is compared.
%   |_ 'marker_fsd': Motion intensity marker of the FSD algorithm.
%   |_ 'T_fsd': Decision signal of the FSD algorithm to which the threshold
%               is compared.
%   |_ 'T_fsd_expanded': Decision signal of the FSD algorithm to which the
%                        threshold is compared. This signal has the same 
%                        length of the input signal. 
%   |_ 'marker_ltsd': Motion intensity marker of the LTSD algorithm.
%   |_ 'T_ltsd': Decision signal of the LTSD algorithm to which the 
%                threshold is compared.
%   |_ 'T_ltsd_expanded': Decision signal of the LTSD algorithm to which 
%                         the threshold is compared. This signal has the 
%                         same length of the input signal. 
%   |_ 'marker_amd': Motion intensity marker of the AMD algorithm.
%   |_ 'T_amd': Decision signal of the AMD algorithm to which the
%                threshold is compared.
%   |_ 'marker_amvd': Motion intensity marker of the AMVD algorithm.
%   |_ 'T_amvd': Decision signal of the AMVD algorithm to which the
%                threshold is compared.
%   |_ 'marker_ared': Motion intensity marker of the ARED algorithm.
%   |_ 'T_ared': Decision signal of the ARED algorithm to which the
%                threshold is compared.
%   |_ 'marker_shod': Motion intensity marker of the SHOD algorithm.
%   |_ 'T_shod': Decision signal of the SHOD algorithm to which the
%                threshold is compared.
%

% 1) Get the decision signal of the MBGT algorithm and the marker.
T_mbgt = mbgt(input_signal, lwin_mbgt);
marker_mbgt = zeros(1, length(T_mbgt));
for i = 1 : length(T_mbgt)
    if T_mbgt(i) > threshold_mbgt;
       marker_mbgt(i) = 1;
    end
end

% 2) Get the decision signal of the MBCUSUM algorithm and the marker.
T_mbcusum = mbcusum(input_signal, lwin_mbcusum, lambda);
marker_mbcusum = zeros(1, length(T_mbcusum));
for i = 1 : length(T_mbcusum)
    if T_mbcusum(i) > threshold_mbcusum;
       marker_mbcusum(i) = 1;
    end
end

% 3) Get the decision signal of the FSD algorithm and the marker.
[V_fsd, T_fsd] = fsd(input_signal, lwin_fsd, shift_fsd, 512, ...
    threshold_fsd);
[marker_fsd, T_fsd_expanded] = compEstMark(V_fsd, T_fsd, input_signal, ...
    lwin_fsd, shift_fsd);

% 4) Get the decision signal of the LTSD algorithm and the marker.
[V_ltsd, T_ltsd] = ltsd(input_signal, lwin_ltsd, shift_ltsd, 512, ...
    threshold_ltsd);
[marker_ltsd, T_ltsd_expanded] = compEstMark(V_ltsd, T_ltsd, ...
    input_signal, lwin_ltsd, shift_ltsd);

% 5) Get the decision signal of the ARED algorithm and the marker.
T_ared = ared(gx, gy, gz, lwin_ared);
marker_ared = zeros(1, length(T_ared));
for i = 1 : length(T_ared)
    if T_ared(i) > threshold_ared;
       marker_ared(i) = 1;
    end
end

% 6) Get the decision signal of the AMD algorithm and the marker.
T_amd = amd2(ax, ay, az, lwin_amd);
marker_amd = zeros(1, length(T_amd));
for i = 1 : length(T_amd)
    if T_ared(i) > threshold_amd;
       marker_amd(i) = 1;
    end
end

% 7) Get the decision signal of the AMVD algorithm and the marker.
T_amvd = amvd(ax, ay, az, lwin_amvd);
marker_amvd = zeros(1, length(T_amvd));
for i = 1 : length(T_amvd)
    if T_ared(i) > threshold_amvd;
       marker_amvd(i) = 1;
    end
end

% 8) Get the decision signal of the SHOD algorithm and the marker.
T_shod = shod(ax, ay, az, gx, gy, gz, lwin_shod);
marker_shod = zeros(1, length(T_shod));
for i = 1 : length(T_shod)
    if T_ared(i) > threshold_shod;
       marker_shod(i) = 1;
    end
end

end
% END OF BUILD_INT_MARKERS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = mbgt(signal, lwin)

% FUNCTION MBGT applies the Memory Based Graph Theoretic Detector to 
% inertial signals in order to get a profile of the intensity of the
% motion.
%
% - INPUT PARAMETERS: 
%   |_ 'signal': Signal containing the inertial information. This signal
%                could be the magnitude of the acceleration, angular 
%                velocity, magnetic field or a linear combination of them.
%   |_ 'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_ 'T': Decision function.
%
% For more information about the MBGT algorithm see: Daniel Nikovski and 
% Ankur Jain, "Memory-Based Algorithms for Abrupt Change Detection in 
% Sensor Data Streams",Industrial Informatics, 2007 5th IEEE International 
% Conference on , vol.1, no., pp.547-552, 23-27 June 2007 .
%

L = length(signal);
for z = 1 : L - lwin 
    frame = signal(z : z + lwin);
    N = length(frame);
    Cprev = zeros(1, N + 1);
    maxim = 0;
    for i = 1 : N - 1
        beta = zeros(1, N + 1);
        for j = i + 1 : N + i - 1
            if N - i < N - j + i + 1
                beta(N - j + i + 1) = beta(N - j + i + 2) + ...
                    sqrt((frame(N - i) - frame(N - j + i + 1)) ^ 2);
            end
        end
        C = Cprev + beta;
        Cprev = C;
        frame_max = max(C);
        if frame_max > maxim
            maxim = frame_max;
        end
    end
    T(z + lwin) = maxim;
end

end
% END OF MBGT FUNCTION
% 
% % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function T = mbcusum(signal, lwin, lambda)

% FUNCTION MBCUSUM applies the Memory Based CUSUM Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'signal': Signal containing the inertial information. This signal
%                could be the magnitude of the acceleration, angular 
%                velocity, magnetic field or a linear combination of them.
%   |_ 'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_ 'T': Decision function.
%
% For more information about the FSD algorithm see.
%

L = length(signal);
for z = 1 : L - lwin 
    frame = signal(z : z + lwin - 1);
    N = length(frame);
    vprev = zeros(1, N + 1);
    mu = zeros(1, N + 1);
    v = zeros(1, N + 1);
    s = zeros(1, N + 1);
    maxim = 0;
    for i = N - 1 : -1 : 1
        for j = N : -1 : i + 1
            sum_weights = 0;
            if i < j
                for k = j : N
                    value = 1 / (N * (2 * lambda ^ 2 * pi) ^ (1 /2 )) ...
                        * exp(-1 / 2 * (abs(frame(j) - frame(k)) / ...
                        lambda) ^ 2);
                    sum_weights = sum_weights + value;
                end
            end
            mu(j) = sum_weights;
            v(j) = vprev(j) + 1 / (N * (2 * lambda ^ 2 * pi) ^ (1 / 2)) ...
                * exp(-1 / 2 * (abs(frame(j) - frame(i)) / lambda) ^ 2);
            vprev = v;
            s(j) = s(j + 1) + log(mu(j)) - log(v(j)) + log(j - i) - ...
                log(N - j + 1);
        end
        frame_max = max(s);
        if frame_max > maxim
            maxim = frame_max;
        end
    end
     T(z) = maxim;
end

end
% END OF MBCUSUM FUNCTION
% 
% % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function [V, T] = fsd(x, LWIN, SHIFT, NFFT, Threshold)

% FUNCTION FSD applies the Framed Spectrum Detector to inertial signals in 
% order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'x': Signal containing the inertial information. This signal could 
%           be the magnitude of the acceleration, angular velocity, 
%           magnetic field or a linear combination of them.
%   |_ 'LWIN': Length of sliding window.
%   |_ 'SHIFT': Overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%   |_ 'Threshold': Decision threshold.
% 
% - OUTPUT PARAMETERS:
%   |_ 'V': Reduced activity marker.
%   |_ 'T': Decision function.
%
% For more information about the FSD algorithm check the following 
% publication: 'Detection of (In)activity Periods in Human Body Motion 
% Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%

% Computation of the number of frames.
Nframes = floor((length(x) - LWIN) / SHIFT) + 1;

% Computation of spectrum.
SP = spectrum(x, LWIN, SHIFT, NFFT);

% Initialization of the reduced marker.
V = zeros(Nframes, 1);

% Noise averaging initizalization time (frames)
FI = 1;

% Noise estimation.
NE = 0;
alfa = 0.98;

% Initialization of the decision function.
T = zeros(Nframes, 1);

for frame = 1 : Nframes
       
    if (frame <= FI)
       NE = (1 - 1 /frame) * NE + 1 / frame * SP(:, frame);
       V(frame) = 0;
    else
       T(frame) = 10 * log10(1 / (NFFT / 2) * sum(SP(:, frame) .^ 2 ./ ...
           NE .^ 2));
       if (T(frame) > Threshold)
           V(frame) = 1; 
       else
           V(frame) = 0;
           NE = alfa * NE + (1 - alfa) * SP(:, frame);
       end
    end
end

end
% END OF FSD FUNCTION

% % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [y1, y2] = compEstMark(V, T, x, LWIN, SHIFT)

% FUNCTION COMPESTMARK transforms the FSD and LTSD decision signals into a
% marker of the same length of the input signal.
%
% - INPUT PARAMETERS:
%   |_'V': Decision function computed with LTSD or FSD.
%   |_'T': Marker of reduced length.
%   |_'x': Signal used as input to compute the figure of merit with LTSD or
%          FSD, i.e. acceleration or angular rate magnitude, sum of 
%          magnitudes or product of magnitudes.
%   |_'LWIN': Window length used to compute the figure of merit with LTSD 
%             or FSD.
%   |_'SHIFT': Overlapping used to compute the figure of merit with LTSD or
%              FSD.
% 
% - OUTPUT PARAMETERS:
%   |_'y1': Intensity binary marker (extended length).
%   |_'y2': Decision function (extended length).
%

y1 = zeros(1, length(x));
y2 = zeros(1, length(x));

N = min(floor((length(x) - LWIN) / SHIFT) + 1, length(V));

for i = 1 : N
    y2((i - 1) * SHIFT + (1 : LWIN)) = T(i);
   if (V(i) == 0)
          y1((i - 1) * SHIFT + (1 : LWIN)) = 0; 
   else
          y1((i - 1) * SHIFT + (1 : LWIN)) = 1; 
   end
end

end
% END OF COMPESTMARK FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function [V, T] = ltsd(x, LWIN, SHIFT, NFFT, threshold)

% FUNCTION LTSD applies the Long Term Spectrum Detector to inertial signals
% in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'x': Signal containing the inertial information. This signal could 
%           be the magnitude of the acceleration, angular velocity, 
%           magnetic field or a linear combination of them.
%   |_ 'LWIN': Length of sliding window.
%   |_ 'SHIFT': Overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%   |_ 'threshold': Decision threshold.
% 
% - OUTPUT PARAMETERS:
%   |_ 'V': Reduced activity marker.
%   |_ 'T': Decision function.
%
% For more information about the LTSD algorithm check the following 
% publication: 'Detection of (In)activity Periods in Human Body Motion 
% Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%

% No-delay backward long-term spectral envelope.
N = 8;

Nframes = floor((length(x) - LWIN) / SHIFT) + 1;
SP = spectrum(x, LWIN, SHIFT, NFFT);
V = zeros(Nframes, 1);

% Noise averaging initizalization time (frames).
FI = 1;

% Noise estimation.
NE = 0;
alfa = 0.98;

% Decision function.
T = zeros(Nframes, 1);
for frame = 1 : Nframes
       
    if (frame <= FI)
       NE = (1 - 1 / frame) * NE + 1 / frame * SP(:, frame);
       V(frame) = 0;
    else
       N1 = min(frame - 1, N);
       LTSE = max(SP(:, frame - N1 : frame)')'; 
          
       T(frame) = 10 * log10(1 / (NFFT / 2) * sum(LTSE .^ 2 ./ NE .^ 2));
       if (T(frame) > threshold)
           V(frame) = 1; 
       else
           V(frame) = 0;
           NE = alfa * NE + (1 - alfa) * SP(:, frame);
       end
    end
end

end
% END OF LTSD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function T = amd2(ax, ay, az, lwin)

% FUNCTION AMD2 applies the Acceleration Magnitude Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the AMD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
%

% Set noise variance of accelerometer and gyroscope.
var_a = 0.009;

% Compute AMD algorithm.
for i = 1 : length(ax) - lwin     
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);
    N = length(ax_frame);
    mod_aF = sqrt(ax_frame .^ 2 + ay_frame .^ 2 + az_frame .^ 2);
    T(i : i + lwin) = 1 / (N * var_a) * (mod_aF - 1) .^ 2;  
end

end
% END OF AMD2 FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function T = amvd(ax, ay, az, lwin)

% FUNCTION AMVD applies the Acceleration Magnitude Variance Detector to 
% inertial signals in order to get a profile of the intensity of the 
% motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the AMVD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
%

% Set variance of accelerometer noise.
var_a = 0.009;

% Compute AMVD algorithm.
for i = 1 : length(ax) - lwin   
    
    % Build the signal window.
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);

    % Compute the average acceleration values of the window.
    av_x = mean(ax_frame);
    av_y = mean(ay_frame);
    av_z = mean(az_frame);
    N = length(ax_frame);

    % Subtract the average acceleration to each sample of the window.
    axC = ax_frame - av_x;
    ayC = ay_frame - av_y;
    azC = az_frame - av_z;
    
    % Compute the magnitude of the modified acceleration.
    mod_aC = sqrt(axC .^ 2 + ayC .^ 2 + azC .^ 2);

    % Compute the detection signal.
    T(i : i + lwin) = 1 / (N * var_a) * mod_aC .^ 2;   
end

end
% END OF AMVD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function T = ared(gx, gy, gz, lwin)

% FUNCTION ARED applies the Angular Rate Energy Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'gx': Angular rate signal (X axis).
%   |_'gy': Angular rate signal (Y axis).
%   |_'gz': Angular rate signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the ARED algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
%

% Compute ARED algorithm.
for i = 1 : length(gx) - lwin   
    gx_frame = gx(i : i + lwin);
    gy_frame = gy(i : i + lwin);
    gz_frame = gz(i : i + lwin);
    N = length(gx_frame);
    mod_gF = sqrt(gx_frame .^ 2 + gy_frame .^ 2 + gz_frame .^ 2);
    T(i : i + lwin) = 1 / N *(mod_gF) .^ 2;  
end

end
% END OF ARED FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function T = shod(ax, ay, az, gx, gy, gz, lwin)

% FUNCTION SHOD applies the Stance Hypothesis Optimal Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'gx': Angular rate signal (X axis).
%   |_'gy': Angular rate signal (Y axis).
%   |_'gz': Angular rate signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the SHOD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
%

% Set noise variance of accelerometer and gyroscope.
var_a = 0.009;
var_g = 5;

% Compute SHOD algorithm.
for i = 1 : length(ax) - lwin   
    
    gx_frame = gx(i : i + lwin);
    gy_frame = gy(i : i + lwin);
    gz_frame = gz(i : i + lwin);
    
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);
    N = length(ax_frame);
    
    av_x = mean(ax_frame);
    av_y = mean(ay_frame);
    av_z = mean(az_frame);
    mod_av_a = sqrt(av_x ^ 2 + av_y ^ 2 + av_z ^ 2);
    term_ax = ax_frame - 1 * (av_x / mod_av_a);
    term_ay = ay_frame - 1 * (av_y / mod_av_a);
    term_az = az_frame - 1 * (av_z / mod_av_a);
    mod_a = sqrt(term_ax .^ 2 + term_ay .^ 2 + term_az .^ 2) .^ 2;
    
    mod_g = sqrt(gx_frame .^ 2 + gy_frame .^ 2 + gz_frame .^ 2) .^ 2;
    
    T(i : i + lwin) = 1 / N * (1 / var_a * mod_a + 1 / var_g * mod_g);  
    
end

end
% END OF SHOD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function SP = spectrum(x, LWIN, SHIFT, NFFT)

% FUNCTION SPECTRUM computes the spectrum of a given signal. 
%
% - INPUT PARAMETERS:
%   |_ 'x': Input signal. 
%   |_ 'LWIN': Length of the sliding window.
%   |_ 'SHIFT': Length of the overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%
% - OUTPUT PARAMETERS:
%   |_ 'SP': Spectrum of the signal.
%

Nframes = floor((length(x) - LWIN) / SHIFT) + 1;

SP = zeros(NFFT / 2, Nframes);

for frame = 1 : Nframes
    r = (frame - 1) * SHIFT + [1 : LWIN];
    s_w = hamming(LWIN) .* x(r);   
    magnitude = abs(fft(s_w, NFFT));
    SP(:, frame) = magnitude(1 : NFFT / 2);
end

end
% END OF SPECTRUM FUNCTION

% % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 
function [x1, x2, v_output, K_output] = fusionGKF(gyro, obs, freq, ...
    alpha1, alpha2, beta1, beta2, marker)

% FUNCTION FUSIONGKF Applies a Gated Kalman Filter sensor fusion approach 
% to estimate the orientation of a body using the acceleration, magnetic 
% field and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'gyro': vector containing the gyroscope signal for a determined 
%               axis. 
%   |_ 'obs': vector containing the angle values obtained via the linear
%             acceleration relations.
%   |_ 'freq': sampling frecuency. Must be real positive.
%   |_ 'alpha1': Weighing factor which multiplies the variance of the
%                measurement noise when the intensity of the motion is low. 
%                It is a parameter which tunes the filter.
%   |_ 'alpha2': Weighing factor which multiplies the variance of the
%                measurement noise when the intensity of the motion is 
%                high. It is a parameter which tunes the filter.
%   |_ 'beta1': Weighing factor which multiplies the variance of the 
%               process noise when the intensity of the motion is low. It 
%               is a parameter which tunes the filter.
%   |_ 'beta2': Weighing factor which multiplies the variance of the 
%               process noise when the intensity of the motion is high. It 
%               is a parameter which tunes the filter.  
%
% - OUTPUT PARAMETERS:
%   |_ 'x1': estimated orientation angle. (First element of the state
%            vector).
%   |_ 'x2': estimated gyroscope bias. (Second element of the state
%            vector).
%   |_ 'v_output': Vector containing the state vector corrected 'a
%                  posteriori' by the observation.
%   |_ 'K_output': Vector containing the value of the Kalman gain for each
%                  time instant.
%

% 1) Definition of variables.
% -------------------------------------------------------------------------
obsVar = var(obs);
gyroVar = var(gyro);
len = length(gyro);

% Sampling period.
dt = 1 / freq;

% Initialization of covariance matrix.
P = [1 0; 0 1];

% Set measurement variance.
Rk1 = alpha1 * obsVar;
Rk2 = alpha2 * obsVar;

% Assuming process noise is white for each component, covariance matrix of
% process noise must be a diagonal matrix, with noise variance for each
% component in each position of the diagonal.
Q1 = beta1 * [gyroVar 0; 0 obsVar];
Q2 = beta2 * [gyroVar 0; 0 obsVar];

% Set state transition matrix.
A = [1 -dt; 0 1];

% Set measurement matrix. 
C = [1 0];

% Initialize state vector.
X = [0; 0];
x1 = zeros(len, 1);
x2 = zeros(len, 1);

% 2) Body of the Kalman Filter.
% -------------------------------------------------------------------------
for i = 1 : len
    
    if marker(i) == 0
        Rk = Rk1;
        Q = Q1;
    end
    if marker(i) == 1
        Rk = Rk2;
        Q = Q2;
    end
    
    % Prediction phase.
    X(1) = X(1) + (gyro(i) - X(2))*dt;
    X(2) = X(2);
    P = A * P * A' + Q;

    % Update phase.
    v = obs(i) - C * X;
    v_output(i) = v;
    Sk = C * P * C' + Rk;
    K = (P * C') / Sk;
    K_output(i, :) = K';
    X = X + K * v;
    P = P - K * Sk * K';
    x1(i) = X(1);
    x2(i) = X(2);   
end

end
% END OF FUSIONGKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\