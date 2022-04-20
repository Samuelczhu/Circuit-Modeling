%% Q5 Circuit with Noise
% In Q5, a current souce In in parallel with R3 is added to model the
% thermal noise generated in the resistor R3. A capacitor Cn in parallel
% with the resistor is also added to bandwidth limit the noise.
% TODO: is this time domain simulation? do we need to simulate all different
% input cases?
simTime = 1;  % Simulate for 1 second
simSteps = 1000;  % Use 1000 steps
% Calculate the deltaT
deltaT = simTime/simSteps;  % second/step
% Declare the time vector
vecTime = linspace(0,1,simSteps);

% Declare the capacitor Cn
Cn = 0.00001;
magIn = 0.001;  % Magnitude for In

% Declare the input for a step from 0 to 1 at 0.03s
vecInputV = zeros(1, simSteps);
vecInputV(0.03*simSteps:simSteps) = 1;
% Generate the vector for noise current
% TODO: verify that this is Guassian distribution for the random number
% vecIn = normrnd(magIn,10,[1,5]);  % parameter: mean, std, [dimension]
vecIn = magIn*randn(1, simSteps);
% Hold the output vector
vecOutputV = zeros(1, simSteps);

% Declare the vectors 
vectorV = zeros(10, 1);  % solution vector: [N1, N2, N3, N4, N5, I1, IL, I3, I4, In]
vectorF = zeros(10, 1);  % F vector: F(1) = VIN, F(10) = In

% Declare the C matrix
matrixC = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           C, -C, 0, 0, 0, 0, 0, 0, 0, 0;
          -C, C, 0, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, -L, 0, 0, 0;
           0, 0, Cn, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% Declare the G matrix
matrixG = [1,     0,        0,    0,          0, 0,  0,      0, 0, 0;
          1/R1, -1/R1,      0,    0,          0, 1,  0,      0, 0, 0;
          -1/R1, 1/R1+1/R2, 0,    0,          0, 0,  1,      0, 0, 0;
            0,    1,       -1,    0,          0, 0,  0,      0, 0, 0;
            0,    0,        0,    0,          0, 0, -1,      1, 0, 1;
            0,    0,    -1/R3,    0,          0, 0,  0,      1, 0, 0;
            0,    0,        0, 1/R4,      -1/R4, 0,  0,      0, 1, 0;
            0,    0,        0,    1,          0, 0,  0, -alpha, 0, 0;
            0,    0,        0, -1/R4, 1/R4+1/RO, 0,  0,      0, 0, 0;
            0,    0,        0,    0,          0, 0,  0,      0, 0, 1];

% Construct the A matrix
matrixA = matrixC/deltaT + matrixG;

% Loop through the simulation
for iSim = 1:simSteps
    % Update the F vector for Vin and In
    vectorF(1) = vecInputV(iSim);
    vectorF(10) = vecIn(iSim);
    % Update the V vector
    vectorV = matrixA^-1 * (matrixC * vectorV / deltaT + vectorF);
    % Save the output voltage
    vecOutputV(iSim) = vectorV(5);
end

% Plot of completed transient simulation for step input
figure(14)
% Time domain plot
subplot(1,2,1)
plot(vecTime, vecInputV, "-b.")  % Vin versus time
hold on
plot(vecTime, vecOutputV, "-r.")  % Vo versus time
hold off
title("Transient simulation for a step input with noise")
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vin versus time", "Vo versus time")
grid on
% Frequency domain plot (fft)
subplot(1,2,2)
% TODO: verify fft!
fftVin = 20*log10(abs(fftshift(fft(vecInputV)))); % Input fft in dB
plot(fftVin, "-b.")  % Plot the input fft
hold on 
fftVo = 20*log10(abs(fftshift(fft(vecOutputV)))); % Output fft in dB
plot(fftVo, "-r.")  % Plot the output fft
hold off
title("FFT of the step input with noise")
xlabel("Frequency (Hz)")
ylabel("V (dB)")
legend("Vin versus time", "Vo versus time")
grid on
snapnow