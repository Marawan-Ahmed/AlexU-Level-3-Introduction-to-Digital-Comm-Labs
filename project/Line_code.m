%% Setup

% Generate random bits
num_bits = 10;

% Set the pulse duration
pulse_duration = 1;  % Duration of each pulse in seconds

% set redundancy to display abrupt changes 
redundancy = 100;

% Display the bit sequence
bits = generate_random_bits(num_bits);
%bits = [0,1,0,1,0,1,0,1,0,1];
%disp ("generated bit sequence:"); disp(bits);

% Generate the time vector
num_samples = length(bits) * redundancy;
time = linspace(0, pulse_duration * length(bits), num_samples); 

% Sampling time
Ts = 1 / redundancy ; % As n bits repesent n seconds

% Plot the random generated bit sequence
% Generate the voltage vector like abrupt changes
bits_volt = repelem(bits, redundancy); 
% About repelem(A,B) function
% A: is the input array to be repeated
% B: is the number of times each element in A should be repeated.

figure;
plot(time, bits_volt);
xlabel('Time');
ylabel('Voltage');
title('The random generated bit sequence');
ylim([-0.5, 1.5]);  % Set y-axis limits to show only 0 and 1 values
saveas(gcf, 'input_bits.png');

%% Line codes Modulation

% Modulate the bits using different line codes
waveform_nrz = nrz_modulation(bits);    %polar
waveform_nrzi = nrzi_modulation(bits);  %polar
waveform_rz = rz_modulation(bits);      %polar RZ
waveform_ami = ami_modulation(bits);    %Bipolar NRZ
waveform_manchester = manchester_modulation(bits);
waveform_mlt3 = mlt3_modulation(bits);
 %
% Plot different line codes modulation
figure;

waveform_nrz_volt = repelem(waveform_nrz, redundancy); 
subplot(6, 1, 1);
plot(time, waveform_nrz_volt);
title('NRZ Modulation');

waveform_nrzi_volt = repelem(waveform_nrzi, redundancy); 
subplot(6, 1, 2);
plot(time, waveform_nrzi_volt);
title('NRZI Modulation');

waveform_rz_volt = repelem(waveform_rz, redundancy/2); % 50 repetitions to fit with time vector
subplot(6, 1, 3);
plot(time, waveform_rz_volt);
title('RZ Modulation');

ylabel('Voltage');

waveform_ami_volt = repelem(waveform_ami, redundancy); 
subplot(6, 1, 4);
plot(time, waveform_ami_volt);
title('AMI Modulation');

waveform_manchester_voltage = repelem(waveform_manchester, redundancy/2); % 50 repetitions to fit with time vector
subplot(6, 1, 5);
plot(time, waveform_manchester_voltage);
title('Manchester Modulation');

waveform_mlt3_volt = repelem(waveform_mlt3, redundancy);
subplot(6, 1, 6);
plot(time, waveform_mlt3_volt);
title('MLT-3 Modulation');

xlabel('Time');

saveas(gcf, 'Modulations.png');
 %}
%% Calculate power spectral density

%{
About periodogram()

periodogram is a non-parametric function that calculates the PSD directly 
from the samples of the signal. 
As input signal is divided into overlapping segments (based on the window) 
and for each segment FFT is applied to obtain the spectrum of the segment 
then we square magnitude of the spectrum to obtain the power spectrum for 
each segment.
all segments are averaged together to obtain an estimate of the PSD.

takes 5 parameters:
1-waveform to estimate its PSD
2-The window taken for each segment to apply fft on it (default is rect)
3-Number of FFT points (the default value is the next power of 2 greater
  than the length of the input signal.)
4- Sampling frequency 
5- center the output (optional)

returns the power and frequency vectors
%}

[P_nrz, f_nrz] = periodogram(waveform_nrz,[],[],1/Ts,'centered');
[P_nrzi, f_nrzi] = periodogram(waveform_nrzi,[],[],1/Ts,'centered');
[P_rz, f_rz] = periodogram(waveform_rz,[],[],1/Ts,'centered');
[P_ami, f_ami] = periodogram(waveform_ami,[],[],1/Ts,'centered');
[P_manchester, f_manchester] = periodogram(waveform_manchester,[],[],1/Ts,'centered');
[P_mlt3, f_mlt3] = periodogram(waveform_mlt3,[],[],1/Ts,'centered');

% Plot power spectral density
figure;

subplot(6,1,1); 
plot(f_nrz,10*log10(P_nrz)); 
title('Non-return to zero');

subplot(6,1,2); 
plot(f_nrzi,10*log10(P_nrzi)); 
title('Non-return to zero inverted');

subplot(6,1,3); 
plot(f_rz,10*log10(P_rz)); 
title('Return to zero');

ylabel('Power (dB)');

subplot(6,1,4); 
plot(f_ami,10*log10(P_ami)); 
title('Alternate mark inversion');

subplot(6,1,5); 
plot(f_manchester,10*log10(P_manchester)); 
title('Manchester coding');

subplot(6,1,6); 
plot(f_mlt3,10*log10(P_mlt3)); 
title('Multi-level transmission3');

xlabel('Frequency (Hz)');

saveas(gcf, 'PSD.png');

%% Functions
% Function to generate random bits
function bits = generate_random_bits(num_bits)
    bits = randi([0, 1], 1, num_bits);
end

% Function for Non-Return to Zero (NRZ) modulation
function waveform = nrz_modulation(bits)
    waveform = 2 * bits - 1;  % Map 0 to -1 and 1 to +1
end

% Function for Non-Return to Zero Inverted (NRZI) modulation
function waveform = nrzi_modulation(bits)
    waveform = ones(size(bits));
    for i = 2:length(bits)
        waveform(i) = waveform(i-1) * (-1)^(bits(i));
    end
end

% Function for Return to Zero (RZ) modulation
function waveform = rz_modulation(bits)
    waveform = zeros(1, length(bits) * 2);
    for i = 1:length(bits)
        if bits(i) == 1
            waveform(2*i-1) = 1;
            %waveform(2*i) = 0;
        else
            waveform(2*i-1) = -1;
            %waveform(2*i) = 0;
        end
    end
end

% Function for Alternative Mark Inversion (AMI) modulation (Bipolar NRZ)
function waveform = ami_modulation(bits)
    waveform = zeros(size(bits));
    last_polarity = 0;
    for i = 1:length(bits)
        if bits(i) == 1
            waveform(i) = (-1)^last_polarity;
            last_polarity = ~last_polarity;
        end
    end
end

% Function for Manchester coding modulation
function waveform = manchester_modulation(bits)
    waveform = zeros(1, length(bits) * 2);
    for i = 1:length(bits)
        if bits(i) == 1
            waveform(2*i-1) = 1;
            waveform(2*i) = -1;
        else
            waveform(2*i-1) = -1;
            waveform(2*i) = 1;
        end
    end
end

% Function for Multilevel Transmission 3 (MLT-3) modulation
function waveform = mlt3_modulation(bits)
    waveform = zeros(1, length(bits));
    level = 0;  % Initial level
    past_level = -1;

    for i = 1:length(bits)
        if bits(i) == 0
            waveform(i) = level;  % No transition, maintain current level
        else
            if (level == 1)
                waveform(i) = 0;
                level = 0;
                past_level = 1;

            elseif (level == -1)
                waveform(i) = 0;
                level = 0;
                past_level = -1;

            elseif (level == 0)
                if (past_level == -1)
                    waveform(i) = 1;
                    level = 1;

                elseif (past_level == 1)
                    waveform(i) = -1;
                    level = -1;

                end
                past_level = 0;
            end
        end
    end
end

