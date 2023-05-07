
% Generate random bits
num_bits = 100;
bits = generate_random_bits(num_bits);

% Modulate the bits using different line codes
waveform_nrz = nrz_modulation(bits);
waveform_nrzi = nrzi_modulation(bits);
waveform_rz = rz_modulation(bits);
waveform_ami = ami_modulation(bits);
waveform_manchester = manchester_modulation(bits);
%waveform_mlt3 = mlt3_modulation(bits);


% Plot the line code modulations
time = (0 : num_bits-1) * 2;  % Assuming each bit lasts for 2 units of time

figure;
stem(time, bits, 'filled', 'MarkerSize', 5);
xlabel('Time');
ylabel('Voltage');
title('The random generated bit sequence');
ylim([-0.5, 1.5]);  % Set y-axis limits to show only 0 and 1 values

figure;
subplot(6, 1, 1);
plot(time, waveform_nrz);
xlabel('Time');
ylabel('Voltage');
title('NRZ Modulation');

subplot(6, 1, 2);
plot(time, waveform_nrzi);
xlabel('Time');
ylabel('Voltage');
title('NRZI Modulation');

resampled_waveform_rz = interp1(linspace(0, 200, 200), waveform_rz, time);
subplot(6, 1, 3);
plot(time, resampled_waveform_rz);
xlabel('Time');
ylabel('Voltage');
title('RZ Modulation');

subplot(6, 1, 4);
plot(time, waveform_ami);
xlabel('Time');
ylabel('Voltage');
title('AMI Modulation');

resampled_waveform_manchester = interp1(linspace(0, 200, 200), waveform_manchester, time);
subplot(6, 1, 5);
plot(time, resampled_waveform_manchester);
xlabel('Time');
ylabel('Voltage');
title('Manchester Modulation');

%subplot(6, 1, 6);
%plot(time, waveform_mlt3);
%xlabel('Time');
%ylabel('Voltage');
%title('MLT-3 Modulation');


%axis off;

% Adjust the spacing between subplots
%subplot_spacing = 0.04;
%subplot_position = get(gcf, 'Position');
%subplot_position(4) = subplot_position(4) - subplot_spacing;
%set(gcf, 'Position', subplot_position);


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
    waveform = repelem(bits, 2);
    waveform(waveform == 0) = -1;  % Map 0 to -1
end

% Function for Alternative Mark Inversion (AMI) modulation
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
        if bits(i) == 0
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
    waveform = zeros(1, length(bits) * 2)';
    current_level = 1;
    for i = 1:length(bits)
        if bits(i) == 0
            waveform(2*i-1:2*i) = waveform(2*i-2);
        else
            waveform(2*i-1:2*i) = -waveform(2*i-2) * current_level;
            current_level = -current_level;
        end
    end
end
