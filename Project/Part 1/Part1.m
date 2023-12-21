clear all
close all

%% Generate pulses
pulse_width = 1/100000;
start_time = -5*pulse_width;
end_time = 5*pulse_width;
fs = 10e6; %use high sampling frequency to approximate an ideal pulse
t = start_time:1/fs:end_time;

rect_pulse1 = zeros(size(t)); % initialize the y vector
rect_pulse1(t >= 0 & t <= pulse_width) = 1; % assign the value 2 to the pulse interval

rect_pulse2 = zeros(size(t)); % initialize the y vector
rect_pulse2(t >= 2 * pulse_width & t <= 3 * pulse_width) = 1; % assign the value 2 to the pulse interval

figure
plot(t,rect_pulse1) % plot the pulse
xlabel('x')
ylabel('y')
title('Rectangular pulse of height 2 and width 0.5 at x = 1')
hold on
plot(t,rect_pulse2) % plot the pulse
hold off
%% getting the frequency spectrum of the sugnal
f = linspace(-fs/2,fs/2,length(t));

% Get the Fourier transform of the pulse
FT_rect_pulse1 = fftshift(fft (double (rect_pulse1)));
FT_rect_pulse2 = fftshift(fft (double (rect_pulse2)));

% Plot the magnitude and phase spectra
figure
subplot(2,2,1)
plot(f,abs(FT_rect_pulse1))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the first pulse')

subplot(2,2,2)
plot(f,angle(FT_rect_pulse1))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the first pulse')

subplot(2,2,3)
plot(f,abs(FT_rect_pulse2))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the second pulse')

subplot(2,2,4)
plot(f,angle(FT_rect_pulse2))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the second pulse')

%% Generate a filter and apply it
% Assign the cutoff frequency b
bandwidth = 100e3; % you can change this value accordingly

% Create a lowpass filter function
LPF = zeros(size(f)); % initialize the y vector
LPF(f >= -bandwidth & f <= bandwidth) = 1; % assign the value 2 to the pulse interval

% Plot the magnitude and phase spectra
figure
plot(f,LPF)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the first pulse')


% Apply the lowpass filter to the fft output
FT_filtered_pulse1 = FT_rect_pulse1.*LPF; % Filtered fft output
FT_filtered_pulse2 = FT_rect_pulse2.*LPF; % Filtered fft output


% Plot the magnitude and phase spectra
figure
subplot(2,2,1)
plot(f,abs(FT_filtered_pulse1))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the first pulse')

subplot(2,2,2)
plot(f,angle(FT_filtered_pulse1))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the first pulse')

subplot(2,2,3)
plot(f,abs(FT_filtered_pulse2))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the second pulse')

subplot(2,2,4)
plot(f,angle(FT_filtered_pulse2))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the second pulse')

% Transform the filtered signal back to the time domain
filtered_pulse1 = ifft(ifftshift(FT_filtered_pulse1),length(t)); % Filtered signal
filtered_pulse2 = ifft(ifftshift(FT_filtered_pulse2),length(t)); % Filtered signal

figure
plot(t,filtered_pulse1) % plot the pulse
xlabel('x')
ylabel('y')
title('Rectangular pulse of height 2 and width 0.5 at x = 1')
hold on
plot(t,filtered_pulse2) % plot the pulse
hold off


%% Mtigating the ISI problem using raised cosine
% Parameters
beta = 1; % Rolloff factor


% Generate the raised cosine filter coefficients
g = zeros(size(t)); % Initialize the filter coefficients
% ind1 = find(t == pulse_width/(2*beta)); % Find the indices of the singularities
% ind2 = find(t == -pulse_width/(2*beta));
g(t == pulse_width/(2*beta)| t == -pulse_width/(2*beta)) = (pi/4*pulse_width).*(sinc(1/(2*beta))); % Assign the values at the singularities
g(t ~= pulse_width/(2*beta)| t ~= -pulse_width/(2*beta)) = (1/pulse_width).*(sinc(t/pulse_width)).*(cos(pi * beta * t / pulse_width)) ./ (1 - (2*beta*t/pulse_width).^2); % Assign the values at the non-singularities
g = g/max(g); % Normalize the filter coefficients

figure
plot(t,g) % plot the pulse
xlabel('x')
ylabel('y')
title('RC')

G = fftshift(fft (double (g)));

% Generate the data
%% getting the frequency spectrum of the sugnal
f = linspace(-fs/2,fs/2,length(t));

% Get the Fourier transform of the pulse
FT_rect_pulse1 = fftshift(fft (double (rect_pulse1)));
FT_rect_pulse2 = fftshift(fft (double (rect_pulse2)));

% Apply the Pulse shping to the fft output
FT_RC_pulse1 = FT_rect_pulse1.*G; % Filtered fft output
FT_RC_pulse2 = FT_rect_pulse2.*G; % Filtered fft output


% Plot the magnitude and phase spectra
figure
subplot(3,2,1)
plot(f,abs(FT_RC_pulse1))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the first pulse')

subplot(3,2,2)
plot(f,angle(FT_RC_pulse1))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the first pulse')

subplot(3,2,3)
plot(f,abs(FT_RC_pulse2))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the second pulse')

subplot(3,2,4)
plot(f,angle(FT_RC_pulse2))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the second pulse')

subplot(3,2,5)
plot(f,abs(G))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of the second pulse')

subplot(3,2,6)
plot(f,angle(G))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Phase spectrum of the second pulse')

% Transform the filtered signal back to the time domain
RC_pulse1 = ifft(ifftshift(FT_RC_pulse1),length(t)); % Filtered signal
RC_pulse2 = ifft(ifftshift(FT_RC_pulse2),length(t)); % Filtered signal

figure
plot(t,RC_pulse1) % plot the pulse
xlabel('x')
ylabel('y')
title('RC')
hold on
plot(t,RC_pulse2) % plot the pulse
hold off



