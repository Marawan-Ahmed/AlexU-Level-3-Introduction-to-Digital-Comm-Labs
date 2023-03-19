clear
clc

% reconstruction from oversampling
t=0:0.001:1;% time signal
y=2*cos(2*pi*5*t);

[B,A] = butter(3,1000/100000,'low' ); % butter fly filter

zero_added_signal=zeros(1,length(y)*10); %zero padding like dsp

for i=1:length(y)
    zero_added_signal(i*10)=y(i); %the new signal (padded signal)
end

zero_added_signal(1:9)=[]; 

t=linspace(0,1,length(zero_added_signal)); 
filtered_signal = filter(B,A,zero_added_signal); 

% Quantize the filtered signal to 8 bits
bits = 8;
quantized_signal = quantize(filtered_signal, linspace(-2, 2, 2^bits));

plot(t, quantized_signal,'r' )
title(sprintf('Quantized signal (%d bits)', bits));
xlabel('time')
ylabel('quantized signal')


% reconstruction from minimum sampling
t = 0:0.1:1;
y = 2*cos(2*pi*5*t);

[B,A] = butter(10,0.1,'low');
zero_added_signal = zeros(1,length(y)*10);
for i = 1:length(y)
    zero_added_signal(i*10) = y(i);
end
zero_added_signal(1:9) = [];

t = linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);

% Quantization
bits = 8; % number of bits
levels = 2^bits; % number of quantization levels
max_val = max(filtered_signal);
min_val = min(filtered_signal);
step = (max_val - min_val) / levels;
quantized_signal = round((filtered_signal - min_val) / step) * step + min_val;

% Plotting
figure
plot(t, filtered_signal, 'r')
hold on
plot(t, quantized_signal, 'b')
hold off
xlabel('Time')
ylabel('Signal')
legend('Filtered Signal', 'Quantized Signal')

s = fft(filtered_signal);
s = fftshift(s); 
fs = 100; % sampling frequency
freq = linspace(-fs/2,fs/2,length(s)); % frequency vector

figure
plot(freq, abs(s))
xlabel('Frequency')
ylabel('Magnitude of minimum sampled signals')

% reconstruction from undersampling sampling
t=0:0.2:1;
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.2,'low' );

zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];

t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);

% Quantize the signal using 8 bits
quantized_signal = quantize(filtered_signal, 8);

% Plot the quantized signal
plot(t, quantized_signal, 'r')

s=fft(quantized_signal);
s=fftshift(s);
fs=100;

freq=linspace(-fs/2,fs/2,length(s));

% Plot the magnitude spectrum of the quantized signal
figure
plot(freq,abs(s))
xlabel('freq')
ylabel('magnitude of under sampled signals')


