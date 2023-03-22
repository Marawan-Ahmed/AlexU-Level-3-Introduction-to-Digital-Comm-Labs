% this part below must be copied to your m file and complete the %
clear
clc

% reconstruction from oversampling
t=0:0.001:1;% time signal
y=2*cos(2*pi*5*t);

[B,A] = butter(3,1000/100000,'low' ); % butter fly filter 
%a Butterworth low-pass filter is designed using the BUTTER function. 
%The filter has a order of 3 and a cutoff frequency of 1000 Hz (normalized with respect to the sampling frequency of 100 kHz). 
%The filter coefficients are stored in the vectors B and A.
zero_added_signal=zeros(1,length(y)*10); %zero padding like dsp

for i=1:length(y)
    zero_added_signal(i*10)=y(i); %the new signal (padded signal)
end

zero_added_signal(1:9)=[]; %To display the signal with a higher sampling rate without changing its spectrum,
				%the first 9 samples of the zero-added signal are removed

% Adding zeros enhances the signal display and don't change the
%spectrum,it changes sampling freq. only (just higher resolution)

t=linspace(0,1,length(zero_added_signal)); %new time vector - 9
filtered_signal = filter(B,A,zero_added_signal); % pass the filter coefficients A B  and the over sampled signal
plot(t,filtered_signal,'r' )
% XLABEL('time')
% YLABEL('oversampled signals')

%The filter function is a general-purpose function for applying digital filters to signals. 
%It can be used with any filter design method, not just Butterworth filters. 
%The BUTTER function, on the other hand, is a specific function for designing Butterworth filters.


%% construction from minimum sampling
figure
t=0:0.1:1; % replace ?? with the suitable number 
%considering the Nyquist rate in the signal to avoid aliasing. In this case, the signal has a frequency of 5 Hz, so the Nyquist rate is 10 Hz.
%Therefore, the step size should be no larger than 0.1 seconds

y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.1,'low' ); %cutoff = 0.1 (normalized wrt sampling frequency)
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
plot(t,filtered_signal,'r' )
% XLABEL('time')
% YLABEL('minimum sampled signals')

s=fft(filtered_signal);
s=fftshift(s); 
fs=100; % why 100?? Write your comments in the m file
%If the original signal has a frequency of 10 Hz and it was resampled by a factor of 10,
% then the resulting signal would have a sampling frequency of 100 Hz


freq=linspace(-fs/2,fs/2,length(s)); %length of spectrum
figure
plot(freq,abs(s)) 
% XLABEL('freq')
% YLABEL('magnitude of minimum sampled signals')
%% construction from undersampling sampling
figure
t=0:0.2:1;
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.2,'low' );
% complete this part as shown in the construction from minimum sampling
%and do the necessary changes , you have to do low pass filtering and %

zero_added_signal=zeros(1,length(y)*10);

for i=1:length(y)
    zero_added_signal(i*10)=y(i);
end

zero_added_signal(1:9)=[];
% Adding zeros enhances the signal display and don't change the
%spectrum,it changes sampling freq. only

t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
plot(t,filtered_signal,'r' )

s=fft(filtered_signal);
s=fftshift(s);
fs=100;


freq=linspace(-fs/2,fs/2,length(s));

figure
plot(freq,abs(s))
 XLABEL('freq')
 YLABEL('magnitude of under sampled signals')