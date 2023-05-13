close all;
clear;

% Simulation parameters
n = 1e4; % number of bits 
SNR = 0:2:30; % signal-to-noise ratio range (in dB)
m = 20; % number of samples that represents waveform
T = 20; % sampling instant
s1 = ones(1,m); % rectangular signal with amplitude 1
s2 = zeros(1,m); % zero signal
matched_filter = (s1 - s2); % The fliplr function is used to flip the waveform so that it can be used as a filter.
matched_filter = matched_filter(end:-1:1); % The fliplr function is used to flip the waveform so that it can be used as a filter.

% Generate random binary data vector
data = randi([0 1],1,n) > 0.5;

% Represent each bit with proper waveform
waveform = [];
for i = 1:length(data)
    if data(i) == 0
        waveform = [waveform s2]; % concatenate zero signal for bit value 0
    else
        waveform = [waveform s1]; % concatenate rectangular pulse for bit value 1
    end
end

% Apply noise to samples
signal_power = sum(waveform.^2)/n;
for i = 1:length(SNR)
    noise_power = signal_power/(10^(SNR(i)/10));
    noise_signal = sqrt(noise_power) * randn(size(waveform));
    Rx_signal = waveform + noise_signal;
    
    % Decide whether the Rx_sequence is '1' or '0' by comparing the samples with threshold
    simple_Rx = simple_receiver(Rx_signal,m,mean(Rx_signal),T-1);
    % Perform convolution process in the receiver
    y = filter(matched_filter,1,Rx_signal);

    % Sample the output of the matched filter;
    mf_Rx = simple_receiver(y,m,mean(y),T-1);
    
    collerator_Rx = collerator_receiver(Rx_signal,m,s1,s2,T-1);
    % Compare the original bits with the detected bits and calculate number of errors
    BER_simple_Rx(i) = ComputeBER(data,simple_Rx);
    BER_mf(i) = ComputeBER(data,mf_Rx);
    BER_collerator(i) = ComputeBER(data,collerator_Rx);

end

% Plot the BER curve against SNR (use semilogy)
semilogy(SNR,BER_simple_Rx,'-o')
xlabel('Signal-to-Noise Ratio (dB)')
ylabel('Bit Error Rate')
% xlim([0 SNR(end)])
% ylim([0 0.5])
title('Matched filter vs. simple Rx')
hold on

semilogy(SNR,BER_mf,'-o')

hold off
figure;
% Plot the BER curve against SNR (use semilogy)
semilogy(SNR,BER_simple_Rx,'-o')
xlabel('Signal-to-Noise Ratio (dB)')
ylabel('Bit Error Rate')
% xlim([0 SNR(end)])
% ylim([0 0.5])
title('Collerator vs. simple Rx')
hold on

semilogy(SNR,BER_collerator,'-o')

hold off
function BER = ComputeBER(bit_seq,rec_bit_seq)
    % Use bitxor to compute the bit-wise XOR of x and y
    BER= sum(bitxor(bit_seq,rec_bit_seq)) / length(bit_seq);
end

function rec_bit_seq = simple_receiver(rec_sample_seq,m,v_th,T)
    number_bits = length(rec_sample_seq)/m;
    rec_bit_seq = zeros(1 ,number_bits);

    for i = 1:number_bits
        if rec_sample_seq((i-1)*m+ T +1) >= v_th
            rec_bit_seq(i) = 1;
        end
    end
end
function rec_bit_seq = collerator_receiver(rec_sample_seq,m,s1,s2,T)
    number_bits = length(rec_sample_seq)/m;
    for i = 1:number_bits
        temp = rec_sample_seq((i-1)*m+1:i*m);
        temp = temp.*(s1-s2);
        temp = sum(temp);
        rec_bit_seq(i) = temp;
    end
    
    v_th = mean(rec_bit_seq);
    rec_bit_seq = rec_bit_seq > v_th;
end