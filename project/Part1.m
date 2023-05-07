% Simulation parameters
n = 1e3; % number of bits
SNR = 0:2:30; % signal-to-noise ratio range (in dB)
m = 20; % number of samples that represents waveform
T = 20; % sampling instant
s1 = ones(1,m); % rectangular signal with amplitude 1
s2 = zeros(1,m); % zero signal

% Generate random binary data vector
data = randi([0 1],1,n);

% Represent each bit with proper waveform
% waveform = repmat(data,m,1);
% waveform = reshape(waveform',[],1);
waveform = [];
for i = 1:length(data)
    if data(i) == 0
        waveform = [waveform s2]; % concatenate zero signal for bit value 0
    else
        waveform = [waveform s1]; % concatenate rectangular pulse for bit value 1
    end
end

% Apply noise to samples
for i = 1:length(SNR)
    snr = SNR(i);
    noise_power = sum(data.^2)/10^(snr/10);
    noise = sqrt(noise_power)*randn(size(waveform));
    Rx_signal = waveform + noise;
    
    % Decide whether the Rx_sequence is '1' or '0' by comparing the samples with threshold
    threshold = sum(s1)/2;    
    simple_Rx = simple_receiver(Rx_signal,T);
    
    % Perform convolution process in the receiver

    
    % Sample the output of the matched filter

    
    % Compare the original bits with the detected bits and calculate number of errors
    num_errors(i) = ComputeBER(data,simple_Rx);
end

% Save the probability of error of each SNR in matrix, BER
BER = num_errors/n;

% Plot the BER curve against SNR (use semilogy)
semilogy(SNR,BER,'-o')
xlabel('Signal-to-Noise Ratio (dB)')
ylabel('Bit Error Rate')
title('Bit Error Rate vs. Signal-to-Noise Ratio')



function BER = ComputeBER(bit_seq,rec_bit_seq)
    number_error_bits = 0;
    number_bits = length(bit_seq);

    for i = 1:number_bits
        if (bit_seq(i) ~= rec_bit_seq(i))
            number_error_bits = number_error_bits + 1;
        end
    end
    BER = number_error_bits / number_bits;
end

function rec_bit_seq = simple_receiver(rec_sample_seq,T)
    number_bits = length(rec_sample_seq)/T;
    rec_bit_seq = zeros(1 ,number_bits);

    for i = 1:number_bits
        temp = 0;
        for j = 1:T
           temp = temp + rec_sample_seq((i-1)*T+j);
        end
        if temp >= T/2
            rec_bit_seq(i) = 1;
        end
    end
end