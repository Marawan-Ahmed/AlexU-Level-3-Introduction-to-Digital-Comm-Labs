function rec_bit_seq = DecodeBitsFromSamples(rec_sample_seq,case_type,fs)
%
% Inputs:
%   rec_sample_seq: The input sample sequence to the channel
%   case_type:      The sampling frequency used to generate the sample sequence
%   fs:             The bit flipping probability
% Outputs:
%   rec_sample_seq: The sequence of sample sequence after passing through the channel
%
% This function takes the sample sequence after passing through the
% channel, and decodes from it the sequence of bits based on the considered
% case and the sampling frequence

if (nargin <= 2)
    fs = 1;
end

switch case_type
    
    case 'part_1'
        %%% WRITE YOUR CODE FOR PART 1 HERE [DONE]
        rec_bit_seq = rec_sample_seq;
        %%%
    case 'part_2'
        %%% WRITE YOUR CODE FOR PART 2 HERE [DONE]
        number_bits = length(rec_sample_seq)/fs;
        rec_bit_seq = zeros(1 ,number_bits);
        
        for i = 1:number_bits
            temp = 0;
            for j = 1:fs
               temp = temp + rec_sample_seq((i-1)*fs+j);
            end
            if temp >= fs/2
                rec_bit_seq(i) = 1;
            end
        end
        %%%
    case 'part_3'
        %%% WRITE YOUR CODE FOR PART 3 HERE [DONE]
        
        % number of received bits
        number_bits = length(rec_sample_seq)/fs;
        
        % number of bit copies in each redundancy block
        number_copies = 5; 
        
        % decoded sequence (after error correction)
        rec_bit_seq = zeros(1 ,number_bits);
        
        % probability of 1 in each redundancy block
        prob_of_one = 0;
        
        % index of the current bit in the redundancy block being processed
        copy_index = 0;
        
        
        %%%
        % The transmitter sends a sequence of "fs" ones in order to send 
        % a "1" or "fs" bits of zeros to send a "0". For the received 
        % fs bits to be interpreted as a "1", the number of successfully 
        % identified ones in the fs-long sequence must be greater than 
        % the number of zeros. Part_3, differs from part_1 as it adds an
        % extra error correction technique usimg redudancy blocks. To send 
        % a "1" the transmitter sends a number of copies of the "1", 5 
        % copies for instance, and the receiver identifies the received
        % block as a "1" only if the probability of "1" in the redundancy
        % block is greater than the probability of "0".
        %%%
   
        for i = 1:number_bits
            
            temp = 0;
            copy_index = copy_index + 1;
            
            % count the number of ones in the fs-long sequence 
            for j = 1:fs
                
               temp = temp + rec_sample_seq((i-1)*fs+j);
            end
            
            if temp >= fs/2
       
                rec_bit_seq(i) = 1;
                
                % calculate the prob. of one in the current redundancy
                % block to be used in error correction (
                prob_of_one = prob_of_one + (1/number_copies);
            end

            %%% correct errors by comparing the redundant bits:
            % ie: if the first receieved five bits are "01101"
            % then the intended value was most probably "1"
            % as the probability of "1" is higher than "0".
            % thus, the bits are corrected to be "11111"

            if(copy_index == number_copies)
                
                if(prob_of_one >= 0.5)
                    
                    % correct the previous copies to be "1"
                    while(copy_index >= 1)
                        
                        rec_bit_seq(i - (copy_index-1))= 1;
                        copy_index = copy_index - 1;
                    end
                    
                    % reset the prob_of_one for the next redundancy block
                    prob_of_one = 0;
                end
                
                if(prob_of_one < 0.5)
                    % correct the previous copies to be "0"
                    while(copy_index >= 1)
           
                        rec_bit_seq(i - (copy_index-1))= 0;
                        copy_index = copy_index - 1;
                    end
                    
                    % reset the prob_of_one for the next redundancy block
                    prob_of_one = 0;
                end
            end
        end
    end
    
        %%%
end