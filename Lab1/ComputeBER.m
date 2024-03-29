function BER = ComputeBER(bit_seq,rec_bit_seq)
%
% Inputs:
%   bit_seq:     The input bit sequence
%   rec_bit_seq: The output bit sequence
% Outputs:
%   BER:         Computed BER
%
% This function takes the input and output bit sequences and computes the
% BER

%%% WRITE YOUR CODE HERE
number_error_bits = 0;
number_bits = length(bit_seq);

for i = 1:number_bits
    if (bit_seq(i) ~= rec_bit_seq(i))
        number_error_bits = number_error_bits + 1;
    end
end
BER = number_error_bits / number_bits;
%%%
