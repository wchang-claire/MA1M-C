function [BER]= test_channel(mod_input,Nbps_input,EbN0_input)
%% input parameters
%   clear; clc; close all;
%   mod_input = 'qam'; % ['pam','qam','qpsk','bpsk']
%   Nbps_input = 4; %[1,4,6,2,1]
%   EbN0_input = 8;

%% modulation parameters
Modu.mod = mod_input; 
Modu.Nbps = Nbps_input; % Number of Bits Per Symbol = length of data packet (changing with types of modulation)
Modu.fsymb = 5e6; % Symbolrate [Symb/s]
Modu.bps = Modu.fsymb*Modu.Nbps; % Bitrate [bits/s]
Modu.Tsymb = 1/Modu.fsymb; % Period of a symbol [s/Symb]

%% RCC filter
RCC.fcutoff = 1e6; % Cutoff frequency [1MHz]
RCC.beta = 0.3; % Roll-off factor [0.25]
RCC.taps = 151; % Number of points of the RCC
RCC.M = 2; % Upsampling Factor in order to satisfy the ISI Nyquist Criterion
RCC.fs = RCC.M*Modu.fsymb; % Sampling Frequency [Hz]

%% LDPC
LDPC.M = 128; % Number of rows
LDPC.N = 256; % Number of columns
LDPC.method = 0; % Method for distributing non-zero element: 0 for Evencol, 1 for Evenboth
LDPC.noCycle = 1; % Length-4 cycle
LDPC.onePerCol = 3; % Number of ones per column
LDPC.strategy = 0; % Strategy for finding the next non-zero diagonal elements
LDPC.iteration = 10; % Numbers of maximum iterations

%% generate a random binary bitstream
tx_wide = 48*2;
tx_bit = randint(LDPC.M,tx_wide);

%% encoding
H=makeLdpc(LDPC.M, LDPC.N, LDPC.method, LDPC.noCycle, LDPC.onePerCol); % Generation of H
[c, newH]=makeParityChk(tx_bit, H, LDPC.strategy); % Encoding
tx_bin_ini = [c;tx_bit]; % Reshape the coded bits to a bitstream
tx_bin = tx_bin_ini(:);
tx_len = length(tx_bin); % Length of bitstream

% encoding test
% tx_test = zeros(LDPC.N - LDPC.M, tx_wide);
% tx_bin_test = reshape(tx_bin, 256, 40);
%  for i=1:tx_wide % Start decoding row by row
%     tx_test(:,i) = tx_bin_test(LDPC.M+1:end,i);
%  end
%% modulation of bitstream
if mod(tx_len,Modu.Nbps) ~= 0 % ZERO PADDING -> Allows to avoid to get an error if the Bitstream length is not a multiple of Nbps
    tx_bin = [tx_bin; zeros(Modu.Nbps - mod(tx_len,Modu.Nbps),1)];
    tx_len = tx_len + Modu.Nbps - mod(tx_len,Modu.Nbps);
end
tx_symb = mapping(tx_bin, Modu.Nbps, Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps

%% upsampling
tx_symb_UP = upsample(tx_symb,RCC.M); % Upsample

%% filtering
RCC.h = RRCFDesign(RCC.beta, RCC.taps, RCC.fs, Modu.Tsymb); % Impulse Response of filter
tx_symb_UP_F = conv(tx_symb_UP,RCC.h); % convolution of signal and filter

%% channel awgn
% In this simulation, we directly use the baseband signal instead of
% narrowband in reality, but the energy is not the same
Es_BB = (trapz(abs(tx_symb_UP_F).^2))*(1/(RCC.fs)); % Basedband Signal Energy
Es_BP = Es_BB/2; % Bandpass Signal Energy
Eb = Es_BP/tx_len; % Average energy per bit
EbN0_dB = EbN0_input; % Energy per Bit to Noise PSD ratio [dB]
EbN0_ratio = 10^(EbN0_dB/10); % Energy per Bit to Noise PSD ratio
N0 = Eb/EbN0_ratio; % Noise PSD
NoisePower = 2*N0*RCC.fs; % Noise Power
noise_AWGN = sqrt(NoisePower/2)*(randn(length(tx_symb_UP_F), 1) + 1i*randn(length(tx_symb_UP_F), 1)); % Noise generation
rx_symb_UP_F = tx_symb_UP_F + noise_AWGN; % Add the noise

%% filtering
rx_symb_UP = conv(rx_symb_UP_F,RCC.h); % Apply the filter
rx_symb_UP = rx_symb_UP(RCC.taps: length(rx_symb_UP)-RCC.taps+1); % Removes the irrelevant values

%% downsampling
rx_symb = downsample(rx_symb_UP,RCC.M);

%% demodulation
rx_bit = demapping(rx_symb,Modu.Nbps,Modu.mod); 

%% decoding
 bin_width = length(rx_bit)/LDPC.N; % Width of coded bits
 rx_bit_ldpc = reshape(rx_bit,LDPC.N,bin_width); % Input of the decoder must be of same length of LDPC.N
 rx_bin = zeros(LDPC.N,bin_width); % Initialization of decoded bits
 rx_decoded_bin = zeros(LDPC.N - LDPC.M, bin_width); 
 %rx_test = zeros(LDPC.N - LDPC.M, bin_width); % decoding test

 for i=1:bin_width % Start decoding row by row
    coded_bit = rx_bit_ldpc(:,i);
    [rx_bin(:,i)] = LDPC_HARDdecoder(coded_bit, newH, LDPC.iteration);
    rx_decoded_bin(:,i) = rx_bin(LDPC.M+1:end,i);
    %rx_test(:,i) = rx_bit_ldpc(LDPC.M+1:end,i); % decoding test
 end
fprintf('Message decoded.\n');
rx_decoded_bin = rx_decoded_bin(:); % received a bitstream

%% output_plot
%SNR = 10*log10(EbN0_input*Modu.bps/RCC.fs);
BER = 1 - sum(rx_decoded_bin == tx_bit(:))/(tx_len/2); % Bit Error Ratio
%BER_test = 1 - sum(tx_bit(:) == rx_test(:))/(tx_len/2); % Bit Error Ratio
%BER_test2 = 1 - sum(rx_bin(:) == rx_bit)/tx_len; % Bit Error Ratio
%SER = 1 - sum( bi2de(fliplr(reshape(rx_bin,Modu.Nbps,tx_len/Modu.Nbps)')) == bi2de(fliplr(reshape(tx_bin,Modu.Nbps,tx_len/Modu.Nbps)'))) / (tx_len/Modu.Nbps) % Symbol Error Ratio

%h=scatterplot(tx_symb)
%scatterplot(rx_symb,[],[],'rx',h);
%title('BPSK constellation with channel AWGN')
%title('64QAM constellation')
%xlabel('Real Axis')
%ylabel('Imaginary Axis')
end