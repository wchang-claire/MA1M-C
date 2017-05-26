%% input parameters
clear; clc; close all;
mod_input = 'pam'; % ['pam','qam','qpsk','bpsk']
Nbps_input = 1; %[1,4,6,2,1]
EbN0_input = 10;

%% modulation parameters
Modu.mod = mod_input; 
Modu.Nbps = Nbps_input; % Number of Bits Per Symbol = length of data packet (changing with types of modulation)
Modu.fsymb = 5e6; % Symbolrate [Symb/s]
Modu.bps = Modu.fsymb*Modu.Nbps; % Bitrate [bits/s]
Modu.Tsymb = 1/Modu.fsymb; % Period of a symbol [s/Symb]

%% RCC filter
RCC.fcutoff = 1e6; % Cutoff frequency [1MHz]
RCC.beta = 0.3; % Roll-off factor [0.25]
RCC.taps = 33; % Number of points of the RCC
RCC.M = 2; % Upsampling Factor in order to satisfy the ISI Nyquist Criterion
RCC.fs = RCC.M*Modu.fsymb; % Sampling Frequency [Hz]

%% generate a random binary bitstream
tx_len = 1024; % Length of the bitstream
tx_bin = randi([0 1], tx_len,1); % Bitstream - Generate binary sequence (0 and 1 are equiprobable)

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
rx_bin = demapping(rx_symb,Modu.Nbps,Modu.mod); 

%% output_plot
%SNR = 10*log10(EbN0_input*Modu.bps/RCC.fs);
BER = 1 - sum(rx_bin == tx_bin)/tx_len; % Bit Error Ratio
%SER = 1 - sum( bi2de(fliplr(reshape(rx_bin,Modu.Nbps,tx_len/Modu.Nbps)')) == bi2de(fliplr(reshape(tx_bin,Modu.Nbps,tx_len/Modu.Nbps)'))) / (tx_len/Modu.Nbps) % Symbol Error Ratio

%h=scatterplot(tx_symb)
%scatterplot(rx_symb,[],[],'rx',h);
%title('BPSK constellation with channel AWGN')
%title('64QAM constellation')
%xlabel('Real Axis')
%ylabel('Imaginary Axis')
end
