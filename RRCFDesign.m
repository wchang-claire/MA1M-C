function [RRC_h, RRC_sp, t_axis, f_axis] = RRCFDesign(Beta, taps, fs, Tsymb)
%% impulse response (using ifft)
fmax = (1/taps)*fs*(taps-1)/2; % Maximal frequency on the axis 
f_axis = linspace(-fmax, fmax, taps); % Frequency axis
t_axis = (-(taps-1)/2:(taps-1)/2)./(2*fmax); % Time axis
RRC_FT = zeros(1,taps); % Initialization of Fourier Transform array

% Fourier Transform of the RRC
for i=1:taps
    if abs(f_axis(i))<=((1-Beta)/(2*Tsymb))
        RRC_FT(i) = Tsymb;
    elseif abs(f_axis(i))>((1-Beta)/(2*Tsymb)) && abs(f_axis(i))<=((1+Beta)/(2*Tsymb))
        RRC_FT(i) = Tsymb/2*(1+cos(pi*Tsymb/Beta*(abs(f_axis(i))-(1-Beta)/(2*Tsymb))));
    else
        RRC_FT(i) = 0;
    end
end

RRC_sp = sqrt(RRC_FT); % Spectrum of the filter
RRC_h = fftshift(ifft(ifftshift(RRC_sp), 'symmetric')); % Impulse response
normCoeff = sqrt(max(abs(conv(RRC_h,RRC_h)))); % Normalize coefficient
RRC_h=RRC_h/normCoeff; % Normalize -> the impulse response will then be equal to 1 @ t=0 ( not the spectrum )
end