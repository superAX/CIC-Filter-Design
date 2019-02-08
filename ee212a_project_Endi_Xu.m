% EE212A Digital Filter Term Project
% Endi Xu 005030030
% Decimation Filter 

% Set parameters
OSR = 128;                 % Oversample ratio is 128
M = 1;                     % Differential Delay, choose M = 1
N = 8;                     % Number of Stages, choose N = 8
Fin = 1024e3;              % Input sample rate frequency
Fout = Fin/OSR;            % Output sample rate frequency   
Fpass = 2e3;               % Pass Frequency
Fstop = 3e3;               % Stop Frequency  
Bq = 10;                   % Bits used for quantization  

% Transfer function of the CIC filter  
syms z;
hcic = 0;
for i = 0 : OSR*M-1
    hcic = hcic + z^-i;
end
hcic = (hcic/(OSR*M))^N;
[num,den] = numden(hcic);    
hnCIC = double(sym2poly(num));             % Get the coefficients of numerator
hdCIC = double(sym2poly(den));             % Get the coefficients of denominator  
hnCIC = hnCIC/hdCIC(1);                    % Transfer function of CIC

% CIC Compensation FIR Filter (Single rate at 8kHZ)
L = 30;                                                                % Order of filter taps 
fp = 0:1e-4:(Fpass+0.1e3)/Fout;                                        % Passband frequency
fs = Fstop/Fout:1e-4:0.5;                                              % Stopband frequency
f = [fp fs] * 2;                                                       % Normalized frequency
f(end) = 1;
Mp = ones(1, length(fp));                                              % Passband response
Mp(2:end) = abs(M*OSR*sin(pi*fp(2:end)/OSR)./sin(pi*M*fp(2:end))).^N;  % Inverse Sinc
Mf = [Mp zeros(1, length(fs))];
hciccomp = fir2(L,f,Mf);   
hciccomp = hciccomp/sum(hciccomp);                                     % transfer function of CIC compensator filter 
hciccompQ = quant(hciccomp,1/(2^(Bq-1)));                              % transfer function of CIC COMP with quantization
kCICCOMPQ = sum(hciccompQ);
hciccompQ = hciccomp/kCICCOMPQ;
    
fvtool(hnCIC, 1, hciccomp, 1);
title('CIC filter and CIC compensatentor filter');
legend('CIC filter', 'CIC compensatentor filter');

% Combination
hciccompup = upsample(hciccomp, OSR);                                  % To make the two filters run at same rate which is 1024KHZ      
hciccompupQ = upsample(hciccompQ, OSR); 
hcombine =  conv(hnCIC, hciccompup); 
hcombineQ = conv(hnCIC, hciccompupQ);
fvtool(hcombine, 1, hcombineQ, 1);                                     % Generate plots for magnitude response
title('Combine Filter Magnitude response');
legend("CIC filter", "CIC fitler with Quantization");
figure(1);
subplot(121);
zplane(hcombine, 1);                                                   % Poles and Zeros without quantization
title('Poles and Zeros with quantization');
subplot(122);
zplane(hcombineQ, 1);                                                  % Poles and Zeros with quantization
title('Poles and Zeros without quantization');

% Noise Calculation
Bmaxcic = ceil(N*log2(OSR*M));                                         % Maximum bit for CIC then no error
strl = sprintf('Bmax for CIC filter is %d', Bmaxcic);
disp(strl);
% Quantization 1 to 16 bits for CIC compensator filter
for Bq = 10 : 1 : 17
    noiseCICCOMPVar = sqrt(L*(1/(2^Bq))^2/12);                                    % Noise one each node
    interNoiseVar = sqrt((1/(2^Bmaxcic) - 1/(2^Bq))^2/12);                        % Noise caused by 56-bit to 8-bit   

    t=0:1/Fin:5;                                                            % sample rate is 1024khz 
    inputSignal = sin(2*pi*1e3*t);                                          % Generate input signal
                                       
    outputCIC = filter(hnCIC, 1, inputSignal);                              % CIC filter
    outputDownsample = downsample(outputCIC, OSR);                          % Downsample OSR = 128
    output = filter(hciccompQ, 1, outputDownsample);                        % CIC compensator filter
    interNoise = normrnd(0, interNoiseVar, size(output));                   
    noiseOutput = filter(hciccompQ, 1, interNoise) + normrnd(0, noiseCICCOMPVar, size(output));
    SNR = snr(output, noiseOutput);
    str1 = sprintf('When Bq = %d, the SNR is %f db', Bq, SNR);
    disp(str1);
end    

% Power Calculation
syms p;
NumAdderCIC = N;
NumRegisterCIC = N;
NumMultiplierCOMP = L+1;
NumAdderCOMP = L;
NumRegisterCOMP = L;
power = NumAdderCIC*4*Bmaxcic*(Fin/Fout)*p + NumRegisterCIC*5*Bmaxcic*(Fin/Fout)*p ...
        + NumAdderCIC*4*Bmaxcic*p+NumRegisterCIC*5*Bmaxcic*p ...
        +NumAdderCOMP*4*Bq*p+NumMultiplierCOMP*2*Bq^2*p+NumRegisterCOMP*5*Bq*p;
fprintf('The power consumed is %s\n', power);

% Testing
tin=0:1/Fin:0.5;                                                             % input sample rate is 1024khz 
tout = 0:1/Fout:0.5;                                                         % output sample rate is 8khz 
inputSignal = 0;
for inputF = 5e2 : 5e2 : 1e4
    inputSignal = inputSignal + sin(2*pi*inputF*tin);
end

outputCIC = filter(hnCIC, 1, inputSignal);                               % CIC filter
outputDownsample = downsample(outputCIC, OSR);                           % Downsample OSR = 128
output = filter(hciccompQ, 1, outputDownsample);                         % CIC compensator filter

figure(2);
subplot(211);
plot(tin, inputSignal);
title("Input signal with addition of sinwave from 500HZ to 10000HZ");
xlabel("time(s)");
ylabel("amplitude");
subplot(212);
plot(tout, output);
title("Output Signal");
xlabel("time(s)");
ylabel("amplitude");