clc;
clear all;
close all;
snr_arr = 0:30;
ber_arr = zeros(length(snr_arr),1);
n_packet=1000;


RayleighChannelFlag = false;

%specify path delays and gains using ITU pedestrian B channel parameters
fs = 1e6;   
%fs = 3.84e6;                                     % Hz
pathDelays = [0 200 800 1200 2300 3700]*1e-9;    % sec
avgPathGains = [0 -0.9 -4.9 -8 -7.8 -23.9];      % dB
fD = 0;                                         % Hz

channel = 0;
if RayleighChannelFlag
	%https://www.mathworks.com/help/comm/ref/comm.rayleighchannel-system-object.html
    channel = comm.RayleighChannel('SampleRate',fs, 'PathDelays',pathDelays, 'AveragePathGains',avgPathGains, 'MaximumDopplerShift',fD);
else
	%https://www.mathworks.com/help/comm/ref/comm.ricianchannel-system-object.html
    channel = comm.RicianChannel('SampleRate',fs, 'PathDelays',pathDelays, 'AveragePathGains',avgPathGains,'KFactor',10, 'MaximumDopplerShift',fD);
end
% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 4096;          % APEP length in bytes
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
Rs = wlanSampleRate(cfgVHT);       % Sampling rate



lstf = wlanLSTF(cfgVHT);
lltf = wlanLLTF(cfgVHT);
lsig = wlanLSIG(cfgVHT);
nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields
vhtsiga = wlanVHTSIGA(cfgVHT);
vhtstf = wlanVHTSTF(cfgVHT);
vhtltf = wlanVHTLTF(cfgVHT);
vhtsigb = wlanVHTSIGB(cfgVHT);
preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];


for snr = snr_arr
    numErr = 0;
    disp(['SNR = ' num2str(snr)]);
    parfor symb = 1:n_packet %peredaem packet_num paketov

        txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
        data = wlanVHTData(txPSDU,cfgVHT);

        % A VHT waveform is constructed by prepending the non-HT and VHT
        % preamble fields with data
        txWaveform = [preamble;data]; % Transmit VHT PPDU

        rxWaveform = awgn(txWaveform,snr,0);
        rxWaveform = channel(rxWaveform);

        indField = wlanFieldIndices(cfgVHT);
        indLLTF = indField.LLTF(1):indField.LLTF(2);
        demodLLTF = wlanLLTFDemodulate(rxWaveform(indLLTF),cfgVHT);
        % Estimate noise power in VHT fields
        nVar = helperNoiseEstimate(demodLLTF,cfgVHT.ChannelBandwidth,cfgVHT.NumSpaceTimeStreams);

        indVHTLTF = indField.VHTLTF(1):indField.VHTLTF(2);
        demodVHTLTF = wlanVHTLTFDemodulate(rxWaveform(indVHTLTF,:),cfgVHT);
        chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);

        indData = indField.VHTData(1):indField.VHTData(2);

        % Recover the bits and equalized symbols in the VHT Data field using the
        % channel estimates from VHT-LTF
        [rxPSDU,~,eqSym] = wlanVHTDataRecover(rxWaveform(indData,:),chanEstVHTLTF,nVar,cfgVHT);

        % Compare transmit and receive PSDU bits
        [number,ratio] = biterr(txPSDU,rxPSDU);
        numErr = numErr + ratio;
    end
    ber_arr(snr+1) = numErr/n_packet; % veroyatnost ochibki dlya kaghdogo packeta
end

txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
data = wlanVHTData(txPSDU,cfgVHT);
txWaveform = [preamble;data]; % Transmit VHT PPDU

txWaveformRe = real(txWaveform);
txWaveformIm = imag(txWaveform);

ber_arr_awgn = ber_arr;

save('wifi802_11_ac_awgn.mat', 'snr_arr', 'ber_arr_awgn')

semilogy(snr_arr, ber_arr, '-'); 
grid on;
xlabel('osh,dB'); 
ylabel('veroyatnost oshibki'); 
legend('WiFi 802.11ac');
