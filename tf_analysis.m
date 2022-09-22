[tones, Fs] = helperDTMFToneGenerator();
p = audioplayer(tones, Fs, 16)
play(p)

N = numel(tones);
t = (0:N-1)/Fs;
subplot(2,1,1)
plot(1e3*t,tones)
xlabel('Time(ms)')
ylabel('Amplitude')
title('DTMF Signal')
subplot(2,1,2)
pspectrum(tones,Fs,'Leakage',1,'FrequencyLimits',[650,1500])

env = envelope(tones,80,'rms');
pulsewidth(env,Fs)
title('PUlse Width of RMS Envelope')

f = [meanfreq(tones,Fs,[700 800]), ...
     meanfreq(tones,Fs,[800 900]), ...
     meanfreq(tones,Fs,[900 1000]), ...
     meanfreq(tones,Fs,[1300 1400])];
round(f)

pspectrum(tones,Fs,'spectrogram','Leakage',1,'OverlapPercent',0, ...
    'MinThreshold',-10,'FrequencyLimits',[650, 1500]);