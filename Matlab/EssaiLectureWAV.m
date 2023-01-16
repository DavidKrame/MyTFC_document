[y,fs]=audioread('LETTRE.wav');
samples=[1,5*fs];
[yc,fs]=audioread('LETTRE.wav',samples);
sound(yc,fs);
stem (yc);