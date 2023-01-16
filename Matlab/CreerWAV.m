% filename='entree.wav';
fs=44100;
filename='D:\entree.wav';
% audiowrite(filename,myrecordings,44100)
[myrecordings,fs]=audioread(filename);