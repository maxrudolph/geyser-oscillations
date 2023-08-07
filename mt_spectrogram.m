function [all_Pxx,frequency,times] = mt_spectrogram(pp,window_length,overlap,Fs)
%% use Rob's multitaper code to make a spectrogram

% [Pxx,Exx,pX,f]=mt_cspek_phs(P(2,t1:t2),Fs,1);
% figure()
% plot(f,Pxx);
% window_length = 100*Fs;
% overlap = window_length*1/2;
% pp = P(2,t1:t2);

N1 = length(pp);

ind=1;
start_indices = 1:(window_length-overlap):(N1-window_length-1);
start = start_indices(ind);
signal = pp(start:(start+window_length-1));
%window = hamming(window_length)';
window = ones(size(signal));

[Pxx,Exx,pX,frequency]=mt_phs(signal.*window,Fs,1);
all_Pxx = zeros(length(Pxx),length(start_indices));
all_Pxx(:,ind) = Pxx;

parfor ind = 2:length(start_indices)
    start = start_indices(ind);
    signal = pp(start:(start+window_length-1));
    [Pxx,Exx,pX,frequency]=mt_phs(window.*signal,Fs,1);
    if ind==1
       all_Pxx(:,ind) = Pxx;
    else
       all_Pxx(:,ind) = Pxx; 
    end
end

times = 1/Fs*(start_indices-1+window_length/2);