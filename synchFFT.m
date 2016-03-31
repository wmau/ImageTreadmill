function [f,P1] = synchFFT(FT,inds)
%
%
%

%%
    Fs = 20; 
    L = length(inds); 
    T = 1/Fs; 
    t = (0:L-1)*T;
    
    synch = sum(FT(:,inds));
    Y = fft(synch);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
     
    f = Fs*(0:L/2)/L; 
end