function E = generateSingleJONSWAP(Hs,w0,gammaVal,w)
% Using Thor Fossen's wavepsec function

E = wavespec(7,[Hs,w0,gammaVal],w,0);

end