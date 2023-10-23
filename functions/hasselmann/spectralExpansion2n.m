function [p_s_2n] = spectralExpansion2n(f_v_r,n)

p_s_2n = (f_v_r.^n)./factorial(n);

p_s_2n = fftshift(fft2(p_s_2n));

end