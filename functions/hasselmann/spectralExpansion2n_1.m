function [p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_negative,f_v_r,n)
p_s_2n_1 = (1i.*(f_rv_r-f_rv_r_negative).*f_v_r.^(n-1))./factorial(n-1); % Update when I know what -r means

p_s_2n_1 = fftshift(fft2(p_s_2n_1));
end