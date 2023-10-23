function [p_s_2n_2] = spectralExpansion2n_2(f_r_r,f_v_r,f_rv_r,f_rv_r_zero,f_rv_r_negative,n)
if(n-2 < 0)
    p_s_2n_2_val = (1./factorial(n-1)).*f_r_r.*f_v_r.^(n-1) + (0).*(f_rv_r - f_rv_r_zero).*(f_rv_r_negative-f_rv_r_zero); %.* f'(r)^(n-2)
    p_s_2n_2 = fft2(p_s_2n_2_val);
    return;
end
p_s_2n_2_val = (1./factorial(n-1)).*f_r_r.*f_v_r.^(n-1) + (1./factorial(n-2)).*(f_rv_r - f_rv_r_zero).*(f_rv_r_negative-f_rv_r_zero); %.* f'(r)^(n-2)

p_s_2n_2 = fftshift(fft2(p_s_2n_2_val));
end