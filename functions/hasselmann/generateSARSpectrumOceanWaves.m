function generatedSARSpectrum = generateSARSpectrumOceanWaves(k,k_y,k_x,E_k,metadata,incidenceAngle,r,nonLinOrder,w)

func = helperFunctions;
look = func.getLook(metadata);
look = func.look(look);
polarisation = func.getPolarisation(metadata);
beta = func.getBeta(metadata);

% Co and autocovariance functions
[f_v_r] = orbitalVelocityCovariance(k,k_y,E_k,metadata,r,incidenceAngle);
[f_r_r] = rarImageIntensityAutocovariance(k,k_y,E_k,metadata,r,incidenceAngle);
[f_rv_r] = rarImageIntensityCovariance(k,k_y,E_k,metadata,r,incidenceAngle);
[p_s_coeff,beta,xi_sqr] = quasilinearCoeff(k,k_y,k_x,E_k,metadata,incidenceAngle);

% Spectral Expansions
[p_s_2n] = spectralExpansion2n(f_v_r,nonLinOrder);
f_rv_r_inv = rarImageIntensityCovariance(k,k_y,E_k,metadata,-1.*r,incidenceAngle);
[p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_inv,f_v_r,nonLinOrder);
f_rv_r_zero = rarImageIntensityCovariance(k,k_y,E_k,metadata,0.*r,incidenceAngle);
[p_s_2n_2] = spectralExpansion2n_2(f_r_r,f_v_r,f_rv_r,f_rv_r_zero,f_rv_r_inv,1);

% Full spectral expansion
P_s = zeros(512,512);

for nonLinOrder = min(nonLinOrder):max(nonLinOrder)
    for m = 2*nonLinOrder-2:2*nonLinOrder
        fprintf('m = %d, ',m)
        fprintf('n = %d \n',nonLinOrder)
        coeff = ((k_x.*beta).^m);
        if (m == 2*nonLinOrder)
            fprintf('P_(n,2n) used \n')
            P_s = P_s + coeff.*p_s_2n;
        end
        if (m == 2*nonLinOrder-1)
            fprintf('P_(n,2n-1) used \n')
            P_s = P_s + coeff.*p_s_2n_1;
        end        
        if (m == 2*nonLinOrder-2)
            fprintf('P_(n,2n-2) used\n')
            P_s = P_s + coeff.*p_s_2n_2;
        end
        
    end
end
P_s = p_s_coeff.*P_s;
generatedSARSpectrum = P_s;
end