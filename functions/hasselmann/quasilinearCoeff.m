function [p_s_ql,beta,xi_sqr,k_x] = quasilinearCoeff(k,k_y,k_x,waveSpectrum,SARmetadata,th)

func = helperFunctions;

beta = func.getBeta(SARmetadata);
look = func.getLook(SARmetadata);

look = func.look(look);
k_new = k;
k_l = func.kl(look,k_y);
omega = func.omega(k_new);

waveSpectrum = func.resize(waveSpectrum,th);

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k = func.resize(Tv_k(:,2:end-1),waveSpectrum);

for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

for i=2:length(k_x)
    dk_x(i) = k_x(i)-k_x(i-1);
end
dk_x = dk_x(2:end);
dk_x = abs(mean(dk_x));

for i=2:length(k_y)
    dk_y(i) = k_y(i)-k_y(i-1);
end
dk_y = dk_y(2:end);
dk_y = abs(mean(dk_y));

xi_sqr = beta.^2.*cumtrapz(cumtrapz(waveSpectrum.*abs(Tv_k).^2).*dk_x).*dk_y;
xi_sqr = func.resize(xi_sqr(:,2:end),waveSpectrum);
xi_sqr = beta.^2.*trapz(trapz(waveSpectrum.*abs(Tv_k).^2).*dk_x).*dk_y;
xi = sqrt(xi_sqr);
test_k_neg = -(k_x).^2;
k_x_cutoff = xi^(-1);
negative = 0;

% Cutoff k_x
if mean(k_x) < 0
    k_x_cutoff = -k_x_cutoff;
    negative = 1;
end    
% Initialize the index variable
index = 0;
% Iterate through the array to find the first index greater than the threshold
if negative
    k_x_cutoff_index = k_x > k_x_cutoff;
    k_x = k_x.*k_x_cutoff_index;
else    
    k_x_cutoff_index = k_x < k_x_cutoff;
    k_x = k_x.*k_x_cutoff_index;
end

p_s_ql = exp(-(k_x).^2*xi_sqr);
end