function nValidation(k,endDepth,step)
% Plots different values of n at different depths (starting at 0m, and
% advancing to the endDepth at the defined step)

% Calculate n at different depths
i = 0.001;
n = 0.5*(1+(2.*k.*i)./sinh(2.*k.*i)); % sinh in radians (eq. 5.4.32 in holthuisjen)
n_vals(1) = n(60);
d_vals(1) = 0;
i = step;
j = 2;
while i <= endDepth
    n = 0.5*(1+(2.*k.*i)./sinh(2.*k.*i)); % sinh in radians (eq. 5.4.32 in holthuisjen)
    n_vals(j) = n(60);
    d_vals(j) = i;
    i = i + step;
    j = j+1;
end
% Plot range of d over n
figure;
plot(d_vals,n_vals);
grid on;
xlabel('Depth (m)'), ylabel('Dispersion relation, n'), title('Comparison of dispersion relation, n, over varying depths');
xticks(d_vals);
matlab2tikz('../plots/nValidation.tex');

end