function twoDWaveSpectrumPlot(E,wOrK,contourOrSurf,w,theta,k_x,k_y)
% Freq-direction = 0
if wOrK
    E_k = E;
    if contourOrSurf
        % Contour of 2D wave spectrum
        figure;
        contour(k_x,k_y,E_k)
        grid on;
        grid minor;
        yline(0);
        xline(0);
        xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
        title('Contour plot of E(k_x,k_y)')
        %set(gca,'XTick',-pi/2:pi/6:pi/2) 
        %set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
        %set(gca,'YTick',-pi/2:pi/6:pi/2) 
        %set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
        set(gca,'FontSize',12)
        %matlab2tikz('../plots/Valcontour_E_k_inv.tex');
    else
        %% 3D Spectrum
        figure;
        h=surf(k_x,k_y,E_k);
        grid on;
        xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex'),zlabel('$E(k_{x},k_{y})$','interpreter','latex');
        set(h,'LineStyle','none');
        c = colorbar;
        %c.Label.String = barStr;
        hL = ylabel(c,'[m^2/rad/Hz]');
        %set(hL,'Rotation',0);
        %set(gca,'FontSize',32)
        title('Surface plot of E(k_x,k_y)')
        %set(gca,'XTick',-pi/2:pi/6:pi/2) 
        %set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
        %set(gca,'YTick',-pi/2:pi/6:pi/2) 
        %set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
        %matlab2tikz('../plots/surf_E_k.tex');
    end   
else
    % Define for plotting
    [t,r] = meshgrid(theta,w);
    [x,y] = pol2cart(t,r);
    if contourOrSurf
        % Contour of 2D wave spectrum
        figure;
        contour(x,y,E)
        yline(0);
        xline(0);
        grid on;
        xlabel('\omega [rad/s]'), ylabel('\omega [rad/s]');
        title('Contour plot of E(\omega,\theta)')
        set(gca,'XTick',-2*pi:pi/2:2*pi) 
        set(gca,'XTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
        set(gca,'YTick',-2*pi:pi/2:2*pi) 
        set(gca,'YTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
        set(gca,'FontSize',12)
        %matlab2tikz('../plots/contour_E_omega.tex');
    else
        %% 3D Spectrum
        figure;
        h=surf(x,y,E);
        grid on;
        xlabel('\omega [rad/s]'), ylabel('\omega [rad/s]');
        set(gca,'XTick',-2*pi:pi/2:2*pi) 
        set(gca,'XTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
        set(gca,'YTick',-2*pi:pi/2:2*pi) 
        set(gca,'YTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
        %set(gca,'FontSize',22)
        set(h,'LineStyle','none');
        title('Surface plot of E(\omega,\theta)')
        c = colorbar;
        %c.Label.String = barStr;
        hL = ylabel(c,'[m^2/rad/Hz]');
        %set(hL,'Rotation',0);
        % figure(6);
        % h = surf(w,theta,E_cleaned);
        % set(h,'LineStyle','none');
        % xlabel('\omega (rad/s)')
        % ylabel('\theta (rad)')
        zlabel('E(\omega,\theta)')
        %matlab2tikz('../plots/surf_E_omega.tex');
    end 


end    

end