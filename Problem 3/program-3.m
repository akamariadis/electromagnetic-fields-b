function program
    clear; 
    close all; 
    clf; 
    clc; 

    ratio_ab = linspace(0, 0.95, 200); 
    mu_opt_r = (1 + ratio_ab.^2) ./ (1 - ratio_ab.^2);

    figure(1); 
    set(gcf, 'Color', 'w');
    plot(ratio_ab, mu_opt_r, 'LineWidth', 2.5, 'Color', 'b');

    grid on;
    xlabel('Λόγος Ακτίνων (a/b)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Βέλτιστη Σχετική Διαπερατότητα', 'FontSize', 12, 'FontWeight', 'bold');
    title('Βέλτιστη Διαπερατότητα Φλοιού', 'FontSize', 14);

    xlim([0 1]);
    ylim([1 10]);

    a = 1;
    mu0 = 4*pi*1e-7; 
    H0 = 1;

    b_partC = 2 * a;
    ratio = a / b_partC;
     
    mur2_opt = (1 + ratio^2) / (1 - ratio^2);
     
    cases_mu = [1e-5, 1, mur2_opt, 2*mur2_opt];
    titles_C = {
        '\mu_{r2} \approx 0', ...
        '\mu_{r2} = 1', ...
        ['Optimal \mu_{r2} = ' num2str(mur2_opt, '%.2f')], ...
        ['2 \times Optimal \mu_{r2} = ' num2str(2*mur2_opt, '%.2f')]
    };

    offset_C = 1; 

    for i = 1:length(cases_mu)
        figure(i + offset_C);
        plot_field(a, b_partC, cases_mu(i), H0, titles_C{i});
    end

    ratios_D = [1.1, 2.1, 3.1, 4.1];
    start_fig = length(cases_mu) + 1 + offset_C; 
     
    for i = 1:length(ratios_D)
        ba_ratio = ratios_D(i);
        b_curr = ba_ratio * a;
        inv_ratio = a / b_curr;
        mur2_curr = (1 + inv_ratio^2) / (1 - inv_ratio^2);
        fig_title = sprintf('b/a=%.1f, Optimal \\mu_{r2}=%.2f', ba_ratio, mur2_curr);
        figure(start_fig + i);
        plot_field(a, b_curr, mur2_curr, H0, fig_title);
    end
end

function plot_field(a, b, mur2, H0, plot_title)
    mu0 = 4*pi*1e-7;
    mu2 = mur2 * mu0;
     
    A = [1, 1/a^2, 0; -b, 1/b, -1/b; mu2, mu2/b^2, -mu0/b^2];
    RHS = [0; -H0*b; mu0*H0];
     
    X = A \ RHS;
    C2 = X(1); D2 = X(2); D3 = X(3);
     
    L = 2.5 * b;
    N = 300;
    x = linspace(-L, L, N);
    y = linspace(-L, L, N);
    [XX, YY] = meshgrid(x, y);
    [PHI, RR] = cart2pol(XX, YY);
     
    Bx = zeros(size(XX));
    By = zeros(size(YY));
     
    mask2 = (RR >= a) & (RR < b);
    Br_2 = mu2 * (C2 + D2 ./ RR(mask2).^2) .* cos(PHI(mask2));
    Bphi_2 = -mu2 * (C2 - D2 ./ RR(mask2).^2) .* sin(PHI(mask2));
     
    Bx(mask2) = Br_2 .* cos(PHI(mask2)) - Bphi_2 .* sin(PHI(mask2));
    By(mask2) = Br_2 .* sin(PHI(mask2)) + Bphi_2 .* cos(PHI(mask2));
     
    mask3 = RR >= b;
    Br_3 = mu0 * (H0 + D3 ./ RR(mask3).^2) .* cos(PHI(mask3));
    Bphi_3 = -mu0 * (H0 - D3 ./ RR(mask3).^2) .* sin(PHI(mask3));
     
    Bx(mask3) = Br_3 .* cos(PHI(mask3)) - Bphi_3 .* sin(PHI(mask3));
    By(mask3) = Br_3 .* sin(PHI(mask3)) + Bphi_3 .* cos(PHI(mask3));
     
    B_mag = sqrt(Bx.^2 + By.^2);
    
    pcolor(XX, YY, B_mag); 
    shading interp; 
    colormap(jet); 
    caxis([0 2.5*mu0*H0]);
    colorbar; 
    hold on;
     
    streamslice(XX, YY, Bx, By, 2, 'method', 'cubic');
     
    viscircles([0 0], a, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
    viscircles([0 0], b, 'Color', 'k', 'LineWidth', 2);
     
    axis equal;
    xlim([-L L]); ylim([-L L]);
    title(plot_title, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x'); ylabel('y');
    box on;
end
