clear;
close all;
clc;

gmin = 0.1; gmax = 40; Ng = 400;
gammas = linspace(gmin, gmax, Ng);

gmax_vals = nan(size(gammas));

xs = 10.^linspace(-6, 2, 4000);

for i = 1:numel(gammas)
    g = gammas(i);
    g_of_x  = @(x) x .* (1 + g ./ (1 + x.^2));
    dg_dx   = @(x) 1 + g ./ (1 + x.^2) - (2*g .* x.^2) ./ (1 + x.^2).^2;
    dvals = dg_dx(xs);
    idx = find( dvals(1:end-1) .* dvals(2:end) < 0 );
    rootsCrit = [];
    for k = 1:numel(idx)
        a = xs(idx(k)); b = xs(idx(k)+1);
        try
            rt = fzero(@(t) dg_dx(t), [a b]);   % fzero είναι core MATLAB
            if rt > 0
                if isempty(rootsCrit) || all(abs(rootsCrit-rt) > 1e-8)
                    rootsCrit(end+1) = rt; %#ok<SAGROW>
                end
            end
        catch
        end
    end
    rootsCrit = sort(rootsCrit);

    if ~isempty(rootsCrit)
        gvals = g_of_x(rootsCrit);
        gmax_vals(i) = max(gvals);
    else
        gmax_vals(i) = NaN;  % μονοτονική περίπτωση (π.χ. γ<=8)
    end
end

mask = (gammas > 8) & ~isnan(gmax_vals);
G  = gammas(mask);
Yb = (gmax_vals(mask)).^2;   % όριο: κE0^2 = g_max(γ)^2

figure('Color','w'); hold on; box on; grid on;

if ~isempty(G)
    fill([G, fliplr(G)], [zeros(size(Yb)), fliplr(Yb)], ...
         [0.95 0.85 0.55], 'EdgeColor','none', 'FaceAlpha',0.35);
    plot(G, Yb, 'LineWidth', 2);
end

if ~isempty(Yb)
    ymax = max(Yb)*1.05;
else
    ymax = 1;
end
plot([8 8], [0 ymax], '--', 'LineWidth', 1.2);

xlabel('\gamma','Interpreter','tex');
ylabel('\kappa E_0^2','Interpreter','tex');
title('Μέρος (Α): Περιοχή διευστάθειας στο επίπεδο (γ, κE0^2)', ...
      'Interpreter','tex');
xlim([gmin gmax]);
ylim([0 ymax]);
