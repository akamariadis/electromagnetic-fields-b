% ΕΡΩΤΗΜΑ Α

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

% ΕΡΩΤΗΜΑ Β

gamma = 20;
kappa = 20;
a = 1;

% --- Ορισμοί συναρτήσεων ---
fE   = @(E) E .* (1 + gamma ./ (1 + kappa .* E.^2));
dfdE = @(E) 1 + gamma ./ (1 + kappa .* E.^2) ...
            - (2*gamma*kappa .* E.^2) ./ (1 + kappa .* E.^2).^2;
Escan = logspace(-6, 1, 4000);
dval  = dfdE(Escan);

idx = find(dval(1:end-1) .* dval(2:end) < 0);
rootsCrit = [];
for k = 1:numel(idx)
    a_br = Escan(idx(k)); b_br = Escan(idx(k)+1);
    try
        rt = fzero(@(x) dfdE(x), [a_br b_br]);  % fzero = core MATLAB
        if rt > 0
            if isempty(rootsCrit) || all(abs(rootsCrit - rt) > 1e-8)
                rootsCrit(end+1) = rt; %#ok<SAGROW>
            end
        end
    catch
    end
end
rootsCrit = sort(rootsCrit);

if numel(rootsCrit) < 2
    error('Δεν βρέθηκαν δύο κρίσιμα σημεία. Ελέγξτε τις παραμέτρους/σάρωση.');
end

fvals = fE(rootsCrit);
f_local_min = min(fvals);
f_local_max = max(fvals);

ra = linspace(0.01, 0.999, 600);
E0_min = f_local_min .* (ra.^2);
E0_max = f_local_max .* (ra.^2);

E0sq_min = E0_min.^2;
E0sq_max = E0_max.^2;

figure('Color','w'); hold on; box on; grid on;

fill([ra, fliplr(ra)], [E0sq_min, fliplr(E0sq_max)], ...
     [0.95 0.85 0.55], 'EdgeColor','none', 'FaceAlpha',0.35);

plot(ra, E0sq_min, '--k', 'LineWidth', 1.3);
plot(ra, E0sq_max, '--k', 'LineWidth', 1.3);

xlabel('r / a', 'Interpreter','tex');
ylabel('E0^2 (V^2/m^2)', 'Interpreter','tex');
title('Μέρος (Β): Περιοχή διευστάθειας σε (E0^2,\ r/a) για γ=20 και κ=20', ...
      'Interpreter','tex');
xlim([0 1]);
ylim([0 max(E0sq_max)*1.05]);

% ΕΡΩΤΗΜΑ Γ

clear;
close all;
clc;

gamma = 20;
kappa = 20;
a     = 1;
E0min = 0.0; 
E0max = 3.0; 
NE0   = 800;
Emax  = 12.0;
Ngrid = 6000;

r_over_a_list = [0.5, 0.6, 0.75, 0.95];
fE = @(E) E .* (1 + gamma ./ (1 + kappa .* E.^2));
E0_vals = linspace(E0min, E0max, NE0);
for idxPos = 1:numel(r_over_a_list)
    ra = r_over_a_list(idxPos);
    S_vals = E0_vals .* (a/ra)^2;    % S = E0*(a/r)^2
    allRoots = cell(size(S_vals));
    for i = 1:numel(S_vals)
        S = S_vals(i);
        Es   = linspace(1e-12, Emax, Ngrid);
        vals = fE(Es) - S;
        iz      = find(abs(vals) < 1e-14);
        signchg = find(vals(1:end-1) .* vals(2:end) < 0);
        rootsList = [];
        if ~isempty(iz)
            rootsList = [rootsList, Es(iz)]; %#ok<AGROW>
        end
        for k = 1:numel(signchg)
            iL = signchg(k);
            a_br = Es(iL);
            b_br = Es(iL+1);
            try
                rt = fzero(@(x) fE(x) - S, [a_br b_br]);
                if rt > 0
                    rootsList(end+1) = rt; %#ok<AGROW>
                end
            catch
            end
        end
        if isempty(rootsList)
            allRoots{i} = [];
        else
            rootsList = sort(rootsList(:).');
            tol = 1e-8;
            keep = true(size(rootsList));
            for j = 2:numel(rootsList)
                if abs(rootsList(j) - rootsList(j-1)) < tol
                    keep(j) = false;
                end
            end
            allRoots{i} = rootsList(keep);
        end
    end
    inc = nan(size(S_vals));  prev = NaN;
    for i = 1:numel(S_vals)
        rts = allRoots{i};
        if isempty(rts)
            inc(i) = NaN; prev = NaN;
        else
            if isnan(prev)
                choice = rts(1);
            else
                [~,ix] = min(abs(rts - prev));
                choice = rts(ix);
            end
            inc(i) = choice; prev = choice;
        end
    end
    dec = nan(size(S_vals));  prev = NaN;
    for i = numel(S_vals):-1:1
        rts = allRoots{i};
        if isempty(rts)
            dec(i) = NaN; prev = NaN;
        else
            if isnan(prev)
                choice = rts(end);
            else
                [~,ix] = min(abs(rts - prev));
                choice = rts(ix);
            end
            dec(i) = choice; prev = choice;
        end
    end
    figure('Color','w'); hold on; box on; grid on;
    for i = 1:numel(E0_vals)
        rts = allRoots{i};
        if ~isempty(rts)
            plot(E0_vals(i)+0*rts, rts, '.k', 'MarkerSize', 7);
        end
    end
    plot(E0_vals, inc, '-',  'LineWidth', 2.0);
    plot(E0_vals, dec, '--', 'LineWidth', 2.0);
    xlabel('E_0 (V/m)');
    ylabel('E(r) (V/m)');
    title(sprintf('Μέρος (\\Gamma): E(r) vs E_0 στο r=%.2f a', ra), 'Interpreter','tex');
    xlim([E0min E0max]);
    ymax = max([3, nanmax([inc(:); dec(:)])*1.1]);
    if ~isfinite(ymax), ymax = 3; end
    ylim([0 ymax]);
    legend({'όλες οι ρίζες','ανιούσα','κατιούσα'}, 'Location','northwest');
end
