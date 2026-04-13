close all;
clear;
clf;
clc;

a = 1.0;
b = 2.0*a;
E0 = 1.0;

Nx = 301; Ny = 601;
xE = linspace(-3*a, 3*a, Nx);
yE = linspace(-6*a, 6*a, Ny);
[xgE,ygE] = meshgrid(xE,yE);

theta = 0.0;
eps_list = [2,4,10,Inf];
Nmax = 40;

for k = 1:numel(eps_list)
    epsr = eps_list(k);
    if isinf(epsr)
        xi = 1.0; denom_inside = Inf;
    else
        xi = (epsr-1)/(epsr+1); denom_inside = epsr+1;
    end
    S = (a^2)/(b^2) * (pi^2/12);
    E0x = E0*cos(theta); E0y = E0*sin(theta);
    denom_x = 1 + xi*S; if abs(denom_x)<1e-12, denom_x = 1e-12*(2*(denom_x>=0)-1); end
    denom_y = 1 - xi*S; if abs(denom_y)<1e-12, denom_y = 1e-12*(2*(denom_y>=0)-1); end
    Ax = (xi/denom_x)*E0x; Ay = (xi/denom_y)*E0y;
    if isinf(denom_inside)
        Bx = 0; By = 0;
    else
        Elocx = E0x - S*Ax; Elocy = E0y + S*Ay;
        Bx = -2/denom_inside * Elocx; By = -2/denom_inside * Elocy;
    end

    V = -E0x.*xgE - E0y.*ygE;
    for n = -Nmax:Nmax
        yc = 2*b*n;
        dx = xgE;
        dy = ygE - yc;
        R2 = dx.^2 + dy.^2;
        inside = (R2 < a^2);
        mask = (~inside) & (R2>0);
        Vout = zeros(size(V));
        Vout(mask) = a^2 .* ( Ax.*dx(mask)./R2(mask) + Ay.*dy(mask)./R2(mask) );
        Vin = Bx.*dx + By.*dy;
        V = V + Vin.*inside + Vout.*(~inside);
    end

    dy_ = yE(2)-yE(1); dx_ = xE(2)-xE(1);
    [Vy,Vx] = gradient(V,dy_,dx_);
    Ex = -Vx; Ey = -Vy;

    Xn = xgE./a; Yn = ygE./a; rba = b/a;

    figure; % V
    vabs = max(abs(V(:))); if vabs==0, vabs=1; end
    contourf(Xn, Yn, V, linspace(-vabs,vabs,31), 'LineStyle','none'); colorbar; hold on;
    ymin = min(Yn(:)); ymax = max(Yn(:));
    nmin = floor(ymin/(2*rba)) - 1; nmax = ceil(ymax/(2*rba)) + 1;
    t = linspace(0,2*pi,200);
    for nn = nmin:nmax
        ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-');
    end
    if isinf(epsr), tag='PEC'; else, tag=sprintf('eps%g',epsr); end
    title(['(E) V, ' tag ', b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;

    figure; % E
    step = max(1, floor(max(size(Xn))/40));
    quiver(Xn(1:step:end,1:step:end), Yn(1:step:end,1:step:end), ...
           Ex(1:step:end,1:step:end), Ey(1:step:end,1:step:end), 0.9); hold on;
    for nn = nmin:nmax
        ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-');
    end
    title(['(E) E, ' tag ', b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;
end

epsr = 4; thetas = [10 30 45 60]*pi/180;

for theta = thetas
    xi = (epsr-1)/(epsr+1);
    S = (a^2)/(b^2) * (pi^2/12);
    E0x = E0*cos(theta); E0y = E0*sin(theta);
    denom_x = 1 + xi*S; if abs(denom_x)<1e-12, denom_x = 1e-12*(2*(denom_x>=0)-1); end
    denom_y = 1 - xi*S; if abs(denom_y)<1e-12, denom_y = 1e-12*(2*(denom_y>=0)-1); end
    Ax = (xi/denom_x)*E0x; Ay = (xi/denom_y)*E0y;
    Elocx = E0x - S*Ax; Elocy = E0y + S*Ay;
    Bx = -2/(epsr+1)*Elocx; By = -2/(epsr+1)*Elocy;

    V = -E0x.*xgE - E0y.*ygE;
    for n = -Nmax:Nmax
        yc = 2*b*n;
        dx = xgE; dy = ygE - yc; R2 = dx.^2 + dy.^2;
        inside = (R2 < a^2);
        mask = (~inside) & (R2>0);
        Vout = zeros(size(V));
        Vout(mask) = a^2 .* ( Ax.*dx(mask)./R2(mask) + Ay.*dy(mask)./R2(mask) );
        Vin = Bx.*dx + By.*dy;
        V = V + Vin.*inside + Vout.*(~inside);
    end

    dy_ = yE(2)-yE(1); dx_ = xE(2)-xE(1);
    [Vy,Vx] = gradient(V,dy_,dx_); Ex = -Vx; Ey = -Vy;

    Xn = xgE./a; Yn = ygE./a; rba = b/a;
    ymin = min(Yn(:)); ymax = max(Yn(:));
    nmin = floor(ymin/(2*rba)) - 1; nmax = ceil(ymax/(2*rba)) + 1; t = linspace(0,2*pi,200);

    figure; % V
    vabs = max(abs(V(:))); if vabs==0, vabs=1; end
    contourf(Xn, Yn, V, linspace(-vabs,vabs,31), 'LineStyle','none'); colorbar; hold on;
    for nn = nmin:nmax, ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-'); end
    title(['(Z) V, eps=4, b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;

    figure; % E
    step = max(1, floor(max(size(Xn))/40));
    quiver(Xn(1:step:end,1:step:end), Yn(1:step:end,1:step:end), ...
           Ex(1:step:end,1:step:end), Ey(1:step:end,1:step:end), 0.9); hold on;
    for nn = nmin:nmax, ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-'); end
    title(['(Z) E, eps=4, b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;
end

epsr = 10; theta = pi/4; rlist = [1.1 1.5 2.5 4.0];
NmaxH = 60;

for r = rlist
    b2 = r*a;
    xH = linspace(-3*a, 3*a, Nx);
    yH = linspace(-8*b2, 8*b2, Ny);
    [xgH,ygH] = meshgrid(xH,yH);

    xi = (epsr-1)/(epsr+1);
    S = (a^2)/(b2^2) * (pi^2/12);
    E0x = E0*cos(theta); E0y = E0*sin(theta);
    denom_x = 1 + xi*S; if abs(denom_x)<1e-12, denom_x = 1e-12*(2*(denom_x>=0)-1); end
    denom_y = 1 - xi*S; if abs(denom_y)<1e-12, denom_y = 1e-12*(2*(denom_y>=0)-1); end
    Ax = (xi/denom_x)*E0x; Ay = (xi/denom_y)*E0y;
    Elocx = E0x - S*Ax; Elocy = E0y + S*Ay;
    Bx = -2/(epsr+1)*Elocx; By = -2/(epsr+1)*Elocy;

    V = -E0x.*xgH - E0y.*ygH;
    for n = -NmaxH:NmaxH
        yc = 2*b2*n;
        dx = xgH; dy = ygH - yc; R2 = dx.^2 + dy.^2;
        inside = (R2 < a^2);
        mask = (~inside) & (R2>0);
        Vout = zeros(size(V));
        Vout(mask) = a^2 .* ( Ax.*dx(mask)./R2(mask) + Ay.*dy(mask)./R2(mask) );
        Vin = Bx.*dx + By.*dy;
        V = V + Vin.*inside + Vout.*(~inside);
    end

    dy_ = yH(2)-yH(1); dx_ = xH(2)-xH(1);
    [Vy,Vx] = gradient(V,dy_,dx_); Ex = -Vx; Ey = -Vy;

    Xn = xgH./a; Yn = ygH./a; rba = b2/a;
    ymin = min(Yn(:)); ymax = max(Yn(:));
    nmin = floor(ymin/(2*rba)) - 1; nmax = ceil(ymax/(2*rba)) + 1; t = linspace(0,2*pi,200);

    figure; % V
    vabs = max(abs(V(:))); if vabs==0, vabs=1; end
    contourf(Xn, Yn, V, linspace(-vabs,vabs,31), 'LineStyle','none'); colorbar; hold on;
    for nn = nmin:nmax, ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-'); end
    title(['(H) V, eps=10, b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;

    figure; % E
    step = max(1, floor(max(size(Xn))/40));
    quiver(Xn(1:step:end,1:step:end), Yn(1:step:end,1:step:end), ...
           Ex(1:step:end,1:step:end), Ey(1:step:end,1:step:end), 0.9); hold on;
    for nn = nmin:nmax, ycN = 2*rba*nn; plot(cos(t), ycN + sin(t), 'k-'); end
    title(['(H) E, eps=10, b/a=' num2str(rba) ', theta=' num2str(theta) ' rad']);
    xlabel('x/a'); ylabel('y/a'); axis equal tight;
end
