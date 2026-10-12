% MAKE_FIGURES  Figures from results.mat (run sim_main first).
%   figures/1_bore : field throughout the coil interior
%   figures/2_dsv  : 40 mm dia x 40 mm long cylinder at the centre only
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'results.mat'));
s = S.s; d = S.d; F = S.F; Bs = S.B; G = S.G; fil = S.fil; lines = S.lines;
dirA = fullfile(here, 'figures', '1_bore'); dirB = fullfile(here, 'figures', '2_dsv');
if ~exist(dirA,'dir'), mkdir(dirA); end
if ~exist(dirB,'dir'), mkdir(dirB); end
set(groot, 'defaultAxesFontSize', 10, 'defaultTextColor', [.04 .04 .04], ...
    'defaultAxesXColor', [.32 .32 .31], 'defaultAxesYColor', [.32 .32 .31], ...
    'defaultAxesZColor', [.32 .32 .31], 'defaultFigureColor', 'w');

% palette
cBx = hex('2a78d6'); cBy = hex('eb6834'); cBz = hex('1baf7a');
cF  = hex('2a78d6'); cB  = hex('eb6834'); cLay = [cF; cB]; ink2 = hex('52514e'); muted = [.6 .6 .58];
DIV = divmap(); SEQ = seqmap();
mm = 1e3; mT = 1e3;
I0 = max(G.I);

% ---- full fields on the bore grid (NaN outside r = 31 mm) -------------------
toFull = @(v, M) fill(M, v);
Bx = toFull(F.grid(:,1), s.in); By = toFull(F.grid(:,2), s.in); Bz = toFull(F.grid(:,3), s.in);
Bm = sqrt(Bx.^2 + By.^2 + Bz.^2);
xv = squeeze(s.X(:,1,1))*mm; yv = squeeze(s.Y(1,:,1))'*mm; zv = squeeze(s.Z(1,1,:))*mm;
ix0 = find(abs(xv) < 1e-9); iy0 = find(abs(yv) < 1e-9); iz0 = find(abs(zv) < 1e-9);
zCoil = [min(G.toZ([G.p1(:,2);G.p2(:,2)])) max(G.toZ([G.p1(:,2);G.p2(:,2)]))]*mm;
rC = G.rMid*mm;

% gradient at isocentre from the on-axis FEM profile
tz = lines.z(:,3)*mm; sel = abs(tz) <= 5;
p = polyfit(tz(sel), F.lz(sel,1)*mT, 1);  G0 = p(1)*1e3;              % mT/m
pb = polyfit(tz(sel), Bs.lz(sel,1)*mT, 1); G0bs = pb(1)*1e3;

% ============================ PART 1: whole bore ===============================
% 01 flat PCB, both layers
f = figure('Position', [50 50 1500 620]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
lay = {'F.Cu (outside of roll)', 'B.Cu (inside of roll, toward former)'};
for L = 1:2
    ax = nexttile; hold on; axis equal; box on
    m = G.layer == L;
    X = [G.p1(m,1) G.p2(m,1) nan(nnz(m),1)]'*mm; Y = [G.p1(m,2) G.p2(m,2) nan(nnz(m),1)]'*mm;
    plot(X(:), Y(:), '-', 'Color', cLay(L,:), 'LineWidth', 0.8);
    % current direction arrows (every ~6 mm of track)
    mid = (G.p1(m,:) + G.p2(m,:))/2*mm; dv = (G.p2(m,:) - G.p1(m,:));
    len = vecnorm(dv,2,2)*mm; keep = mod(cumsum(len), 6) < len & len > 0.3;
    u = dv(keep,:)./vecnorm(dv(keep,:),2,2);
    quiver(mid(keep,1), mid(keep,2), u(:,1)*2, u(:,2)*2, 0, 'Color', [.1 .1 .1], ...
        'LineWidth', 0.6, 'MaxHeadSize', 2);
    plot(G.vias(:,1)*mm, G.vias(:,2)*mm, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [.1 .1 .1], 'LineWidth', 1.2);
    jp = [G.J1; G.J2]*mm;
    rectangle('Position', [jp(L,:)-1.5 3 3], 'FaceColor', [.85 .85 .83], 'EdgeColor', 'k');
    jl = {'J1+  (3.8 A in)', 'J2-  (3.8 A out)'};
    text(jp(L,1)+2.5, jp(L,2), jl{L}, 'FontWeight', 'bold');
    set(ax, 'YDir', 'reverse'); xlim([5 218]); ylim([5 147]);
    xlabel('PCB x (mm)  \rightarrow wraps around former'); ylabel('PCB y (mm)  \leftarrow +z (bore)');
    title(lay{L});
end
title(tl, sprintf('Z coil PCB as parsed: 2 layers in series, %d tracks, %d vias, I = %.1f A (arrows = current direction)', ...
    numel(G.I), size(G.vias,1), I0), 'FontWeight', 'bold');
save_fig(f, dirA, '01_pcb_layers');

% 02 rolled coil in 3-D
f = figure('Position', [50 50 1100 850]); hold on; axis equal; grid on; box on
for L = [2 1]
    m = fil.layer == L;
    X = [fil.A(m,1) fil.B(m,1) nan(nnz(m),1)]'*mm; Y = [fil.A(m,2) fil.B(m,2) nan(nnz(m),1)]'*mm;
    Z = [fil.A(m,3) fil.B(m,3) nan(nnz(m),1)]'*mm;
    h(L) = plot3(X(:), Y(:), Z(:), '-', 'Color', cLay(L,:), 'LineWidth', 0.7); %#ok<SAGROW>
end
draw_cyl(20, [-20 20], hex('1baf7a'), 0.25);
quiver3(-55, 0, 0, 30, 0, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.6);
text(-58, 0, 6, 'B_0 (+x)', 'FontWeight', 'bold');
plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
legend(h([1 2]), {'F.Cu (outer layer)', 'B.Cu (inner layer)'}, 'Location', 'northeast');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm, bore axis)');
view(35, 18); title('Rolled coil on 65 mm former, 40 mm DSV in green', 'FontWeight', 'bold');
save_fig(f, dirA, '02_coil_3d');

% 03 FEM mesh, xz cut
f = figure('Position', [50 50 1300 600]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nd = S.msh.nodes*mm; k = abs(nd(2,:)) < 0.6;
nexttile; plot(nd(3,k), nd(1,k), '.', 'Color', ink2, 'MarkerSize', 2); axis equal; box on
xlabel('z (mm)'); ylabel('x (mm)'); title('Mesh nodes near y = 0 (full domain)');
nexttile; plot(nd(3,k), nd(1,k), '.', 'Color', ink2, 'MarkerSize', 3); axis equal; box on
xlim([-80 80]); ylim([-36 36]); xlabel('z (mm)'); ylabel('x (mm)');
title('Zoom: copper shells at r \approx 32.6 mm, refined DSV cell');
title(tl, sprintf('PDE Toolbox mesh: %d nodes, %d tetrahedra (solve %.0f s)', ...
    size(S.msh.nodes,2), size(S.msh.elems,2), S.solveTime), 'FontWeight', 'bold');
save_fig(f, dirA, '03_fem_mesh');

% 04 profiles along the three axes
f = figure('Position', [50 50 1400 900]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
prof = {lines.z(:,3), F.lz, Bs.lz, 'z (mm), along bore axis (x = y = 0)'; ...
        lines.x(:,1), F.lx, Bs.lx, 'x (mm), along B_0 direction (y = z = 0)'; ...
        lines.y(:,2), F.ly, Bs.ly, 'y (mm) (x = z = 0)'};
for r = 1:3
    nexttile; hold on; box on; grid on
    t = prof{r,1}*mm; cols = [cBx; cBy; cBz]; nm = {'B_x','B_y','B_z'};
    for c = 1:3
        hh(c) = plot(t, prof{r,2}(:,c)*mT, '-', 'Color', cols(c,:), 'LineWidth', 2); %#ok<SAGROW>
        plot(t, prof{r,3}(:,c)*mT, '--', 'Color', cols(c,:)*0.6, 'LineWidth', 1);
    end
    hh(4) = plot(nan, nan, 'k--', 'LineWidth', 1);
    if r == 1
        xline(zCoil, ':', 'Color', muted); xline([-20 20], '-', 'Color', cBz, 'Alpha', 0.4);
        text(zCoil(2), 0, ' coil end', 'Color', ink2); text(20, 0, ' DSV', 'Color', cBz);
    else
        xline([-rC rC], ':', 'Color', muted); xline([-20 20], '-', 'Color', cBz, 'Alpha', 0.4);
    end
    xlabel(prof{r,4}); ylabel('B (mT)');
    if r == 1, legend(hh, [strcat(nm, ' FEM'), {'Biot-Savart (dashed)'}], 'Location', 'eastoutside'); end
end
title(tl, sprintf('Field along the axes at %.1f A.  G_z = dB_x/dz at centre = %.2f mT/m (FEM), %.2f mT/m (Biot-Savart)', ...
    I0, G0, G0bs), 'FontWeight', 'bold');
save_fig(f, dirA, '04_axis_profiles');

% 05 longitudinal slices
f = figure('Position', [50 50 1500 900]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
Rg = hypot(s.X, s.Y)*mm;
lim = prctile(abs(Bx(Rg <= 26)), 99.5)*mT;   % colour range of the bore interior (wall reaches ~2.7 mT)
sl = {squeeze(Bx(:,iy0,:))*mT, 'B_x in xz plane (y = 0)', 'x (mm)', DIV, lim*[-1 1];
      squeeze(Bx(ix0,:,:))*mT, 'B_x in yz plane (x = 0)', 'y (mm)', DIV, lim*[-1 1];
      squeeze(Bm(:,iy0,:))*mT, '|B| in xz plane (y = 0), with B_x-B_z field lines', 'x (mm)', SEQ, [0 lim];
      squeeze(Bm(ix0,:,:))*mT, '|B| in yz plane (x = 0)', 'y (mm)', SEQ, [0 lim]};
for k = 1:4
    ax = nexttile; imagesc(zv, xv, sl{k,1}, 'AlphaData', ~isnan(sl{k,1})); axis xy equal tight
    colormap(ax, sl{k,4}); clim(sl{k,5}); cb = colorbar; cb.Label.String = 'mT';
    hold on; coil_outline_long(zCoil, rC);
    rectangle('Position', [-20 -20 40 40], 'EdgeColor', cBz, 'LineWidth', 1.2);
    if k == 3
        streamslice(zv, xv, squeeze(Bz(:,iy0,:))*mT, squeeze(Bx(:,iy0,:))*mT, 1.2);
        set(findobj(ax, 'Type', 'line', '-not', 'Color', cBz), 'Color', [.2 .2 .2]);
        coil_outline_long(zCoil, rC);
    end
    xlabel('z (mm)'); ylabel(sl{k,3}); title(sl{k,2});
end
title(tl, 'Field in longitudinal planes through the bore (coil at r = 32.6 mm dashed, DSV boxed)', 'FontWeight', 'bold');
save_fig(f, dirA, '05_slices_longitudinal');

% 06 transverse slices
zs = [-60 -40 -20 0 20 40 60];
f = figure('Position', [50 50 1700 520]);
tl = tiledlayout(1, numel(zs), 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(zs)
    ax = nexttile; iz = find(abs(zv - zs(k)) < 1e-9);
    Mz = squeeze(Bx(:,:,iz))'*mT;
    imagesc(xv, yv, Mz, 'AlphaData', ~isnan(Mz)); axis xy equal tight; hold on
    colormap(ax, DIV); clim(lim*[-1 1]);
    th = linspace(0, 2*pi, 200); plot(rC*cos(th), rC*sin(th), '--', 'Color', ink2);
    if abs(zs(k)) <= 20, plot(20*cos(th), 20*sin(th), '-', 'Color', cBz, 'LineWidth', 1.2); end
    title(sprintf('z = %+d mm', zs(k))); xlabel('x (mm)'); if k == 1, ylabel('y (mm)'); end
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'B_x (mT)';
title(tl, 'B_x in transverse planes (B_0 along +x)', 'FontWeight', 'bold');
save_fig(f, dirA, '06_slices_transverse');

% 07 3-D volume slices
[Ym, Xm, Zm, Bxm] = mg(s, Bx);
f = figure('Position', [50 50 1100 850]);
hs = slice(Xm*mm, Ym*mm, Zm*mm, Bxm*mT, 0, 0, [-50 0 50]);
set(hs, 'EdgeColor', 'none', 'FaceAlpha', 0.9); hold on
colormap(DIV); clim(lim*[-1 1]); cb = colorbar; cb.Label.String = 'B_x (mT)';
light_coil(fil); axis equal; box on; grid on
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); view(40, 20);
title('B_x on slice planes x = 0, y = 0, z = -50 / 0 / +50 mm (coil copper faint)', 'FontWeight', 'bold');
save_fig(f, dirA, '07_volume_slices_3d');

% 08 isosurfaces
f = figure('Position', [50 50 1100 850]); hold on
lv = round(lim*[0.2 0.45 0.7 0.95]*20)/20;
Bxt = Bxm; Bxt(hypot(Xm, Ym)*mm > 28) = nan;
for v = [-fliplr(lv) lv]
    pp = patch(isosurface(Xm*mm, Ym*mm, Zm*mm, Bxt*mT, v));
    c = interp1(linspace(-lim, lim, size(DIV,1)), DIV, v);
    set(pp, 'FaceColor', c, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end
light_coil(fil); camlight; lighting gouraud; axis equal; box on; grid on
colormap(DIV); clim(lim*[-1 1]); cb = colorbar; cb.Label.String = 'B_x (mT)';
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); view(40, 20);
title(sprintf('Iso-surfaces of B_x at \\pm%s mT (r < 28 mm)', strjoin(string(lv), ', ')), 'FontWeight', 'bold');
save_fig(f, dirA, '08_isosurfaces_3d');

% 09 local gradient dBx/dz
BxB = fill(s.in, Bs.grid(:,1));
[~, ~, gF] = gradient(Bx, 1e-3, 1e-3, 2e-3);   % 3rd output = d/d(dim 3) = d/dz
[~, ~, gB] = gradient(BxB, 1e-3, 1e-3, 2e-3);
gF = gF*mT; gB = gB*mT;                          % T/m -> mT/m
f = figure('Position', [50 50 1500 900]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
gm = {gF, 'FEM'; gB, 'Biot-Savart'};
for r = 1:2, for k = 1:2
    ax = nexttile;
    if k == 1, M = squeeze(gm{r,1}(:,iy0,:)); lab = 'x (mm)'; ttl = 'xz plane (y = 0)';
    else,      M = squeeze(gm{r,1}(ix0,:,:)); lab = 'y (mm)'; ttl = 'yz plane (x = 0)'; end
    imagesc(zv, xv, M, 'AlphaData', ~isnan(M)); axis xy equal tight; hold on
    colormap(ax, DIV); clim(G0 + abs(G0)*[-1 1]); cb = colorbar; cb.Label.String = 'mT/m';
    contour(zv, xv, M/G0, [0.95 1.05], 'k', 'LineWidth', 1);
    coil_outline_long(zCoil, rC); rectangle('Position', [-20 -20 40 40], 'EdgeColor', cBz, 'LineWidth', 1.2);
    xlabel('z (mm)'); ylabel(lab); title([gm{r,2} ', ' ttl]);
end, end
title(tl, sprintf(['Local gradient dB_x/dz; grey = centre value %.2f mT/m, red = stronger, blue = weaker/reversed; ' ...
    'black = \\pm5%%.  FEM shows element-scale noise'], G0), 'FontWeight', 'bold');
save_fig(f, dirA, '09_gradient_map');

% 10 FEM vs Biot-Savart
f = figure('Position', [50 50 1300 520]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on; grid on
rr = hypot(s.X(s.in), s.Y(s.in))*mm; inner = rr < 26;
plot(Bs.grid(inner,1)*mT, F.grid(inner,1)*mT, '.', 'Color', cBx, 'MarkerSize', 3);
plot(lim*[-1 1], lim*[-1 1], 'k-');
xlabel('Biot-Savart B_x (mT)'); ylabel('FEM B_x (mT)'); axis equal tight
title('All grid points with r < 26 mm');
nexttile; hold on; box on; grid on
e = vecnorm(F.grid - Bs.grid, 2, 2)./max(abs(Bs.grid(:,1)));
ed = vecnorm(F.dsv - Bs.dsv, 2, 2)./max(abs(Bs.dsv(:,1)));
histogram(100*e(inner), 60, 'FaceColor', cBx, 'EdgeColor', 'w');
histogram(100*ed, 60, 'FaceColor', cBy, 'EdgeColor', 'w');
legend({'bore, r < 26 mm', 'DSV'}); xlabel('|B_{FEM} - B_{BS}| / max|B_x|  (%)'); ylabel('points');
title(sprintf('Difference: DSV median %.2f%%, max %.2f%%', 100*median(ed), 100*max(ed)));
title(tl, 'Cross-check: PDE Toolbox FEM vs direct Biot-Savart of every track', 'FontWeight', 'bold');
save_fig(f, dirA, '10_fem_vs_biot_savart');

% 11 contributions per winding quadrant and layer
f = figure('Position', [50 50 1400 560]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
qn = {'z>0, x>0', 'z>0, x<0', 'z<0, x>0', 'z<0, x<0'};
Gq = zeros(4, 2);
for q = 1:4, for L = 1:2
    pq = polyfit(tz(sel), Bs.parts(sel,1,q,L)*mT, 1); Gq(q,L) = pq(1)*1e3;
end, end
nexttile; b = bar(categorical(qn, qn), Gq, 'grouped'); b(1).FaceColor = cF; b(2).FaceColor = cB;
ylabel('contribution to dB_x/dz at centre (mT/m)'); grid on
legend({'F.Cu', 'B.Cu'}, 'Location', 'best');
title(sprintf('All 8 contributions add with the same sign (sum %.2f mT/m)', sum(Gq(:))));
nexttile; hold on; box on; grid on
hl = plot(tz, squeeze(sum(Bs.parts(:,1,:,1),3))*mT, '-', 'Color', cF, 'LineWidth', 2);
hl(2) = plot(tz, squeeze(sum(Bs.parts(:,1,:,2),3))*mT, '-', 'Color', cB, 'LineWidth', 2);
hl(3) = plot(tz, Bs.lz(:,1)*mT, 'k-', 'LineWidth', 1);
legend(hl, {'F.Cu only', 'B.Cu only', 'total'}, 'Location', 'northwest');
xlabel('z (mm)'); ylabel('B_x on axis (mT)'); xline(zCoil, ':', 'Color', muted);
title('On-axis B_x from each layer');
title(tl, 'Design check: every quadrant and both layers drive the gradient the same way', 'FontWeight', 'bold');
save_fig(f, dirA, '11_layer_quadrant_contributions');

% ============================ PART 2: DSV only ================================
Dx = fill(d.in, F.dsv(:,1)); Dy = fill(d.in, F.dsv(:,2)); Dz = fill(d.in, F.dsv(:,3));
dx = squeeze(d.X(:,1,1))*mm; dy = squeeze(d.Y(1,:,1))'*mm; dz = squeeze(d.Z(1,1,:))*mm;
c0 = find(abs(dx) < 1e-9);
[DYm, DXm, DZm, Dxm] = mg(d, Dx);
dlim = max(abs(Dx(:)), [], 'omitnan')*mT;
NL = 100*(Dx*mT - G0*1e-3*d.Z*mm)/(abs(G0)*1e-3*20);  % % of ideal field at z = 20 mm
NL(~d.in) = nan;
[~, ~, ~, NLm] = mg(d, NL);

% 01 slices (DSV only)
f = figure('Position', [50 50 1100 850]);
hs = slice(DXm*mm, DYm*mm, DZm*mm, Dxm*mT, 0, 0, [-20 0 20]);
set(hs, 'EdgeColor', 'none'); hold on
colormap(DIV); clim(dlim*[-1 1]); cb = colorbar; cb.Label.String = 'B_x (mT)';
draw_cyl(20, [-20 20], [.3 .3 .3], 0); axis equal; box on; grid on
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); view(40, 20);
title('DSV (40 mm dia x 40 mm): B_x on planes x = 0, y = 0, z = -20 / 0 / +20 mm', 'FontWeight', 'bold');
save_fig(f, dirB, '01_dsv_slices_3d');

% 02 field on the DSV surface
f = figure('Position', [50 50 1500 650]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
[TH, ZZ] = meshgrid(linspace(0, 2*pi, 181), linspace(-20, 20, 81));
XS = 20*cos(TH); YS = 20*sin(TH);
bS = interpolateLike(S, XS, YS, ZZ);
for k = 1:2
    ax = nexttile; hold on
    if k == 1, V = reshape(bS(:,1), size(XS))*mT; ttl = 'B_x on the DSV surface (FEM)'; cm = DIV; cl = dlim*[-1 1]; cbl = 'B_x (mT)';
    else
        bb = biot_savart(fil.A, fil.B, fil.I, [XS(:) YS(:) ZZ(:)]/mm);
        V = 100*(reshape(bb(:,1), size(XS))*mT - G0bs*1e-3*ZZ)/(abs(G0bs)*1e-3*20);
        ttl = 'Gradient non-linearity on the DSV surface (Biot-Savart)'; cm = DIV; cl = max(abs(V(:)))*[-1 1]; cbl = '% of G_0 \times 20 mm';
    end
    surf(XS, YS, ZZ, V, 'EdgeColor', 'none');
    [Xc, Yc] = meshgrid(-20:0.5:20); Rc = hypot(Xc, Yc); Xc(Rc > 20) = nan;
    for zc = [-20 20]
        if k == 1
            bc = interpolateLike(S, Xc, Yc, zc + 0*Xc); Vc = reshape(bc(:,1), size(Xc))*mT;
        else
            ok = ~isnan(Xc); bc = nan(numel(Xc), 3);
            bc(ok,:) = biot_savart(fil.A, fil.B, fil.I, [Xc(ok) Yc(ok) zc+0*Xc(ok)]/mm);
            Vc = 100*(reshape(bc(:,1), size(Xc))*mT - G0bs*1e-3*zc)/(abs(G0bs)*1e-3*20);
        end
        surf(Xc, Yc, zc + 0*Xc, Vc, 'EdgeColor', 'none');
    end
    colormap(ax, cm); clim(cl); cb = colorbar; cb.Label.String = cbl;
    axis equal; view(40, 20); box on; grid on; xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
    title(ttl);
end
title(tl, 'DSV boundary (side + end caps)', 'FontWeight', 'bold');
save_fig(f, dirB, '02_dsv_surface');

% 03 isosurfaces in the DSV
f = figure('Position', [50 50 1000 850]); hold on
lv = linspace(-dlim, dlim, 9); lv = lv(2:end-1);
for v = lv
    pp = patch(isosurface(DXm*mm, DYm*mm, DZm*mm, Dxm*mT, v));
    set(pp, 'FaceColor', interp1(linspace(-dlim, dlim, size(DIV,1)), DIV, v), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
end
draw_cyl(20, [-20 20], [.3 .3 .3], 0); camlight; lighting gouraud
colormap(DIV); clim(dlim*[-1 1]); cb = colorbar; cb.Label.String = 'B_x (mT)';
axis equal; box on; grid on; xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); view(40, 15);
title(sprintf('Iso-B_x surfaces in the DSV every %.2f mT (flat & evenly spaced = linear gradient)', lv(2)-lv(1)), 'FontWeight', 'bold');
save_fig(f, dirB, '03_dsv_isosurfaces');

% 04 planes inside the DSV
f = figure('Position', [50 50 1600 480]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
P = {squeeze(Dx(:,c0,:))*mT, dz, dx, 'z (mm)', 'x (mm)', 'B_x, xz plane (y = 0)', DIV, dlim*[-1 1], 'mT';
     squeeze(Dx(c0,:,:))*mT, dz, dy, 'z (mm)', 'y (mm)', 'B_x, yz plane (x = 0)', DIV, dlim*[-1 1], 'mT';
     squeeze(Dx(:,:,end))'*mT, dx, dy, 'x (mm)', 'y (mm)', 'B_x, end cap z = +20 mm', DIV, dlim*[-1 1], 'mT';
     squeeze(NL(:,c0,:)), dz, dx, 'z (mm)', 'x (mm)', 'Non-linearity, xz plane', DIV, [], '%';
     squeeze(NL(c0,:,:)), dz, dy, 'z (mm)', 'y (mm)', 'Non-linearity, yz plane', DIV, [], '%';
     squeeze(NL(:,:,end))', dx, dy, 'x (mm)', 'y (mm)', 'Non-linearity, z = +20 mm', DIV, [], '%'};
nlim = max(abs(NL(:)), [], 'omitnan');
for k = 1:3
    ax = nexttile; M = P{k,1};
    imagesc(P{k,2}, P{k,3}, M, 'AlphaData', ~isnan(M)); axis xy equal tight; hold on
    colormap(ax, P{k,7}); if isempty(P{k,8}), clim(nlim*[-1 1]); else, clim(P{k,8}); end
    cb = colorbar; cb.Label.String = P{k,9};
    if k >= 4, contour(P{k,2}, P{k,3}, M, [-5 -2 2 5], 'k', 'ShowText', 'on'); end
    xlabel(P{k,4}); ylabel(P{k,5}); title(P{k,6});
end
title(tl, sprintf('B_x in DSV planes (FEM), centre gradient G_0 = %.2f mT/m', G0), 'FontWeight', 'bold');
save_fig(f, dirB, '04_dsv_planes');

% 05 profiles inside the DSV
f = figure('Position', [50 50 1400 560]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
offs = [0 0; 10 0; 20 0; 0 10; 0 20; 14 14];
oc = [hex('0b0b0b'); cBx; hex('6250d6'); cBy; hex('e34948'); cBz];
ax1 = nexttile; hold on; box on; grid on; ax2 = nexttile; hold on; box on; grid on;
for k = 1:size(offs,1)
    i1 = find(abs(dx - offs(k,1)) < 1e-9); i2 = find(abs(dy - offs(k,2)) < 1e-9);
    v = squeeze(Dx(i1, i2, :))*mT;
    plot(ax1, dz, v, '-', 'Color', oc(k,:), 'LineWidth', 1.8);
    plot(ax2, dz(2:end-1), (v(3:end) - v(1:end-2))/2*1e3, '-', 'Color', oc(k,:), 'LineWidth', 1.8);
end
plot(ax1, dz, G0*1e-3*dz, 'k:', 'LineWidth', 1);
lg = compose('x = %d, y = %d mm', offs);
legend(ax1, [lg; {'ideal G_0 z'}], 'Location', 'northwest');
xlabel(ax1, 'z (mm)'); ylabel(ax1, 'B_x (mT)'); title(ax1, 'B_x along z at several (x, y)');
xlabel(ax2, 'z (mm)'); ylabel(ax2, 'dB_x/dz (mT/m)'); title(ax2, 'Local gradient along z');
yline(ax2, G0*[0.95 1.05], ':', 'Color', muted);
title(tl, 'Inside the DSV', 'FontWeight', 'bold');
save_fig(f, dirB, '05_dsv_profiles');

% 06 components transverse to B0 (concomitant) in the DSV
f = figure('Position', [50 50 1500 520]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cm3 = {squeeze(Dy(:,c0,:))*mT, 'B_y, xz plane'; squeeze(Dz(:,c0,:))*mT, 'B_z, xz plane'; squeeze(Dz(c0,:,:))*mT, 'B_z, yz plane'};
tl3 = max(abs([Dy(:); Dz(:)]), [], 'omitnan')*mT;
for k = 1:3
    ax = nexttile; imagesc(dz, dx, cm3{k,1}, 'AlphaData', ~isnan(cm3{k,1})); axis xy equal tight
    colormap(ax, DIV); clim(tl3*[-1 1]); cb = colorbar; cb.Label.String = 'mT';
    xlabel('z (mm)'); ylabel([char('x'*(k<3) + 'y'*(k==3)) ' (mm)']); title(cm3{k,2});
end
title(tl, sprintf('Field components perpendicular to B_0 in the DSV (max %.3f mT, vs max B_x %.3f mT)', tl3, dlim), 'FontWeight', 'bold');
save_fig(f, dirB, '06_dsv_transverse_components');

% 07 linearity from the exact Biot-Savart field (FEM element noise is a few %)
BXb = fill(d.in, Bs.dsv(:,1));
NLb = 100*(BXb*mT - G0bs*1e-3*d.Z*mm)/(abs(G0bs)*1e-3*20); NLb(~d.in) = nan;
nlimb = max(abs(NLb(:)), [], 'omitnan');
f = figure('Position', [50 50 1600 820]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
Q = {squeeze(NL(:,c0,:)), 'FEM, xz plane'; squeeze(NL(c0,:,:)), 'FEM, yz plane'; squeeze(NL(:,:,end))', 'FEM, z = +20 mm';
     squeeze(NLb(:,c0,:)), 'Biot-Savart, xz plane'; squeeze(NLb(c0,:,:)), 'Biot-Savart, yz plane'; squeeze(NLb(:,:,end))', 'Biot-Savart, z = +20 mm'};
ql = max(nlim, nlimb);
for k = 1:6
    ax = nexttile; M = Q{k,1};
    if mod(k-1,3) == 2, xa = dx; ya = dy; xl = 'x (mm)'; yl = 'y (mm)';
    else, xa = dz; ya = dx; xl = 'z (mm)'; yl = [char('x'*(mod(k-1,3)==0) + 'y'*(mod(k-1,3)==1)) ' (mm)']; end
    imagesc(xa, ya, M, 'AlphaData', ~isnan(M)); axis xy equal tight; hold on
    colormap(ax, DIV); clim(ql*[-1 1]); cb = colorbar; cb.Label.String = '%';
    contour(xa, ya, M, [-10 -5 -2 2 5 10], 'k', 'ShowText', 'on');
    xlabel(xl); ylabel(yl); title(Q{k,2});
end
title(tl, sprintf(['Gradient non-linearity in the DSV, %% of G_0 \\times 20 mm.  FEM: max %.1f%%, RMS %.1f%%   |   ' ...
    'Biot-Savart: max %.1f%%, RMS %.1f%%'], nlim, rms(NL(d.in)), nlimb, rms(NLb(d.in))), 'FontWeight', 'bold');
save_fig(f, dirB, '07_dsv_linearity_fem_vs_biot_savart');

% ---- numbers -----------------------------------------------------------------
fid = fopen(fullfile(here, 'figures', 'summary.txt'), 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];
pr('Current: %.2f A (J1 F.Cu in -> J2 B.Cu out)\n', I0);
pr('Gradient dBx/dz at centre: FEM %.3f mT/m, Biot-Savart %.3f mT/m\n', G0, G0bs);
pr('Efficiency: %.3f mT/m/A (FEM), %.3f mT/m/A (BS)\n', G0/I0, G0bs/I0);
Gfit = (d.Z(d.in)*mm) \ (Dx(d.in)*mT) * 1e3;
pr('Least-squares gradient over DSV: %.3f mT/m\n', Gfit);
pr('Bx at DSV ends (axis, z = +/-20 mm): %+.4f / %+.4f mT\n', interp1(tz, F.lz(:,1), [20 -20])*mT);
pr('DSV non-linearity vs centre gradient: max %.2f %%, RMS %.2f %%\n', nlim, rms(NL(d.in)));
pr('DSV non-linearity (Biot-Savart): max %.2f %%, RMS %.2f %%\n', nlimb, rms(NLb(d.in)));
pr('Least-squares gradient over DSV (Biot-Savart): %.3f mT/m\n', (d.Z(d.in)*mm) \ (Bs.dsv(:,1)*mT) * 1e3);
pr('Max |By|, |Bz| in DSV: FEM %.4f, %.4f mT; Biot-Savart %.4f, %.4f mT\n', max(abs(Dy(:)))*mT, ...
    max(abs(Dz(:)))*mT, max(abs(Bs.dsv(:,2)))*mT, max(abs(Bs.dsv(:,3)))*mT);
pr('FEM vs Biot-Savart in DSV: median %.3f %%, max %.3f %% of max|Bx|\n', 100*median(ed), 100*max(ed));
pr('Quadrant/layer contributions (mT/m), rows %s, cols F.Cu B.Cu:\n', strjoin(qn, ' | '));
pr('  %+.3f %+.3f\n', Gq');
fclose(fid);
close all

% ============================== helpers ======================================
function A = fill(M, v), A = nan(size(M)); A(M) = v; end
function c = hex(h), c = sscanf(h, '%2x%2x%2x')'/255; end
function M = divmap()
k = [hex('1c4f8f'); hex('2a78d6'); hex('f0efec'); hex('e34948'); hex('8f2222')];
M = interp1([0 .3 .5 .7 1], k, linspace(0, 1, 256));
end
function M = seqmap()
k = [hex('f7f7f5'); hex('a9c8f0'); hex('2a78d6'); hex('123a6b')];
M = interp1([0 .3 .7 1], k, linspace(0, 1, 256));
end
function [Ym, Xm, Zm, Vm] = mg(g, V)
Xm = permute(g.X, [2 1 3]); Ym = permute(g.Y, [2 1 3]); Zm = permute(g.Z, [2 1 3]);
Vm = permute(V, [2 1 3]);
end
function save_fig(f, dirn, name)
exportgraphics(f, fullfile(dirn, [name '.png']), 'Resolution', 150);
end
function coil_outline_long(zc, rc)
plot(zc, [rc rc], '--', 'Color', [.32 .32 .31], 'LineWidth', 1.2);
plot(zc, -[rc rc], '--', 'Color', [.32 .32 .31], 'LineWidth', 1.2);
end
function draw_cyl(r, zr, c, alpha)
[th, zz] = meshgrid(linspace(0, 2*pi, 73), zr);
if alpha > 0
    surf(r*cos(th), r*sin(th), zz, 'FaceColor', c, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end
t = linspace(0, 2*pi, 145);
for z = zr, plot3(r*cos(t), r*sin(t), z + 0*t, '-', 'Color', c*0.8, 'LineWidth', 1.2); end
for a = 0:pi/2:3*pi/2, plot3(r*cos(a)*[1 1], r*sin(a)*[1 1], zr, '-', 'Color', c*0.8, 'LineWidth', 1); end
end
function light_coil(fil)
k = 1:3:size(fil.A,1);
X = [fil.A(k,1) fil.B(k,1) nan(numel(k),1)]'*1e3; Y = [fil.A(k,2) fil.B(k,2) nan(numel(k),1)]'*1e3;
Z = [fil.A(k,3) fil.B(k,3) nan(numel(k),1)]'*1e3;
plot3(X(:), Y(:), Z(:), '-', 'Color', [.55 .55 .53 .35], 'LineWidth', 0.4);
end
function b = interpolateLike(S, X, Y, Z)
% DSV-surface values: linear interpolation of the stored FEM bore grid
% (1 x 1 x 2 mm, valid out to r = 31 mm so the DSV boundary is well inside it)
persistent Fi
if isempty(Fi)
    for c = 1:3
        V = nan(size(S.s.X)); V(S.s.in) = S.F.grid(:,c);
        Fi{c} = griddedInterpolant(S.s.X*1e3, S.s.Y*1e3, S.s.Z*1e3, V, 'linear', 'none'); %#ok<AGROW>
    end
end
b = [Fi{1}(X(:), Y(:), Z(:)) Fi{2}(X(:), Y(:), Z(:)) Fi{3}(X(:), Y(:), Z(:))];
end
