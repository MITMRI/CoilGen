% SIM_MAIN  Solve the rolled Z-coil with PDE Toolbox and sample B on grids.
% Results -> results.mat (FEM field + independent Biot-Savart reference).
here = fileparts(mfilename('fullpath'));
out = run_fem();
R = out.R; G = out.G;

% ---- sample grids (m) ---------------------------------------------------------
mm = 1e-3;
% whole interior of the coil
[s.X, s.Y, s.Z] = ndgrid((-31:1:31)*mm, (-31:1:31)*mm, (-80:2:80)*mm);
s.in = hypot(s.X, s.Y) <= 31*mm;
% DSV: 40 mm dia x 40 mm long cylinder
[d.X, d.Y, d.Z] = ndgrid((-20:1:20)*mm, (-20:1:20)*mm, (-20:1:20)*mm);
d.in = hypot(d.X, d.Y) <= 20*mm + 1e-9 & abs(d.Z) <= 20*mm + 1e-9;
% lines
t = (-100:0.5:100)'*mm;  ta = (-31:0.25:31)'*mm;
lines.z  = [0*t 0*t t];                 % bore axis
lines.x  = [ta 0*ta 0*ta];
lines.y  = [0*ta ta 0*ta];

fem = @(P) femB(R, P);
F.grid = fem([s.X(s.in) s.Y(s.in) s.Z(s.in)]);
F.dsv  = fem([d.X(d.in) d.Y(d.in) d.Z(d.in)]);
F.lz = fem(lines.z); F.lx = fem(lines.x); F.ly = fem(lines.y);

% ---- Biot-Savart reference (actual copper radii, track centre lines) -----------
fil = coil_filaments(G);
bs = @(P) biot_savart(fil.A, fil.B, fil.I, P);
tic;
B.grid = bs([s.X(s.in) s.Y(s.in) s.Z(s.in)]);
B.dsv  = bs([d.X(d.in) d.Y(d.in) d.Z(d.in)]);
B.lz = bs(lines.z); B.lx = bs(lines.x); B.ly = bs(lines.y);
% per group / per layer contributions on the bore axis (design check)
grp = zeros(size(fil.I));
tz = G.toZ((G.p1(:,2)+G.p2(:,2))/2); tx = G.toTheta((G.p1(:,1)+G.p2(:,1))/2);
grp(:) = 1 + (tz(fil.track) < 0)*2 + (cos(tx(fil.track)) < 0);   % 4 quadrants
B.parts = zeros(numel(t), 3, 4, 2);
for gI = 1:4, for L = 1:2
    m = grp == gI & fil.layer == L;
    B.parts(:,:,gI,L) = biot_savart(fil.A(m,:), fil.B(m,:), fil.I(m), lines.z);
end, end
fprintf('Biot-Savart: %.1f s\n', toc);

msh.nodes = R.Mesh.Nodes; msh.elems = R.Mesh.Elements(1:4,:);
cells = out.cells; radii = out.radii; zRange = out.zRange; solveTime = out.solveTime;
save(fullfile(here, 'results.mat'), 's', 'd', 'lines', 'F', 'B', 'G', 'fil', ...
    'cells', 'radii', 'zRange', 'solveTime', 'msh', '-v7.3');
fprintf('saved\n');

function Bq = femB(R, P)
b = interpolateMagneticFlux(R, P(:,1), P(:,2), P(:,3));
Bq = [b.Bx(:) b.By(:) b.Bz(:)];
end
