function out = run_fem(opts)
%RUN_FEM  3-D magnetostatic FEM (PDE Toolbox) of the rolled 2-layer Z coil.
%
%   Geometry (m):  big air cylinder  >  coil assembly
%     coil assembly = inner air (r < rIn) + copper shell rIn..rOut
%     inner air contains a 40 mm dia x 40 mm DSV cell
%   The shell carries J = J(F.Cu) + J(B.Cu): both layers' real track
%   patterns (current spread across each track's width), each mapped with
%   its own arc-length scaling. The two layers are 60 um apart, far below any
%   mesh size, so they share one 1 mm shell; the sum is also the only
%   divergence-free combination (current hops layers at the vias).
%   The PDE Toolbox 3-D magnetostatic solver needs ~7 GB at ~15k nodes, so
%   the sheet current is smoothed to the shell mesh size (sigma) before
%   sampling; this only removes track-scale ripple that cannot reach the bore.
arguments
    opts.pcb = fullfile(fileparts(mfilename('fullpath')), 'Z_coils.kicad_pcb')
    opts.I0 = 3.8
    opts.tShell = 1e-3
    opts.Hshell = 5e-3
    opts.Hdsv = 5e-3
    opts.Hin = 8e-3
    opts.Hmax = 40e-3
    opts.meshOnly = false
end

C = parse_kicad_coil(opts.pcb, opts.I0);
G = coil_geometry(C);
S = current_sheet(G, 0.25e-3, opts.Hshell/2);

% ---- geometry -------------------------------------------------------------
rIn = G.rMid - opts.tShell/2; rOut = G.rMid + opts.tShell/2;
zAll = G.toZ([G.p1(:,2); G.p2(:,2)]);
zLo = min(zAll) - 3e-3; zHi = max(zAll) + 3e-3;
Rair = 0.12; Hair = 0.40;

gCoil = translate(multicylinder([rIn rOut], zHi - zLo), [0 0 zLo]);
gDSV  = translate(multicylinder(0.020, 0.040), [0 0 -0.020]);
gCoil = addCell(gCoil, gDSV);
g = addCell(translate(multicylinder(Rair, Hair), [0 0 -Hair/2]), gCoil);
g = fegeometry(g);

model = femodel(AnalysisType = "magnetostatic", Geometry = g);
model.VacuumPermeability = 4*pi*1e-7;
model.MaterialProperties = materialProperties(RelativePermeability = 1);

% identify cells with a coarse probe mesh
probe = generateMesh(model, Hmax = 0.03);
cellInfo = zeros(g.NumCells, 4);              % min r, min z, max z, max r
for c = 1:g.NumCells
    e = findElements(probe.Geometry.Mesh, 'region', 'Cell', c);
    nd = probe.Geometry.Mesh.Nodes(:, unique(probe.Geometry.Mesh.Elements(1:4, e)));
    r = hypot(nd(1,:), nd(2,:));
    cellInfo(c,:) = [min(r), min(nd(3,:)), max(nd(3,:)), max(r)];
end
[~, cAir] = max(cellInfo(:,3) - cellInfo(:,2));
cCu  = find(abs(cellInfo(:,4) - rOut) < 1e-4);   % classify by outer radius
cDSV = find(abs(cellInfo(:,4) - 0.02) < 1e-4);
cIn  = setdiff(1:g.NumCells, [cAir cCu(:)' cDSV(:)']);
assert(isscalar(cCu) && isscalar(cDSV) && isscalar(cIn), 'cell id failed');

% outer boundary faces
fOut = unique([nearestFace(g, [0 0 Hair/2]), nearestFace(g, [0 0 -Hair/2]), ...
               nearestFace(g, [Rair 0 0])]);
model.FaceBC(fOut) = faceBC(MagneticPotential = [0;0;0]);

% ---- current density in the copper shell (sum of both layers) -------------
xmin = S.xg(1);
FX = griddedInterpolant({S.yg, S.xg}, sum(S.Kx, 3), 'linear', 'none');
FY = griddedInterpolant({S.yg, S.xg}, sum(S.Ky, 3), 'linear', 'none');
t = opts.tShell;
model.CellLoad(cCu) = cellLoad(CurrentDensity = @(loc, ~) shellJ(loc, FX, FY, G, t, xmin));

fSh  = cellFaces(g, cCu);
fDSV = cellFaces(g, cDSV);
fIn  = setdiff(cellFaces(g, cIn), [fSh fDSV]);
model = generateMesh(model, Hmax = opts.Hmax, Hgrad = 1.4, ...
    Hface = {fSh, opts.Hshell, fDSV, opts.Hdsv, fIn, opts.Hin});
msh = model.Geometry.Mesh;
fprintf('mesh: %d nodes, %d elements (%s)\n', size(msh.Nodes,2), size(msh.Elements,2), msh.GeometricOrder);
out.cells = struct('air',cAir,'cu',cCu,'dsv',cDSV,'in',cIn);
out.radii = [rIn rOut]; out.zRange = [zLo zHi];
out.G = G; out.opts = opts;
if opts.meshOnly, out.model = model; return; end

tic; R = solve(model); out.solveTime = toc;
fprintf('solved in %.1f s\n', out.solveTime);
out.R = R; out.model = model;
end

function J = shellJ(loc, FX, FY, G, t, xmin)
% flat-PCB (x,y) for each quadrature point; arc length uses the flex neutral
% plane, J_theta scaled by r/rMid so that div J = 0 in the shell.
th = atan2(loc.y, loc.x);
r  = hypot(loc.x, loc.y);
xp = G.x0 + (th + G.th0)*G.rMid;
circ = 2*pi*G.rMid;
xp = xmin + mod(xp - xmin, circ);
yp = G.y0 - loc.z;
kx = FX(yp, xp); ky = FY(yp, xp);
kx(isnan(kx)) = 0; ky(isnan(ky)) = 0;
Jt = kx.*r/(G.rMid*t);                 % along +theta (PCB +x)
Jz = -ky/t;                            % PCB +y is -z
J = [-Jt.*sin(th); Jt.*cos(th); Jz];
end
