function C = parse_kicad_coil(pcbFile, I0)
%PARSE_KICAD_COIL  Read copper of a 2-layer KiCad flex coil and solve the DC
%current in every track piece for a current I0 injected at J1 (F.Cu) and
%extracted at J2 (B.Cu).
%
%   Conductors: (segment ...) tracks, copper fp_lines inside the NetTie
%   footprints (they bridge the turn-to-turn nets), and (via ...) between
%   layers. Current directions/magnitudes come from a resistive network solve,
%   so segment start/end order in the file does not matter and any parallel
%   paths are split correctly.
%
%   Output C (all lengths in mm, PCB coordinates):
%     C.p1, C.p2   [N x 2] segment endpoints, oriented along positive current
%     C.layer      [N x 1] 1 = F.Cu, 2 = B.Cu
%     C.width      [N x 1] track width
%     C.I          [N x 1] current (A, >= 0 after orientation)
%     C.vias       [M x 2] via positions, C.viaI their current (F->B positive)
%     C.J1, C.J2   pad positions

txt = fileread(pcbFile);
num = '(-?[\d.]+)';

% ---- top-level tracks ------------------------------------------------------
pat = ['\(segment\s*\(start ' num ' ' num '\)\s*\(end ' num ' ' num ...
       '\)\s*\(width ' num '\)\s*\(layer "([FB])\.Cu"\)'];
t = regexp(txt, pat, 'tokens');
t = vertcat(t{:});
seg = str2double(t(:,1:5));
segLayer = 1 + strcmp(t(:,6), 'B');

% ---- footprints: net-tie copper lines and the J1/J2 pads ------------------
fpStart = regexp(txt, '\n\t\(footprint ');
fpEnd   = [fpStart(2:end), numel(txt)];
tie = zeros(0,5); tieLayer = zeros(0,1);
J1 = []; J2 = [];
for k = 1:numel(fpStart)
    b = txt(fpStart(k):fpEnd(k));
    at = str2double(regexp(b, ['\n\t\t\(at ' num ' ' num '( -?[\d.]+)?\)'], 'tokens', 'once'));
    if numel(at) == 3 && ~isnan(at(3)) && at(3) ~= 0
        error('Rotated footprint not supported');
    end
    if contains(b, '"Reference" "J1')
        J1 = at(1:2);
    elseif contains(b, '"Reference" "J2')
        J2 = at(1:2);
    end
    lt = regexp(b, ['\(fp_line\s*\(start ' num ' ' num '\)\s*\(end ' num ' ' num ...
        '\)\s*\(stroke\s*\(width ' num '\)[^)]*\)\s*\)\s*\(layer "([FB])\.Cu"\)'], 'tokens');
    for j = 1:numel(lt)
        v = str2double(lt{j}(1:5));
        tie(end+1,:) = [v(1:2)+at(1:2), v(3:4)+at(1:2), v(5)]; %#ok<AGROW>
        tieLayer(end+1,1) = 1 + strcmp(lt{j}{6}, 'B');         %#ok<AGROW>
    end
end

% ---- vias -------------------------------------------------------------------
vt = regexp(txt, ['\(via\s*\(at ' num ' ' num '\)'], 'tokens');
vias = str2double(vertcat(vt{:}));

allSeg   = [seg; tie];
allLayer = [segLayer; tieLayer];
nS = size(allSeg,1);

% ---- node graph -------------------------------------------------------------
% node key = rounded (x,y) + layer; vias join the two layers at one (x,y)
tol = 5e-3;   % mm; net-tie ends are rounded to 0.1 um in the file
P = [allSeg(:,1:2), allLayer; allSeg(:,3:4), allLayer; ...
     vias, ones(size(vias,1),1); vias, 2*ones(size(vias,1),1); ...
     J1, 1; J2, 2];
[~, ~, id] = uniquetol([P(:,1:2), 1e3*P(:,3)], tol, 'ByRows', true, 'DataScale', 1);
nN = max(id);
a = id(1:nS); b = id(nS+1:2*nS);
nV = size(vias,1);
va = id(2*nS+(1:nV)); vb = id(2*nS+nV+(1:nV));
nJ1 = id(end-1); nJ2 = id(end);

% conductances: copper 35 um, sigma(Cu) ; via barrel ~ 0.3 mm drill, 20 um plating
sigma = 5.8e4;            % S/mm
tCu   = 0.035;            % mm
L  = max(hypot(allSeg(:,3)-allSeg(:,1), allSeg(:,4)-allSeg(:,2)), 1e-4);
g  = sigma*tCu*allSeg(:,5)./L;
gv = sigma*pi*0.3*0.02/0.095*ones(nV,1);

ea = [a; va]; eb = [b; vb]; ge = [g; gv];
G = sparse([ea; eb; ea; eb], [eb; ea; ea; eb], [-ge; -ge; ge; ge], nN, nN);
rhs = zeros(nN,1); rhs(nJ1) = I0;
% only copper connected to J1/J2 can carry current; ground J2
comp = conncomp(graph(ea, eb, [], nN));
if comp(nJ1) ~= comp(nJ2), error('J1 and J2 are not connected'); end
live = comp(:) == comp(nJ1);
keep = find(live & (1:nN)' ~= nJ2);
V = zeros(nN,1);
V(keep) = G(keep,keep) \ rhs(keep);
% floating islands: set potential equal at both ends -> zero current
V(~live) = 0;

Iseg = g .* (V(a) - V(b));      % current a -> b
Ivia = gv .* (V(va) - V(vb));   % F -> B

% orient segments with positive current
flip = Iseg < 0;
p1 = allSeg(:,1:2); p2 = allSeg(:,3:4);
tmp = p1(flip,:); p1(flip,:) = p2(flip,:); p2(flip,:) = tmp;

% connectivity diagnostics
deg = accumarray([ea; eb], 1, [nN 1]);
C = struct('p1',p1, 'p2',p2, 'layer',allLayer, 'width',allSeg(:,5), ...
    'I',abs(Iseg), 'vias',vias, 'viaI',Ivia, 'J1',J1, 'J2',J2, ...
    'R', V(nJ1)/I0, 'nodeDegree',deg, 'nTie',size(tie,1), 'nTrack',size(seg,1));
end
