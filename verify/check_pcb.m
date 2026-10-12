function R = check_pcb(pcbFile)
%CHECK_PCB  Connectivity + clearance verification of the coil copper.
%   - every copper node must have exactly 2 connections (series path), except
%     J1 and J2 which have 1; anything else is a stub / branch / open
%   - all copper must be in the J1-J2 network (no floating islands)
%   - copper-to-copper clearance between different nets on each layer
%     (tracks, vias, J pads, net-tie lines), including across the wrap seam

txt = fileread(pcbFile);
num = '(-?[\d.]+)';
t = regexp(txt, ['\(segment\s*\(start ' num ' ' num '\)\s*\(end ' num ' ' num ...
    '\)\s*\(width ' num '\)\s*\(layer "([FB])\.Cu"\)\s*\(net (\d+)\)'], 'tokens');
t = vertcat(t{:});
S = str2double(t(:,[1:5 7]));                 % x1 y1 x2 y2 w net
S(:,7) = 1 + strcmp(t(:,6),'B');              % layer

% net ties: copper line + the two nets it bridges
fpStart = regexp(txt, '\n\t\(footprint ');
fpEnd = [fpStart(2:end), numel(txt)];
T = zeros(0,7); tiePairs = zeros(0,2); tieXY = zeros(0,3);
for k = 1:numel(fpStart)
    b = txt(fpStart(k):fpEnd(k));
    if ~contains(b, 'NetTie'), continue; end
    at = str2double(regexp(b, ['\n\t\t\(at ' num ' ' num], 'tokens', 'once'));
    lt = str2double(regexp(b, ['\(fp_line\s*\(start ' num ' ' num '\)\s*\(end ' num ' ' num ...
        '\)\s*\(stroke\s*\(width ' num], 'tokens', 'once'));
    L = 1 + contains(regexp(b, '\(fp_line.*?\(layer "(.)\.Cu"', 'tokens', 'once'), 'B');
    pn = regexp(b, '\(pad "\d" smd \w+.*?\(net (\d+)', 'tokens');
    pn = str2double([pn{:}]);
    T(end+1,:) = [lt(1:2)+at, lt(3:4)+at, lt(5), pn(1), L]; %#ok<AGROW>
    tiePairs(end+1,:) = sort(pn(1:2));                        %#ok<AGROW>
    tieXY(end+1,:) = [at L];                                  %#ok<AGROW>
end
vt = regexp(txt, ['\(via\s*\(at ' num ' ' num '\)\s*\(size ' num '\).*?\(net (\d+)\)'], 'tokens');
V = str2double(vertcat(vt{:}));               % x y size net

% ---------------- connectivity ----------------------------------------------
C = parse_kicad_coil(pcbFile, 1);
R.degreeHist = accumarray(C.nodeDegree+1, 1)';  % index 1 -> degree 0
R.nDeg1 = nnz(C.nodeDegree == 1);
R.nDegGt2 = nnz(C.nodeDegree > 2);
R.nDead = nnz(C.I < 1e-6);                      % copper carrying no current
R.minI = min(C.I); R.maxI = max(C.I);
R.resistance = C.R;

% ---------------- clearance ---------------------------------------------------
h = 0.05;                                       % mm sampling along track
P = zeros(0,5);                                 % x y halfwidth net layer
allS = [S; T];
for i = 1:size(allS,1)
    a = allS(i,1:2); b = allS(i,3:4);
    n = max(2, ceil(norm(b-a)/h)+1);
    s = linspace(0,1,n)';
    P = [P; a + s.*(b-a), repmat([allS(i,5)/2, allS(i,6), allS(i,7)], n, 1)]; %#ok<AGROW>
end
for i = 1:size(V,1)       % via annulus on both layers
    P = [P; V(i,1:2) V(i,3)/2 V(i,4) 1; V(i,1:2) V(i,3)/2 V(i,4) 2]; %#ok<AGROW>
end
[gx, gy] = meshgrid(-1.5:0.1:1.5);              % J1 (F, net 1) and J2 (B, net 62) 3x3 pads
P = [P; 12.5+gx(:) 12.5+gy(:) zeros(numel(gx),1) ones(numel(gx),1) ones(numel(gx),1)];
P = [P; 12.5+gx(:) 12.5+gy(:) zeros(numel(gx),1) 62*ones(numel(gx),1) 2*ones(numel(gx),1)];

% wrap seam: board spans x = 10 .. 213.7035; after rolling, the two board
% edges meet. Copy points near the far edge to the near side.
rMid = 32.5 + 0.0975;                           % neutral plane of the 0.195 mm flex
circ = 2*pi*rMid;
near = P(:,1) > 213.7035 - 10;
Pw = P(near,:); Pw(:,1) = Pw(:,1) - circ;
Pall = [P; Pw];

R.minClear = inf(2,1); R.worst = cell(2,1);
viol = zeros(0,7);
K = 300;
for L = 1:2
    Q = Pall(Pall(:,5)==L,:);
    nets = unique(Q(:,4));
    best = inf; bw = [];
    for n = nets'
        A = Q(Q(:,4)==n,:); B = Q(Q(:,4)~=n,:);
        [idx, d] = knnsearch(B(:,1:2), A(:,1:2), 'K', K);
        c = d - A(:,3) - reshape(B(idx,3), size(idx));
        [cm, j] = min(c, [], 2);
        jj = idx(sub2ind(size(idx), (1:size(idx,1))', j));
        % allowed contact: the two nets bridged by a net tie, near the tie
        pr = sort([repmat(n,numel(jj),1), B(jj,4)], 2);
        isTie = ismember(pr, tiePairs, 'rows');
        dTie = min(pdist2(A(:,1:2), tieXY(tieXY(:,3)==L,1:2)), [], 2);
        ok = isTie & dTie < 2.5;
        % J1/J2 lead nets touch the pads they connect to
        cm(ok) = inf;
        bad = cm < 0.1;
        viol = [viol; A(bad,1:2), repmat(n,nnz(bad),1), B(jj(bad),4), cm(bad), repmat(L,nnz(bad),1), B(jj(bad),1)]; %#ok<AGROW>
        [m, k] = min(cm);
        if m < best
            best = m; bw = [A(k,1:2), n, B(jj(k),1:2), B(jj(k),4)];
        end
    end
    R.minClear(L) = best; R.worst{L} = bw;
end
R.violations = viol;     % x y net otherNet clearance layer otherX
R.circumference = circ; R.seamGap = circ - (213.7035 - 10);
end
