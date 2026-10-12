function F = coil_filaments(G, maxLen)
%COIL_FILAMENTS  Straight PCB tracks -> chains of short 3-D chords on the
%rolled cylinder (each layer at its own copper radius). Track centre lines.
if nargin < 2, maxLen = 0.5e-3; end
L = vecnorm(G.p2 - G.p1, 2, 2);
n = max(1, ceil(L/maxLen));
idx = repelem((1:numel(L))', n);
k = cumsum(ones(size(idx))) - repelem(cumsum(n) - n, n);    % 1..n within track
s0 = (k-1)./n(idx); s1 = k./n(idx);
a = G.p1(idx,:) + s0.*(G.p2(idx,:) - G.p1(idx,:));
b = G.p1(idx,:) + s1.*(G.p2(idx,:) - G.p1(idx,:));
r = G.rB*ones(size(idx)); r(G.layer(idx)==1) = G.rF;
F.A = to3d(a, r); F.B = to3d(b, r);
F.I = G.I(idx); F.layer = G.layer(idx); F.track = idx;

    function P = to3d(p, r)
        th = G.toTheta(p(:,1));
        P = [r.*cos(th), r.*sin(th), G.toZ(p(:,2))];
    end
end
