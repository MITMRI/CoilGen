function S = current_sheet(G, h, sigma)
%CURRENT_SHEET  Surface current density K (A/m) of each copper layer on a
%regular grid in flat-PCB coordinates (x = arc length, y = bore axis).
%   Current is spread uniformly across each track's actual width, then
%   lightly Gaussian-smoothed (sigma) so it can be sampled on the FE mesh.
if nargin < 2, h = 0.2e-3; end
if nargin < 3, sigma = 0.4e-3; end
pad = 3e-3;
lo = min([G.p1; G.p2]) - pad; hi = max([G.p1; G.p2]) + pad;
S.xg = lo(1):h:hi(1); S.yg = lo(2):h:hi(2);
nx = numel(S.xg); ny = numel(S.yg);
S.Kx = zeros(ny, nx, 2); S.Ky = zeros(ny, nx, 2);
ds = h/3;                                    % sub-sample spacing
for L = 1:2
    m = find(G.layer == L);
    X = cell(numel(m),1); Y = X; W = X;
    for k = 1:numel(m)
        i = m(k);
        a = G.p1(i,:); d = G.p2(i,:) - a; len = norm(d);
        if len == 0, continue; end
        t = d/len; nrm = [-t(2) t(1)]; w = G.width(i);
        su = ((1:ceil(len/ds)) - 0.5)/ceil(len/ds);
        sv = ((1:ceil(w/ds)) - 0.5)/ceil(w/ds) - 0.5;
        [uu, vv] = meshgrid(su*len, sv*w);
        X{k} = a(1) + uu(:)*t(1) + vv(:)*nrm(1);
        Y{k} = a(2) + uu(:)*t(2) + vv(:)*nrm(2);
        % each sub-sample carries I*t * (dl*dw)/w  [A*m]
        wt = G.I(i)*(len/numel(su))*(w/numel(sv))/w;
        W{k} = repmat(wt*t, numel(uu), 1);
    end
    X = vertcat(X{:}); Y = vertcat(Y{:}); W = vertcat(W{:});
    ix = round((X - S.xg(1))/h) + 1; iy = round((Y - S.yg(1))/h) + 1;
    S.Kx(:,:,L) = accumarray([iy ix], W(:,1), [ny nx]) / h^2;
    S.Ky(:,:,L) = accumarray([iy ix], W(:,2), [ny nx]) / h^2;
    S.Kx(:,:,L) = imgaussfilt(S.Kx(:,:,L), sigma/h);
    S.Ky(:,:,L) = imgaussfilt(S.Ky(:,:,L), sigma/h);
end
S.h = h;
end
