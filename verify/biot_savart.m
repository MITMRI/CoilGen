function B = biot_savart(A, Bp, I, P)
%BIOT_SAVART  Exact field of straight finite current filaments A->Bp
%(N x 3, m, current I) at points P (M x 3). Returns B (M x 3) in tesla.
mu0 = 4*pi*1e-7;
B = zeros(size(P,1), 3);
L = Bp - A; Lh = L ./ vecnorm(L,2,2);
chunk = max(1, floor(3e6 / size(A,1)));   % ~25 MB per M x N temp
for i0 = 1:chunk:size(P,1)
    ii = i0:min(i0+chunk-1, size(P,1));
    px = P(ii,1); py = P(ii,2); pz = P(ii,3);
    % r1 = P - A, r2 = P - B  (M x N)
    r1x = px - A(:,1)'; r1y = py - A(:,2)'; r1z = pz - A(:,3)';
    r2x = px - Bp(:,1)'; r2y = py - Bp(:,2)'; r2z = pz - Bp(:,3)';
    n1 = sqrt(r1x.^2 + r1y.^2 + r1z.^2); n2 = sqrt(r2x.^2 + r2y.^2 + r2z.^2);
    % L x r1
    cx = L(:,2)'.*r1z - L(:,3)'.*r1y;
    cy = L(:,3)'.*r1x - L(:,1)'.*r1z;
    cz = L(:,1)'.*r1y - L(:,2)'.*r1x;
    c2 = cx.^2 + cy.^2 + cz.^2;
    c2(c2 < 1e-24) = inf;                      % on the filament axis
    f = (I' * mu0/(4*pi)) .* ...
        ((Lh(:,1)'.*r1x + Lh(:,2)'.*r1y + Lh(:,3)'.*r1z)./n1 - ...
         (Lh(:,1)'.*r2x + Lh(:,2)'.*r2y + Lh(:,3)'.*r2z)./n2) ...
        .* vecnorm(L,2,2)' ./ c2;
    B(ii,:) = [sum(f.*cx,2), sum(f.*cy,2), sum(f.*cz,2)];
end
end
