function G = coil_geometry(C)
%COIL_GEOMETRY  Map the flat PCB onto the 65 mm former (all lengths in m).
%
%   PCB x wraps around the cylinder, PCB y is the bore axis.  Mapping uses
%   the CoilGen header in the board file (x0, y0, th0) and agrees with all
%   the silkscreen: +X at PCB x~62, -X at PCB x~164, "+Z into bore -->"
%   pointing toward the J1/J2 tab end, and F.Cu facing outward:
%       theta = (x - x0)/rMid - th0,   z = y0 - y
%   Magnet frame: B0 along +x, bore along z.
%
%   Stack-up rolled onto a 65.0 mm OD former, B.Cu toward the former:
%     former | B.Mask 50 | B.Cu 35 | PI 25 | F.Cu 35 | F.Mask 50 (um)

G.x0 = 11.25e-3; G.y0 = 81.5e-3; G.th0 = 1.61695017;
G.rFormer = 32.5e-3;
G.rB = G.rFormer + 50e-6 + 17.5e-6;          % B.Cu centre
G.rF = G.rB + 60e-6;                         % F.Cu centre
G.rMid = G.rFormer + 97.5e-6;                % neutral plane of the flex
G.I = C.I; G.layer = C.layer; G.width = C.width*1e-3;
G.p1 = C.p1*1e-3; G.p2 = C.p2*1e-3;
G.vias = C.vias*1e-3; G.J1 = C.J1*1e-3; G.J2 = C.J2*1e-3;
G.toTheta = @(x) (x - G.x0)/G.rMid - G.th0;
G.toZ = @(y) G.y0 - y;
end
