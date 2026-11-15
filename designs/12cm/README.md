# 12 cm Halbach gradient set: wound coils on printed formers

Three gradient coils, each one layer of 14 AWG magnet wire laid in grooves on the
**outside** of a 3D printed cylindrical former, nested inside each other. The insides of the
formers are plain plastic. Imaging volume: a 120 mm diameter sphere.

## The design

Magnet frame: **x = B0 (Halbach field direction), z = bore axis**, y completes a right-handed
frame. All lengths in mm, copper at 20 °C unless noted.

| Coil | Former ID / OD | Wire-centre dia. | Length | Turns | Wire | mT/m/A | R20 | R80 | L | mT/m per √W |
|---|---|---|---|---|---|---|---|---|---|---|
| Z (bore, inner) | 132 / 140 | 138.3 | 330 | 38 | 29.1 m | 1.94 | 0.24 Ω | 0.30 Ω | 153 µH | 3.96 |
| Y (middle) | 142 / 150 | 148.3 | 324 | 28 | 29.5 m | 3.06 | 0.24 Ω | 0.30 Ω | 153 µH | 6.19 |
| X (B0, outer) | 152 / 160 | 158.3 | 318 | 30 | 37.1 m | 3.25 | 0.31 Ω | 0.38 Ω | 237 µH | 5.86 |

Wire: 14 AWG heavy build, 1.628 mm bare and at most 1.732 mm over the enamel (NEMA MW1000).
Neighbouring wires are at least 0.8 mm apart, a 2.53 mm pitch at the densest place.
Turns are CoilGen contour levels; wire lengths exclude the leads from the terminals to the amplifier.
For example, 10 mT/m takes (10/3.96)² = 6.4 W and 5.2 A on Z.

Maximum nonlinearity of the real wire, `max|Bz − G·r − c| / (G·R)` over a sphere of radius R, in %:

| Sphere diameter | 60 | 80 | 100 | 110 | 120 mm |
|---|---|---|---|---|---|
| Z | 0.19 | 0.50 | 1.20 | 1.82 | 3.61 |
| Y | 0.07 | 0.19 | 0.46 | 0.69 | 1.01 |
| X | 0.06 | 0.11 | 0.23 | 0.35 | 0.59 |

The worst error is at the very edge of the sphere. There it is only 9 mm from the wire, and one
turn's own field is a large part of the gradient field. At 10 mT/m and 45 mT, Z's worst edge error
(22 µT, 2.2 mm of image shift) is about the same as 500 ppm of B0 (22.5 µT, 2.25 mm). Everything
else is several times smaller. It is a fixed distortion that can be corrected in reconstruction
from these wire paths.

## Files

| File | What it is |
|---|---|
| `wire_paths/{Z,Y,X}_wire_path.csv` | **Wire centreline for CAD**, magnet frame, mm (see below) |
| `design_summary.txt` | Per coil: former, groove and slot depths, crossings, terminals, z extent, R, L, efficiency |
| `results/{Z,Y,X}_coil.mat` | Full result (`res`): CoilGen layout, wire path in the CoilGen frame, field check, spacing, R, L, AC resistance |
| `figures/*_unrolled.png` | Each former unrolled at true scale, crossings (leads) in black, terminals in yellow, with a 1:1 zoom of the densest 50×30 mm (wire drawn at its real diameter) |
| `figures/nested_coils.png` | The three coils nested, with the 120 mm sphere |
| `figures/Z_crossing_closeup.png` | One crossing of Z in 3D and cut open, wire at its real diameter |
| `figures/*_contours_2D.png`, `*_contours_3D.png` | CoilGen's own plots of the stream function and contours |
| `figures/winding_summary.txt` | Turns, crossings, bends tighter than 5 mm radius |
| `all_layouts_0p8mm.txt` | Every wire layout tried at 0.8 mm (the chosen three are in it) |
| `length_study_1mm.txt` | Best gradient per watt per linearity limit over coil length (1 mm spacing, 14 and 16 AWG) |

### Wire path CSV

One row per point along the wire, in order. The first and last rows are the two terminals.
`design_summary.txt` gives where they are and how far apart.

| Column | Meaning |
|---|---|
| `x_mm, y_mm, z_mm` | Wire centre, magnet frame (x = B0, z = bore) |
| `r_mm`, `theta_deg` | Radius, and angle from +x toward +y (0–360) |
| `centre_depth_mm` | Wire centre below the former's outer surface: 0.87 normally (wire flush with the surface), 2.25 on a lead, −0.4+0.87 = 0.47 on a lifted turn |
| `crossing` | 1 on a lead: wire that passes **under** other turns |
| `lifted_mm` | How far a turn is raised above its normal position where it passes over a lead (0–0.4, ramping over 3 mm) |

### Crossings

A gradient lobe needs all its turns to circulate the same way, so its inner end is enclosed by its
own turns. With one layer of wire, that end has to pass under the turns to get out: every lobe needs
at least one crossing. (A single closed curve on a plane or cylinder can go around a point only once.)
Z has 7 crossings, Y 15 (8 of them pass under only 1–2 turns where a lobe's inner loop exits), X 7.

They are built so that both wires stay within the 4 mm walls (`stack_crossings.m`):

- The lead lies in an **open slot**, laid in from outside before the turns: no tunnels, nothing threaded.
  Lead centre 2.25 mm below the surface, slot 3.11 mm deep, **0.89 mm of plastic** left inside.
- The lead touches the turns above it, enamel on enamel (0.05 mm nominal gap). The enamel is rated for
  kV; neighbouring wires differ by a few volts.
- Only the short pieces of turn over a lead sit **0.4 mm proud** of the surface, ramping back to flush
  over 3 mm. That uses 0.4 mm of the 1 mm gap to the next former. X has no former outside it, so
  whatever is outside 80 mm needs 0.4 mm clearance at those spots (`lifted_mm` > 0).
- The extra heat where wires stack is negligible: each turn passes over a lead for about 2 mm.
  At 10 mT/m on Z (5.2 A) the wire makes 0.22 W/m, so each overlap adds about 0.5 mW.

Field, resistance and inductance were evaluated again after this step (the table above is final).

### Terminals

**All terminals are on the outermost turn at the +z end of the coil** (Z at θ ≈ 0°, Y at θ ≈ 45°,
X at θ ≈ 185°), 3–7 mm in from the former end. The leads leave axially without crossing any
turns. CoilGen's wire is one series loop opened at a single cut, and for Y and X that cut was mid-former.
`move_wire_terminals.m` closed it with a straight jumper (5.9 mm on Y, 3.2 mm on X, at least 2.5 mm from
any other turn) and opened the wire at the outermost +z point instead, with a 3.2–3.6 mm gap between the
terminals. Field, resistance and inductance were evaluated again afterwards and changed by less than 0.1 %.
The yellow dots in `figures/*_unrolled.png` mark the terminals. Run the two leads as a twisted pair.

## How the design was chosen

Code is in `Examples/` and `sub_functions/` on branch `12cm`. CoilGen works in a frame with
B0 = z_cg and bore = x_cg. The magnet frame is x = z_cg, y = −y_cg, z = x_cg.

1. **Inputs.** Formers at radii 66–70 (Z), 71–75 (Y) and 76–80 mm (X); maximum lengths 330, 324 and 318 mm.
   The wire sits in a groove flush with the outside, so its centre is at the outer radius minus half a wire diameter.
   The goal was the most gradient per watt with as good linearity as possible over a 120 mm sphere.
   Z goes on the inner former because it is the least efficient axis.
2. **Ideal current patterns over length** (`halbach_12cm_sf_sweep.m`, `sweep_sf_tradeoff.m`).
   For each axis and length (140–340 mm in 20 mm steps, plus 318, 324 and 330), CoilGen builds the
   cylinder mesh and the field and resistance matrices for a 1 T/m target on 2500 points on the 120 mm sphere.
   The stream function is then solved for 36 power weights at once, by SVD of a power-weighted
   least-squares fit (power x'Rx rather than CoilGen's |Rx|²). Each solution is checked on 4500 different points.
   About 20 s per length; two MATLAB processes with 6 threads each, under 8 GB.
3. **Gradient per watt for round wire.** The number of turns is set by the densest place, where
   neighbouring wires are one pitch apart. That gives:
   efficiency = 1/Δψ, wire length = ∫|K|dA/Δψ, and gradient per √W = 1/√(Ω/m · Δψ · ∫|K|dA), with Δψ = K_max · pitch.
   Results (`length_study_1mm.txt`):
   - **Longer is better per watt.** With the sphere this close to the coil, short coils need large end
     currents. Z improves up to the maximum length; Y and X level off from about 280–300 mm.
   - **14 AWG beats 16 AWG** by 15–18 % per watt on every axis, because more of the groove is copper.
4. **Real wire layouts** (`halbach_12cm_wound_set.m`, `design_wound_coil.m`). For each candidate, the
   stream function with the most gradient per watt within a linearity limit (0.05–3 %) is passed to
   CoilGen. It cuts the contours, joins them into one wire and moves the crossings 2.03 mm into the former.
   The wire is then checked with exact Biot–Savart, for spacing (no two parts closer than one wire
   diameter), and for inductance and resistance.
   - The ideal patterns can reach 0.05 % over the sphere, but the real wire can't. At the sphere's edge,
     linearity is limited by the number of turns, so a slightly looser pattern with more turns does better.
     For example, Z's 0.05 % pattern gave 17 % with real wire, while the 3 % pattern gave 4.4 %
     (3.6 % after the crossing step below).
   - The plastic between wires went from 1 mm to 0.8 mm, which raised gradient per watt by 4–6 %.
5. **Chosen:** Z 330 mm with the 3 % pattern; Y 324 mm with the 1 % pattern (the 0.2 % one had wire
   collisions); X 318 mm with the 0.2 % pattern. A Z alternative with the 2 % pattern gets 0.87 %
   instead of 1.20 % at 100 mm for 7 % less gradient per √W (`Z_L330_AWG14_lin200` in `all_layouts_0p8mm.txt`).
6. **Export** (`export_halbach_12cm_designs.m`): terminals moved to the +z end (`move_wire_terminals.m`),
   then the crossings rebuilt with touching wires and lifted turns (`stack_crossings.m`), then the field,
   spacing, R and L evaluated again. After this step no two wire parts are closer than one wire
   diameter (1.72–1.75 mm centre to centre where they touch at crossings).

## Regenerating

From `Examples/`, in MATLAB:

```
run('halbach_12cm_sf_sweep.m')                                   % ideal patterns (SF_JOB="Z,0.330" for one)
setenv('WOUND_CLEARANCE_MM','0.8'); run('halbach_12cm_wound_set.m') % layouts (WOUND_JOB="Z,330,AWG14,3.00")
run('halbach_12cm_wound_summary.m')                              % linearity vs sphere size, collisions
run('export_halbach_12cm_designs.m')                             % this folder (terminals, crossings)
run('view_halbach_12cm_coils.m')                                 % figures in this folder
run('view_halbach_12cm_crossing.m')                              % crossing close-up
```

Intermediate results go to `Wound_gradient_coils/halbach_12cm/` (not in git, about 200 MB).
Results that already exist are not recomputed.
