# 12 cm Halbach RF coil: variable-pitch litz solenoid on a printed former

72 mT Halbach, B0 along x, bore along z → proton frequency **3.0656 MHz**. A solenoid along the bore
gives B1 along z, which is already perpendicular to B0, so no birdcage or saddle is needed.

Magnet frame as in `designs/12cm`: **x = B0, z = bore**, theta measured from +x toward +y. Lengths in mm.

The imaging volume is a **120 mm sphere that touches the former's inner wall** (ID 120). So nothing can go
inside the bore, and the edge of the sphere at the coil's midplane is only 3.75 mm from the turns.

## The design (18 turns)

| | |
|---|---|
| Former | ID 120, OD 130, printed; fits inside the Z gradient former (ID 132) |
| Wire | 1000 × 48 AWG litz, 2.5 OD, in a groove on the outside, flush (groove 2.5 deep, ~2.6 wide), 2.5 of plastic under it |
| Turns | 18, ring centres at z = ±70.20, ±66.90, ±63.60, ±60.30, ±57.00, ±53.70, ±30.30, ±16.78, ±5.05 (`turns.csv`) |
| Pitch | 6 turns at the 3.3 minimum (0.8 plastic between grooves) at each end; 23.4, 13.5, 11.7, 10.1 towards the middle |
| Length | 142.9 over the litz, leads to ±74.2 |
| Litz | 7.22 m including leads |
| Return | copper foil under the seam at theta = 90° |
| L / R / Q | 26.5 µH, 0.349 Ω at f0 (litz dc 0.164), unloaded Q ≈ 1460 |
| B1 | 95.8 µT per A at the centre (perpendicular to B0), 229 µT per √W of coil loss |
| Tuning | 102 pF around the loop, in **4 equal shares of 407 pF**: **C1, C2, C3 in the winding (2 × ~204 pF 1206 C0G in parallel, stacked)**, and the port's tune/match network; 127 V peak per A peak across each |

B1 (the field perpendicular to B0) of the real wire path, including the jogs and the return:

- **96.7 % of the 120 mm sphere's volume is within ±5 % of the centre value, and 99.8 % is within ±10 %.**
- Worst case over a sphere, `(max−min)/(max+min)`:

| Sphere diameter | 80 | 100 | 110 | 116 | 120 |
|---|---|---|---|---|---|
| spread | 0.68 % | 2.3 % | 4.1 % | 7.0 % | 13 % |

The worst points are in a thin shell against the wall, between turns. There each turn's own field
dominates, and moving turns can't fix that. For RF this matters much less than for gradients. A B1 that is
10 % off turns a 90° pulse into 99°, which loses about 1.2 % of the signal.

## Why 18 turns

`sweep.txt` gives the optimum turn positions, optimised over the 120 mm sphere, for 4–30 turns and lengths of 100–160 mm.
`candidates.txt` evaluates the real paths:

| Turns | Length | Litz | B1 per √W | worst over 110 / 120 mm | volume of 120 mm sphere within ±5 % / ±10 % |
|---|---|---|---|---|---|
| 12 | 121 | 4.8 m | 209 µT | 8.4 / 22 % | 89 / 98.0 % |
| 16 | 135 | 6.4 m | 226 µT | 5.7 / 15 % | 94 / 99.7 % |
| **18** | **143** | **7.2 m** | **230 µT** | **4.1 / 13 %** | **97 / 99.8 %** |
| 22 | 145 | 8.8 m | 253 µT | 5.1 / 9 % | 92 / 100 % |
| 28 | 160 | 11.2 m | 270 µT | 3.2 / 5.4 % | 99.9 / 100 % |

18 turns gave the best uniformity for the litz it uses. Other turn counts: `N=12; export_rf_solenoid`.

The turns cannot be evenly spaced. The optimum puts 6 turns at the minimum pitch at each end, which
replaces the field a longer coil would add, and spreads the rest out in the middle.

## How it's wound

See `figures/unrolled.png` (the former unrolled at true scale, centred on the seam) and `figures/B1_maps.png`
(B1 in the x-z and y-z planes, with the litz drawn to scale and the foil drawn edge-on).

- **Rings and jogs.** Each turn is a ring at constant z. The step to the next turn happens in a 30° "jog" (33 mm
  of arc), raised-cosine in z and centred on the **seam at theta = 90°** (+y). All jogs are at the seam, so the
  grooves elsewhere are plain circles. Neighbouring grooves are at least 3.3 apart everywhere, jogs included.
  Winding goes counter-clockwise seen from +z, from the +z end (litz start at theta 105°, z 74.2) to the −z end
  (litz end at theta 75°, z −74.2).
- **Return.** A solenoid's two ends are at opposite ends of the coil, so one has to come back. It can't go in the
  bore, because the sphere touches the wall. A 2.5 mm litz can't pass under or over the turns in a 5 mm wall
  either. The return is therefore **copper foil**, 6 wide, 2 oz (≥ 70 µm), insulated with Kapton. It sits in a slot
  whose floor is 2.65 below the OD (0.15 under the groove floor), so the jogs over it stay flush, with 2.35 of
  plastic left inside. Its path: a tab at z −74.2 from the litz end (theta 75°) to theta 90°, axially at theta 90°
  under the jogs up to z +74.2, then a tab to theta 99.6° at the port. It is 175 long and about 9 mΩ.
  Theta = 90° is where its near field points along B0, which does nothing in RF. Under the seam it also cancels the
  jogs' axial current.
- **Port** at the +z end, z 74.2, between the litz start (theta 105°) and the foil tab (theta 99.6°), 6 mm apart.
  The tune/match board connects here.
- **Capacitors in the winding** (`caps.csv`). **C1** is at theta 322°, z +57.0; **C2** at theta 168°, z −5.05;
  **C3** at theta 16°, z −57.0. They split the loop into 4 chunks of 1.85 m each (port→C1→C2→C3→port), equal to
  within 0.1 %. C1 and C3 sit in the end clusters, where the pitch is 3.3. There the two 1206 must be **stacked**
  (one soldered on top of the other, 1.6 wide, about 1.9 high for 0.8–0.95 mm thick parts). They then sit in the
  groove itself, widened to about 1.7 over 8 along the wire. The litz ends are soldered onto the stack's end
  terminals.

## Capacitors and voltages

The loop needs 102 pF in series in total. It is split in 4 equal shares so every point sees the same voltage:
C1, C2, C3 and the port network (tune + match, about 407 pF equivalent). 407 pF is the value for an unloaded coil in
free space. Use for example **2 × 200 pF** (or 220 + 180) C0G/NP0 1206 per cap, put the trimming on the match board,
and tune on the bench. The gradient coils 1 mm outside the former (there is no RF shield) will lower L and shift the
tuning.

Voltage across each cap: 127 V peak per A peak. The two caps in a stack share the current, half each.
A 100 µs hard 90° pulse needs 58.7 µT rotating, 117 µT linear, so 1.2 A peak, about **160 V per cap** and
0.26 W of coil loss. A 50 µs pulse needs about 310 V. Use caps rated **≥ 500 V**, or 1 kV for margin.

## How B1 is computed

The litz centreline is cut into straight pieces about 1.1 mm long (1° of arc). The foil is modelled as 5 parallel
filaments across its width. The field is summed with the exact finite-segment Biot-Savart formula (`biot_savart.m`).
This is quasi-static, which is valid because the wavelength (98 m) is far longer than the 7 m of litz and the caps
keep the current nearly uniform. B1 is `√(By² + Bz²)`; Bx is along B0 and does nothing in RF.
The turn positions are optimised with exact ring fields (elliptic integrals, `loop_field_rz.m`).

## Not modelled

The Q (~1460) and B1 per √W are coil-only: litz with proximity loss (Bessel model, `wire_ac_resistance.m`)
plus the foil. Cap ESR, solder joints, the sample, and coupling to the gradient coils and magnet will all lower Q.
Expect a few hundred in practice.

## Files

| File | What it is |
|---|---|
| `litz_path.csv` | Litz centreline, one row per point from the port around the coil (`part`: lead / ring / jog), then the foil centreline (`return`), with depth below the OD |
| `turns.csv` | Ring z positions |
| `caps.csv` | Positions and values of C1–C3 |
| `design_summary.txt` | Everything above, as numbers |
| `figures/unrolled.png`, `figures/B1_maps.png` | Unrolled former; B1 deviation maps in the x-z and y-z planes |
| `sweep.txt`, `candidates.txt` | Turn-count/length sweep (ring model), real-path candidates |
| `rf_solenoid_N18.mat` | Full result: path, return, caps, field, L, R |

Code: `Examples/rf_solenoid_sweep.m`, `rf_solenoid_candidates.m`, `export_rf_solenoid.m`. Sub-functions:
`rf_12cm_spec.m` (all inputs), `optimize_solenoid_turns.m` (fmincon minimax over the sphere, exact ring
fields from `loop_field_rz.m`), `build_rf_solenoid_path.m`, `evaluate_rf_solenoid.m`.
