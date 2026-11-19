# 12 cm Halbach RF coil: variable-pitch litz solenoid on a printed former

72 mT Halbach, B0 along x, bore along z → proton frequency **3.0656 MHz**. A solenoid along the bore
gives B1 along z, which is already perpendicular to B0, so no birdcage or saddle is needed.

Magnet frame as in `designs/12cm`: **x = B0, z = bore**, theta measured from +x toward +y. Lengths in mm.

## The design (12 turns)

| | |
|---|---|
| Former | ID 120, OD 130, printed; fits inside the Z gradient former (ID 132) |
| Wire | 1000 × 48 AWG litz, 2.5 OD, in a groove on the outside, flush (groove 2.5 deep, ~2.6 wide), 2.5 of plastic under it |
| Turns | 12, ring centres at z = ±66.17, ±62.87, ±59.57, ±56.27, ±27.14, ±8.58 (`turns.csv`) |
| Pitch | 3.3 (minimum: 0.8 plastic between grooves) in the 4-turn end clusters, 17–29 in the middle |
| Length | 134.8 over the litz, leads to ±70.2 |
| Litz | 4.82 m including leads |
| L / R / Q | 12.9 µH, 0.226 Ω at f0 (litz dc 0.110), unloaded Q ≈ 1100 |
| B1 | 64.2 µT per A at the centre (perpendicular to B0), 191 µT per √W of coil loss |
| Return | the same litz, inside the bore at theta = 90°, half buried in the inner wall (centre r = 60) |
| Tuning | 208 pF around the loop, in 3 equal shares: **C1, C2 in the winding = 625 pF each (2 × ~313 pF 1206 C0G in parallel)**, the 3rd share is the port's tune/match network; 83 V peak per A peak across each |

B1 uniformity of the real wire path (litz, jogs and return), `(max−min)/(max+min)` over a sphere:

| Sphere diameter | 60 | 80 | 100 | 110 | 120 |
|---|---|---|---|---|---|
| spread | 0.46 % | 0.97 % | 2.5 % | 4.6 % | — |

The edge of the 120 mm sphere is 3.75 mm from the turns, so there each turn's own field dominates. That
can't be fixed by moving turns. The return also runs along its surface. The optimisation therefore targets a 110 mm sphere.

## Why 12 turns

`sweep.txt` gives the optimum turn positions for 4–30 turns and lengths of 100–160 mm. `candidates.txt`
evaluates the real paths. The goal was the least litz that still gives a uniform field:

| Turns | Length | Litz | B1 per √W | spread, 100 / 110 mm sphere |
|---|---|---|---|---|
| 8 | 114 | 3.3 m | 172 µT | 6.2 / 9.8 % |
| 10 | 128 | 4.1 m | 183 µT | 4.2 / 6.6 % |
| **12** | **135** | **4.9 m** | **191 µT** | **2.5 / 4.4 %** |
| 16 | 153 | 6.5 m | 206 µT | 1.4 / 2.3 % |
| 19 | 160 | 7.7 m | 219 µT | 0.8 / 1.5 % |

Uniformity stops improving quickly past about 12 turns. Going from 8 to 19 turns more than doubles the litz
for only 27 % more B1 per √W. With 13 turns or fewer the best coil is shorter than 160 mm, because the optimizer
doesn't use the extra length. Other turn counts: `N=10; export_rf_solenoid`.

The turns cannot be evenly spaced. The optimum puts 4 turns at the minimum pitch at each end, which
replaces the field a longer coil would add, and spreads the rest out in the middle.

## How it's wound

See `figures/unrolled.png` (the former unrolled at true scale, centred on the seam) and `figures/B1_maps.png`
(B1 in the x-z and y-z planes, with the litz drawn to scale).

- **Rings and jogs.** Each turn is a ring at constant z. The step to the next turn happens in a 30° "jog" (33 mm
  of arc), raised-cosine in z and centred on the **seam at theta = 90°** (+y). All jogs are at the seam, so the
  grooves elsewhere are plain circles. Neighbouring grooves are at least 3.3 apart everywhere, jogs included.
  Winding goes counter-clockwise seen from +z, from the +z end (litz start at theta 105°, z 70.2) to the −z end
  (litz end at theta 75°, z −70.2).
- **Return inside the bore.** A solenoid's two ends are at opposite ends of the coil, so one has to come back. The
  litz goes through a hole at the −z end (theta 75°, z −70.2) and runs along z = −70.2 to theta 90°. From there it
  runs axially up the inner wall to z +70.2, then along z = +70.2 to theta 99.6°, and out through a hole next to the
  litz start. It lies in a 1.25 deep groove in the inner wall (centre r = 60): it sticks 1.25 into the bore and
  leaves 1.25 of plastic to the outer grooves. It sits at theta = 90° because there, close to it, its field points
  along x = B0, which does nothing in RF. Under the seam it also cancels the jogs' axial current. At theta = 0 the same
  return gives 5.5 % instead of 4.6 % over the 110 mm sphere. A copper foil under the seam
  (`'return_type','strip'` in `build_rf_solenoid_path`) gives 4.4 %.
- **Port** at the +z end, z 70.2, between the litz start (theta 105°) and the return coming out (theta 99.6°),
  5.4 mm apart. The tune/match board connects here.
- **Capacitors in the winding** (`caps.csv`): **C1** at theta 147°, z +27.1, and **C2** at theta 68°, z −27.1.
  Each is two 1206 side by side in a pocket 3.6 wide × 8 long along the wire, about 2.5 deep, with the litz ends
  soldered onto the cap terminals. These are the only rings with room for a 3.6 wide pocket, because the end
  clusters don't have it. That is why the three chunks are only roughly equal (port→C1 1.67, C1→C2 1.53,
  C2→port 1.80 m). Stacking the two caps on top of each other would fit in the groove width instead.

## Capacitors and voltages

The loop needs 208 pF in series in total. It is split in 3 equal shares so the three points see the same voltage:
C1 and C2 in the winding, and the port network (tune + match, about 625 pF equivalent). If the port network had
no reactance, C1 and C2 would each have to be 417 pF and would see 125 V per A instead of 83 V per A.
625 pF is for an unloaded coil in free space. Use for example **330 pF + 300 pF** C0G/NP0 1206 for C1 and C2,
put the trimming on the match board, and tune on the bench. The gradient coils 1 mm outside the former
(there is no RF shield) will lower L and shift the tuning.

Voltage across C1 or C2: 83 V peak per A peak. The two caps in parallel share the current, half each. A 100 µs hard 90° pulse needs 58.7 µT rotating, 117 µT
linear, so 1.8 A peak, about **150 V across C1 and C2** and 0.4 W of coil loss. A 50 µs pulse needs about 300 V. Use
caps rated **≥ 500 V**, or 1 kV for margin.

## Not modelled

The Q (~1100) and B1 per √W are coil-only: litz with proximity loss (Bessel model, `wire_ac_resistance.m`)
plus the foil. Cap ESR, solder joints, the sample, and coupling to the gradient coils will all lower Q.
Expect a few hundred in practice.

## Files

| File | What it is |
|---|---|
| `litz_path.csv` | Litz centreline, one row per point from the port around the coil and back through the return (`part`: lead / ring / jog / return), with depth below the OD |
| `turns.csv` | Ring z positions |
| `caps.csv` | Positions and values of C1 and C2 |
| `design_summary.txt` | Everything above, as numbers |
| `figures/unrolled.png`, `figures/B1_maps.png` | Unrolled former; B1 deviation maps in the x-z and y-z planes |
| `sweep.txt`, `candidates.txt` | Turn-count/length sweep (ring model), real-path candidates |
| `rf_solenoid_N12.mat` | Full result: path, return, caps, field, L, R |

Code: `Examples/rf_solenoid_sweep.m`, `rf_solenoid_candidates.m`, `export_rf_solenoid.m`. Sub-functions:
`rf_12cm_spec.m` (all inputs), `optimize_solenoid_turns.m` (fmincon minimax over the sphere, exact ring
fields from `loop_field_rz.m`), `build_rf_solenoid_path.m`, `evaluate_rf_solenoid.m`.
