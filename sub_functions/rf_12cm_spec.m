function spec=rf_12cm_spec()
%Specification of the 12 cm Halbach RF solenoid (see designs/rf_12cm/README.md)
spec.B0=0.072;                          %T, along x
spec.f0=42.577478e6*spec.B0;            %Hz, 1H Larmor frequency
spec.former_id=0.120; spec.former_od=0.130;
spec.litz_od=0.0025;                    %1000 x 48 AWG litz, outer diameter
spec.litz_strands=1000;
spec.litz_strand_radius=0.0315e-3/2;    %48 AWG bare
spec.litz_lay_factor=1.03;
spec.a=spec.former_od/2-spec.litz_od/2; %turn radius: litz flush with the outer surface
spec.pitch_min=spec.litz_od+0.0008;     %0.8 mm plastic between turns
spec.roi_radius=0.060;                  %optimise B1 over the 120 mm imaging sphere (touches the inner wall)
spec.num_breaks=4;                      %equal shares of the loop: the port network + 3 caps in the winding (2 x 1206 in parallel each)
end
