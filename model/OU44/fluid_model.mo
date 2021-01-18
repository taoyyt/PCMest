within OU44;
model fluid_model
  Buildings.Fluid.MixingVolumes.MixingVolume vol
    annotation (Placement(transformation(extent={{-40,-12},{-20,8}})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary
    annotation (Placement(transformation(extent={{-72,-34},{-52,-14}})));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end fluid_model;
