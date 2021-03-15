within OU44.Components;
model CO2test
  CO2 cO2_1 annotation (Placement(transformation(extent={{-6,-8},{14,12}})));
  Modelica.Blocks.Sources.Step step(height=300, startTime=3600)
    annotation (Placement(transformation(extent={{-78,18},{-58,38}})));
  Modelica.Blocks.Sources.Step step1(height=10, startTime=1800)
    annotation (Placement(transformation(extent={{-76,-32},{-56,-12}})));
equation
  connect(cO2_1.Vve, step.y) annotation (Line(points={{-6.4,8},{-30,8},{-30,28},
          {-57,28}}, color={0,0,127}));
  connect(step1.y, cO2_1.persons) annotation (Line(points={{-55,-22},{-30,-22},
          {-30,-4},{-6.4,-4}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end CO2test;
