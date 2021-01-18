within OU44;
model test_ou44
  ZoneR2C2_SJH zoneR2C2_SJH
    annotation (Placement(transformation(extent={{-14,-16},{30,28}})));
  Modelica.Blocks.Sources.Constant Tvestp(k=21)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Modelica.Blocks.Sources.Pulse solrad(
    amplitude=400,
    width=33,
    period=86400,
    startTime=28800)
    annotation (Placement(transformation(extent={{-80,58},{-60,78}})));
  Modelica.Blocks.Sources.Pulse Tout(
    width=33,
    period=86400,
    startTime=28800,
    offset=3,
    amplitude=7)
    annotation (Placement(transformation(extent={{-80,22},{-60,42}})));
  Modelica.Blocks.Sources.Pulse verate(
    period=86400,
    startTime=28800,
    offset=0,
    width=20,
    amplitude=0.4)
    annotation (Placement(transformation(extent={{-80,-72},{-60,-52}})));
  Modelica.Blocks.Sources.Pulse qrad(
    period=86400,
    startTime=28800,
    offset=0,
    width=20,
    amplitude=0)
    annotation (Placement(transformation(extent={{-52,-94},{-32,-74}})));
  Modelica.Blocks.Sources.Pulse occ(
    width=33,
    period=86400,
    startTime=28800,
    amplitude=30,
    offset=0)
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
equation
  connect(Tvestp.y, zoneR2C2_SJH.Tvestp) annotation (Line(points={{-59,-30},{
          -36,-30},{-36,1.4},{-14.6,1.4}}, color={0,0,127}));
  connect(solrad.y, zoneR2C2_SJH.solrad) annotation (Line(points={{-59,68},{-38,
          68},{-38,24},{-14,24}}, color={0,0,127}));
  connect(Tout.y, zoneR2C2_SJH.Tout) annotation (Line(points={{-59,32},{-36,32},
          {-36,15.8},{-14,15.8}}, color={0,0,127}));
  connect(verate.y, zoneR2C2_SJH.verate) annotation (Line(points={{-59,-62},{
          -36,-62},{-36,-1.8},{-14.6,-1.8}}, color={0,0,127}));
  connect(qrad.y, zoneR2C2_SJH.qrad) annotation (Line(points={{-31,-84},{-22,
          -84},{-22,-11.2},{-14,-11.2}}, color={0,0,127}));
  connect(occ.y, zoneR2C2_SJH.occ) annotation (Line(points={{-59,0},{-36,0},{
          -36,7.6},{-14,7.6}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=604800));
end test_ou44;
