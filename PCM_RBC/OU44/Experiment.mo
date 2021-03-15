within OU44;
model Experiment
  Modelica.Blocks.Sources.Constant const1(k=50)
    annotation (Placement(transformation(extent={{-70,-26},{-50,-6}})));
  Modelica.Blocks.Sources.Constant const2(k=21)
    annotation (Placement(transformation(extent={{-70,-56},{-50,-36}})));
  ZoneR4C3 zoneR4C3_1(
    CO2pp=0.041526774,
    shgc=0.726811677,
    Cair=61313.01535,
    Cin=100.3197209,
    Cwall=292280.9472,
    RExt=0.004107916,
    Ri=9.996425197,
    Rin=0.481974341,
    Rinf=9.684514317,
    Tve=21.28598152,
    Vinf=2432.499992,
    maxVent=1200,
    occheff=0,
    Vi=486.5,
    CO2n=461.4400024,
    maxHeat=2689)
    annotation (Placement(transformation(extent={{-14,-6},{12,18}})));
  Modelica.Blocks.Sources.Constant const5(k=20)
    annotation (Placement(transformation(extent={{-72,8},{-52,28}})));
  Modelica.Blocks.Sources.Constant const6(k=20)
    annotation (Placement(transformation(extent={{-70,36},{-50,56}})));
  Modelica.Blocks.Sources.Constant const7(k=200)
    annotation (Placement(transformation(extent={{-72,68},{-52,88}})));
equation
  connect(const7.y, zoneR4C3_1.solrad) annotation (Line(points={{-51,78},{-34,
          78},{-34,16.8},{-14.8,16.8}}, color={0,0,127}));
  connect(const6.y, zoneR4C3_1.Tout) annotation (Line(points={{-49,46},{-36,46},
          {-36,14.1},{-14.9,14.1}}, color={0,0,127}));
  connect(const5.y, zoneR4C3_1.Occ) annotation (Line(points={{-51,18},{-38,18},
          {-38,9.8},{-15,9.8}}, color={0,0,127}));
  connect(const1.y, zoneR4C3_1.dpos) annotation (Line(points={{-49,-16},{-32,
          -16},{-32,0.5},{-14.9,0.5}}, color={0,0,127}));
  connect(const2.y, zoneR4C3_1.Tsetp) annotation (Line(points={{-49,-46},{-30,
          -46},{-30,-2.5},{-14.9,-2.5}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Experiment;
