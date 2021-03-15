within OU44.Components;
model AirMix_new
  "Calculates heat flow based on indoor temperature and ventilation airflow and temperature"

  Modelica.Blocks.Sources.Constant h_to_s(k=1/3600)
    annotation (Placement(transformation(extent={{-64,58},{-44,78}})));
  Modelica.Blocks.Sources.Constant AirRhoCp(k=1.2*1005)
    annotation (Placement(transformation(extent={{-64,26},{-44,46}})));
  Modelica.Blocks.Math.MultiProduct multiProduct(      significantDigits=8, nu=4)
    annotation (Placement(transformation(extent={{-10,22},{14,46}})));
  Modelica.Blocks.Interfaces.RealInput Tve
    annotation (Placement(transformation(extent={{-128,14},{-88,54}})));
  Modelica.Blocks.Interfaces.RealInput Ti annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-108})));
  Modelica.Blocks.Interfaces.RealInput Vve "[m3/h]"
    annotation (Placement(transformation(extent={{-130,-70},{-90,-30}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{48,-10},{68,10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  Modelica.Blocks.Math.Add add1(k1=+1, k2=-1)
    annotation (Placement(transformation(extent={{-54,-10},{-34,10}})));
equation
  connect(h_to_s.y, multiProduct.u[1]) annotation (Line(points={{-43,68},{-22,
          68},{-22,40.3},{-10,40.3}},
                                    color={0,0,127}));
  connect(AirRhoCp.y, multiProduct.u[2]) annotation (Line(points={{-43,36},{-43,
          36.1},{-10,36.1}},   color={0,0,127}));
  connect(Vve, multiProduct.u[3]) annotation (Line(points={{-110,-50},{-18,-50},
          {-18,31.9},{-10,31.9}},   color={0,0,127}));
  connect(multiProduct.y, prescribedHeatFlow.Q_flow) annotation (Line(points=
          {{16.04,34},{32,34},{32,0},{48,0}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, port_b)
    annotation (Line(points={{68,0},{74,0},{100,0}}, color={191,0,0}));
  connect(Tve, add1.u1) annotation (Line(points={{-108,34},{-84,34},{-84,6},{
          -56,6}}, color={0,0,127}));
  connect(add1.y, multiProduct.u[4]) annotation (Line(points={{-33,0},{-22,0},{
          -22,27.7},{-10,27.7}}, color={0,0,127}));
  connect(add1.u2, Ti) annotation (Line(points={{-56,-6},{-62,-6},{-62,-40},{0,
          -40},{0,-108}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          lineThickness=0.5),
        Line(points={{-68,46},{-62,28},{-56,46}}, color={28,108,200}),
        Line(points={{-78,46},{-78,26}}, color={28,108,200}),
        Line(points={{-74,46},{-82,46}}, color={28,108,200}),
        Line(points={{-84,-38},{-78,-56},{-72,-38}}, color={28,108,200}),
        Line(points={{-68,-38},{-62,-56},{-56,-38}}, color={28,108,200}),
        Line(points={{-8,-62},{-8,-82}}, color={28,108,200}),
        Line(points={{-4,-62},{-12,-62}}, color={28,108,200}),
        Line(points={{2,-62},{2,-82}}, color={28,108,200}),
        Line(points={{0,-50},{0,0},{90,0}}, color={28,108,200}),
        Line(points={{-50,36},{-40,36},{-20,0},{0,0}}, color={28,108,200}),
        Line(points={{-50,-48},{-40,-48},{-20,0}}, color={28,108,200}),
        Ellipse(
          extent={{-4,4},{4,-4}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}));
end AirMix_new;
