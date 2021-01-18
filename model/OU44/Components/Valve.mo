within OU44.Components;
model Valve
  "Generic linear valve. The model translates the valve position [%] to heat or air flow."

  Modelica.Blocks.Interfaces.RealInput inp "Actuator position [%]"
    annotation (Placement(transformation(extent={{-126,-20},{-86,20}})));
  Modelica.Blocks.Interfaces.RealOutput out "Heat [W] or airflow [m/h]"
    annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  Modelica.Blocks.Sources.Constant capacityMax(k=capacity)
    "\"Maximum heat [W] or airflow [m3/h]\""
    annotation (Placement(transformation(extent={{-64,-36},{-44,-16}})));

  Modelica.Blocks.Sources.Constant FromPercentsToFraction(k=0.01)
    annotation (Placement(transformation(extent={{-64,16},{-44,36}})));
  parameter Real capacity=100 "Maximum heat [W] or airflow [m3/h]";
  Modelica.Blocks.Math.MultiProduct multiProduct(nu=3)
    annotation (Placement(transformation(extent={{-6,-16},{26,16}})));
equation
  connect(multiProduct.u[1], FromPercentsToFraction.y) annotation (Line(
        points={{-6,7.46667},{-24,7.46667},{-24,26},{-43,26}}, color={0,0,127}));
  connect(inp, multiProduct.u[2]) annotation (Line(points={{-106,0},{-6,0},{
          -6,8.88178e-016}}, color={0,0,127}));
  connect(multiProduct.u[3], capacityMax.y) annotation (Line(points={{-6,
          -7.46667},{-24,-7.46667},{-24,-26},{-43,-26}}, color={0,0,127}));
  connect(multiProduct.y, out)
    annotation (Line(points={{28.72,0},{106,0}},         color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          lineThickness=0.5), Line(
          points={{80,80},{-80,-80}},
          color={28,108,200},
          thickness=0.5)}));
end Valve;
