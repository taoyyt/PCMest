within OU44.Components;
model FanFit "Fan electric energy consumption curve fit"

  parameter Real a "Parameter in y=ax^2+bx+c";
  parameter Real b "Parameter in y=ax^2+bx+c";
  parameter Real c "Parameter in y=ax^2+bx+c";

  Modelica.Blocks.Interfaces.RealInput VFR
    annotation (Placement(transformation(extent={{-128,-20},{-88,20}})));
  Modelica.Blocks.Interfaces.RealOutput qel
    annotation (Placement(transformation(extent={{96,-10},{116,10}})));
equation
  qel = a * VFR^2 + b * VFR + c;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{2,28},{56,-26}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.None),
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,127}),
        Line(
          points={{-68,40},{-50,-16},{-12,-46},{68,-52}},
          color={28,108,200},
          smooth=Smooth.Bezier),
        Line(points={{-80,50},{-80,-70},{80,-70}}, color={28,108,200}),
        Line(points={{-80,50},{-82,42}}, color={28,108,200}),
        Line(points={{-80,50},{-78,42}}, color={28,108,200}),
        Line(points={{80,-70},{74,-68}}, color={28,108,200}),
        Line(points={{80,-70},{74,-72}}, color={28,108,200}),
        Polygon(
          points={{12,22},{12,-20},{56,0},{12,22}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}),                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end FanFit;
