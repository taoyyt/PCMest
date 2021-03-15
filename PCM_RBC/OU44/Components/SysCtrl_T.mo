within OU44.Components;
model SysCtrl_T

  Modelica.Blocks.Interfaces.RealInput Tstp
    annotation (Placement(transformation(extent={{-128,40},{-88,80}})));
  Modelica.Blocks.Interfaces.RealInput Tm annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-40,110})));
  Modelica.Blocks.Continuous.LimPID PID(
    yMax=100,
    yMin=0,
    y_start=0,
    xi_start=0,
    xd_start=0,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    controllerType=Modelica.Blocks.Types.SimpleController.PID,
    k=1,
    Ti=30,
    Td=60)  annotation (Placement(transformation(extent={{-30,50},{-10,70}})));
  Modelica.Blocks.Interfaces.RealOutput hPos
    annotation (Placement(transformation(extent={{96,50},{116,70}})));
equation
  connect(Tstp, PID.u_s)
    annotation (Line(points={{-108,60},{-32,60}}, color={0,0,127}));
  connect(Tm, PID.u_m) annotation (Line(points={{-40,110},{-40,40},{-20,40},{-20,
          48}}, color={0,0,127}));
  connect(PID.y, hPos)
    annotation (Line(points={{-9,60},{106,60}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={236,236,236},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,60},{-40,60}}, color={28,108,200}),
        Rectangle(extent={{-20,80},{60,40}}, lineColor={28,108,200}),
        Line(points={{-40,88}}, color={28,108,200}),
        Line(points={{-40,100},{-40,60}}, color={28,108,200}),
        Line(points={{-40,60},{-20,60}}, color={28,108,200}),
        Line(points={{60,60},{96,60}}, color={28,108,200}),
        Line(points={{28,0}}, color={28,108,200}),
        Polygon(
          points={{4,0},{-4,4},{-4,-4},{4,0}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          origin={-10,74},
          rotation=90),
        Line(points={{-10,70},{-10,42}},  color={0,0,0}),
        Line(points={{-18,50},{54,50}},   color={0,0,0}),
        Line(points={{-10,50},{-10,60},{0,68},{40,68}},    color={0,0,127}),
        Polygon(
          points={{58,50},{50,54},{50,46},{58,50}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-6,30},{46,8}},
          lineColor={28,108,200},
          textString="T")}),   Diagram(coordinateSystem(preserveAspectRatio=false)));
end SysCtrl_T;
