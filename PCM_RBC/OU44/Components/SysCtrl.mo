within OU44.Components;
model SysCtrl

  Modelica.Blocks.Interfaces.RealInput Tstp
    annotation (Placement(transformation(extent={{-128,40},{-88,80}})));
  Modelica.Blocks.Interfaces.RealInput CO2stp
    annotation (Placement(transformation(extent={{-128,-60},{-88,-20}})));
  Modelica.Blocks.Interfaces.RealInput CO2m annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,-108})));
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
    k=k,
    Td=Td,
    initType=Modelica.Blocks.Types.InitPID.NoInit,
    Ti=Ti,
    controllerType=Modelica.Blocks.Types.SimpleController.PI)
            annotation (Placement(transformation(extent={{-30,50},{-10,70}})));
  Modelica.Blocks.Continuous.LimPID PID1(
    yMax=0,
    yMin=-100,
    limitsAtInit=true,
    y_start=0,
    xi_start=0,
    xd_start=0,
    k=kCO2,
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    initType=Modelica.Blocks.Types.InitPID.NoInit)
    annotation (Placement(transformation(extent={{-30,-50},{-10,-30}})));
  Modelica.Blocks.Interfaces.RealOutput hPos
    annotation (Placement(transformation(extent={{96,50},{116,70}})));
  Modelica.Blocks.Interfaces.RealOutput vPos
    annotation (Placement(transformation(extent={{96,-50},{116,-30}})));
  Modelica.Blocks.Math.Abs abs1
    annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
  parameter Real k=1 "Gain of controller";
  parameter Real kCO2=1 "Gain of controller";
  parameter Modelica.SIunits.Time Td=0.1
    "Time constant of Derivative block";
  parameter Modelica.SIunits.Time Ti=0.5
    "Time constant of Integrator block";
equation
  connect(Tstp, PID.u_s)
    annotation (Line(points={{-108,60},{-32,60}}, color={0,0,127}));
  connect(Tm, PID.u_m) annotation (Line(points={{-40,110},{-40,40},{-20,40},{-20,
          48}}, color={0,0,127}));
  connect(PID.y, hPos)
    annotation (Line(points={{-9,60},{106,60}}, color={0,0,127}));
  connect(CO2stp, PID1.u_s)
    annotation (Line(points={{-108,-40},{-32,-40}}, color={0,0,127}));
  connect(CO2m, PID1.u_m) annotation (Line(points={{-40,-108},{-40,-70},{-20,-70},
          {-20,-52}}, color={0,0,127}));
  connect(PID1.y, abs1.u)
    annotation (Line(points={{-9,-40},{18,-40}}, color={0,0,127}));
  connect(abs1.y, vPos)
    annotation (Line(points={{41,-40},{106,-40}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={28,108,200},
          fillColor={236,236,236},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,60},{-40,60}}, color={28,108,200}),
        Line(points={{-40,-98},{-40,-40}}, color={28,108,200}),
        Rectangle(extent={{-20,80},{60,40}}, lineColor={28,108,200}),
        Line(points={{-40,88}}, color={28,108,200}),
        Line(points={{-40,100},{-40,60}}, color={28,108,200}),
        Line(points={{-40,60},{-20,60}}, color={28,108,200}),
        Line(points={{60,60},{96,60}}, color={28,108,200}),
        Line(points={{28,0}}, color={28,108,200}),
        Rectangle(extent={{-20,-20},{60,-60}}, lineColor={28,108,200}),
        Line(points={{-100,-40},{-20,-40}}, color={28,108,200}),
        Line(points={{60,-40},{96,-40}}, color={28,108,200}),
        Line(points={{-10,-24},{-10,-58}},color={0,0,0}),
        Line(points={{-18,-50},{54,-50}}, color={0,0,0}),
        Polygon(
          points={{58,-50},{50,-46},{50,-54},{58,-50}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-10,-50},{-10,-40},{0,-32},{40,-32}},color={0,0,127}),
        Polygon(
          points={{4,0},{-4,4},{-4,-4},{4,0}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          origin={-10,-26},
          rotation=90),
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
          textString="T"),
        Text(
          extent={{-16,-70},{58,-92}},
          lineColor={28,108,200},
          textString="CO2")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
end SysCtrl;
