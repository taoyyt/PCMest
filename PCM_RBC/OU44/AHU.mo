within OU44;
model AHU

  replaceable package Air = Buildings.Media.Air;

  Modelica.Blocks.Interfaces.RealInput Vin
    annotation (Placement(transformation(extent={{-128,-60},{-88,-20}})));
  Modelica.Blocks.Interfaces.RealInput Vex
    annotation (Placement(transformation(extent={{128,20},{88,60}})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary(use_m_flow_in=true,
    redeclare package Medium = Air,
    use_T_in=true,
    nPorts=1)
    annotation (Placement(transformation(extent={{46,22},{26,42}})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary1(          use_m_flow_in=true,
    redeclare package Medium = Air,
    use_T_in=true,
    nPorts=1)
    annotation (Placement(transformation(extent={{-40,-36},{-20,-16}})));
  Buildings.Fluid.HeatExchangers.ConstantEffectiveness hex(
    redeclare package Medium1 = Air,
    redeclare package Medium2 = Air,
    m1_flow_nominal=0.1,
    m2_flow_nominal=0.1,
    dp1_nominal=50,
    dp2_nominal=50,
    eps=eps,
    allowFlowReversal1=false,
    allowFlowReversal2=false)
    annotation (Placement(transformation(extent={{10,-10},{-10,10}})));
  Buildings.Fluid.HeatExchangers.Heater_T hea(
    redeclare package Medium = Air,
    m_flow_nominal=0.1,
    dp_nominal=50,
    allowFlowReversal=false,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState)
    annotation (Placement(transformation(extent={{50,-36},{70,-16}})));
  Components.FanFit fanFit(
    a=a,
    b=b,
    c=c)
    annotation (Placement(transformation(extent={{-42,-64},{-22,-44}})));
  Modelica.Blocks.Interfaces.RealInput Tsp
    annotation (Placement(transformation(extent={{-128,-102},{-88,-62}})));
  Modelica.Blocks.Math.Gain gain(k=1.2)
    annotation (Placement(transformation(extent={{-76,-28},{-56,-8}})));
  Buildings.Fluid.Sources.FixedBoundary bou(nPorts=1, redeclare package Medium =
        Air)
    annotation (Placement(transformation(extent={{98,-36},{78,-16}})));
  Buildings.Fluid.Sources.FixedBoundary bou1(nPorts=1, redeclare package Medium =
        Air)
    annotation (Placement(transformation(extent={{-88,16},{-68,36}})));
  Modelica.Blocks.Math.Gain gain1(k=1.2)
    annotation (Placement(transformation(extent={{76,30},{56,50}})));
  Modelica.Blocks.Interfaces.RealInput Tin
    annotation (Placement(transformation(extent={{-128,-20},{-88,20}})));
  Modelica.Blocks.Interfaces.RealInput Tex
    annotation (Placement(transformation(extent={{128,-16},{88,24}})));
  Modelica.Blocks.Interfaces.RealOutput qh
    annotation (Placement(transformation(extent={{96,-66},{116,-46}})));
  Modelica.Blocks.Interfaces.RealOutput qel
    annotation (Placement(transformation(extent={{96,-90},{116,-70}})));
  parameter Real a=1 "Parameter in y=ax^2+bx+c";
  parameter Real b=1 "Parameter in y=ax^2+bx+c";
  parameter Real c=1 "Parameter in y=ax^2+bx+c";
  parameter Modelica.SIunits.Efficiency eps=0.8 "Heat exchanger effectiveness";
  Buildings.Fluid.Sensors.TemperatureTwoPort Thxex(redeclare package Medium =
        Air, m_flow_nominal=0.1)
    annotation (Placement(transformation(extent={{-28,16},{-48,36}})));
  Buildings.Fluid.Sensors.TemperatureTwoPort Thxin(redeclare package Medium =
        Air, m_flow_nominal=0.1)
    annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
equation
  connect(Tsp, hea.TSet) annotation (Line(points={{-108,-82},{44,-82},{44,-18},
          {48,-18}},color={0,0,127}));
  connect(Vin, fanFit.VFR) annotation (Line(points={{-108,-40},{-82,-40},{-82,-54},
          {-42.8,-54}}, color={0,0,127}));
  connect(Vin, gain.u) annotation (Line(points={{-108,-40},{-82,-40},{-82,-18},
          {-78,-18}},color={0,0,127}));
  connect(boundary1.m_flow_in, gain.y)
    annotation (Line(points={{-42,-18},{-55,-18}}, color={0,0,127}));
  connect(hea.port_b, bou.ports[1])
    annotation (Line(points={{70,-26},{78,-26}}, color={0,127,255}));
  connect(Vex, gain1.u)
    annotation (Line(points={{108,40},{78,40}}, color={0,0,127}));
  connect(gain1.y, boundary.m_flow_in)
    annotation (Line(points={{55,40},{48,40}}, color={0,0,127}));
  connect(Tin, boundary1.T_in) annotation (Line(points={{-108,0},{-50,0},{-50,
          -22},{-42,-22}}, color={0,0,127}));
  connect(Tex, boundary.T_in) annotation (Line(points={{108,4},{54,4},{54,36},{
          48,36}}, color={0,0,127}));
  connect(hea.Q_flow, qh) annotation (Line(points={{71,-18},{74,-18},{74,-56},{
          106,-56}}, color={0,0,127}));
  connect(fanFit.qel, qel) annotation (Line(points={{-21.4,-54},{64,-54},{64,
          -80},{106,-80}}, color={0,0,127}));
  connect(Thxex.port_b, bou1.ports[1])
    annotation (Line(points={{-48,26},{-68,26}}, color={0,127,255}));
  connect(Thxin.port_b, hea.port_a)
    annotation (Line(points={{38,-26},{50,-26}}, color={0,127,255}));
  connect(boundary1.ports[1], hex.port_a2) annotation (Line(points={{-20,-26},{
          -16,-26},{-16,-6},{-10,-6}}, color={0,127,255}));
  connect(hex.port_b2, Thxin.port_a) annotation (Line(points={{10,-6},{14,-6},{
          14,-26},{18,-26}}, color={0,127,255}));
  connect(boundary.ports[1], hex.port_a1) annotation (Line(points={{26,32},{20,
          32},{20,6},{10,6}}, color={0,127,255}));
  connect(hex.port_b1, Thxex.port_a) annotation (Line(points={{-10,6},{-18,6},{
          -18,26},{-28,26}}, color={0,127,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end AHU;
