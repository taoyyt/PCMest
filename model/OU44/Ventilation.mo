within OU44;
model Ventilation

  Buildings.Fluid.HeatExchangers.ConstantEffectiveness hex(
    dp1_nominal=50,
    dp2_nominal=50,
    show_T=false,
    eps=hxEta,
    m1_flow_nominal=0.1,
    m2_flow_nominal=0.1,
    allowFlowReversal1=false,
    allowFlowReversal2=false,
    redeclare package Medium1 = Buildings.Media.Specialized.Air.PerfectGas,
    redeclare package Medium2 = Buildings.Media.Specialized.Air.PerfectGas)
    annotation (Placement(transformation(extent={{-70,-6},{-50,14}})));

  Buildings.Fluid.HeatExchangers.PrescribedOutlet heat_coil(
    dp_nominal=50,
    QMin_flow=0,
    m_flow_nominal=0.1,
    T_start(displayUnit="degC") = 294.15,
    allowFlowReversal=false,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas,
    show_T=true,
    tau=60,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    QMax_flow=hCoilCapacity,
    use_X_wSet=false)
    annotation (Placement(transformation(extent={{112,0},{132,20}})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary(
    use_T_in=true,
    use_m_flow_in=true,
    nPorts=1,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas)
    annotation (Placement(transformation(extent={{-94,0},{-74,20}})));
  Buildings.Fluid.Sources.FixedBoundary bou(
                                      nPorts=1, redeclare package Medium =
        Buildings.Media.Specialized.Air.PerfectGas)
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
        rotation=-90,
        origin={-70,-34})));
  Modelica.Blocks.Interfaces.RealInput TAmb "[degC]"
    annotation (Placement(transformation(extent={{-228,-34},{-188,6}})));
  Modelica.Blocks.Interfaces.RealInput VFresh "Ventilation airflow rate [m3/h]"
    annotation (Placement(transformation(extent={{-228,10},{-188,50}})));
  Modelica.Blocks.Interfaces.RealOutput QFlow "Heat added to the air [J]"
    annotation (Placement(transformation(extent={{196,122},{216,142}})));
  parameter Modelica.SIunits.Efficiency hxEta=0.7
    "Heat exchanger effectiveness [-]";
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-174,-24},{-154,-4}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
    annotation (Placement(transformation(extent={{-144,132},{-124,152}})));
  Modelica.Blocks.Interfaces.RealOutput qEl
    "Fan electric energy consumption [W]"
    annotation (Placement(transformation(extent={{194,-110},{214,-90}})));
  Modelica.Blocks.Interfaces.RealOutput QEl
    "Total fan electric energy consumption [J]"
    annotation (Placement(transformation(extent={{194,-152},{214,-132}})));
  Modelica.Blocks.Continuous.Integrator integrator(initType=Modelica.Blocks.Types.Init.InitialState,
      y_start=0)
    annotation (Placement(transformation(extent={{140,-152},{160,-132}})));
  Modelica.Blocks.Interfaces.RealOutput qFlow "Heat added to the air [W]"
    annotation (Placement(transformation(extent={{196,166},{216,186}})));
  Buildings.Fluid.Sensors.TemperatureTwoPort T_after_hex(
    m_flow_nominal=1,
    initType=Modelica.Blocks.Types.Init.SteadyState,
    allowFlowReversal=false,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas,
    transferHeat=false,
    tau=15)
    annotation (Placement(transformation(extent={{66,0},{86,20}})));
  Modelica.Blocks.Interfaces.RealInput TIndoor "[degC]"
    annotation (Placement(transformation(extent={{-228,-72},{-188,-32}})));
  Buildings.Fluid.HeatExchangers.PrescribedOutlet IndoorAir(
    m_flow_nominal=1,
    dp_nominal=50,
    allowFlowReversal=false,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas,
    tau=15,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    use_X_wSet=false)
    annotation (Placement(transformation(extent={{90,-36},{70,-16}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
    annotation (Placement(transformation(extent={{-174,-62},{-154,-42}})));
  Modelica.Blocks.Interfaces.RealOutput TAfterHX
    "Air temperature after leaving HX [degC]" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={76,206})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={76,110})));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold=0)
    annotation (Placement(transformation(extent={{-136,-152},{-116,-132}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal(            realFalse=0, realTrue=
       qFan)
    annotation (Placement(transformation(extent={{78,-152},{98,-132}})));
  parameter Real qFan=1000 "Constant fan electricity consumption rate [W]";
  Modelica.Blocks.Math.Gain h_coil_efficiency(k=1/coilEta)
                                                         annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={144,52})));
  parameter Real coilEta=0.6 "Heating coil efficiency";
  Modelica.Blocks.Math.Gain m3h_to_kgs(k=1.2/3600)
    annotation (Placement(transformation(extent={{-180,20},{-160,40}})));
  Modelica.Blocks.Continuous.Integrator integrator1(
                                                   initType=Modelica.Blocks.Types.Init.InitialState,
      y_start=0)
    annotation (Placement(transformation(extent={{160,122},{180,142}})));
  parameter Modelica.SIunits.HeatFlowRate hCoilCapacity=15000
    "Heating capacity of the heating coil";
  Modelica.Blocks.Math.BooleanToReal booleanToReal1(           realFalse=0, realTrue=
       1)
    annotation (Placement(transformation(extent={{20,-130},{40,-110}})));
  Modelica.Blocks.Math.Product product annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={138,98})));
  Buildings.Fluid.HeatExchangers.PrescribedOutlet correction(
    m_flow_nominal=1,
    dp_nominal=50,
    allowFlowReversal=false,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas,
    tau=15,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    use_X_wSet=false)
    annotation (Placement(transformation(extent={{26,0},{46,20}})));
  Modelica.Blocks.Logical.Switch switch1
    annotation (Placement(transformation(extent={{-12,40},{8,60}})));
  Buildings.Fluid.Sensors.TemperatureTwoPort T_after_hex1(
    m_flow_nominal=1,
    initType=Modelica.Blocks.Types.Init.SteadyState,
    allowFlowReversal=false,
    redeclare package Medium = Buildings.Media.Specialized.Air.PerfectGas,
    transferHeat=false,
    tau=15)
    annotation (Placement(transformation(extent={{-34,0},{-14,20}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=basicFlow)
    annotation (Placement(transformation(extent={{-172,52},{-152,72}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-138,26},{-118,46}})));
  parameter Modelica.Blocks.Interfaces.RealOutput basicFlow=0.1
    "Additional flow in the main ventilation loop";

  Modelica.Blocks.Interfaces.RealInput Tvestp "Ventilation airflow rate [m3/h]"
    annotation (Placement(transformation(extent={{-226,122},{-186,162}})));
equation
  connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-153,
          -14},{-110,-14},{-110,14},{-96,14}},
                                   color={0,0,127}));
  connect(QEl, integrator.y)
    annotation (Line(points={{204,-142},{204,-142},{161,-142}},
                                                           color={0,0,127}));
  connect(hex.port_b2, bou.ports[1])
    annotation (Line(points={{-70,-2},{-70,-2},{-70,-24}}, color={0,127,255}));
  connect(TIndoor, toKelvin1.Celsius) annotation (Line(points={{-208,-52},{-208,
          -52},{-176,-52}},    color={0,0,127}));
  connect(T_after_hex.T, fromKelvin.Kelvin)
    annotation (Line(points={{76,21},{76,98}},         color={0,0,127}));
  connect(greaterThreshold.y, booleanToReal.u) annotation (Line(points={{-115,
          -142},{-115,-142},{76,-142}},
                                color={255,0,255}));
  connect(booleanToReal.y, qEl) annotation (Line(points={{99,-142},{128,
          -142},{128,-100},{204,-100}},
                           color={0,0,127}));
  connect(booleanToReal.y, integrator.u)
    annotation (Line(points={{99,-142},{99,-142},{138,-142}},
                                                           color={0,0,127}));
  connect(qEl, qEl)
    annotation (Line(points={{204,-100},{204,-100}}, color={0,0,127}));
  connect(VFresh,m3h_to_kgs. u)
    annotation (Line(points={{-208,30},{-182,30}},           color={0,0,127}));
  connect(TAmb, toKelvin.Celsius) annotation (Line(points={{-208,-14},{-192,-14},
          {-176,-14}}, color={0,0,127}));
  connect(boundary.ports[1], hex.port_a1)
    annotation (Line(points={{-74,10},{-70,10}}, color={0,127,255}));
  connect(m3h_to_kgs.y, greaterThreshold.u) annotation (Line(points={{-159,30},
          {-148,30},{-148,-142},{-138,-142}},
                                            color={0,0,127}));
  connect(TAfterHX, TAfterHX)
    annotation (Line(points={{76,206},{76,201},{76,206}}, color={0,0,127}));
  connect(heat_coil.Q_flow, h_coil_efficiency.u) annotation (Line(points={{133,18},
          {144,18},{144,40}},         color={0,0,127}));
  connect(fromKelvin.Celsius, TAfterHX) annotation (Line(points={{76,121},{
          76,121},{76,206}}, color={0,0,127}));
  connect(QFlow, integrator1.y) annotation (Line(points={{206,132},{194,132},
          {181,132}}, color={0,0,127}));
  connect(heat_coil.port_b, IndoorAir.port_a) annotation (Line(points={{132,
          10},{140,10},{140,-26},{90,-26}}, color={0,127,255}));
  connect(IndoorAir.port_b, hex.port_a2) annotation (Line(points={{70,-26},
          {-50,-26},{-50,-2}}, color={0,127,255}));
  connect(toKelvin2.Kelvin, heat_coil.TSet) annotation (Line(points={{-123,
          142},{102,142},{102,18},{110,18}}, color={0,0,127}));
  connect(greaterThreshold.y, booleanToReal1.u) annotation (Line(points={{
          -115,-142},{6,-142},{6,-120},{18,-120}}, color={255,0,255}));
  connect(h_coil_efficiency.y, product.u2) annotation (Line(points={{144,63},
          {144,63},{144,86}}, color={0,0,127}));
  connect(product.y, integrator1.u) annotation (Line(points={{138,109},{138,
          109},{138,132},{158,132}}, color={0,0,127}));
  connect(product.y, qFlow) annotation (Line(points={{138,109},{138,109},{
          138,176},{206,176}}, color={0,0,127}));
  connect(booleanToReal1.y, product.u1) annotation (Line(points={{41,-120},
          {52,-120},{52,78},{132,78},{132,86}}, color={0,0,127}));
  connect(toKelvin.Kelvin, switch1.u3) annotation (Line(points={{-153,-14},
          {-110,-14},{-110,42},{-14,42}}, color={0,0,127}));
  connect(greaterThreshold.y, switch1.u2) annotation (Line(points={{-115,
          -142},{-40,-142},{-40,50},{-14,50}}, color={255,0,255}));
  connect(heat_coil.port_a, T_after_hex.port_b)
    annotation (Line(points={{112,10},{86,10}}, color={0,127,255}));
  connect(T_after_hex.port_a, correction.port_b)
    annotation (Line(points={{66,10},{46,10}}, color={0,127,255}));
  connect(hex.port_b1, T_after_hex1.port_a)
    annotation (Line(points={{-50,10},{-34,10}}, color={0,127,255}));
  connect(correction.port_a, T_after_hex1.port_b)
    annotation (Line(points={{26,10},{-14,10}}, color={0,127,255}));
  connect(T_after_hex1.T, switch1.u1) annotation (Line(points={{-24,21},{
          -24,21},{-24,58},{-14,58}}, color={0,0,127}));
  connect(toKelvin1.Kelvin, IndoorAir.TSet) annotation (Line(points={{-153,
          -52},{104,-52},{104,-18},{92,-18}}, color={0,0,127}));
  connect(m3h_to_kgs.y, add.u2)
    annotation (Line(points={{-159,30},{-140,30}}, color={0,0,127}));
  connect(add.u1, realExpression.y) annotation (Line(points={{-140,42},{
          -146,42},{-146,62},{-151,62}}, color={0,0,127}));
  connect(add.y, boundary.m_flow_in) annotation (Line(points={{-117,36},{-104,
          36},{-104,18},{-96,18}},      color={0,0,127}));
  connect(switch1.y, correction.TSet) annotation (Line(points={{9,50},{18,50},
          {18,18},{24,18}},     color={0,0,127}));
  connect(toKelvin2.Celsius, Tvestp) annotation (Line(points={{-146,142},{-156,142},
          {-206,142}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,
        extent={{-200,-200},{200,200}},
        initialScale=0.1),                                      graphics={
        Rectangle(extent={{-200,200},{200,-200}}, lineColor={28,108,200}),
        Ellipse(extent={{0,40},{80,-40}}, lineColor={28,108,200}),
        Polygon(
          points={{20,34},{20,-34},{80,0},{20,34}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(
          points={{0,40},{-20,40},{-60,-40},{-80,-40}},
          color={255,0,0},
          smooth=Smooth.Bezier,
          thickness=1),
        Line(
          points={{-80,40},{-60,40},{-20,-40},{0,-40}},
          color={28,108,200},
          smooth=Smooth.Bezier,
          thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false,
        extent={{-200,-200},{200,200}},
        initialScale=0.1)));
end Ventilation;
