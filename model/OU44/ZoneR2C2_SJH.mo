within OU44;
model ZoneR2C2_SJH "Simple zone"

  package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
    "Medium model";

  Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
    annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
        iconTransformation(extent={{-234,166},{-206,194}})));
  Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
    annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
        iconTransformation(extent={{-234,84},{-206,112}})));
  Modelica.Blocks.Interfaces.RealInput verate
    "Ventilation airflow rate fraction [0-1]"
                                      annotation (Placement(transformation(
          extent={{-240,-92},{-212,-64}}), iconTransformation(extent={{-240,-92},
            {-212,-64}})));
  Modelica.Blocks.Interfaces.RealOutput T "[K]"
    annotation (Placement(transformation(extent={{214,110},{234,130}}),
        iconTransformation(extent={{214,110},{234,130}})));

  parameter Real TairInit=20 "Initial temperature of indoor air [degC]";
  parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
  parameter Real RExt=3.0 "External wall thermal resistance";
  parameter Real RInt=0.05 "Internal wall thermal resistance";
  parameter Real tmass=25 "Zone thermal mass factor [-]";
  parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
  parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
  parameter Real CO2n=435.45 "CO2 Neutral Level";
  parameter Real CO2pp(unit="m3/h")=0.0398 "CO2 generation per person [m3/h]";
  parameter Real maxVent=4800 "Maximum ventilation flowrate [m3/h]";
  parameter Real maxHeat=2689 "Heating power of radiators [W]";

  Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
    annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{62,110},{82,130}})));
  Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
    annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
  Modelica.Blocks.Math.Product hmltp
    annotation (Placement(transformation(extent={{-94,18},{-74,38}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
    annotation (Placement(transformation(extent={{-14,18},{6,38}})));
  Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
    annotation (Placement(transformation(extent={{-160,-6},{-136,18}})));
  Modelica.Blocks.Interfaces.RealInput qrad "Radiator power fraction [0-1]"
                                                                 annotation (
      Placement(transformation(extent={{-240,-150},{-212,-122}}),
        iconTransformation(extent={{-234,-186},{-206,-158}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
  Modelica.Blocks.Interfaces.RealOutput Qrad "Heat supplied by radiators [J]"
                                     annotation (Placement(transformation(
          extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
            {234,-116}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(extent={{-102,-146},{-82,-126}})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
    annotation (Placement(transformation(extent={{182,110},{202,130}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-200,-56},{-180,-36}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
    annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
  parameter Real Vi=486.5 "Air volume [m3]";
  parameter Real occheff=1.2 "Occupant heat generation effectiveness";
  Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
    annotation (Placement(transformation(extent={{-126,-6},{-106,14}})));
  parameter Real Vinf=2292 "Infiltration rate";
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                  C=tmass*1.2*1005*Vi, T(start=TairInit + 273.15))
    annotation (Placement(transformation(extent={{36,116},{56,136}})));
  Components.AirMix airMix
    annotation (Placement(transformation(extent={{-10,-52},{10,-32}})));
  Modelica.Blocks.Interfaces.RealInput Tvestp
    "Ventilation air setpoint temperature [degC]" annotation (Placement(
        transformation(extent={{-240,-60},{-212,-32}}), iconTransformation(
          extent={{-240,-60},{-212,-32}})));
  Components.CO2 cO2_1(
    Vi=Vi,
    CO2PpmInitial=CO2Init,
    CO2PerPerson=CO2pp,
    CO2Neutral=CO2n)
    annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
  Modelica.Blocks.Interfaces.RealOutput CO2 "[ppm]" annotation (Placement(
        transformation(extent={{214,-8},{234,12}}), iconTransformation(
          extent={{214,-8},{234,12}})));
  Modelica.Blocks.Continuous.Integrator integrator
    annotation (Placement(transformation(extent={{104,-144},{124,-124}})));
  Modelica.Blocks.Sources.Constant const1(k=Vinf)
    annotation (Placement(transformation(extent={{-20,74},{-40,94}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));
  Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
    annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
        iconTransformation(extent={{-234,2},{-206,30}})));
  Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
    annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
        iconTransformation(extent={{214,-76},{234,-56}})));
  Modelica.Blocks.Continuous.Integrator integrator1(k=1/3600)
    annotation (Placement(transformation(extent={{100,-88},{120,-68}})));
  Modelica.Blocks.Math.Max max1
    annotation (Placement(transformation(extent={{-134,-56},{-114,-36}})));
  Modelica.Blocks.Math.Gain ventilation(k=maxVent)
    annotation (Placement(transformation(extent={{-164,-88},{-144,-68}})));
  Modelica.Blocks.Math.Gain heating(k=maxHeat)
    annotation (Placement(transformation(extent={{-160,-146},{-140,-126}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
        *Vi, T(fixed=true, start=293.15))
    annotation (Placement(transformation(extent={{108,158},{128,178}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
    annotation (Placement(transformation(extent={{76,148},{96,168}})));
equation
  connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                      color={0,0,127}));
  connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
          166},{-158,166}},            color={0,0,127}));
  connect(hmltp.y, occHeatGain.Q_flow)
    annotation (Line(points={{-73,28},{-73,28},{-14,28}}, color={0,0,127}));
  connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
        points={{-82,-136},{-68,-136},{-68,-120},{-24,-120}},  color={191,0,0}));
  connect(prescribedTemperature.port, re.port_a)
    annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
  connect(temperatureSensor.T, fromKelvin.Kelvin)
    annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
  connect(T, fromKelvin.Celsius)
    annotation (Line(points={{224,120},{210,120},{203,120}},
                                                   color={0,0,127}));
  connect(toKelvin1.Kelvin, prescribedTemperature.T)
    annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
  connect(toKelvin1.Celsius,Tout)
    annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
  connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-96,22},{-102,22},
          {-102,4},{-105,4}},     color={0,0,127}));
  connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-134.8,6},
          {-132,6},{-132,4},{-128,4}},  color={0,0,127}));
  connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
          {90,120},{90,-14},{-168,-14},{-168,6},{-162.4,6}},
        color={0,0,127}));
  connect(re.port_b, air.port)
    annotation (Line(points={{-108,90},{-80,90},{-80,116},{46,116}},
                                                 color={191,0,0}));
  connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,166},
          {-80,166},{-80,116},{46,116}}, color={191,0,0}));
  connect(temperatureSensor.port, air.port)
    annotation (Line(points={{62,120},{62,116},{46,116}},
                                                  color={191,0,0}));
  connect(occHeatGain.port, air.port)
    annotation (Line(points={{6,28},{46,28},{46,116}},   color={191,0,0}));
  connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},{
          46,-120},{46,116}},                color={191,0,0}));
  connect(Tvestp, toKelvin.Celsius) annotation (Line(points={{-226,-46},{-226,
          -46},{-202,-46}}, color={0,0,127}));
  connect(cO2_1.CO2, CO2)
    annotation (Line(points={{10.6,2},{118,2},{224,2}},
                                                color={0,0,127}));
  connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
          -130},{-14,-134},{102,-134}},                color={0,0,127}));
  connect(Qrad, integrator.y)
    annotation (Line(points={{226,-134},{125,-134}}, color={0,0,127}));
  connect(cO2_1.Vve, add.y)
    annotation (Line(points={{-10.4,8},{-23,8}}, color={0,0,127}));
  connect(const1.y, add.u1) annotation (Line(points={{-41,84},{-60,84},{-60,14},
          {-46,14}},     color={0,0,127}));
  connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-96,34}},
                color={0,0,127}));
  connect(occ, cO2_1.persons) annotation (Line(points={{-226,34},{-116,34},
          {-116,46},{-70,46},{-70,-4},{-10.4,-4}}, color={0,0,127}));
  connect(integrator1.y, vetot)
    annotation (Line(points={{121,-78},{226,-78}}, color={0,0,127}));
  connect(toKelvin.Kelvin, max1.u2) annotation (Line(points={{-179,-46},{-179,
          -52},{-136,-52}}, color={0,0,127}));
  connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-179,90},{-174,
          90},{-174,-40},{-136,-40}}, color={0,0,127}));
  connect(airMix.port_b, air.port)
    annotation (Line(points={{10,-42},{46,-42},{46,116}}, color={191,0,0}));
  connect(verate, ventilation.u) annotation (Line(points={{-226,-78},{-196,
          -78},{-166,-78}}, color={0,0,127}));
  connect(ventilation.y, add.u2) annotation (Line(points={{-143,-78},{-96,-78},
          {-96,2},{-46,2}}, color={0,0,127}));
  connect(ventilation.y, integrator1.u)
    annotation (Line(points={{-143,-78},{98,-78}}, color={0,0,127}));
  connect(qrad, heating.u) annotation (Line(points={{-226,-136},{-194,-136},{
          -162,-136}}, color={0,0,127}));
  connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{
          -102,-136},{-120,-136},{-139,-136}}, color={0,0,127}));
  connect(re1.port_b, IntMass.port)
    annotation (Line(points={{96,158},{110,158},{118,158}}, color={191,0,0}));
  connect(re1.port_a, air.port) annotation (Line(points={{76,158},{62,158},{62,
          116},{46,116}},
                     color={191,0,0}));
  connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{90,
          120},{90,-60},{0,-60},{0,-52.8},{0,-52.8}}, color={0,0,127}));
  connect(ventilation.y, airMix.Vve) annotation (Line(points={{-143,-78},{-52,
          -78},{-52,-47},{-11,-47}}, color={0,0,127}));
  connect(max1.y, airMix.Tve) annotation (Line(points={{-113,-46},{-66,-46},{
          -66,-38.6},{-10.8,-38.6}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
            -220},{220,220}},
        initialScale=0.1), graphics={
        Rectangle(
          extent={{-208,138},{-100,66}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-208,60},{16,-20}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-202,136},{-104,122}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          textString="External walls (cond. and inf.)"),
        Text(
          extent={{-204,60},{-118,46}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          horizontalAlignment=TextAlignment.Left,
          textString="Occ. heat gains and CO2"),
        Rectangle(
          extent={{-208,212},{-124,144}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-204,212},{-150,202}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          horizontalAlignment=TextAlignment.Left,
          textString="Solar heat gains"),
        Rectangle(
          extent={{-208,-26},{16,-96}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-30,-82},{10,-92}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          textString="Ventilation"),
        Rectangle(
          extent={{-208,-104},{16,-160}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-208,-106},{-170,-116}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          textString="Heating"),
        Rectangle(
          extent={{-74,116},{16,66}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-70,114},{0,104}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          horizontalAlignment=TextAlignment.Left,
          textString="Infiltration (only CO2)")}),
    experiment(Tolerance=1e-006),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
            220}},
        initialScale=0.1), graphics={
        Rectangle(
          extent={{-220,-220},{220,220}},
          lineColor={95,95,95},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-202,202},{200,-202}},
          pattern=LinePattern.None,
          lineColor={117,148,176},
          fillColor={170,213,255},
          fillPattern=FillPattern.Sphere),
        Rectangle(
          extent={{-96,102},{96,-100}},
          lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid)}));
end ZoneR2C2_SJH;
