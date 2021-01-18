within OU44;
model ZoneExp "Simple zone"

  package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
    "Medium model";

  Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
    annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
        iconTransformation(extent={{-234,166},{-206,194}})));
  Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
    annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
        iconTransformation(extent={{-234,84},{-206,112}})));
  Modelica.Blocks.Interfaces.RealOutput T "[K]"
    annotation (Placement(transformation(extent={{214,110},{234,130}}),
        iconTransformation(extent={{214,110},{234,130}})));

  parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
  parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
  parameter Real RExt=1.0 "External wall thermal resistance";
  parameter Real RInt=1.0 "Internal wall thermal resistance";
  parameter Real tmass=5 "Zone thermal mass factor [-]";
  parameter Real imass=10 "Zone internal thermal mass factor [-]";
  parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
  parameter Real CO2n=400 "CO2 Neutral Level";
  parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
  parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
  parameter Real maxHeat=5000 "Heating power of radiators [W]";

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
    annotation (Placement(transformation(extent={{-160,-8},{-136,16}})));
  Modelica.Blocks.Interfaces.RealInput Tstp
    "Indoor temperature setpoint [degC]"                         annotation (
      Placement(transformation(extent={{-240,-150},{-212,-122}}),
        iconTransformation(extent={{-234,-186},{-206,-158}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
  Modelica.Blocks.Interfaces.RealOutput Qrad
    "Heat supplied by radiators [Wh]"
                                     annotation (Placement(transformation(
          extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
            {234,-116}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
    annotation (Placement(transformation(extent={{182,110},{202,130}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-166,-62},{-146,-42}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
    annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
  parameter Real Vi=300 "Air volume [m3]";
  parameter Real occheff=1. "Occupant heat generation effectiveness";
  Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
    annotation (Placement(transformation(extent={{-126,-6},{-106,14}})));
  parameter Real Vinf=50 "Infiltration rate";
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(T(fixed=true,
        start=TAirInit + 273.15), C=tmass*1.2*1005*Vi)
    annotation (Placement(transformation(extent={{36,120},{56,140}})));
  Components.AirMix airMix
    annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
  Components.CO2 cO2_1(
    Vi=Vi,
    CO2PpmInitial=CO2Init,
    CO2PerPerson=CO2pp,
    CO2Neutral=CO2n)
    annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
  Modelica.Blocks.Interfaces.RealOutput CO2 "[ppm]" annotation (Placement(
        transformation(extent={{214,-8},{234,12}}), iconTransformation(
          extent={{214,-8},{234,12}})));
  Modelica.Blocks.Continuous.Integrator integrator(k=1/3600)
    annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
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
    annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
  Modelica.Blocks.Math.Max max1
    annotation (Placement(transformation(extent={{-134,-56},{-114,-36}})));
  Modelica.Blocks.Math.Gain ventilation(k=maxVent)
    annotation (Placement(transformation(extent={{-88,-88},{-68,-68}})));
  Modelica.Blocks.Math.Gain heating(k=maxHeat)
    annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
        *Vi, T(fixed=true, start=293.15))
    annotation (Placement(transformation(extent={{108,156},{128,176}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
    annotation (Placement(transformation(extent={{76,148},{96,168}})));
  Modelica.Blocks.Continuous.LimPID PID_temp(
    yMax=1,
    yMin=0,
    initType=Modelica.Blocks.Types.InitPID.InitialState,
    k=0.5,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    Ti=600)
    annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
    annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
  Modelica.Blocks.Interfaces.RealOutput verate
    "Ventilation airflow rate [m3/h]"
    annotation (Placement(transformation(extent={{216,-58},{236,-38}}),
        iconTransformation(extent={{214,-76},{234,-56}})));
  Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
    annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
        iconTransformation(extent={{214,-136},{234,-116}})));
  Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-124,-78})));
  Modelica.Blocks.Interfaces.RealInput dpos "VAV damper position [%]"
    annotation (Placement(transformation(extent={{-242,-92},{-214,-64}}),
        iconTransformation(extent={{-240,-60},{-212,-32}})));
  Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
    annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
        iconTransformation(extent={{214,-136},{234,-116}})));
  Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={160,-204})));
  Modelica.Blocks.Sources.Constant const2(k=Tve)
    annotation (Placement(transformation(extent={{-202,-62},{-182,-42}})));
  parameter Real Tve=21 "Ventilation air temperature";
equation
  connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                      color={0,0,127}));
  connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
          166},{-158,166}},            color={0,0,127}));
  connect(hmltp.y, occHeatGain.Q_flow)
    annotation (Line(points={{-73,28},{-73,28},{-14,28}}, color={0,0,127}));
  connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
        points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
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
  connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-134.8,
          4},{-128,4}},                 color={0,0,127}));
  connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
          {90,120},{90,-14},{-168,-14},{-168,4},{-162.4,4}},
        color={0,0,127}));
  connect(re.port_b, air.port)
    annotation (Line(points={{-108,90},{-108,90},{-80,90},{-80,120},{46,120}},
                                                 color={191,0,0}));
  connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
          166},{-80,166},{-80,120},{46,120}},
                                         color={191,0,0}));
  connect(temperatureSensor.port, air.port)
    annotation (Line(points={{62,120},{62,120},{46,120}},
                                                  color={191,0,0}));
  connect(occHeatGain.port, air.port)
    annotation (Line(points={{6,28},{46,28},{46,120}},   color={191,0,0}));
  connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
          90,120},{90,-22},{90,-60},{0,-60},{0,-50.8}},          color={0,0,127}));
  connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
          {46,-120},{46,120}},               color={191,0,0}));
  connect(cO2_1.CO2, CO2)
    annotation (Line(points={{10.6,2},{118,2},{224,2}},
                                                color={0,0,127}));
  connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
          -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
  connect(Qrad, integrator.y)
    annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                     color={0,0,127}));
  connect(cO2_1.Vve, add.y)
    annotation (Line(points={{-10.4,8},{-23,8}}, color={0,0,127}));
  connect(const1.y, add.u1) annotation (Line(points={{-41,84},{-60,84},{-60,14},
          {-46,14}},     color={0,0,127}));
  connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-96,34}},
                color={0,0,127}));
  connect(occ, cO2_1.persons) annotation (Line(points={{-226,34},{-116,34},
          {-116,46},{-70,46},{-70,-4},{-10.4,-4}}, color={0,0,127}));
  connect(integrator1.y, vetot)
    annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                   color={0,0,127}));
  connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-179,90},{-174,
          90},{-174,-40},{-136,-40}}, color={0,0,127}));
  connect(max1.y, airMix.Tve) annotation (Line(points={{-113,-46},{-36,-46},{
          -36,-36.6},{-10.8,-36.6}}, color={0,0,127}));
  connect(airMix.port_b, air.port)
    annotation (Line(points={{10,-40},{46,-40},{46,120}}, color={191,0,0}));
  connect(ventilation.y, add.u2) annotation (Line(points={{-67,-78},{-60,-78},
          {-60,2},{-46,2}}, color={0,0,127}));
  connect(ventilation.y, airMix.Vve) annotation (Line(points={{-67,-78},{-28,
          -78},{-28,-45},{-11,-45}}, color={0,0,127}));
  connect(ventilation.y, integrator1.u)
    annotation (Line(points={{-67,-78},{-67,-78},{148,-78}},
                                                   color={0,0,127}));
  connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
          -136},{-54,-136},{-61,-136}},        color={0,0,127}));
  connect(re1.port_b, IntMass.port)
    annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
  connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
          120},{46,120}},
                     color={191,0,0}));
  connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
          -136},{-84,-136}}, color={0,0,127}));
  connect(Tstp, PID_temp.u_s)
    annotation (Line(points={{-226,-136},{-190,-136},{-152,-136}},
                                                       color={0,0,127}));
  connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
          -180},{-140,-180},{-140,-148}}, color={0,0,127}));
  connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{
          -48,-180},{90,-180},{90,120},{82,120}}, color={0,0,127}));
  connect(ventilation.y, verate) annotation (Line(points={{-67,-78},{126.5,
          -78},{126.5,-48},{226,-48}}, color={0,0,127}));
  connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
          -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
  connect(ventilation.u, scale.y)
    annotation (Line(points={{-90,-78},{-113,-78}}, color={0,0,127}));
  connect(dpos, scale.u)
    annotation (Line(points={{-228,-78},{-136,-78}}, color={0,0,127}));
  connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
          {-106,-204},{148,-204}}, color={0,0,127}));
  connect(scale1.y, vpos)
    annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
  connect(max1.u2, toKelvin.Kelvin)
    annotation (Line(points={{-136,-52},{-145,-52}}, color={0,0,127}));
  connect(toKelvin.Celsius, const2.y)
    annotation (Line(points={{-168,-52},{-181,-52}}, color={0,0,127}));
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
          extent={{-208,214},{-124,146}},
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
          extent={{-208,-28},{16,-98}},
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
end ZoneExp;
