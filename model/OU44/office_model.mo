within OU44;
package office_model
  model office_zone "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
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
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";


    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
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
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-86,-78},{-66,-58}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
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
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-238,-80},{-210,-52}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,120},{
            90,120},{90,-60},{0,-60},{0,-48.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix.port_b, air.port)
      annotation (Line(points={{10,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-68},{-28,
            -68},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-190,-136},{-152,-136}},
                                                         color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-186,34},{-186,32},{-146,32}},
                                                     color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-224,-66},{-160,-66},{-160,
            -68},{-88,-68}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
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
            textString="Heating")}),
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
  end office_zone;

  model office_zone_nowindow "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
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
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
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
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,18},{-212,46}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=0)
      annotation (Placement(transformation(extent={{-86,-78},{-66,-58}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=100)
      annotation (Placement(transformation(extent={{-174,-146},{-154,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-236,-82},{-208,-54}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
    Modelica.Blocks.Math.Gain scale(k=0.01)
      annotation (Placement(transformation(extent={{-126,-146},{-106,-126}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,120},{
            90,120},{90,-60},{0,-60},{0,-48.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix.port_b, air.port)
      annotation (Line(points={{10,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-68},{-28,
            -68},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-176,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-164,-180},{-164,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,32},{-146,32}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-222,-68},{-88,-68}},
                             color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    connect(PID_temp.y, scale.u)
      annotation (Line(points={{-153,-136},{-128,-136}}, color={0,0,127}));
    connect(scale.y, heating.u)
      annotation (Line(points={{-105,-136},{-84,-136}}, color={0,0,127}));
    connect(vpos, vpos)
      annotation (Line(points={{226,-200},{226,-200}}, color={0,0,127}));
    connect(PID_temp.y, vpos) annotation (Line(points={{-153,-136},{-146,-136},
            {-146,-200},{226,-200}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
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
            textString="Heating")}),
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
  end office_zone_nowindow;

  model office_zone_MPC "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
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
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{68,110},{88,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,
              -110},{234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{68,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix.Vve) annotation (Line(points={{-234,-70},{-80,-70},{
            -80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-22,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,
            -36},{-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix.Tve) annotation (Line(points={{-107.3,-35},
            {-84.05,-35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-12,
            -48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad) annotation (Line(points={{-12,-48},{
            -12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{88,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{88,120},{
            98,120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
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
  end office_zone_MPC;

  model office_zone_MPC_est_para "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=2.86 "External wall thermal resistance";
    parameter Real RInt=0.102 "Internal wall thermal resistance";
    parameter Real tmass=40 "Zone thermal mass factor [-]";
    parameter Real imass=199 "Zone internal thermal mass factor [-]";
    parameter Real shgc=2.82 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2
          /3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{68,110},{88,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=171.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            293.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(
          2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{68,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix.Vve) annotation (Line(points={{-234,-70},{-80,-70},{-80,-43},
            {-61,-43}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-22,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,-36},
            {-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix.Tve) annotation (Line(points={{-107.3,-35},{-84.05,
            -35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{88,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{88,120},{98,
            120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
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
  end office_zone_MPC_est_para;

  model office_zone_MPC_est_para_KA "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]" annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}), iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,146},{96,166}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-16,-48},{4,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,166},
            {-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}}, color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,156},{118,156}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,156},{60,156},{60,120},
            {46,120}}, color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix_new.Vve) annotation (Line(points={{-234,-70},{-80,-70},
            {-80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix_new.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-16,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{4,-38},{46,-38},{46,120}},  color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,-36},
            {-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix_new.Tve) annotation (Line(points={{-107.3,
            -35},{-84.05,-35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-6,-48},{-6,-58},{142,-58}},   color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-6,-48},{-6,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{92,
            120},{98,120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA;

  model office_zone_nowindow_update "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt= 5 "External wall thermal resistance";
    parameter Real RInt=13 "Internal wall thermal resistance";
    parameter Real Cair=5e7
                           "Zone thermal mass factor [-]";
    //parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real Cin=4e6 "Zone internal thermal mass factor [-]";
    //parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt)
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,108},{84,128}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
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
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                                         T(fixed=false, start=
            TAirInit + 273.15), C=Cair)
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,18},{-212,46}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=0)
      annotation (Placement(transformation(extent={{-86,-78},{-66,-58}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(
               T(fixed=false, start=TAirInit + 273.15), C=Cin)
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt)
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=100)
      annotation (Placement(transformation(extent={{-174,-146},{-154,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-236,-82},{-208,-54}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
    Modelica.Blocks.Math.Gain scale(k=0.01)
      annotation (Placement(transformation(extent={{-126,-146},{-106,-126}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,118},{134,118},{134,120},{180,120}},
                                                               color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,118},{56,118},{56,120},{46,120}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,118},{90,
            118},{90,-60},{0,-60},{0,-48.8}},                      color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix.port_b, air.port)
      annotation (Line(points={{10,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-68},{-28,
            -68},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-176,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-164,-180},{-164,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,118},{84,118}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,32},{-146,32}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-222,-68},{-88,-68}},
                             color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    connect(PID_temp.y, scale.u)
      annotation (Line(points={{-153,-136},{-128,-136}}, color={0,0,127}));
    connect(scale.y, heating.u)
      annotation (Line(points={{-105,-136},{-84,-136}}, color={0,0,127}));
    connect(vpos, vpos)
      annotation (Line(points={{226,-200},{226,-200}}, color={0,0,127}));
    connect(PID_temp.y, vpos) annotation (Line(points={{-153,-136},{-146,-136},
            {-146,-200},{226,-200}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
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
            textString="Heating")}),
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
  end office_zone_nowindow_update;

  model office_zone_nowindow_update1 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
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
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
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
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{-80,90},{-60,110}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,18},{-212,46}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=0)
      annotation (Placement(transformation(extent={{-86,-78},{-66,-58}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{26,90},{46,70}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-34,80},{-14,100}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=100)
      annotation (Placement(transformation(extent={{-174,-146},{-154,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-236,-82},{-208,-54}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
    Modelica.Blocks.Math.Gain scale(k=0.01)
      annotation (Placement(transformation(extent={{-126,-146},{-106,-126}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-70,90}},color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,120},{
            90,120},{90,-60},{0,-60},{0,-48.8}},                   color={0,0,127}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-68},{-28,
            -68},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-176,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-164,-180},{-164,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,32},{-146,32}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-222,-68},{-88,-68}},
                             color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    connect(PID_temp.y, scale.u)
      annotation (Line(points={{-153,-136},{-128,-136}}, color={0,0,127}));
    connect(scale.y, heating.u)
      annotation (Line(points={{-105,-136},{-84,-136}}, color={0,0,127}));
    connect(vpos, vpos)
      annotation (Line(points={{226,-200},{226,-200}}, color={0,0,127}));
    connect(PID_temp.y, vpos) annotation (Line(points={{-153,-136},{-146,-136},
            {-146,-200},{226,-200}}, color={0,0,127}));
    connect(air.port, re1.port_a)
      annotation (Line(points={{-70,90},{-34,90}}, color={191,0,0}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{-14,90},{36,90}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, IntMass.port)
      annotation (Line(points={{-138,166},{36,166},{36,90}}, color={191,0,0}));
    connect(IntMass.port, temperatureSensor.port) annotation (Line(points={{36,
            90},{50,90},{50,120},{64,120}}, color={191,0,0}));
    connect(occHeatGain.port, IntMass.port) annotation (Line(points={{8,34},{22,
            34},{22,90},{36,90}}, color={191,0,0}));
    connect(airMix.port_b, IntMass.port) annotation (Line(points={{10,-38},{22,
            -38},{22,90},{36,90}}, color={191,0,0}));
    connect(heatFlowSensor.port_b, IntMass.port) annotation (Line(points={{-4,
            -120},{22,-120},{22,90},{36,90}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
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
            textString="Heating")}),
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
  end office_zone_nowindow_update1;

  model office_zone_MPC_est_para_KA_cooling_demand "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120
                                "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,122},{56,142}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,146},{96,166}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-76},{-218,-48}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-84,-72},{-64,-52}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=0,
      yMin=-1)
      annotation (Placement(transformation(extent={{-182,-68},{-170,-56}})));
    Modelica.Blocks.Math.Gain gain3(k=10000)
      annotation (Placement(transformation(extent={{-7,-7},{7,7}},
          rotation=0,
          origin={-119,-61})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,122},{46,122}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,122},{46,122}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,122},{46,122}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,122}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,156},{118,156}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,156},{60,156},{
            60,122},{46,122}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,122}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(fromKelvin.Celsius, PID_temp.u_m) annotation (Line(points={{203,120},
            {204,120},{204,-90},{-176,-90},{-176,-69.2}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-232,-62},{-183.2,-62}}, color={0,0,127}));
    connect(occHeatGain1.port, heatFlowSensor1.port_a) annotation (Line(points=
            {{-64,-62},{-44,-62},{-44,-38},{-22,-38}}, color={191,0,0}));
    connect(PID_temp.y, gain3.u) annotation (Line(points={{-169.4,-62},{-148,
            -62},{-148,-61},{-127.4,-61}}, color={0,0,127}));
    connect(gain3.y, occHeatGain1.Q_flow) annotation (Line(points={{-111.3,-61},
            {-98.65,-61},{-98.65,-62},{-84,-62}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-22},{-130,-38}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_cooling_demand;

  package Cpcm_mean_mpc

  end Cpcm_mean_mpc;

  model office_zone_MPC_est_para_KA_copy "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,84},{-108,104}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi,
      der_T(fixed=false),
      T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{34,122},{54,142}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi,
      T(fixed=false, start=293.15),
      der_T(fixed=false))
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-20,-48},{0,-28}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,94},{-128,94}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,94},{-80,94},{-80,122},{44,122}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,122},{44,122}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,122},{44,122}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{44,34},{44,122}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,122},{44,122}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix_new.Vve) annotation (Line(points={{-234,-70},{-80,-70},
            {-80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix_new.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-20,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{0,-38},{44,-38},{44,122}},  color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-10,-48},{-10,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{92,120},
            {98,120},{98,-70},{-50,-70},{-50,-48.8}},      color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},
            {-172,90},{-172,-34},{-60.8,-34},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-10,-48},
            {-10,-58},{142,-58}},      color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-26},{16,-96}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventilation cooling
")}), experiment(Tolerance=1e-11),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_copy;

  model office_zone_MPC_est_para_KA_copy_q "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,124},{56,144}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput q "heat flow rate" annotation (
        Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-102,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,124},{46,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{46,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{46,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{108,158},{108,156},{118,156}},
                                                              color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,124},{46,124}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-114,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-91,34},{-12,34}},  color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,124}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(q, occHeatGain1.Q_flow)
      annotation (Line(points={{-234,-70},{-88,-70}}, color={0,0,127}));
    connect(occHeatGain1.port, heatFlowSensor1.port_a) annotation (Line(points=
            {{-68,-70},{-46,-70},{-46,-38},{-22,-38}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-11),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_copy_q;

  model office_zone_MPC_est_para_KA_copy1 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            293.15))
      annotation (Placement(transformation(extent={{34,124},{54,144}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
    Modelica.Blocks.Math.Add add(k1=+1, k2=-1)
      annotation (Placement(transformation(extent={{-98,-10},{-78,10}})));
    Modelica.Blocks.Math.Gain gain1(k=1/3600)
      annotation (Placement(transformation(extent={{-56,-6},{-44,6}})));
    Modelica.Blocks.Math.Gain gain2(k=1.2*1005)
      annotation (Placement(transformation(extent={{-34,-6},{-22,6}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow VentHeatGain
      annotation (Placement(transformation(extent={{-66,-50},{-46,-30}})));
    Modelica.Blocks.Math.Product product
      annotation (Placement(transformation(extent={{-98,-48},{-78,-28}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,124},{44,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{44,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{44,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{44,34},{44,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,124},{44,124}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{0,-40},{44,-40},{44,124}},  color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-10,-50},{-10,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-10,-50},
            {-10,-58},{142,-58}},      color={0,0,127}));
    connect(toKelvin1.Kelvin, add.u1) annotation (Line(points={{-179,90},{-170,
            90},{-170,6},{-100,6}}, color={0,0,127}));
    connect(temperatureSensor.T, add.u2) annotation (Line(points={{92,120},{96,
            120},{96,-74},{-112,-74},{-112,-6},{-100,-6}}, color={0,0,127}));
    connect(add.y, gain1.u)
      annotation (Line(points={{-77,0},{-57.2,0}}, color={0,0,127}));
    connect(gain1.y, gain2.u)
      annotation (Line(points={{-43.4,0},{-35.2,0}}, color={0,0,127}));
    connect(VentHeatGain.port, heatFlowSensor1.port_a) annotation (Line(points={{-46,-40},
            {-20,-40}},                                color={191,0,0}));
    connect(product.u1, gain2.y) annotation (Line(points={{-100,-32},{-100,-12},
            {-14,-12},{-14,0},{-21.4,0}}, color={0,0,127}));
    connect(Vair, product.u2) annotation (Line(points={{-234,-70},{-116,-70},{
            -116,-44},{-100,-44}}, color={0,0,127}));
    connect(product.y, VentHeatGain.Q_flow) annotation (Line(points={{-77,-38},
            {-72,-38},{-72,-40},{-66,-40}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-26},{16,-96}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventilation cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_copy1;

  model office_zone_MPC_est_para_KA_ou44_copy "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=8 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=0.0918)
      annotation (Placement(transformation(extent={{-126,82},{-106,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=4)
      annotation (Placement(transformation(extent={{-196,158},{-176,178}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(C=3e7, T(fixed=
            false, start=293.15))
      annotation (Placement(transformation(extent={{36,122},{56,142}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=4e7, T(
          fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{110,154},{130,174}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=0.004)
      annotation (Placement(transformation(extent={{76,146},{96,166}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-218,166},{
            -218,168},{-198,168}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            168},{-166,168},{-166,166},{-158,166}},
                                         color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-126,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-106,92},{-80,92},{-80,122},{46,122}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,122},{46,122}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,122},{46,122}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,122}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,156},{108,156},{108,154},{120,154}},
                                                              color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,156},{60,156},{
            60,122},{46,122}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix_new.Vve) annotation (Line(points={{-234,-70},{-80,-70},
            {-80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix_new.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-22,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,122}}, color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,-36},
            {-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix_new.Tve) annotation (Line(points={{-107.3,
            -35},{-84.05,-35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{92,
            120},{98,120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_ou44_copy;

  model office_zone_MPC_est_para_KA_copy_for_TEST "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=15 "Zone thermal mass factor [-]";
    parameter Real imass=15 "Zone internal thermal mass factor [-]";
    parameter Real shgc=8 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(4*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,84},{-108,104}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi,
      der_T(fixed=false),
      T(start=298.15))
      annotation (Placement(transformation(extent={{34,124},{54,144}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi,
      T(fixed=false, start=293.15),
      der_T(fixed=false))
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(2*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-20,-48},{0,-28}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,94},{-128,94}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,94},{-80,94},{-80,124},{44,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{44,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{44,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{44,34},{44,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,124},
            {44,124}}, color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix_new.Vve) annotation (Line(points={{-234,-70},{-80,-70},{-80,
            -43},{-61,-43}},      color={0,0,127}));
    connect(airMix_new.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-20,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{0,-38},{44,-38},{44,124}},  color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-10,-48},{-10,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{92,120},
            {98,120},{98,-70},{-50,-70},{-50,-48.8}},      color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -172,90},{-172,-34},{-60.8,-34},{-60.8,-34.6}},  color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-10,-48},
            {-10,-58},{142,-58}},      color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-26},{16,-96}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventilation cooling
")}), experiment(StopTime=604800, Tolerance=1e-11),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
            fillPattern=FillPattern.Solid)}),
      __Dymola_experimentFlags(
        Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
        Evaluate=false,
        OutputCPUtime=false,
        OutputFlatModelica=false));
  end office_zone_MPC_est_para_KA_copy_for_TEST;

  model test_zone
    office_zone_MPC_est_para_KA_copy_for_TEST
      office_zone_MPC_est_para_KA_copy_for_TEST1
      annotation (Placement(transformation(extent={{20,-22},{62,22}})));
    Modelica.Blocks.Sources.Pulse Vair(
      offset=0,
      period=57600,
      startTime=28800,
      width=50,
      amplitude=0)
      annotation (Placement(transformation(extent={{-30,-86},{-10,-66}})));
    Modelica.Blocks.Sources.Pulse Tout(
      startTime=28800,
      period=57600,
      width=50,
      offset=10,
      amplitude=12)
      annotation (Placement(transformation(extent={{-94,2},{-74,22}})));
    Modelica.Blocks.Sources.Pulse Solrad(
      width=50,
      period=57600,
      startTime=28000,
      amplitude=600)
      annotation (Placement(transformation(extent={{-56,-60},{-36,-40}})));
    Modelica.Blocks.Sources.Pulse occ(
      width=50,
      startTime=28800,
      period=57600,
      amplitude=0)
      annotation (Placement(transformation(extent={{-60,-22},{-40,-2}})));
  equation
    connect(occ.y, office_zone_MPC_est_para_KA_copy_for_TEST1.occ) annotation (
        Line(points={{-39,-12},{6,-12},{6,5.8},{19.4,5.8}}, color={0,0,127}));
    connect(Tout.y, office_zone_MPC_est_para_KA_copy_for_TEST1.Tout)
      annotation (Line(points={{-73,12},{-30,12},{-30,10.4},{19.4,10.4}}, color=
           {0,0,127}));
    connect(Vair.y, office_zone_MPC_est_para_KA_copy_for_TEST1.Vair)
      annotation (Line(points={{-9,-76},{14,-76},{14,-4.6},{19.4,-4.6}}, color=
            {0,0,127}));
    connect(Solrad.y, office_zone_MPC_est_para_KA_copy_for_TEST1.solrad)
      annotation (Line(points={{-35,-50},{-14,-50},{-14,-34},{10,-34},{10,16.8},
            {19.2,16.8}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=604800));
  end test_zone;

  model office_zone_MPC_est_para_KA_copy_for_TEST_q "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=15 "Zone thermal mass factor [-]";
    parameter Real imass=15 "Zone internal thermal mass factor [-]";
    parameter Real shgc=8 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(4*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,84},{-108,104}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi,
      der_T(fixed=false),
      T(start=293.15))
      annotation (Placement(transformation(extent={{34,124},{54,144}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi,
      T(fixed=false, start=293.15),
      der_T(fixed=false))
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(2*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput q "heat/cooling power supply"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-20,-48},{0,-28}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,94},{-128,94}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,94},{-80,94},{-80,124},{44,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{44,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{44,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{44,34},{44,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158}},           color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,124},
            {44,124}}, color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{0,-38},{44,-38},{44,124}},  color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-10,-48},{-10,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-10,-48},
            {-10,-58},{142,-58}},      color={0,0,127}));
    connect(q, occHeatGain1.Q_flow)
      annotation (Line(points={{-234,-70},{-88,-70}}, color={0,0,127}));
    connect(occHeatGain1.port, heatFlowSensor1.port_a) annotation (Line(points=
            {{-68,-70},{-26,-70},{-26,-38},{-20,-38}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-26},{16,-96}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventilation cooling
")}), experiment(
        StopTime=604800,
        Tolerance=1e-13,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
            fillPattern=FillPattern.Solid)}),
      __Dymola_experimentFlags(
        Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
        Evaluate=false,
        OutputCPUtime=false,
        OutputFlatModelica=false));
  end office_zone_MPC_est_para_KA_copy_for_TEST_q;

  model TEST4Simulink
    office_zone_MPC_est_para_KA_copy_q office_zone_MPC_est_para_KA_copy_q1
      annotation (Placement(transformation(extent={{-28,-16},{14,28}})));
    Modelica.Blocks.Sources.Constant solrad(k=400)
      annotation (Placement(transformation(extent={{-66,32},{-52,46}})));
    Modelica.Blocks.Sources.Constant Tout(k=15)
      annotation (Placement(transformation(extent={{-66,8},{-52,22}})));
    Modelica.Blocks.Sources.Constant occ(k=0)
      annotation (Placement(transformation(extent={{-66,-14},{-52,0}})));
    Modelica.Blocks.Sources.Constant q(k=-1000)
      annotation (Placement(transformation(extent={{-66,-38},{-52,-24}})));
  equation
    connect(solrad.y, office_zone_MPC_est_para_KA_copy_q1.solrad) annotation (
        Line(points={{-51.3,39},{-38,39},{-38,22.8},{-28.8,22.8}}, color={0,0,
            127}));
    connect(Tout.y, office_zone_MPC_est_para_KA_copy_q1.Tout) annotation (Line(
          points={{-51.3,15},{-40.65,15},{-40.65,16.4},{-28.6,16.4}}, color={0,
            0,127}));
    connect(occ.y, office_zone_MPC_est_para_KA_copy_q1.occ) annotation (Line(
          points={{-51.3,-7},{-40,-7},{-40,11.8},{-28.6,11.8}}, color={0,0,127}));
    connect(q.y, office_zone_MPC_est_para_KA_copy_q1.q) annotation (Line(points=
           {{-51.3,-31},{-36,-31},{-36,1.4},{-28.6,1.4}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(Tolerance=1e-11));
  end TEST4Simulink;

  model office_zone_MPC_est_para_KA_copy_q_demand "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{106,110},{126,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,124},{56,144}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-102,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
    Modelica.Blocks.Sources.Pulse pulse6(
      width=50,
      period=57600,
      amplitude=24,
      startTime=28800)
      annotation (Placement(transformation(extent={{-210,-114},{-196,-100}})));
    Modelica.Blocks.Sources.Pulse pulse5(
      offset=0,
      width=50,
      period=57600,
      amplitude=30,
      startTime=0)
      annotation (Placement(transformation(extent={{-210,-88},{-198,-76}})));
    Modelica.Blocks.Math.Add Tsp
      annotation (Placement(transformation(extent={{-178,-100},{-170,-92}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=0,
      yMin=-1)
      annotation (Placement(transformation(extent={{-156,-102},{-144,-90}})));
    Modelica.Blocks.Math.Gain gain1(k=10000)   annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-110,-70})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{127,
            120}},           color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,124},{46,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{46,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{46,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{108,158},{108,156},{118,156}},
                                                              color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,124},{46,124}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-114,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-91,34},{-12,34}},  color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,124}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{104,120}}, color={0,0,127}));
    connect(occHeatGain1.port, heatFlowSensor1.port_a) annotation (Line(points=
            {{-68,-70},{-46,-70},{-46,-38},{-22,-38}}, color={191,0,0}));
    connect(pulse5.y, Tsp.u1) annotation (Line(points={{-197.4,-82},{-188,-82},
            {-188,-93.6},{-178.8,-93.6}}, color={0,0,127}));
    connect(pulse6.y, Tsp.u2) annotation (Line(points={{-195.3,-107},{-186.65,
            -107},{-186.65,-98.4},{-178.8,-98.4}}, color={0,0,127}));
    connect(fromKelvin.Celsius, PID_temp.u_m) annotation (Line(points={{127,120},
            {132,120},{132,-106},{-150,-106},{-150,-103.2}}, color={0,0,127}));
    connect(Tsp.y, PID_temp.u_s)
      annotation (Line(points={{-169.6,-96},{-157.2,-96}}, color={0,0,127}));
    connect(gain1.y, occHeatGain1.Q_flow)
      annotation (Line(points={{-99,-70},{-88,-70}}, color={0,0,127}));
    connect(PID_temp.y, gain1.u) annotation (Line(points={{-143.4,-96},{-134,
            -96},{-134,-70},{-122,-70}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-11),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_copy_q_demand;

  model office_zone_MPC_est_para_KA_copy_q_2x "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,124},{56,144}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput q "heat flow rate" annotation (
        Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-102,34})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-88,-80},{-68,-60}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor1
      annotation (Placement(transformation(extent={{138,154},{158,174}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{176,154},{196,174}})));
    Modelica.Blocks.Interfaces.RealOutput Tinterior " degCelsius" annotation (
        Placement(transformation(extent={{210,154},{230,174}}),
          iconTransformation(extent={{214,110},{234,130}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,124},{46,124}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,124},{46,124}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,124},{46,124}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,124}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{108,158},{108,156},{118,156}},
                                                              color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,124},{46,124}},
                       color={191,0,0}));
    connect(occ, gain.u)
      annotation (Line(points={{-234,34},{-114,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-91,34},{-12,34}},  color={0,0,127}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,124}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(q, occHeatGain1.Q_flow)
      annotation (Line(points={{-234,-70},{-88,-70}}, color={0,0,127}));
    connect(occHeatGain1.port, heatFlowSensor1.port_a) annotation (Line(points=
            {{-68,-70},{-46,-70},{-46,-38},{-22,-38}}, color={191,0,0}));
    connect(IntMass.port, temperatureSensor1.port) annotation (Line(points={{
            118,156},{132,156},{132,164},{138,164}}, color={191,0,0}));
    connect(temperatureSensor1.T, fromKelvin1.Kelvin)
      annotation (Line(points={{158,164},{174,164}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, Tinterior)
      annotation (Line(points={{197,164},{220,164}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-11),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_copy_q_2x;

  model office_zone_MPC_est_para_KA_4MPC_controller "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=3.0 "External wall thermal resistance";
    parameter Real RInt=0.05 "Internal wall thermal resistance";
    parameter Real tmass=25 "Zone thermal mass factor [-]";
    parameter Real imass=17.89 "Zone internal thermal mass factor [-]";
    parameter Real shgc=7.94 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.039 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 120 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,82},{-108,102}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{72,110},{92,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]" annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}), iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=486.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{36,118},{56,138}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,144},{96,164}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-16,-48},{4,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor1
      annotation (Placement(transformation(extent={{140,146},{160,166}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{174,146},{194,166}})));
    Modelica.Blocks.Interfaces.RealOutput Tinterior " degCelsius" annotation (
        Placement(transformation(extent={{216,146},{236,166}}),
          iconTransformation(extent={{216,56},{236,76}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-134,90},{-134,92},{-128,92}},
                                                     color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,92},{-80,92},{-80,118},{46,118}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,118},{46,118}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{72,120},{60,120},{60,118},{46,118}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,118}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}}, color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,154},{108,154},{108,156},{118,156}},
                                                              color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,154},{60,154},{
            60,118},{46,118}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix_new.Vve) annotation (Line(points={{-234,-70},{-80,-70},
            {-80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix_new.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-16,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{4,-38},{46,-38},{46,118}},  color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,-36},
            {-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix_new.Tve) annotation (Line(points={{-107.3,
            -35},{-84.05,-35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-6,-48},{-6,-58},{142,-58}},   color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-6,-48},{-6,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{92,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{92,
            120},{98,120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    connect(temperatureSensor1.port, IntMass.port)
      annotation (Line(points={{140,156},{118,156}}, color={191,0,0}));
    connect(temperatureSensor1.T, fromKelvin1.Kelvin)
      annotation (Line(points={{160,156},{172,156}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, Tinterior)
      annotation (Line(points={{195,156},{226,156}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{200,220}},
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
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
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
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{200,
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
  end office_zone_MPC_est_para_KA_4MPC_controller;
end office_model;
