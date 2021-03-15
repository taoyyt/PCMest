within OU44.Components;
model T_simple "Thermal submodel"

  parameter Real CExt=10000000
    "External heat capacity per unit area of external walls [J/(K m2)]";
  parameter Real CAir=1200
    "Zone thermal capacity per unit volume [J/(K m3)]";
  parameter Modelica.SIunits.Temp_C TInitial=20
    "Initial temperature of indoor air";
  //parameter Modelica.SIunits.Temp_K TInitialK=TInitial+273.15
  //  "Initial temperature of indoor air";
  parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
  parameter Modelica.SIunits.Volume Vi=100 "Indoor air volume [m3]";
  parameter Real RExt=0.01 "External wall termal resistance (outer layer)";
  parameter Real RInt=0.01 "External wall termal resistance (inner layer)";
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor REx2(R=RInt/3/
        Vi^(2/3))
    annotation (Placement(transformation(extent={{-70,-12},{-40,18}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
    der_T(start=0),
    T(fixed=true, start=TInitial + 273.15),
    C=Vi*CAir)
    annotation (Placement(transformation(extent={{34,2},{106,74}})));
  Modelica.Blocks.Interfaces.RealInput Tamb "[degC]"
    annotation (Placement(transformation(extent={{-350,-104},{-310,-64}})));
  Modelica.Blocks.Interfaces.RealOutput TInt "[degC]"
                                             annotation (Placement(
        transformation(extent={{270,-20},{310,20}}), iconTransformation(extent={{270,-26},
            {320,24}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow solarHeat
    annotation (Placement(transformation(extent={{-166,76},{-130,112}})));
  Modelica.Blocks.Interfaces.RealInput solRad "[W]"
    annotation (Placement(transformation(extent={{-348,74},{-308,114}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatingHeat
    annotation (Placement(transformation(extent={{-74,158},{-38,194}})));
  Modelica.Blocks.Interfaces.RealInput heating "[W]"
    annotation (Placement(transformation(extent={{-348,156},{-308,196}})));
  Modelica.Blocks.Interfaces.RealInput Vve "[m3/h]"
    annotation (Placement(transformation(extent={{-346,-192},{-306,-152}})));
  Modelica.Blocks.Interfaces.RealInput Tve "[degC]"
    annotation (Placement(transformation(extent={{-346,-294},{-306,-254}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature TAmbPrescribed
    annotation (Placement(transformation(extent={{-282,-14},{-252,16}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor tempSensor
    annotation (Placement(transformation(extent={{108,-18},{146,20}})));

  AirMix ventilationHeatLoad annotation (Placement(transformation(extent=
            {{-162,-170},{-64,-102}})));

  Modelica.Blocks.Tables.CombiTable1D MetHeat(table=[20.,100.; 34.,0.])
    annotation (Placement(transformation(extent={{190,186},{134,242}})));
  Modelica.Blocks.Interfaces.RealInput persons "[-]"
    annotation (Placement(transformation(extent={{-348,224},{-308,264}})));
  Modelica.Blocks.Math.Product product1
    annotation (Placement(transformation(extent={{14,-14},{-14,14}},
        rotation=90,
        origin={58,186})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occupantHeat
    annotation (Placement(transformation(
        extent={{21,-21},{-21,21}},
        rotation=90,
        origin={57,107})));

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor CExternal(
    T(displayUnit="K",
      start=TInitial + 273.15,
      fixed=false),
    der_T(fixed=true, start=0),
    C=CExt*3*Vi^(2/3))
    annotation (Placement(transformation(extent={{-164,2},{-106,60}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor REx3(R=RExt/3/
        Vi^(2/3))
    annotation (Placement(transformation(extent={{-226,-14},{-196,16}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin annotation (Placement(
        transformation(
        extent={{-25,-25},{25,25}},
        rotation=90,
        origin={-295,-45})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
    annotation (Placement(transformation(extent={{194,-22},{240,24}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2 annotation (
      Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=0,
        origin={-276,-274})));
  Modelica.Blocks.Math.Gain gain(k=shgc)
    annotation (Placement(transformation(extent={{-258,84},{-238,104}})));
  Modelica.Blocks.Interfaces.RealOutput qair
    "Total heat supplied to indoor air [W]"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={220,286})));
  Modelica.Blocks.Continuous.Integrator integratorHeat(
    y_start=0,
    k=1,
    initType=Modelica.Blocks.Types.Init.InitialState)
    annotation (Placement(transformation(extent={{-186,200},{-166,220}})));
  Modelica.Blocks.Interfaces.RealOutput QHeat
    "Total heat supplied by heating system [J]"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-140,288})));
  Modelica.Blocks.Continuous.Integrator integratorVent(
    y_start=0,
    k=1,
    initType=Modelica.Blocks.Types.Init.InitialState)
    annotation (Placement(transformation(extent={{82,-170},{102,-150}})));
  Modelica.Blocks.Interfaces.RealOutput QVe
    "Total heat supplied by ventilation [J]" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={40,286})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensorVent
    annotation (Placement(transformation(extent={{-30,-146},{-10,-126}})));
  Modelica.Blocks.Continuous.Derivative derivative(
    initType=Modelica.Blocks.Types.Init.SteadyState,
    x_start=293,
    k=1,
    T=0.01)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={220,90})));

  Modelica.Blocks.Math.Gain gain1(k=Vi*1.2*1005) annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={220,132})));
  Modelica.Blocks.Math.Gain gain2(k=occH)        annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={94,214})));
  parameter Real occH=1 "Occ heat gains multiplier";
  Modelica.Blocks.Interfaces.RealOutput qheat
    "Total heat supplied by heating system [W]" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-210,288})));
equation
  connect(air.port, tempSensor.port)
    annotation (Line(points={{70,2},{108,2},{108,1}}, color={191,0,0}));
  connect(Vve, ventilationHeatLoad.Vve) annotation (Line(points={{-326,-172},
          {-274,-172},{-274,-153},{-166.9,-153}}, color={0,0,127}));
  connect(ventilationHeatLoad.Ti, tempSensor.T) annotation (Line(points={{-113,
          -172.72},{-112,-172.72},{-112,-186},{152,-186},{152,1},{146,1}},
        color={0,0,127}));
  connect(TInt, TInt)
    annotation (Line(points={{290,0},{290,0}}, color={0,0,127}));
  connect(heating, heatingHeat.Q_flow) annotation (Line(points={{-328,176},{-180,
          176},{-74,176}}, color={0,0,127}));
  connect(persons, product1.u1) annotation (Line(points={{-328,244},{51,
          244},{51,202.8},{49.6,202.8}},
                                color={0,0,127}));
  connect(product1.y, occupantHeat.Q_flow) annotation (Line(points={{58,
          170.6},{58,128},{57,128}},
                              color={0,0,127}));
  connect(REx2.port_a, CExternal.port)
    annotation (Line(points={{-70,3},{-70,2},{-135,2}},   color={191,0,0}));
  connect(occupantHeat.port,air. port) annotation (Line(points={{57,86},{
          57,86},{57,76},{26,76},{26,2},{70,2}},
                                               color={191,0,0}));
  connect(heatingHeat.port,air. port) annotation (Line(points={{-38,176},{26,
          176},{26,2},{70,2}}, color={191,0,0}));
  connect(Tamb, toKelvin.Celsius) annotation (Line(points={{-330,-84},{-295,
          -84},{-295,-75}},
                         color={0,0,127}));
  connect(toKelvin.Kelvin, TAmbPrescribed.T) annotation (Line(points={{-295,-17.5},
          {-295,0.9},{-285,0.9},{-285,1}}, color={0,0,127}));
  connect(tempSensor.T, fromKelvin.Kelvin)
    annotation (Line(points={{146,1},{146,1},{189.4,1}}, color={0,0,127}));
  connect(TInt, fromKelvin.Celsius)
    annotation (Line(points={{290,0},{268,0},{268,1},{242.3,1}},
                                               color={0,0,127}));
  connect(Tve, toKelvin2.Celsius) annotation (Line(points={{-326,-274},{-300,
          -274}},         color={0,0,127}));
  connect(toKelvin2.Kelvin, ventilationHeatLoad.Tve) annotation (Line(points={{-254,
          -274},{-250,-274},{-250,-124.44},{-165.92,-124.44}},          color=
         {0,0,127}));
  connect(fromKelvin.Celsius, MetHeat.u[1]) annotation (Line(points={{242.3,1},{
          254,1},{254,214},{195.6,214}}, color={0,0,127}));
  connect(solRad, gain.u)
    annotation (Line(points={{-328,94},{-260,94}}, color={0,0,127}));
  connect(solarHeat.Q_flow, gain.y)
    annotation (Line(points={{-166,94},{-237,94}}, color={0,0,127}));
  connect(heating, integratorHeat.u) annotation (Line(points={{-328,176},
          {-210,176},{-210,210},{-188,210}},
                                  color={0,0,127}));
  connect(integratorHeat.y,QHeat)
    annotation (Line(points={{-165,210},{-140,210},{-140,288}},
                                                    color={0,0,127}));
  connect(ventilationHeatLoad.port_b, heatFlowSensorVent.port_a) annotation (
      Line(points={{-64,-136},{-48,-136},{-30,-136}}, color={191,0,0}));
  connect(heatFlowSensorVent.port_b,air. port) annotation (Line(points={{-10,-136},
          {26,-136},{26,2},{70,2}}, color={191,0,0}));
  connect(heatFlowSensorVent.Q_flow, integratorVent.u) annotation (Line(points={{-20,
          -146},{-20,-146},{-20,-160},{80,-160}},      color={0,0,127}));
  connect(fromKelvin.Kelvin, derivative.u) annotation (Line(points={{
          189.4,1},{189.4,24.5},{220,24.5},{220,78}}, color={0,0,127}));
  connect(integratorVent.y, QVe) annotation (Line(points={{103,-160},{114,
          -160},{114,260},{40,260},{40,286}}, color={0,0,127}));
  connect(REx3.port_a, TAmbPrescribed.port)
    annotation (Line(points={{-226,1},{-226,1},{-252,1}}, color={191,0,0}));
  connect(derivative.y, gain1.u)
    annotation (Line(points={{220,101},{220,101},{220,120}}, color={0,0,127}));
  connect(gain1.y, qair)
    annotation (Line(points={{220,143},{220,286}}, color={0,0,127}));
  connect(solarHeat.port,air. port) annotation (Line(points={{-130,94},{-52,94},
          {-52,76},{26,76},{26,2},{70,2}}, color={191,0,0}));
  connect(REx3.port_b, CExternal.port) annotation (Line(points={{-196,1},{-192,1},
          {-192,2},{-135,2}}, color={191,0,0}));
  connect(gain2.y, product1.u2) annotation (Line(points={{83,214},{66,214},
          {66,202.8},{66.4,202.8}}, color={0,0,127}));
  connect(gain2.u, MetHeat.y[1])
    annotation (Line(points={{106,214},{131.2,214}}, color={0,0,127}));
  connect(qheat, heating) annotation (Line(points={{-210,288},{-210,288},
          {-210,176},{-328,176}}, color={0,0,127}));
  connect(REx2.port_b,air. port) annotation (Line(points={{-40,3},{-17,3},{-17,2},
          {70,2}}, color={191,0,0}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-320,
            -300},{280,280}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-320,-300},{280,280}}), graphics={Rectangle(extent={{-320,
              280},{280,-300}},
                          lineColor={28,108,200}),
        Ellipse(
          extent={{168,32},{234,-34}},
          lineColor={0,127,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Sphere),
        Ellipse(
          extent={{8,34},{74,-32}},
          lineColor={0,127,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Sphere),
        Ellipse(
          extent={{-148,36},{-82,-30}},
          lineColor={0,127,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Sphere),
        Ellipse(
          extent={{172,210},{238,144}},
          lineColor={0,127,255},
          fillColor={175,175,175},
          fillPattern=FillPattern.Sphere),
        Rectangle(
          extent={{90,10},{150,-10}},
          lineColor={0,127,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-70,10},{-10,-10}},
          lineColor={0,127,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-30,10},{30,-10}},
          lineColor={0,127,255},
          origin={202,90},
          rotation=90,
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Line(points={{-170,30}}, color={0,127,255}),
        Rectangle(
          extent={{-224,10},{-164,-10}},
          lineColor={0,127,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-42,-184},{124,-264}},
          lineColor={255,0,0},
          pattern=LinePattern.Dash),
        Rectangle(
          extent={{-266,-84},{-100,-164}},
          lineColor={255,0,0},
          pattern=LinePattern.Dash),
        Line(points={{-162,-300},{-162,-228},{-44,-228}}, color={28,108,
              200}),
        Line(points={{-322,-274},{-240,-274},{-240,-166}}, color={28,108,
              200}),
        Line(points={{-312,-172},{-282,-172},{-282,-130},{-266,-130}},
            color={28,108,200}),
        Line(points={{-320,-84},{-290,-84},{-290,0},{-232,0}}, color={28,
              108,200}),
        Line(points={{-14,-184}}, color={0,127,255}),
        Line(points={{126,-228},{200,-228},{200,-46}}, color={255,0,0}),
        Line(points={{-100,-128},{200,-128}}, color={255,0,0}),
        Line(points={{-274,238},{-282,220}}, color={255,0,0}),
        Line(points={{-268,220},{-274,238},{-274,258}}, color={255,0,0}),
        Line(points={{-284,250},{-264,250}}, color={255,0,0}),
        Ellipse(extent={{-280,270},{-268,258}}, lineColor={255,0,0}),
        Ellipse(extent={{-282,108},{-256,82}}, lineColor={255,0,0}),
        Line(points={{-264,110},{-254,132}}, color={255,0,0}),
        Line(points={{-276,110},{-284,132}}, color={255,0,0}),
        Line(points={{-256,100},{-234,110}}, color={255,0,0}),
        Line(points={{-256,88},{-234,82}}, color={255,0,0}),
        Line(points={{-266,80},{-262,60}}, color={255,0,0}),
        Line(points={{-278,84},{-294,66}}, color={255,0,0}),
        Line(points={{-284,94},{-304,92}}, color={255,0,0}),
        Line(points={{-282,104},{-298,114}}, color={255,0,0}),
        Rectangle(extent={{-292,188},{-234,160}}, lineColor={255,0,0}),
        Line(points={{-242,182},{-242,166}}, color={255,0,0}),
        Line(points={{-250,182},{-250,166}}, color={255,0,0}),
        Line(points={{-258,182},{-258,166}}, color={255,0,0}),
        Line(points={{-266,182},{-266,166}}, color={255,0,0}),
        Line(points={{-274,182},{-274,166}}, color={255,0,0}),
        Line(points={{-282,182},{-282,166}}, color={255,0,0}),
        Line(points={{-128,100}}, color={255,0,0}),
        Ellipse(
          extent={{-264,-192},{-214,-242}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(points={{-258,-232},{-220,-232},{-240,-192},{-258,-232}},
            lineColor={28,108,200}),
        Ellipse(
          extent={{-25,27},{25,-27}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          origin={-105,-229},
          rotation=-90),
        Polygon(
          points={{-19,-20},{19,-20},{-1,20},{-19,-20}},
          lineColor={28,108,200},
          origin={-101,-228},
          rotation=-90),
        Line(points={{270,0},{244,0}}, color={0,127,255}),
        Line(
          points={{-222,100},{-96,102},{28,178},{166,178}},
          color={255,0,0},
          smooth=Smooth.Bezier),
        Line(points={{-218,242},{-180,242},{-180,172},{-220,172}}, color=
              {255,0,0}),
        Line(
          points={{-180,206},{-72,206},{168,34}},
          color={255,0,0},
          smooth=Smooth.Bezier),
        Line(
          points={{-234,-102},{-214,-144},{-194,-102}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{-148,-102},{-174,-102},{-174,-142},{-148,-142}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{-174,-122},{-160,-122}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{10,-206},{10,-246}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{22,-246},{22,-206},{42,-246},{42,-206}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{54,-246},{54,-206},{72,-206}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{54,-226},{64,-226}},
          color={0,0,127},
          thickness=0.5)}));
end T_simple;
