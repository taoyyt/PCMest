within ;
model PCM_model

  Buildings.HeatTransfer.Conduction.SingleLayer lay1(material=matPCM,
    A=0.135,
    T_a_start=295.15,
    T_b_start=295.15)
    annotation (Placement(transformation(extent={{18,-12},{38,8}})));
  Buildings.HeatTransfer.Conduction.SingleLayer lay2(material=matPCM,
    A=0.135,
    T_a_start=295.15,
    T_b_start=295.15)
    annotation (Placement(transformation(extent={{42,-12},{62,8}})));
  Buildings.HeatTransfer.Conduction.SingleLayer lay3(material=matPCM,
    A=0.135,
    T_a_start=295.15,
    T_b_start=295.15)
    annotation (Placement(transformation(extent={{66,-12},{86,8}})));
  Buildings.HeatTransfer.Conduction.SingleLayer lay4(material=matPCM,
    A=0.135,
    T_a_start=295.15,
    T_b_start=295.15)
    annotation (Placement(transformation(extent={{94,-12},{114,8}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature prescribedTemperature
    annotation (Placement(transformation(extent={{-48,-12},{-28,8}})));
  parameter Buildings.HeatTransfer.Data.SolidsPCM.Generic matPCM(
    k=0.6,
    LHea=140000,
    c=2000,
    TSol=273.15 + 12,
    TLiq=273.15 + 24,
    x=0.00375,
    d=1450,
    nStaRef=1)
    annotation (Placement(transformation(extent={{-46,30},{-26,50}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{118,-8},{130,4}})));
  Modelica.Thermal.HeatTransfer.Components.Convection convection
    annotation (Placement(transformation(extent={{6,-12},{-14,8}})));
  Modelica.Blocks.Sources.Constant const1(k=2.753)
    annotation (Placement(transformation(extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={-4,28})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-74,-8},{-62,4}})));
  Modelica.Blocks.Sources.CombiTimeTable Tair(
    tableOnFile=true,
    timeScale=60,
    tableName="Tair1",
    fileName="C:/Users/taoy/Desktop/Modelica PCM model/Tair1.txt") annotation (
      Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=-90,
        origin={-80,14})));
equation
  connect(lay1.port_b, lay2.port_a)
    annotation (Line(points={{38,-2},{42,-2}}, color={191,0,0}));
  connect(lay2.port_b, lay3.port_a)
    annotation (Line(points={{62,-2},{66,-2}}, color={191,0,0}));
  connect(lay3.port_b, lay4.port_a)
    annotation (Line(points={{86,-2},{94,-2}}, color={191,0,0}));
  connect(lay4.port_b, temperatureSensor.port)
    annotation (Line(points={{114,-2},{118,-2}}, color={191,0,0}));
  connect(prescribedTemperature.port, convection.fluid)
    annotation (Line(points={{-28,-2},{-14,-2}}, color={191,0,0}));
  connect(convection.solid, lay1.port_a)
    annotation (Line(points={{6,-2},{18,-2}}, color={191,0,0}));
  connect(const1.y, convection.Gc)
    annotation (Line(points={{-4,21.4},{-4,8}}, color={0,0,127}));
  connect(Tair.y[1], toKelvin.Celsius)
    annotation (Line(points={{-80,5.2},{-80,-2},{-75.2,-2}}, color={0,0,127}));
  connect(toKelvin.Kelvin, prescribedTemperature.T)
    annotation (Line(points={{-61.4,-2},{-50,-2}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{140,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{140,
            100}})),
    uses(Buildings(version="6.0.0"), Modelica(version="3.2.2"),
      AixLib(version="0.7.3")));
end PCM_model;
