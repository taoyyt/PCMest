within OU44;
model energyhub

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=C)
    annotation (Placement(transformation(extent={{8,8},{28,28}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor(R=R)
    annotation (Placement(transformation(extent={{-28,-2},{-8,18}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-74,-2},{-54,18}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-76,-48},{-56,-28}})));
  Modelica.Blocks.Interfaces.RealInput Tout
    annotation (Placement(transformation(extent={{-140,-12},{-100,28}})));
  Modelica.Blocks.Interfaces.RealOutput Tzone
    annotation (Placement(transformation(extent={{100,-2},{120,18}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{38,-2},{58,18}})));
  Modelica.Blocks.Interfaces.RealInput qhvac
    annotation (Placement(transformation(extent={{-140,-58},{-100,-18}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
  Modelica.Blocks.Interfaces.RealOutput q
    annotation (Placement(transformation(extent={{100,-62},{120,-42}})));
  parameter Real R=0.01;
  parameter Real C=1e6;

equation
  connect(thermalResistor.port_b, heatCapacitor.port)
    annotation (Line(points={{-8,8},{18,8}}, color={191,0,0}));
  connect(prescribedTemperature.port, thermalResistor.port_a)
    annotation (Line(points={{-54,8},{-28,8}}, color={191,0,0}));
  connect(heatCapacitor.port, temperatureSensor.port)
    annotation (Line(points={{18,8},{38,8}}, color={191,0,0}));
  connect(temperatureSensor.T, Tzone)
    annotation (Line(points={{58,8},{110,8}}, color={0,0,127}));
  connect(Tzone, Tzone)
    annotation (Line(points={{110,8},{110,8}}, color={0,0,127}));
  connect(qhvac, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-120,-38},{-76,-38}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, heatFlowSensor.port_a)
    annotation (Line(points={{-56,-38},{-22,-38}}, color={191,0,0}));
  connect(heatFlowSensor.port_b, heatCapacitor.port)
    annotation (Line(points={{-2,-38},{18,-38},{18,8}}, color={191,0,0}));
  connect(heatFlowSensor.Q_flow, q)
    annotation (Line(points={{-12,-48},{-12,-52},{110,-52}}, color={0,0,127}));
  connect(Tout, prescribedTemperature.T)
    annotation (Line(points={{-120,8},{-76,8}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end energyhub;
