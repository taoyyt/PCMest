within OU44;
model ZoneR4C3_RBC "single zone model"
    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
    "Medium model";

  Modelica.Blocks.Interfaces.RealInput solrad "solar radiation[W]"
    annotation (Placement(transformation(extent={{-156,100},{-140,116}})));
  Modelica.Blocks.Math.Gain solarcoefficient(k=shgc)
    annotation (Placement(transformation(extent={{-120,104},{-112,112}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-92,102},{-80,114}})));
  Modelica.Blocks.Interfaces.RealInput Tout "ambient temperature (deg C)"
    annotation (Placement(transformation(
        extent={{-9,-9},{9,9}},
        rotation=0,
        origin={-149,81})));
  Modelica.Blocks.Interfaces.RealInput occ "Occupancy"
    annotation (Placement(transformation(extent={{-158,30},{-142,46}})));
  Modelica.Blocks.Interfaces.RealInput dpos "damperpositon(0-100)%"
    annotation (Placement(transformation(extent={{-160,-66},{-138,-44}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-132,74},{-118,88}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-100,74},{-86,88}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Re(R=RExt/(3*Vi^(2/3)))
    "thermal resistence of external wall"
    annotation (Placement(transformation(extent={{-76,70},{-58,88}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_air(C=Cair, T(
      fixed=true,
      displayUnit="K",
      start=TairInit + 273.15)) "air heat capacity"
    annotation (Placement(transformation(extent={{4,90},{22,108}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rinternal(R=Rin)
    "thermal resistance of internal wall/furniture"
    annotation (Placement(transformation(extent={{40,102},{52,114}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_internal_wall(C=Cin)
    annotation (Placement(transformation(extent={{60,108},{74,122}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rinfiltration(R=Rinf)
    "thermal resistance of infltration"
    annotation (Placement(transformation(extent={{-28,92},{-18,102}})));
  Modelica.Blocks.Sources.Constant const1(k=Vinf)
    annotation (Placement(transformation(extent={{-8,58},{-18,68}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{34,82},{46,94}})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
    annotation (Placement(transformation(extent={{68,82},{80,94}})));
  Modelica.Blocks.Interfaces.RealOutput T "indoor temperature"
    annotation (Placement(transformation(extent={{120,78},{140,98}})));
  Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,
        0.])
    annotation (Placement(transformation(extent={{-88,-2},{-74,12}})));
  Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
    annotation (Placement(transformation(extent={{-66,0},{-56,10}})));
  Modelica.Blocks.Math.Product hmltp
    annotation (Placement(transformation(extent={{-38,20},{-26,32}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
    annotation (Placement(transformation(extent={{-10,20},{2,32}})));
  Modelica.Blocks.Sources.Constant const2(k=Tve)
    annotation (Placement(transformation(extent={{-144,-34},{-132,-22}})));
  Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
        extent={{5,-5},{-5,5}},
        rotation=180,
        origin={-105,-55})));
  Modelica.Blocks.Math.Gain ventilation(k=maxVent)
    annotation (Placement(transformation(extent={{-72,-62},{-58,-48}})));
  Modelica.Blocks.Continuous.Integrator integrator(k=1/3600)
    annotation (Placement(transformation(extent={{86,-62},{100,-48}})));
  Modelica.Blocks.Interfaces.RealOutput vetot "aggregated ventilation flow[m3]"
    annotation (Placement(transformation(extent={{120,-66},{142,-44}})));
  Modelica.Blocks.Interfaces.RealOutput verate
    "ventilation air supply rate[m3/h]"
    annotation (Placement(transformation(extent={{120,-42},{142,-20}})));
  Modelica.Blocks.Math.Gain heating(k=maxHeat)
    annotation (Placement(transformation(extent={{-72,-94},{-54,-76}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(extent={{-40,-92},{-26,-78}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-18,-92},{-4,-78}})));
  Modelica.Blocks.Interfaces.RealOutput Qrad
    "aggregated energy supplied by radiators[Wh]"
    annotation (Placement(transformation(extent={{120,-90},{142,-68}})));
  Modelica.Blocks.Interfaces.RealOutput qrad "heating power of radiators[W]"
    annotation (Placement(transformation(extent={{120,-108},{142,-86}})));
  Modelica.Blocks.Continuous.Integrator integrator1(k=1/3600)
    annotation (Placement(transformation(extent={{86,-86},{100,-72}})));
  Modelica.Blocks.Math.Gain scale1(k=0.01)
                                          annotation (Placement(transformation(
        extent={{7,-7},{-7,7}},
        rotation=180,
        origin={-103,-85})));
  OU44.Components.CO2
                 cO2_1(
    Vi=Vi,
    CO2PpmInitial=CO2Init,
    CO2PerPerson=CO2pp,
    CO2Neutral=CO2n)
    annotation (Placement(transformation(extent={{2,0},{14,12}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-16,10},{-8,18}})));
  Modelica.Blocks.Interfaces.RealOutput CO2 "[ppm]"
    annotation (Placement(transformation(extent={{120,-14},{142,8}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin    toKelvin1
    annotation (Placement(transformation(extent={{-126,-34},{-114,-22}})));
  Modelica.Blocks.Math.Max max
    annotation (Placement(transformation(extent={{-104,-30},{-94,-20}})));
  OU44.Components.AirMix
                    airMix
    annotation (Placement(transformation(extent={{-6,-28},{8,-14}})));

  parameter Real shgc=1 "solar heat gains coefficient";
  parameter Real Vinf=50 "inflitration rate";
  parameter Real TairInit=20 "Initial temperature of indoor air [deg C]";
  parameter Real CO2n=400 "CO2 neutral level";
  parameter Real CO2pp=0.02 "CO2 generation per person [m3/h]";
  parameter Real maxVent=2000 "Maximum ventilation flowrate[m3/h]";
  parameter Real maxHeat=5000 "Heating power of radiator [W]";
  parameter Real Tve=21 "Ventilation air temperature [deg C]";
  parameter Real Vi=300 "Air volume [m3]";
  parameter Real occheff=1. "Occupant heat generation effectiveness";
  parameter Real CO2Init=400 "Initial CO2 concentration[ppmv]";
  parameter Real RExt=1.0 "External wall thermal resistance";
  parameter Real Cwall=1e6 "Heat capacity of external walls";
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rindoor(R=Ri)
    "Indoor thermal convective/conductive resistance"
    annotation (Placement(transformation(extent={{-30,72},{-16,86}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_wall(C=Cwall)
    "heat capacity of external wall" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={-50,62})));
  parameter Real Rinf=10 "Thermal resistance of inflitration";
  parameter Real Ri=10 "Indoor thermal resistance";
  parameter Real Rin=10 "thermal resistance of indoor funiture";
  parameter Real Cin=1e6 "heat capacity of indoor furniture";
  parameter Real Cair=1e6 "heat capacity of indoor air";
  Modelica.Blocks.Interfaces.RealInput vpos "valver position [0-100]"
    annotation (Placement(transformation(extent={{-162,-96},{-140,-74}})));
equation
  connect(solrad, solarcoefficient.u)
    annotation (Line(points={{-148,108},{-120.8,108}}, color={0,0,127}));
  connect(solarcoefficient.y, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-111.6,108},{-92,108}}, color={0,0,127}));
  connect(Tout, toKelvin.Celsius)
    annotation (Line(points={{-149,81},{-133.4,81}}, color={0,0,127}));
  connect(toKelvin.Kelvin, prescribedTemperature.T)
    annotation (Line(points={{-117.3,81},{-101.4,81}}, color={0,0,127}));
  connect(prescribedTemperature.port, Re.port_a)
    annotation (Line(points={{-86,81},{-82,81},{-82,79},{-76,79}},
                                                 color={191,0,0}));
  connect(HC_air.port, Rinternal.port_a) annotation (Line(points={{13,90},{26,
          90},{26,108},{40,108}},
                              color={191,0,0}));
  connect(Rinternal.port_b, HC_internal_wall.port)
    annotation (Line(points={{52,108},{67,108}}, color={191,0,0}));
  connect(prescribedTemperature.port, Rinfiltration.port_a) annotation (Line(
        points={{-86,81},{-82,81},{-82,97},{-28,97}}, color={191,0,0}));
  connect(Rinfiltration.port_b,HC_air. port) annotation (Line(points={{-18,97},
          {0,97},{0,90},{13,90}},color={191,0,0}));
  connect(HC_air.port, temperatureSensor.port)
    annotation (Line(points={{13,90},{24,90},{24,88},{34,88}},
                                               color={191,0,0}));
  connect(temperatureSensor.T, fromKelvin.Kelvin)
    annotation (Line(points={{46,88},{66.8,88}}, color={0,0,127}));
  connect(T, T) annotation (Line(points={{130,88},{130,88}}, color={0,0,127}));
  connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{46,
          88},{54,88},{54,-10},{-94,-10},{-94,5},{-89.4,5}}, color={0,0,127}));
  connect(MetabolicHeat.y[1], occeffectiv.u)
    annotation (Line(points={{-73.3,5},{-67,5}}, color={0,0,127}));
  connect(occeffectiv.y, hmltp.u2) annotation (Line(points={{-55.5,5},{-42.75,5},
          {-42.75,22.4},{-39.2,22.4}}, color={0,0,127}));
  connect(occ, hmltp.u1) annotation (Line(points={{-150,38},{-92.5,38},{-92.5,
          29.6},{-39.2,29.6}}, color={0,0,127}));
  connect(hmltp.y, occHeatGain.Q_flow)
    annotation (Line(points={{-25.4,26},{-10,26}}, color={0,0,127}));
  connect(occHeatGain.port,HC_air. port)
    annotation (Line(points={{2,26},{13,26},{13,90}}, color={191,0,0}));
  connect(dpos, scale.u)
    annotation (Line(points={{-149,-55},{-111,-55}}, color={0,0,127}));
  connect(scale.y, ventilation.u)
    annotation (Line(points={{-99.5,-55},{-73.4,-55}}, color={0,0,127}));
  connect(ventilation.y, integrator.u)
    annotation (Line(points={{-57.3,-55},{84.6,-55}}, color={0,0,127}));
  connect(integrator.y, vetot)
    annotation (Line(points={{100.7,-55},{131,-55}}, color={0,0,127}));
  connect(ventilation.y, verate) annotation (Line(points={{-57.3,-55},{74,-55},
          {74,-31},{131,-31}}, color={0,0,127}));
  connect(heating.y, prescribedHeatFlow1.Q_flow) annotation (Line(points={{-53.1,
          -85},{-40,-85}},                           color={0,0,127}));
  connect(prescribedHeatFlow1.port, heatFlowSensor.port_a)
    annotation (Line(points={{-26,-85},{-18,-85}}, color={191,0,0}));
  connect(heatFlowSensor.port_b,HC_air. port) annotation (Line(points={{-4,-85},
          {20,-85},{20,90},{13,90}}, color={191,0,0}));
  connect(heatFlowSensor.Q_flow, integrator1.u) annotation (Line(points={{-11,
          -92},{44,-92},{44,-79},{84.6,-79}}, color={0,0,127}));
  connect(integrator1.y, Qrad) annotation (Line(points={{100.7,-79},{110.35,-79},
          {110.35,-79},{131,-79}}, color={0,0,127}));
  connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-11,-92},{88,
          -92},{88,-97},{131,-97}}, color={0,0,127}));
  connect(ventilation.y, add.u2) annotation (Line(points={{-57.3,-55},{-34,-55},
          {-34,11.6},{-16.8,11.6}}, color={0,0,127}));
  connect(const1.y, add.u1) annotation (Line(points={{-18.5,63},{-22,63},{-22,
          16.4},{-16.8,16.4}}, color={0,0,127}));
  connect(add.y, cO2_1.Vve) annotation (Line(points={{-7.6,14},{-2,14},{-2,9.6},
          {1.76,9.6}}, color={0,0,127}));
  connect(occ, cO2_1.persons) annotation (Line(points={{-150,38},{-50,38},{-50,
          2.4},{1.76,2.4}}, color={0,0,127}));
  connect(cO2_1.CO2, CO2) annotation (Line(points={{14.36,6},{68,6},{68,-3},{
          131,-3}}, color={0,0,127}));
  connect(toKelvin1.Kelvin, max.u2)
    annotation (Line(points={{-113.4,-28},{-105,-28}}, color={0,0,127}));
  connect(toKelvin.Kelvin, max.u1) annotation (Line(points={{-117.3,81},{-110,
          81},{-110,-22},{-105,-22}}, color={0,0,127}));
  connect(max.y, airMix.Tve) annotation (Line(points={{-93.5,-25},{-17.75,-25},
          {-17.75,-18.62},{-6.56,-18.62}}, color={0,0,127}));
  connect(ventilation.y, airMix.Vve) annotation (Line(points={{-57.3,-55},{-14,
          -55},{-14,-24.5},{-6.7,-24.5}}, color={0,0,127}));
  connect(HC_air.port, airMix.port_b) annotation (Line(points={{13,90},{20,90},
          {20,-21},{8,-21}}, color={191,0,0}));
  connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{46,88},{54,
          88},{54,-42},{1,-42},{1,-28.56}}, color={0,0,127}));
  connect(Re.port_b, Rindoor.port_a)
    annotation (Line(points={{-58,79},{-30,79}}, color={191,0,0}));
  connect(Rindoor.port_b,HC_air. port) annotation (Line(points={{-16,79},{0,79},
          {0,90},{13,90}}, color={191,0,0}));
  connect(Re.port_b,HC_wall. port)
    annotation (Line(points={{-58,79},{-50,79},{-50,69}}, color={191,0,0}));
  connect(prescribedHeatFlow.port, Rindoor.port_a) annotation (Line(points={{
          -80,108},{-42,108},{-42,79},{-30,79}}, color={191,0,0}));
  connect(fromKelvin.Celsius, T)
    annotation (Line(points={{80.6,88},{130,88}}, color={0,0,127}));
  connect(const2.y, toKelvin1.Celsius)
    annotation (Line(points={{-131.4,-28},{-127.2,-28}}, color={0,0,127}));
  connect(vpos, scale1.u)
    annotation (Line(points={{-151,-85},{-111.4,-85}}, color={0,0,127}));
  connect(scale1.y, heating.u) annotation (Line(points={{-95.3,-85},{-84.65,-85},
          {-84.65,-85},{-73.8,-85}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{120,
            120}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{
            120,120}})));
end ZoneR4C3_RBC;
