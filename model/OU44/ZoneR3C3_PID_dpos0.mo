within OU44;
model ZoneR3C3_PID_dpos0 "single zone model"
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
  Modelica.Blocks.Interfaces.RealInput Tstp_heating
    annotation (Placement(transformation(extent={{-178,-122},{-160,-104}})));
  Modelica.Blocks.Interfaces.RealInput CO2stp
    annotation (Placement(transformation(extent={{-162,-66},{-142,-46}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
    annotation (Placement(transformation(extent={{-132,74},{-118,88}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-100,74},{-86,88}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Re(R=RExt/(3*Vi^(2/3)))
    "thermal resistence of external wall"
    annotation (Placement(transformation(extent={{-78,72},{-60,90}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_air(
                                  C=tmass*1.2*1005*Vi, T(
      displayUnit="degC",
      start=293.15,
      fixed=false))             "air heat capacity"
    annotation (Placement(transformation(extent={{6,88},{24,106}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rinternal(R=Rint/(3*
        Vi^(2/3)))
    "thermal resistance of internal wall/furniture"
    annotation (Placement(transformation(extent={{40,100},{52,112}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_internal_wall(C=imass*
        1.2*1005*Vi, T(fixed=false))
    annotation (Placement(transformation(extent={{60,106},{74,120}})));
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
  Modelica.Blocks.Math.Gain scale(k=0)    annotation (Placement(transformation(
        extent={{5,-5},{-5,5}},
        rotation=180,
        origin={-83,-55})));
  Modelica.Blocks.Math.Gain ventilation(k=maxVent)
    annotation (Placement(transformation(extent={{-72,-62},{-58,-48}})));
  Modelica.Blocks.Continuous.Integrator integrator(k=1)
    annotation (Placement(transformation(extent={{86,-62},{100,-48}})));
  Modelica.Blocks.Interfaces.RealOutput vetot "aggregated ventilation flow[m3]"
    annotation (Placement(transformation(extent={{120,-66},{142,-44}})));
  Modelica.Blocks.Interfaces.RealOutput verate
    "ventilation air supply rate[m3/h]"
    annotation (Placement(transformation(extent={{120,-42},{142,-20}})));
  Modelica.Blocks.Math.Gain heating(k=maxHeat)
    annotation (Placement(transformation(extent={{-64,-90},{-52,-78}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(extent={{-40,-92},{-26,-78}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-18,-92},{-4,-78}})));
  Modelica.Blocks.Interfaces.RealOutput Qrad
    "aggregated energy supplied by radiators[Wh]"
    annotation (Placement(transformation(extent={{126,-112},{148,-90}})));
  Modelica.Blocks.Interfaces.RealOutput qrad "heating power of radiators[W]"
    annotation (Placement(transformation(extent={{124,-130},{146,-108}})));
  Modelica.Blocks.Continuous.Integrator integrator1(k=1)
    annotation (Placement(transformation(extent={{84,-116},{98,-102}})));
  OU44.Components.CO2 co2(
  Vi=Vi,
  CO2PpmInitial=CO2Init,
  CO2PerPerson=CO2pp,
  CO2Neutral=CO2n)
  annotation (Placement(transformation(extent={{2,-2},{14,10}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-16,10},{-8,18}})));
  Modelica.Blocks.Interfaces.RealOutput CO2 "[ppm]"
    annotation (Placement(transformation(extent={{124,6},{146,28}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin    toKelvin1
    annotation (Placement(transformation(extent={{-126,-34},{-114,-22}})));
  Modelica.Blocks.Math.Max max
    annotation (Placement(transformation(extent={{-104,-30},{-94,-20}})));
  OU44.Components.AirMix
                    airMix
    annotation (Placement(transformation(extent={{-6,-28},{8,-14}})));

  parameter Real shgc=0.5 "solar heat gains coefficient";
  parameter Real Vinf=100 "inflitration rate";
  parameter Real TairInit=20 "Initial temperature of indoor air [deg C]";
  parameter Real CO2n=400 "CO2 neutral level";
  parameter Real CO2pp=0.02 "CO2 generation per person [m3/h]";
  parameter Real maxVent=2000 "Maximum ventilation flowrate[m3/h]";
  parameter Real maxHeat=5000 "Heating power of radiator [W]";
  parameter Real Tve=21 "Ventilation air temperature [deg C]";
  parameter Real tmass=5;
  parameter Real imass=10;
  parameter Real wmass=5;
  parameter Real Vi=300 "Air volume [m3]";
  parameter Real occheff=1. "Occupant heat generation effectiveness";
  parameter Real CO2Init=400 "Initial CO2 concentration[ppmv]";
  parameter Real RExt=1.0 "External wall thermal resistance";

  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rindoor(R=Ri/(3*Vi^(2
        /3)))
    "Indoor thermal convective/conductive resistance"
    annotation (Placement(transformation(extent={{-28,74},{-14,88}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HC_wall(C=wmass*1.2*1005
        *Vi)
    "heat capacity of external wall" annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={-50,60})));
  parameter Real Ri=1.0 "Indoor thermal resistance";
  parameter Real Rint=1.0 "thermal resistance of indoor funiture";
  Modelica.Blocks.Math.Gain scale2(k=0.01)
                                          annotation (Placement(transformation(
        extent={{5,-5},{-5,5}},
        rotation=180,
        origin={-89,-83})));
  Modelica.Blocks.Continuous.LimPID PID_temp(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=100,
    Ti=0.5,
    initType=Modelica.Blocks.Types.InitPID.DoNotUse_InitialIntegratorState,
    yMax=100,
    yMin=0)
    annotation (Placement(transformation(extent={{-128,-116},{-120,-108}})));
  Modelica.Blocks.Interfaces.RealInput Tstp_cooling
    annotation (Placement(transformation(extent={{-178,-90},{-160,-72}})));
  Modelica.Blocks.Continuous.LimPID PID_temp1(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=100,
    Ti=0.5,
    yMax=0,
    yMin=-100,
    initType=Modelica.Blocks.Types.InitPID.DoNotUse_InitialIntegratorState)
    annotation (Placement(transformation(extent={{-128,-84},{-120,-76}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
    annotation (Placement(transformation(extent={{-150,-88},{-136,-74}})));
  Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin3
    annotation (Placement(transformation(extent={{-152,-120},{-138,-106}})));
  Modelica.Blocks.Math.Add add1
    annotation (Placement(transformation(extent={{-112,-88},{-102,-78}})));
  Modelica.Blocks.Continuous.LimPID PID_CO2(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    Ti=0.5,
    yMax=0,
    yMin=-100,
    initType=Modelica.Blocks.Types.InitPID.DoNotUse_InitialIntegratorState,
    k=10) annotation (Placement(transformation(extent={{-126,-60},{-118,-52}})));
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
    annotation (Line(points={{-86,81},{-78,81}}, color={191,0,0}));
  connect(HC_air.port, Rinternal.port_a) annotation (Line(points={{15,88},{26,
        88},{26,106},{40,106}},
                              color={191,0,0}));
  connect(Rinternal.port_b, HC_internal_wall.port)
    annotation (Line(points={{52,106},{67,106}}, color={191,0,0}));
  connect(HC_air.port, temperatureSensor.port)
    annotation (Line(points={{15,88},{34,88}}, color={191,0,0}));
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
    annotation (Line(points={{2,26},{15,26},{15,88}}, color={191,0,0}));
  connect(scale.y, ventilation.u)
    annotation (Line(points={{-77.5,-55},{-73.4,-55}}, color={0,0,127}));
  connect(ventilation.y, integrator.u)
    annotation (Line(points={{-57.3,-55},{84.6,-55}}, color={0,0,127}));
  connect(integrator.y, vetot)
    annotation (Line(points={{100.7,-55},{131,-55}}, color={0,0,127}));
  connect(ventilation.y, verate) annotation (Line(points={{-57.3,-55},{74,-55},
          {74,-31},{131,-31}}, color={0,0,127}));
  connect(heating.y, prescribedHeatFlow1.Q_flow) annotation (Line(points={{-51.4,
          -84},{-46,-84},{-46,-85},{-40,-85}},       color={0,0,127}));
  connect(prescribedHeatFlow1.port, heatFlowSensor.port_a)
    annotation (Line(points={{-26,-85},{-18,-85}}, color={191,0,0}));
  connect(heatFlowSensor.port_b,HC_air. port) annotation (Line(points={{-4,-85},
        {20,-85},{20,88},{15,88}},   color={191,0,0}));
  connect(integrator1.y, Qrad) annotation (Line(points={{98.7,-109},{110.35,
          -109},{110.35,-101},{137,-101}},
                                   color={0,0,127}));
  connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-11,-92},{
        -12,-92},{-12,-119},{135,-119}},
                                    color={0,0,127}));
  connect(ventilation.y, add.u2) annotation (Line(points={{-57.3,-55},{-34,-55},
          {-34,11.6},{-16.8,11.6}}, color={0,0,127}));
  connect(const1.y, add.u1) annotation (Line(points={{-18.5,63},{-22,63},{-22,
          16.4},{-16.8,16.4}}, color={0,0,127}));
connect(add.y, co2.Vve) annotation (Line(points={{-7.6,14},{-2,14},{-2,7.6},{
        1.76,7.6}}, color={0,0,127}));
connect(occ, co2.persons) annotation (Line(points={{-150,38},{-50,38},{-50,
        0.4},{1.76,0.4}}, color={0,0,127}));
connect(co2.CO2, CO2) annotation (Line(points={{14.36,4},{68,4},{68,17},{135,
        17}}, color={0,0,127}));
  connect(toKelvin1.Kelvin, max.u2)
    annotation (Line(points={{-113.4,-28},{-105,-28}}, color={0,0,127}));
  connect(toKelvin.Kelvin, max.u1) annotation (Line(points={{-117.3,81},{-110,
          81},{-110,-22},{-105,-22}}, color={0,0,127}));
  connect(max.y, airMix.Tve) annotation (Line(points={{-93.5,-25},{-17.75,-25},
          {-17.75,-18.62},{-6.56,-18.62}}, color={0,0,127}));
  connect(ventilation.y, airMix.Vve) annotation (Line(points={{-57.3,-55},{-14,
          -55},{-14,-24.5},{-6.7,-24.5}}, color={0,0,127}));
  connect(HC_air.port, airMix.port_b) annotation (Line(points={{15,88},{20,88},
        {20,-21},{8,-21}},   color={191,0,0}));
  connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{46,88},{54,88},
          {54,-64},{1,-64},{1,-28.56}},     color={0,0,127}));
  connect(Re.port_b, Rindoor.port_a)
    annotation (Line(points={{-60,81},{-28,81}}, color={191,0,0}));
  connect(Rindoor.port_b,HC_air. port) annotation (Line(points={{-14,81},{0,
        81},{0,88},{15,88}},
                           color={191,0,0}));
  connect(Re.port_b,HC_wall. port)
    annotation (Line(points={{-60,81},{-50,81},{-50,67}}, color={191,0,0}));
  connect(prescribedHeatFlow.port, Rindoor.port_a) annotation (Line(points={{-80,108},
          {-42,108},{-42,81},{-28,81}},          color={191,0,0}));
  connect(fromKelvin.Celsius, T)
    annotation (Line(points={{80.6,88},{130,88}}, color={0,0,127}));
  connect(const2.y, toKelvin1.Celsius)
    annotation (Line(points={{-131.4,-28},{-127.2,-28}}, color={0,0,127}));
connect(scale2.y, heating.u) annotation (Line(points={{-83.5,-83},{-78,-83},{
          -78,-84},{-65.2,-84}},   color={0,0,127}));
connect(heatFlowSensor.Q_flow, integrator1.u) annotation (Line(points={{-11,-92},
          {34,-92},{34,-109},{82.6,-109}},    color={0,0,127}));
  connect(Tstp_cooling, toKelvin2.Celsius) annotation (Line(points={{-169,-81},
          {-160.5,-81},{-160.5,-81},{-151.4,-81}}, color={0,0,127}));
  connect(Tstp_heating, toKelvin3.Celsius) annotation (Line(points={{-169,-113},
          {-160.5,-113},{-160.5,-113},{-153.4,-113}}, color={0,0,127}));
  connect(toKelvin2.Kelvin, PID_temp1.u_s) annotation (Line(points={{-135.3,-81},
          {-132.65,-81},{-132.65,-80},{-128.8,-80}}, color={0,0,127}));
  connect(temperatureSensor.T, PID_temp1.u_m) annotation (Line(points={{46,88},
          {54,88},{54,-102},{-124,-102},{-124,-84.8}}, color={0,0,127}));
  connect(temperatureSensor.T, PID_temp.u_m) annotation (Line(points={{46,88},{
          56,88},{56,-128},{-124,-128},{-124,-116.8}}, color={0,0,127}));
  connect(toKelvin3.Kelvin, PID_temp.u_s) annotation (Line(points={{-137.3,-113},
          {-132.65,-113},{-132.65,-112},{-128.8,-112}}, color={0,0,127}));
  connect(PID_temp1.y, add1.u1)
    annotation (Line(points={{-119.6,-80},{-113,-80}}, color={0,0,127}));
  connect(PID_temp.y, add1.u2) annotation (Line(points={{-119.6,-112},{-116,
          -112},{-116,-86},{-113,-86}}, color={0,0,127}));
  connect(add1.y, scale2.u) annotation (Line(points={{-101.5,-83},{-98.75,-83},
          {-98.75,-83},{-95,-83}}, color={0,0,127}));
  connect(CO2stp, PID_CO2.u_s)
    annotation (Line(points={{-152,-56},{-126.8,-56}}, color={0,0,127}));
  connect(co2.CO2, PID_CO2.u_m) annotation (Line(points={{14.36,4},{32,4},{32,
          -72},{-122,-72},{-122,-60.8}}, color={0,0,127}));
  connect(PID_CO2.y, scale.u) annotation (Line(points={{-117.6,-56},{-104,-56},
          {-104,-55},{-89,-55}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{120,
            120}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{
            120,120}})));
end ZoneR3C3_PID_dpos0;
