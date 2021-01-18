within OU44;
model VeWrapper

  Ventilation ventilation(
    hxEta=hxEta,
    qFan=qFan,
    coilEta=coilEta,
    hCoilCapacity=hCoilCapacity,
    basicFlow=basicFlow)
    annotation (Placement(transformation(extent={{46,-22},{86,18}})));
  Modelica.Blocks.Interfaces.RealInput verate1 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,86},{-108,126}})));
  Modelica.Blocks.Interfaces.RealInput verate2 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,56},{-108,96}})));
  Modelica.Blocks.Interfaces.RealInput verate3 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,20},{-108,60}})));
  Modelica.Blocks.Interfaces.RealInput verate4 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,-14},{-108,26}})));
  Modelica.Blocks.Interfaces.RealInput verate5 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,-48},{-108,-8}})));
  Modelica.Blocks.Interfaces.RealInput verate6 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,-80},{-108,-40}})));
  Modelica.Blocks.Interfaces.RealInput verate7 "fraction [0-1]"
    annotation (Placement(transformation(extent={{-148,-110},{-108,-70}})));
  Modelica.Blocks.Interfaces.RealInput Tindoor annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-20,-148})));
  Modelica.Blocks.Interfaces.RealInput Tvestp annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={12,-148})));
  Modelica.Blocks.Interfaces.RealInput Tout annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={46,-148})));
  Modelica.Blocks.Math.Gain gain1(k=360)
    annotation (Placement(transformation(extent={{-72,96},{-52,116}})));
  Modelica.Blocks.Math.Gain gain2(k=360)
    annotation (Placement(transformation(extent={{-72,66},{-52,86}})));
  Modelica.Blocks.Math.Gain gain3(k=360)
    annotation (Placement(transformation(extent={{-72,30},{-52,50}})));
  Modelica.Blocks.Math.Gain gain4(k=360)
    annotation (Placement(transformation(extent={{-74,-4},{-54,16}})));
  Modelica.Blocks.Math.Gain gain5(k=360)
    annotation (Placement(transformation(extent={{-74,-38},{-54,-18}})));
  Modelica.Blocks.Math.Gain gain6(k=360)
    annotation (Placement(transformation(extent={{-74,-70},{-54,-50}})));
  Modelica.Blocks.Math.Gain gain7(k=360)
    annotation (Placement(transformation(extent={{-74,-100},{-54,-80}})));
  Modelica.Blocks.Math.MultiSum multiSum(nu=7)
    annotation (Placement(transformation(extent={{-6,34},{6,46}})));
  Modelica.Blocks.Interfaces.RealOutput qheat "Heating coil power [W]"
    annotation (Placement(transformation(extent={{138,50},{158,70}})));
  Modelica.Blocks.Interfaces.RealOutput Qheat "Total heating coil energy [J]"
    annotation (Placement(transformation(extent={{136,10},{156,30}})));
  Modelica.Blocks.Interfaces.RealOutput qel "Fan electric power [W]"
    annotation (Placement(transformation(extent={{138,-18},{158,2}})));
  Modelica.Blocks.Interfaces.RealOutput Qel
    "Total fan electric energy consumption [J]"
    annotation (Placement(transformation(extent={{138,-40},{158,-20}})));
  Modelica.Blocks.Interfaces.RealOutput TAfterHX
    "Temperature after passing HX [degC]"
    annotation (Placement(transformation(extent={{136,90},{156,110}})));
  parameter Modelica.SIunits.Efficiency hxEta=0.7
    "Heat exchanger effectiveness [-]";
  parameter Real qFan=1000 "Constant fan electricity consumption rate [W]";
  parameter Real coilEta=0.6 "Heating coil efficiency";
  parameter Modelica.SIunits.HeatFlowRate hCoilCapacity=15000
    "Heating capacity of the heating coil";
  parameter Modelica.Blocks.Interfaces.RealOutput basicFlow=0.1
    "Additional flow in the main ventilation loop";
equation
  connect(Tindoor, ventilation.TIndoor) annotation (Line(points={{-20,-148},{
          -19.6,-148},{-19.6,-7.2},{45.2,-7.2}},
                                           color={0,0,127}));
  connect(Tvestp, ventilation.Tvestp) annotation (Line(points={{12,-148},{12,
          -148},{12,12.2},{45.4,12.2}},
                                  color={0,0,127}));
  connect(verate1, gain1.u) annotation (Line(points={{-128,106},{-128,106},{-74,
          106}}, color={0,0,127}));
  connect(verate2, gain2.u)
    annotation (Line(points={{-128,76},{-128,76},{-74,76}}, color={0,0,127}));
  connect(verate3, gain3.u)
    annotation (Line(points={{-128,40},{-128,40},{-74,40}}, color={0,0,127}));
  connect(verate4, gain4.u)
    annotation (Line(points={{-128,6},{-128,6},{-76,6}}, color={0,0,127}));
  connect(verate5, gain5.u) annotation (Line(points={{-128,-28},{-92,-28},{-76,-28}},
        color={0,0,127}));
  connect(verate6, gain6.u) annotation (Line(points={{-128,-60},{-116,-60},{-76,
          -60}}, color={0,0,127}));
  connect(verate7, gain7.u) annotation (Line(points={{-128,-90},{-128,-90},{-76,
          -90}}, color={0,0,127}));
  connect(Tout, ventilation.TAmb) annotation (Line(points={{46,-148},{46,-148},
          {46,-70},{36,-70},{36,-3.4},{45.2,-3.4}},color={0,0,127}));
  connect(ventilation.VFresh, multiSum.y) annotation (Line(points={{45.2,1},{
          26,1},{26,40},{7.02,40}},
                                 color={0,0,127}));
  connect(gain1.y, multiSum.u[1]) annotation (Line(points={{-51,106},{-30,106},{
          -30,43.6},{-6,43.6}}, color={0,0,127}));
  connect(gain2.y, multiSum.u[2]) annotation (Line(points={{-51,76},{-30,76},{-30,
          42.4},{-6,42.4}}, color={0,0,127}));
  connect(gain3.y, multiSum.u[3])
    annotation (Line(points={{-51,40},{-6,40},{-6,41.2}}, color={0,0,127}));
  connect(gain4.y, multiSum.u[4]) annotation (Line(points={{-53,6},{-30,6},{-30,
          40},{-6,40}}, color={0,0,127}));
  connect(gain5.y, multiSum.u[5]) annotation (Line(points={{-53,-28},{-30,-28},{
          -30,38.8},{-6,38.8}}, color={0,0,127}));
  connect(gain6.y, multiSum.u[6]) annotation (Line(points={{-53,-60},{-30,-60},{
          -30,37.6},{-6,37.6}}, color={0,0,127}));
  connect(gain7.y, multiSum.u[7]) annotation (Line(points={{-53,-90},{-30,-90},{
          -30,36.4},{-6,36.4}}, color={0,0,127}));
  connect(ventilation.qFlow, qheat) annotation (Line(points={{86.6,15.6},{
          114.3,15.6},{114.3,60},{148,60}}, color={0,0,127}));
  connect(ventilation.QFlow,Qheat)  annotation (Line(points={{86.6,11.2},{
          119.3,11.2},{119.3,20},{146,20}},
                                      color={0,0,127}));
  connect(ventilation.qEl, qel) annotation (Line(points={{86.4,-12},{114,-12},
          {114,-8},{148,-8}},
                         color={0,0,127}));
  connect(ventilation.QEl, Qel) annotation (Line(points={{86.4,-16.2},{112.2,
          -16.2},{112.2,-30},{148,-30}},
                                  color={0,0,127}));
  connect(ventilation.TAfterHX, TAfterHX) annotation (Line(points={{73.6,18.6},
          {73.6,100.3},{146,100.3},{146,100}},color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,-140},
            {140,140}})), Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-120,-140},{140,140}})));
end VeWrapper;
