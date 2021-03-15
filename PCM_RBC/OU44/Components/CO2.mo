within OU44.Components;
model CO2 "CO2 balance model"

  Modelica.Blocks.Interfaces.RealInput Vve "[m3/h]"
    annotation (Placement(transformation(extent={{-120,44},{-88,76}}),
        iconTransformation(extent={{-120,44},{-88,76}})));
  Modelica.Blocks.Interfaces.RealOutput CO2
    annotation (Placement(transformation(extent={{96,-10},{116,10}})));
  CO2Eqn balance(Vi=Vi, CO2ppmv_initial=CO2PpmInitial)
    annotation (Placement(transformation(extent={{-22,-25},{26,25}})));
  Modelica.Blocks.Interfaces.RealInput persons
    annotation (Placement(transformation(extent={{-119,-75},{-89,-45}}),
        iconTransformation(extent={{-119,-75},{-89,-45}})));
  Modelica.Blocks.Sources.Constant CO2Generation(k=CO2PerPerson) "[m3/h]"
    annotation (Placement(transformation(extent={{-76,-38},{-56,-18}})));
  parameter Real CO2PerPerson=0.02 "CO2 generation per person [m3/h]";
  parameter Real Vi=100 "Indoor volume [m3]";
  parameter Real CO2PpmInitial=400 "Initial CO2 concentration [ppmv]";
  parameter Real CO2Neutral=400 "CO2 Neutral Level";
  Modelica.Blocks.Sources.Constant CO2NeutralLevel(k=CO2Neutral)
    annotation (Placement(transformation(extent={{-76,6},{-56,26}})));

equation
  connect(persons, balance.persons) annotation (Line(points={{-104,-60},{
          -84,-60},{-84,-6.5},{-22.96,-6.5}}, color={0,0,127}));
  connect(Vve, balance.Vve) annotation (Line(points={{-104,60},{-32,60},{
          -32,18.5},{-22.96,18.5}}, color={0,0,127}));
  connect(balance.CO2output, CO2) annotation (Line(points={{27.44,0},{
          27.44,0},{106,0}}, color={0,0,127}));
  connect(CO2Generation.y, balance.CO2_per_person) annotation (Line(
        points={{-55,-28},{-40,-28},{-40,-18.5},{-22.96,-18.5}}, color={0,
          0,127}));
  connect(CO2NeutralLevel.y, balance.CO2ppmv_s) annotation (Line(points={
          {-55,16},{-40,16},{-40,6.5},{-22.96,6.5}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={
        Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200},
          lineThickness=0.5),
        Line(
          points={{36,60}},
          color={28,108,200},
          thickness=0.5),
        Line(
          points={{-66,14}},
          color={135,135,135},
          thickness=0.5),
        Line(
          points={{82,10},{-18,-30}},
          color={28,108,200},
          thickness=0.5),
        Polygon(
          points={{32,-10},{24,-30},{40,-30},{32,-10}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillPattern=FillPattern.Sphere,
          fillColor={28,108,200}),
        Line(points={{-72,-54},{-52,-54}}, color={255,0,0}),
        Line(points={{-54,-84},{-62,-66},{-62,-46}}, color={255,0,0}),
        Ellipse(extent={{-68,-34},{-56,-46}}, lineColor={255,0,0}),
        Line(points={{-62,-66},{-70,-84}}, color={255,0,0}),
        Line(
          points={{-82,56},{-76,40},{-66,60},{-56,40},{-46,60},{-38,46}},
          color={0,128,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{-78,72},{-72,56},{-62,76},{-52,56},{-42,76},{-34,62}},
          color={0,128,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{22,52},{8,52},{6,50},{6,14},{8,12},{22,12}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{44,52},{28,52},{26,50},{26,14},{28,12},{44,12},{46,14},
              {46,50},{44,52}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{50,16},{52,18},{54,18},{56,16},{56,14},{50,8},{56,8}},
          color={0,0,127},
          thickness=0.5)}));
end CO2;
