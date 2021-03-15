within OU44.Components;
model CO2Eqn "CO2 balance model (equations)"

  Modelica.Blocks.Interfaces.RealInput Vve
    "Ventilation air flow rate [m3/h]"
    annotation (Placement(transformation(extent={{-116,62},{-92,86}})));
  Modelica.Blocks.Interfaces.RealInput CO2ppmv_s
    "CO2 concentration in incoming air [ppmv]"
    annotation (Placement(transformation(extent={{-116,14},{-92,38}})));
  Modelica.Blocks.Interfaces.RealInput persons "Number of persons"
    annotation (Placement(transformation(extent={{-116,-38},{-92,-14}})));

  Modelica.Blocks.Interfaces.RealInput CO2_per_person
    "CO2 generation per person [ppmv/h]"
    annotation (Placement(transformation(extent={{-116,-86},{-92,-62}})));
  Modelica.Blocks.Interfaces.RealOutput CO2output
    annotation (Placement(transformation(extent={{96,-10},{116,10}})));

  parameter Real Vi "Indoor volume [m3]";
  parameter Real CO2ppmv_initial=400 "Initial CO2 concentration [ppmv]";

  // Real VCO2_g_rate "CO2 human generation rate [m3/h]";
  // Real VCO2_s_rate "CO2 supply rate [m3/h]";
  // Real VCO2_ex_rate "C02 extraction rate [m3/h]";
  // Real VCO2_i "CO2 indoor concentration [m3]";
  Real CO2ppmv_i(fixed=true, start=CO2ppmv_initial)
    "CO2 indoor concentration [ppmv]";

equation
  // Short version of the balance (works with EstimationPy!):
  3600 * Vi * der(CO2ppmv_i) = persons * CO2_per_person * 10^6 + Vve * (CO2ppmv_s - CO2ppmv_i);
  CO2output = CO2ppmv_i;

  // Long version (more intuitive) of the balance:
  /*
  der(VCO2_i) = VCO2_g_rate + VCO2_s_rate - VCO2_ex_rate "Transient balance";
  VCO2_g_rate = persons * CO2_per_person / 3600 "Generation rate";
  VCO2_s_rate = Vve * CO2ppmv_s / 10^6 / 3600 "Supply rate";
  VCO2_ex_rate = Vve * CO2ppmv_i / 10^6 / 3600 "Extraction rate";
  VCO2_i = CO2ppmv_i * Vi / 10^6 "Relationship between PPMV and m3";
  CO2output = CO2ppmv_i "Output in m3";
  */

   annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={
        Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
        Line(
          points={{84,-4},{-16,-44}},
          color={28,108,200},
          thickness=0.5),
        Polygon(
          points={{34,-24},{26,-44},{42,-44},{34,-24}},
          lineColor={28,108,200},
          lineThickness=0.5,
          fillPattern=FillPattern.Sphere,
          fillColor={28,108,200}),
        Line(points={{-76,-18},{-56,-18}}, color={255,0,0}),
        Line(points={{-60,-48},{-66,-30},{-66,-10}}, color={255,0,0}),
        Ellipse(extent={{-72,2},{-60,-10}}, lineColor={255,0,0}),
        Line(points={{-66,-30},{-74,-48}}, color={255,0,0}),
        Ellipse(
          extent={{-76,-84},{-74,-86}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-72,-80},{-70,-82}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-64,-78},{-62,-80}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-68,-72},{-66,-74}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-60,-62},{-58,-64}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-74,-64},{-72,-66}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-76,-74},{-74,-76}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-88,24},{-82,8},{-72,28},{-62,8},{-52,28},{-44,14}},
          color={0,255,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{-84,40},{-78,24},{-68,44},{-58,24},{-48,44},{-40,30}},
          color={0,255,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Ellipse(
          extent={{-54,18},{-52,16}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-56,36},{-54,34}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-78,36},{-76,34}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-72,30},{-70,28}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-88,70},{-82,54},{-72,74},{-62,54},{-52,74},{-44,60}},
          color={0,255,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{-86,78},{-80,62},{-70,82},{-60,62},{-50,82},{-42,68}},
          color={0,255,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{-84,86},{-78,70},{-68,90},{-58,70},{-48,90},{-40,76}},
          color={0,255,255},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Ellipse(
          extent={{-64,26},{-62,24}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-74,16},{-72,14}},
          lineColor={135,135,135},
          lineThickness=0.5,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Line(
          points={{24,38},{10,38},{8,36},{8,0},{10,-2},{24,-2}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{46,38},{30,38},{28,36},{28,0},{30,-2},{46,-2},{48,0},{
              48,36},{46,38}},
          color={0,0,127},
          thickness=0.5),
        Line(
          points={{52,2},{54,4},{56,4},{58,2},{58,0},{52,-6},{58,-6}},
          color={0,0,127},
          thickness=0.5)}));
end CO2Eqn;
