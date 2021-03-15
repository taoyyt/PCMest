within OU44;
model test_ou44
  office_model.office_zone_MPC_est_para_KA_4MPC_controller
               office_zone_MPC_est_para_KA_4MPC_controller
    annotation (Placement(transformation(extent={{-14,-16},{30,28}})));
  Modelica.Blocks.Sources.Constant Tvestp(k=20)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Modelica.Blocks.Sources.Pulse solrad(
    amplitude=400,
    width=33,
    period=86400,
    startTime=28800)
    annotation (Placement(transformation(extent={{-80,58},{-60,78}})));
  Modelica.Blocks.Sources.Pulse Tout(
    width=33,
    period=86400,
    startTime=28800,
    offset=3,
    amplitude=7)
    annotation (Placement(transformation(extent={{-80,22},{-60,42}})));
  Modelica.Blocks.Sources.Pulse occ(
    width=33,
    period=86400,
    startTime=28800,
    amplitude=30,
    offset=0)
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Blocks.Sources.Constant Tvestp1(k=500)
    annotation (Placement(transformation(extent={{-80,-70},{-60,-50}})));
equation
  connect(Tout.y, office_zone_MPC_est_para_KA_4MPC_controller.Tout) annotation (
     Line(points={{-59,32},{-36,32},{-36,16.4},{-14.6286,16.4}}, color={0,0,127}));
  connect(solrad.y, office_zone_MPC_est_para_KA_4MPC_controller.Solrad)
    annotation (Line(points={{-59,68},{-36,68},{-36,22.8},{-14.8381,22.8}},
        color={0,0,127}));
  connect(occ.y, office_zone_MPC_est_para_KA_4MPC_controller.Occ) annotation (
      Line(points={{-59,0},{-36,0},{-36,11.8},{-14.6286,11.8}}, color={0,0,127}));
  connect(Tvestp.y, office_zone_MPC_est_para_KA_4MPC_controller.Tsup)
    annotation (Line(points={{-59,-30},{-36,-30},{-36,6.6},{-14.6286,6.6}},
        color={0,0,127}));
  connect(Tvestp1.y, office_zone_MPC_est_para_KA_4MPC_controller.Vair)
    annotation (Line(points={{-59,-60},{-38,-60},{-38,1.4},{-14.6286,1.4}},
        color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=604800));
end test_ou44;
