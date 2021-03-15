within ;
package Cpcm_mean
  model MPC_demo_one_stack_officezone_emulator
    extends Buildings.BaseClasses.BaseIcon;
    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"});

    parameter Real TairInit=20;

    Modelica.Blocks.Math.Gain gain1(k=50)
      annotation (Placement(transformation(extent={{-20,-34},{-4,-18}})));
    Modelica.Blocks.Math.Gain gain2(k=3600)
      annotation (Placement(transformation(extent={{48,-34},{64,-18}})));
    OU44.office_model.office_zone_MPC_est_para_KA
               mPC_demo4_1
      annotation (Placement(transformation(extent={{98,-10},{130,22}})));
    Modelica.Blocks.Interfaces.RealInput solrad "solar radiation [w/m2]"
      annotation (Placement(transformation(extent={{-136,40},{-120,56}})));
    Modelica.Blocks.Interfaces.RealInput Vair "volume flow rate [m3/h]"
      annotation (Placement(transformation(extent={{-136,8},{-120,24}})));
    Modelica.Blocks.Interfaces.RealInput Tout
      "air temperature [celsius degree]"
      annotation (Placement(transformation(extent={{-136,-14},{-120,2}})));
    Modelica.Blocks.Interfaces.RealOutput Tin
      "indoor temperature [celsius degree]"
      annotation (Placement(transformation(extent={{160,2},{180,22}})));
    PCM.PCM_MPC.AHU_PCM_one_stack_copy aHU_PCM_one_stack
      annotation (Placement(transformation(extent={{-70,-6},{-42,22}})));
    Modelica.Blocks.Interfaces.RealOutput Tsup
      "supply ventilation temperature"
      annotation (Placement(transformation(extent={{160,24},{180,44}})));
    Modelica.Blocks.Interfaces.RealOutput Tpcm "PCM temperature"
      annotation (Placement(transformation(extent={{160,78},{180,98}})));
    Modelica.Blocks.Interfaces.RealOutput Qc "cooling capacity"
      annotation (Placement(transformation(extent={{160,-32},{180,-12}})));
    Modelica.Blocks.Interfaces.RealInput occ "occupancy"
      annotation (Placement(transformation(extent={{-136,26},{-120,42}})));
    Modelica.Blocks.Interfaces.RealOutput qc "cooling power"
      annotation (Placement(transformation(extent={{160,-50},{180,-30}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{106,82},{118,94}})));
  equation
    connect(gain1.y, gain2.u)
      annotation (Line(points={{-3.2,-26},{46.4,-26}}, color={0,0,127}));
    connect(gain2.y, mPC_demo4_1.Vair) annotation (Line(points={{64.8,-26},{84,
            -26},{84,2.65455},{97.5429,2.65455}},
                                           color={0,0,127}));
    connect(mPC_demo4_1.Solrad,solrad)  annotation (Line(points={{97.3905,
            18.2182},{80,18.2182},{80,48},{-128,48}},
                                          color={0,0,127}));
    connect(mPC_demo4_1.Tin, Tin) annotation (Line(points={{131.829,14.7273},{
            162,14.7273},{162,12},{170,12}},
                              color={0,0,127}));
    connect(Vair, aHU_PCM_one_stack.VF) annotation (Line(points={{-128,16},{
            -100,16},{-100,10.66},{-72.275,10.66}}, color={0,0,127}));
    connect(Tout, aHU_PCM_one_stack.Te) annotation (Line(points={{-128,-6},{
            -102,-6},{-102,3.66},{-71.925,3.66}}, color={0,0,127}));
    connect(aHU_PCM_one_stack.Tsup, mPC_demo4_1.Tsup) annotation (Line(points={
            {-40.25,9.68},{58,9.68},{58,6.43636},{97.5429,6.43636}}, color={0,0,
            127}));
    connect(aHU_PCM_one_stack.Vair_out, gain1.u) annotation (Line(points={{
            -40.25,3.24},{-30,3.24},{-30,-26},{-21.6,-26}}, color={0,0,127}));
    connect(aHU_PCM_one_stack.Tsup, Tsup) annotation (Line(points={{-40.25,9.68},
            {48,9.68},{48,34},{170,34}}, color={0,0,127}));
    connect(Tout, mPC_demo4_1.Tout) annotation (Line(points={{-128,-6},{-102,-6},
            {-102,-48},{90,-48},{90,13.5636},{97.5429,13.5636}},     color={0,
            0,127}));
    connect(occ, mPC_demo4_1.Occ) annotation (Line(points={{-128,34},{-10,34},{
            -10,40},{94,40},{94,10.2182},{97.5429,10.2182}},  color={0,0,127}));
    connect(mPC_demo4_1.Qrad, Qc) annotation (Line(points={{131.829,-1.27273},{
            156,-1.27273},{156,-22},{170,-22}},  color={0,0,127}));
    connect(mPC_demo4_1.qrad, qc) annotation (Line(points={{131.981,-4.03636},{
            144,-4.03636},{144,-40},{170,-40}},  color={0,0,127}));
    connect(aHU_PCM_one_stack.Tpcm, fromKelvin.Kelvin) annotation (Line(points=
            {{-40.25,15.28},{34,15.28},{34,88},{104.8,88}}, color={0,0,127}));
    connect(fromKelvin.Celsius, Tpcm)
      annotation (Line(points={{118.6,88},{170,88}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,
              -140},{160,160}}), graphics={Rectangle(extent={{-34,122},{416,
                -132}}, lineColor={28,108,200}),
                                        Text(
            extent={{34,-70},{334,-30}},
            lineColor={0,0,255},
            fontSize=20,
            textStyle={TextStyle.Bold},
            textString="ZONE 

")}),                       Diagram(coordinateSystem(preserveAspectRatio=false,
            extent={{-120,-140},{160,160}})),
      experiment(StopTime=432000));
  end MPC_demo_one_stack_officezone_emulator;

  model MPC_demo_one_stack_officezone_controller
    extends Buildings.BaseClasses.BaseIcon;
    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"});

    parameter Real TairInit=20;

    Modelica.Blocks.Math.Gain gain1(k=50)
      annotation (Placement(transformation(extent={{-20,-34},{-4,-18}})));
    Modelica.Blocks.Math.Gain gain2(k=3600)
      annotation (Placement(transformation(extent={{48,-34},{64,-18}})));
    OU44.office_model.office_zone_MPC_est_para_KA
               mPC_demo4_1
      annotation (Placement(transformation(extent={{98,-10},{130,22}})));
    Modelica.Blocks.Interfaces.RealInput solrad "solar radiation [w/m2]"
      annotation (Placement(transformation(extent={{-136,40},{-120,56}})));
    Modelica.Blocks.Interfaces.RealInput Vair "volume flow rate [m3/h]"
      annotation (Placement(transformation(extent={{-136,8},{-120,24}})));
    Modelica.Blocks.Interfaces.RealInput Tout
      "air temperature [celsius degree]"
      annotation (Placement(transformation(extent={{-136,-14},{-120,2}})));
    Modelica.Blocks.Interfaces.RealOutput Tin
      "indoor temperature [celsius degree]"
      annotation (Placement(transformation(extent={{160,2},{180,22}})));
    PCM.PCM_MPC.AHU_PCM_one_stack_cp_mean aHU_PCM_one_stack
      annotation (Placement(transformation(extent={{-70,-10},{-42,18}})));
    Modelica.Blocks.Interfaces.RealOutput Tsup
      "supply ventilation temperature"
      annotation (Placement(transformation(extent={{160,24},{180,44}})));
    Modelica.Blocks.Interfaces.RealOutput Tpcm "pcm temperature"
      annotation (Placement(transformation(extent={{160,78},{180,98}})));
    Modelica.Blocks.Interfaces.RealOutput Qc "cooling capacity"
      annotation (Placement(transformation(extent={{160,-32},{180,-12}})));
    Modelica.Blocks.Interfaces.RealInput occ "occupancy"
      annotation (Placement(transformation(extent={{-136,26},{-120,42}})));
    Modelica.Blocks.Interfaces.RealOutput qc "cooling power"
      annotation (Placement(transformation(extent={{160,-50},{180,-30}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{116,82},{128,94}})));
  equation
    connect(gain1.y, gain2.u)
      annotation (Line(points={{-3.2,-26},{46.4,-26}}, color={0,0,127}));
    connect(gain2.y, mPC_demo4_1.Vair) annotation (Line(points={{64.8,-26},{84,
            -26},{84,2.65455},{97.5429,2.65455}},
                                           color={0,0,127}));
    connect(mPC_demo4_1.Solrad,solrad)  annotation (Line(points={{97.3905,
            18.2182},{80,18.2182},{80,48},{-128,48}},
                                          color={0,0,127}));
    connect(mPC_demo4_1.Tin, Tin) annotation (Line(points={{131.829,14.7273},{
            162,14.7273},{162,12},{170,12}},
                              color={0,0,127}));
    connect(Vair, aHU_PCM_one_stack.VF) annotation (Line(points={{-128,16},{
            -100,16},{-100,6.66},{-72.275,6.66}}, color={0,0,127}));
    connect(Tout, aHU_PCM_one_stack.Te) annotation (Line(points={{-128,-6},{
            -102,-6},{-102,-0.34},{-71.925,-0.34}}, color={0,0,127}));
    connect(aHU_PCM_one_stack.Tsup, mPC_demo4_1.Tsup) annotation (Line(points={
            {-40.25,5.68},{58,5.68},{58,6.43636},{97.5429,6.43636}}, color={0,0,
            127}));
    connect(aHU_PCM_one_stack.Vair_out, gain1.u) annotation (Line(points={{
            -40.25,-0.76},{-30,-0.76},{-30,-26},{-21.6,-26}}, color={0,0,127}));
    connect(aHU_PCM_one_stack.Tsup, Tsup) annotation (Line(points={{-40.25,5.68},
            {48,5.68},{48,34},{170,34}}, color={0,0,127}));
    connect(Tout, mPC_demo4_1.Tout) annotation (Line(points={{-128,-6},{-102,-6},
            {-102,-48},{90,-48},{90,13.5636},{97.5429,13.5636}},     color={0,
            0,127}));
    connect(occ, mPC_demo4_1.Occ) annotation (Line(points={{-128,34},{-10,34},{
            -10,40},{94,40},{94,10.2182},{97.5429,10.2182}},  color={0,0,127}));
    connect(mPC_demo4_1.Qrad, Qc) annotation (Line(points={{131.829,-1.27273},{
            156,-1.27273},{156,-22},{170,-22}},  color={0,0,127}));
    connect(mPC_demo4_1.qrad, qc) annotation (Line(points={{131.981,-4.03636},{
            144,-4.03636},{144,-40},{170,-40}},  color={0,0,127}));
    connect(aHU_PCM_one_stack.Tpcm, fromKelvin.Kelvin) annotation (Line(points=
            {{-40.25,11.28},{36,11.28},{36,88},{114.8,88}}, color={0,0,127}));
    connect(fromKelvin.Celsius, Tpcm)
      annotation (Line(points={{128.6,88},{170,88}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,
              -140},{160,160}}), graphics={Rectangle(extent={{-34,122},{416,
                -132}}, lineColor={28,108,200}),
                                        Text(
            extent={{34,-70},{334,-30}},
            lineColor={0,0,255},
            fontSize=20,
            textStyle={TextStyle.Bold},
            textString="ZONE 

")}),                       Diagram(coordinateSystem(preserveAspectRatio=false,
            extent={{-120,-140},{160,160}})),
      experiment(StopTime=432000));
  end MPC_demo_one_stack_officezone_controller;
  annotation (uses(
      Buildings(version="6.0.0"),
      Modelica(version="3.2.2"),
      OU44(version="2")));
end Cpcm_mean;
