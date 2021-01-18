within ;
package PCM
    extends Buildings.BaseClasses.BaseIcon;

  model HeatCapacitor_PCM "Lumped thermal element storing heat"
    parameter Modelica.SIunits.SpecificHeatCapacity cp_s=2000
      "Specific heat capacity of solid element cp";
    parameter Modelica.SIunits.SpecificHeatCapacity cp_l=2000
    "Specific heat capacity of liqiud element ";
     parameter Modelica.SIunits.SpecificEnthalpy hf_pcm=140000
    "Specific heat capacity of liqiud element ";
    Modelica.SIunits.SpecificHeatCapacity cp_pcm
    "Specific heat capacity of PCM ";

    //parameter Modelica.Siunits.Mass m_pcm=2
    parameter Real m_pcm=2 "Mass(kg)";
    Modelica.SIunits.Temperature T(start=293.15, displayUnit="degC")
      "Temperature of element";
    Modelica.SIunits.TemperatureSlope der_T(start=0)
      "Time derivative of temperature (= der(T))";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation (
        Placement(transformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Blocks.Interfaces.RealInput dHdT
      annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
    Modelica.Blocks.Interfaces.RealInput H "liquid fraction"
      annotation (Placement(transformation(extent={{-140,12},{-100,52}})));

  equation
    T = port.T;
    der_T = der(T);
    cp_pcm = H*(cp_l+dHdT*hf_pcm)+(1-H)*(cp_s+dHdT*hf_pcm);
    m_pcm*cp_pcm*der(T)= port.Q_flow;
    annotation (Diagram(graphics={
          Polygon(
            points={{0,71},{-20,67},{-40,61},{-52,47},{-58,39},{-68,29},{-72,17},
                {-76,3},{-78,-11},{-76,-27},{-76,-39},{-76,-49},{-70,-61},{-64,
                -69},{-48,-73},{-30,-79},{-18,-79},{-2,-81},{8,-85},{22,-85},{
                32,-83},{42,-77},{54,-71},{56,-69},{66,-57},{68,-49},{70,-47},{
                72,-31},{76,-17},{78,-9},{78,7},{74,19},{66,29},{54,37},{44,45},
                {36,61},{26,69},{0,71}},
            lineColor={160,160,164},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-58,39},{-68,29},{-72,17},{-76,3},{-78,-11},{-76,-27},{-76,
                -39},{-76,-49},{-70,-61},{-64,-69},{-48,-73},{-30,-79},{-18,-79},
                {-2,-81},{8,-85},{22,-85},{32,-83},{42,-77},{54,-71},{42,-73},{
                40,-73},{30,-75},{20,-77},{18,-77},{10,-77},{2,-73},{-12,-69},{
                -22,-69},{-30,-67},{-40,-61},{-50,-51},{-56,-39},{-58,-31},{-58,
                -21},{-60,-9},{-60,-1},{-60,11},{-58,21},{-56,23},{-52,31},{-48,
                39},{-44,49},{-40,61},{-58,39}},
            lineColor={0,0,0},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid)}), Icon(graphics={
          Polygon(
            points={{2,79},{-18,75},{-38,69},{-50,55},{-56,47},{-66,37},{-70,25},
                {-74,11},{-76,-3},{-74,-19},{-74,-31},{-74,-41},{-68,-53},{-62,
                -61},{-46,-65},{-28,-71},{-16,-71},{0,-73},{10,-77},{24,-77},{
                34,-75},{44,-69},{56,-63},{58,-61},{68,-49},{70,-41},{72,-39},{
                74,-23},{78,-9},{80,-1},{80,15},{76,27},{68,37},{56,45},{46,53},
                {38,69},{28,77},{2,79}},
            lineColor={160,160,164},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-56,47},{-66,37},{-70,25},{-74,11},{-76,-3},{-74,-19},{-74,
                -31},{-74,-41},{-68,-53},{-62,-61},{-46,-65},{-28,-71},{-16,-71},
                {0,-73},{10,-77},{24,-77},{34,-75},{44,-69},{56,-63},{44,-65},{
                42,-65},{32,-67},{22,-69},{20,-69},{12,-69},{4,-65},{-10,-61},{
                -20,-61},{-28,-59},{-38,-53},{-48,-43},{-54,-31},{-56,-23},{-56,
                -13},{-58,-1},{-58,7},{-58,19},{-56,29},{-54,31},{-50,39},{-46,
                47},{-42,57},{-38,69},{-56,47}},
            lineColor={0,0,0},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid)}));
  end HeatCapacitor_PCM;

  model convection_coefficient

    // Parameters
    parameter Modelica.SIunits.Mass m_PCM=100;

    // Air properties
    parameter Modelica.SIunits.Conductivity k_a=0.02364;
    parameter Modelica.SIunits.Density rho_a=1.288;
    parameter Modelica.SIunits.DynamicViscosity my_a=0.00001729;
    parameter Modelica.SIunits.PrandtlNumber Pr_a=0.7344;

    // Geometry, flow etc. values for heat transfer calculation.
    parameter Modelica.SIunits.Distance L_CSM=0.45; // Length of CSM plate
    parameter Modelica.SIunits.Distance H_c=0.005;  // Height of air channel between CSM
    parameter Modelica.SIunits.Distance W_CSM=0.15; // Width of CSM plate
    //parameter Modelica.SIunits.Distance W_CSM=0.3; // Width of CSM plate
    parameter Modelica.SIunits.NusseltNumber Nu_rest=7.54;

    Modelica.SIunits.Distance L_c; // Characteristic length
    Modelica.SIunits.NusseltNumber Nu_entry; // Nusselt number for entrance region

    Modelica.SIunits.ReynoldsNumber Re;
    Modelica.SIunits.Velocity v_a; // Air velocity
    Real N_channels; // Number of air channels.
    Modelica.SIunits.Area A_PCM; // Heat transfer surface area in module.
    Modelica.SIunits.Distance l; // entrance length

    Modelica.Blocks.Interfaces.RealInput Vh "air volume flow rate (m3/s) "
      annotation (Placement(transformation(extent={{-130,-10},{-100,20}})));
    Modelica.Blocks.Interfaces.RealOutput h "convection coefficient"
      annotation (Placement(transformation(extent={{100,2},{120,22}})));
    Modelica.Blocks.Interfaces.RealOutput Nu "Nusselt number"
      annotation (Placement(transformation(extent={{100,-26},{120,-6}})));


  equation
     // Geometry:
     N_channels=m_PCM/2-1;
     A_PCM=(m_PCM/2-1)*(L_CSM*W_CSM)*2;
     Vh=v_a*(N_channels*H_c*L_CSM);
     //Vh=v_a*(m_PCM/4*H_c*L_CSM);
     L_c=2*H_c;

     // Determining Re, entrance length and Nusselt number:
     Re=(rho_a*v_a*L_c)/my_a;
     l=0.05*Re*Pr_a*L_c;
     Nu_entry=7.54+(0.03*(L_c/W_CSM)*Re*Pr_a)/(1+0.016*(0.01+(L_c/W_CSM)*Re*Pr_a).^(2/3)); // some minor change may need in this equation

     // Calculation of the Nusselt number based on the entrance length.
     if l>W_CSM then
       Nu=Nu_entry; // The air flow is not fully developed in the air channel.
     else
       Nu=Nu_entry*l/W_CSM+Nu_rest*(W_CSM-l)/W_CSM; // Weighted average of the fully developed flow and the entrance length flow.
     end if;

     Nu=h*L_c/k_a;

  //algorithm
    //when Vh<=0 then
      //h:=0;
    //end when;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-60,42},{-38,-66}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Forward),
          Line(points={{-32,28},{50,28}}, color={28,108,200}),
          Line(points={{46,32},{50,28},{46,24}}, color={28,108,200}),
          Line(points={{-32,2},{50,2}}, color={28,108,200}),
          Line(points={{46,6},{50,2},{46,-2}}, color={28,108,200}),
          Line(points={{-32,-24},{50,-24}}, color={28,108,200}),
          Line(points={{46,-20},{50,-24},{46,-28}}, color={28,108,200}),
          Line(points={{-32,-50},{50,-50}}, color={28,108,200}),
          Line(points={{46,-46},{50,-50},{46,-54}}, color={28,108,200})}),
                                                                   Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end convection_coefficient;

  model LiquidFraction_modify_2_2
      // SP21EK PCM hysteresis curve fitting properties:
    parameter Real Bh=2.118;
    parameter Real Bc=18.37;
    parameter Real Tmc=21.56+273.15;
    parameter Real Tmh=22.32+273.15;
    parameter Real vh=2.804;
    parameter Real vc=31.11;

    //Real Ah(stateSelect= StateSelect.always);
    //Real Kc(stateSelect= StateSelect.always);
    Real Ah;
    Real Kc;
    Real Kh=1.0;
    Real Ac=0.0;
    //Real Ah_new;
    //Real Kc_new;

    Boolean heating;
    Real mode;

    Modelica.Blocks.Interfaces.RealInput Tpcm annotation (Placement(
          transformation(extent={{-140,0},{-100,40}}), iconTransformation(extent={
              {-130,6},{-100,36}})));
    Modelica.Blocks.Interfaces.RealInput Tair annotation (Placement(
          transformation(extent={{-140,-44},{-100,-4}}), iconTransformation(
            extent={{-130,-34},{-100,-4}})));
    Modelica.Blocks.Interfaces.RealOutput H annotation (Placement(transformation(
            extent={{100,0},{142,42}}), iconTransformation(extent={{100,14},{130,44}})));
    Modelica.Blocks.Interfaces.RealOutput dHdT annotation (Placement(
          transformation(extent={{100,-46},{142,-4}}), iconTransformation(extent={
              {100,-36},{130,-6}})));

  initial equation
    heating=true;
    Kc=1;
    Ah=0.;
    //dHdT=0.05;


  equation
      if pre(heating) then //"previously heating mode"
        if Tair>=Tpcm then
          H=Ah+(Kh-Ah)/(1+exp(-Bh*(Tpcm-Tmh)))^(1/vh);
          dHdT=-(Bh*exp(Bh*(Tmh-Tpcm))*(Ah-Kh))/(vh*(exp(Bh*(Tmh-Tpcm))+1)^(1/vh+1));

        //Kc_heatingTocooling=Kc;
        //Ah_coolingToheating=Ah;
          mode=1.0;
        //"continue heating mode"
        else
          H=Ah+(Kh-Ah)/(1+exp(-Bh*(Tpcm-Tmh)))^(1/vh);
        //H= Ac+(Kc-Ac)/(1+exp(-Bc*(Tpcm-Tmc)))^(1/vc);
        //if time<=27000. then
          //dHdT=0.05;
          //else
          dHdT=-(Bc*exp(Bc*(Tmc - Tpcm))*(Ac - Kc))/(vc*(exp(Bc*(Tmc - Tpcm)) + 1)^(1/vc + 1));
        //end if;
        //dHdT=(Kc-Ac)*vc*((exp(Bc*(Tmc-Tpcm))+1)^(vc-1))*(exp(Bc*(Tmc-Tpcm)))*(-Bc);
        //Kc_heatingTocooling=(H-Ah)*((1+exp(Bh*Tmh-Bh*Tpcm))^(1/vh))+Ah;
        //Ah_coolingToheating=Ah;
          mode=2.0; //"switch heating to cooling mode"
        end if;
      else   //"previously cooling mode"
        if Tair>=Tpcm then
          H= Ac+(Kc-Ac)/(1+exp(-Bc*(Tpcm-Tmc)))^(1/vc);
        //H=Ah+(Kh-Ah)/(1+exp(-Bh*(Tpcm-Tmh)))^(1/vh);
          dHdT=-(Bh*exp(Bh*(Tmh-Tpcm))*(Ah-Kh))/(vh*(exp(Bh*(Tmh-Tpcm))+1)^(1/vh+1));

        //dHdT=(Kh-Ah)*vh*((exp(Bh*(Tmh-Tpcm))+1)^(vh-1))*(exp(Bh*(Tmh-Tpcm)))*(-Bh);
        //dHdT=(Kc-Ac)*vc*((exp(Bc*(Tmc-Tpcm))+1)^(vc-1))*(exp(Bc*(Tmc-Tpcm)))*(-Bc);
        //Kc_heatingTocooling=Kc;
          mode=3.0; //"switch cooling to heating mode"
        else
          H= Ac+(Kc-Ac)/(1+exp(-Bc*(Tpcm-Tmc)))^(1/vc);
          dHdT=-(Bc*exp(Bc*(Tmc - Tpcm))*(Ac - Kc))/(vc*(exp(Bc*(Tmc - Tpcm)) + 1)^(1/vc + 1));
        //dHdT=(Kc-Ac)*vc*((exp(Bc*(Tmc-Tpcm))+1)^(vc-1))*(exp(Bc*(Tmc-Tpcm)))*(-Bc);
        //Ah_coolingToheating=Ah;
        //Kc_heatingTocooling=Kc;
          mode=4.0; //"continue cooling mode"
        end if;
      end if;

  algorithm
    when pre(heating) and Tair>Tpcm then
    //when {mode>=0.5 and mode<=1.5} then
      heating:=true;
      //dHdT:=-(Bh*exp(Bh*(Tmh-Tpcm))*(Ah-Kh))/(vh*(exp(Bh*(Tmh-Tpcm))+1)^(1/vh+1));
    end when;
    when pre(heating) and Tair<=Tpcm then
    //when {mode>=1.5 and mode<=2.5} then
      heating:=false;
      Kc:=(pre(H) - Ac)*((1 + exp(Bc*Tmc - Bc*Tpcm))^(1/vc)) + Ac;
      //dHdT:=-(Bc*exp(Bc*(Tmc - Tpcm))*(Ac - Kc))/(vc*(exp(Bc*(Tmc - Tpcm)) + 1)^(1/vc + 1));
    end when;
    when not pre(heating) and Tair>Tpcm then
    //when {mode>=2.5 and mode<=3.5} then
      heating:=true;
      Ah:=(pre(H) - Kh/((1 + exp(-Bh*(Tpcm - Tmh)))^(1/vh)))/(1 - 1/((1 + exp(-Bh*(Tpcm - Tmh)))^(1/vh)));
      //dHdT:=-(Bh*exp(Bh*(Tmh-Tpcm))*(Ah-Kh))/(vh*(exp(Bh*(Tmh-Tpcm))+1)^(1/vh+1));
    end when;
    when not pre(heating) and Tair<=Tpcm then
    //when {mode>=3.5 and mode<=4.5} then
      heating:=false;
      //Ah:=(H - Kh/((1 + exp(-Bh*(Tpcm - Tmh)))^(1/vh)))/(1 - 1/((1 + exp(-Bh*(Tpcm - Tmh)))^(1/vh)));
      //dHdT:=-(Bc*exp(Bc*(Tmc - Tpcm))*(Ac - Kc))/(vc*(exp(Bc*(Tmc - Tpcm)) + 1)^(1/vc + 1));
    end when;


   annotation (Line(points={{121,-25},{111.5,-25},{111.5,-25},
           {121,-25}}, color={0,0,127}),
                Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{-50,40},{48,-42}},
            lineColor={170,213,255},
            fillColor={28,108,200},
            fillPattern=FillPattern.Backward), Polygon(
            points={{-20,26},{-20,-30},{32,0},{-20,26}},
            lineColor={170,213,255},
            fillColor={238,46,47},
            fillPattern=FillPattern.Backward,
            lineThickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end LiquidFraction_modify_2_2;

  package PCM_HX
    extends Buildings.BaseClasses.BaseIcon;
    package Air = Buildings.Media.Air;

    model PCM_HX_partial_modify_2_media

      Modelica.SIunits.Area A=0.135
      "heat transfer area ";
      parameter Real Tinit=22.0;

      HeatCapacitor_PCM heatCapacitor(
        hf_pcm=140000,
        der_T(start=0, fixed=false),
        m_pcm=2,
        T(start=288.15))
        annotation (Placement(transformation(extent={{-38,54},{-24,68}})));
      convection_coefficient convection_coefficient1
        annotation (Placement(transformation(extent={{-88,-32},{-70,-12}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=15)
                annotation (Placement(transformation(
            extent={{-5,-5},{5,5}},
            rotation=90,
            origin={-33,23})));
      Modelica.Thermal.HeatTransfer.Components.Convection convection annotation (
          Placement(transformation(
            extent={{-7,7},{7,-7}},
            rotation=-90,
            origin={-33,-23})));
      Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor Tpcm
        annotation (Placement(transformation(extent={{-16,44},{-8,52}})));
      Modelica.Blocks.Math.Gain Area(k=0.135)
        annotation (Placement(transformation(extent={{-56,-28},{-46,-18}})));
      LiquidFraction_modify_2_2 Lf(Tmc=20.65 + 273.15, Tmh=23.9 + 273.15)
        annotation (Placement(transformation(extent={{6,36},{26,56}})));
      Modelica.Blocks.Interfaces.RealInput Vair annotation (Placement(
            transformation(extent={{-124,-32},{-100,-8}}),
                                                         iconTransformation(
              extent={{-122,-30},{-100,-8}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (
          Placement(transformation(extent={{-10,-56},{8,-38}}),
            iconTransformation(extent={{-12,-58},{8,-38}})));
      Modelica.Blocks.Interfaces.RealInput Tair "K"
        annotation (Placement(transformation(extent={{-124,8},{-100,32}}),
            iconTransformation(extent={{-124,8},{-100,32}})));
    equation
      connect(convection.solid, thermalConductor.port_a)
        annotation (Line(points={{-33,-16},{-33,18}},color={191,0,0}));
      connect(thermalConductor.port_b, heatCapacitor.port) annotation (Line(
            points={{-33,28},{-33,54},{-31,54}}, color={191,0,0}));
      connect(heatCapacitor.port, Tpcm.port) annotation (Line(points={{-31,54},
              {-24,54},{-24,48},{-16,48}}, color={191,0,0}));
      connect(convection.Gc, Area.y)
        annotation (Line(points={{-40,-23},{-45.5,-23}}, color={0,0,127}));
      connect(Tpcm.T, Lf.Tpcm)
        annotation (Line(points={{-8,48},{-8,48.1},{4.5,48.1}}, color={0,0,127}));
      connect(Tair, Lf.Tair) annotation (Line(points={{-112,20},{-6,20},{-6,
              44.1},{4.5,44.1}},
                      color={0,0,127}));
      connect(Lf.H, heatCapacitor.H) annotation (Line(points={{27.5,48.9},{30,
              48.9},{30,62},{-48,62},{-48,63.24},{-39.4,63.24}}, color={0,0,127}));
      connect(Lf.dHdT, heatCapacitor.dHdT) annotation (Line(points={{27.5,43.9},
              {34,43.9},{34,66},{-50,66},{-50,59.6},{-39.4,59.6}}, color={0,0,
              127}));
      connect(convection_coefficient1.h, Area.u) annotation (Line(points={{-69.1,
              -20.8},{-60.7,-20.8},{-60.7,-23},{-57,-23}},
                                                     color={0,0,127}));
      connect(convection.fluid, port_a) annotation (Line(points={{-33,-30},{-32,
              -30},{-32,-47},{-1,-47}},
                                   color={191,0,0}));
      connect(Vair, convection_coefficient1.Vh) annotation (Line(points={{-112,
              -20},{-100,-20},{-100,-21.5},{-89.35,-21.5}},
                                                      color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -140,-100},{100,100}}),                             graphics={
            Rectangle(
              extent={{-40,38},{40,-40}},
              lineColor={244,125,35},
              lineThickness=1,
              fillColor={244,125,35},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-20,22},{20,-16}},
              lineColor={0,0,0},
              lineThickness=1,
              fillColor={244,125,35},
              fillPattern=FillPattern.None),
            Text(
              extent={{-16,66},{16,42}},
              lineColor={0,0,0},
              lineThickness=1,
              fillColor={244,125,35},
              fillPattern=FillPattern.None,
              textString="PCM module
",            fontName="Arial Black")}),                             Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{
                100,100}})),
        experiment(StopTime=432000));
    end PCM_HX_partial_modify_2_media;

    model sensitivity_analysis
      PCM_HX_partial_modify_2_media pcm
        annotation (Placement(transformation(extent={{-28,24},{-8,44}})));
      Modelica.Blocks.Sources.RealExpression T_air(y=Tair_1.y[1] + 273.15)
        annotation (Placement(transformation(extent={{-50,32},{-36,44}})));
      Modelica.Blocks.Sources.RealExpression V_air1(y=Vair.k)
        annotation (Placement(transformation(extent={{-50,20},{-36,32}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
        redeclare package Medium = Buildings.Media.Air,
        tau=0.1,
        m_flow_nominal=0.1)
        annotation (Placement(transformation(extent={{-28,-8},{-22,0}})));
      Modelica.Blocks.Math.Gain rho(k=1.288)
        annotation (Placement(transformation(extent={{-80,2},{-74,8}})));
      Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
        annotation (Placement(transformation(extent={{-80,-16},{-74,-10}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_m_flow_in=true,
        redeclare package Medium = Buildings.Media.Air,
        nPorts=1,
        use_T_in=true)
        annotation (Placement(transformation(extent={{-50,-10},{-38,2}})));
      Buildings.Fluid.Sources.Boundary_pT bou(redeclare package Medium =
            Buildings.Media.Air, nPorts=1) annotation (Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=180,
            origin={24,-4})));
      Modelica.Blocks.Sources.CombiTimeTable Tair_1(
        tableOnFile=true,
        timeScale=60,
        tableName="Tair1",
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/Tair1.txt")
        annotation (Placement(transformation(extent={{-100,-18},{-90,-8}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem2(
        redeclare package Medium = Buildings.Media.Air,
        tau=0.1,
        m_flow_nominal=0.1)
        annotation (Placement(transformation(extent={{2,-8},{8,0}})));
      Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
        annotation (Placement(transformation(
            extent={{-4,4},{4,-4}},
            rotation=-90,
            origin={-18,16})));
      Modelica.Blocks.Sources.RealExpression OP(y=pcm.heatCapacitor.T - 273.15)
        annotation (Placement(transformation(extent={{54,-6},{68,6}})));
      Modelica.Blocks.Interfaces.RealOutput Tpcm
        annotation (Placement(transformation(extent={{100,-10},{120,10}})));
      Mixing_volume mixing_volume
        annotation (Placement(transformation(extent={{-18,-10},{-6,2}})));
      Modelica.Blocks.Sources.Constant Vair(k=0.1388) "m3/s"
        annotation (Placement(transformation(extent={{-100,26},{-90,36}})));
      Modelica.Blocks.Sources.Constant Vair1(k=0.1388/50) "m3/s"
        annotation (Placement(transformation(extent={{-100,0},{-90,10}})));
    equation
      connect(T_air.y, pcm.Tair) annotation (Line(points={{-35.3,38},{-25.6667,
              38},{-25.6667,36}},
                          color={0,0,127}));
      connect(V_air1.y, pcm.Vair) annotation (Line(points={{-35.3,26},{-30,26},
              {-30,32.1},{-25.5833,32.1}},
                                        color={0,0,127}));
      connect(rho.y, boundary.m_flow_in) annotation (Line(points={{-73.7,5},{
              -64,5},{-64,0.8},{-51.2,0.8}}, color={0,0,127}));
      connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-73.7,
              -13},{-64,-13},{-64,-1.6},{-51.2,-1.6}}, color={0,0,127}));
      connect(boundary.ports[1], senTem1.port_a)
        annotation (Line(points={{-38,-4},{-28,-4}}, color={0,127,255}));
      connect(Tair_1.y[1], toKelvin.Celsius)
        annotation (Line(points={{-89.5,-13},{-80.6,-13}}, color={0,0,127}));
      connect(senTem2.port_b, bou.ports[1])
        annotation (Line(points={{8,-4},{20,-4}}, color={0,127,255}));
      connect(pcm.port_a, heatFlowSensor.port_a) annotation (Line(points={{-16.5,
              29.2},{-16.5,28},{-18,28},{-18,20}},       color={191,0,0}));
      connect(OP.y, Tpcm)
        annotation (Line(points={{68.7,0},{110,0}}, color={0,0,127}));
      connect(Vair1.y, rho.u) annotation (Line(points={{-89.5,5},{-84.75,5},{
              -84.75,5},{-80.6,5}}, color={0,0,127}));
      connect(senTem1.port_b, mixing_volume.port_a)
        annotation (Line(points={{-22,-4},{-18,-4}},   color={0,127,255}));
      connect(mixing_volume.port_b, senTem2.port_a)
        annotation (Line(points={{-6.6,-4},{2,-4}},
                                                  color={0,127,255}));
      connect(mixing_volume.Heatport_a, heatFlowSensor.port_b) annotation (Line(
            points={{-12.24,2.36},{-12.24,12},{-18,12}},
                                                      color={191,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},
                {100,100}})), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-140,-100},{100,100}})));
    end sensitivity_analysis;

    model para_est
      PCM_HX_partial_modify_2_media pcm
        annotation (Placement(transformation(extent={{-28,22},{-8,42}})));
      Modelica.Blocks.Sources.RealExpression V_air1(y=Vair.k)
        annotation (Placement(transformation(extent={{-50,20},{-36,32}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
        redeclare package Medium = Buildings.Media.Air,
        tau=0.1,
        m_flow_nominal=0.1)
        annotation (Placement(transformation(extent={{-28,-8},{-22,0}})));
      Modelica.Blocks.Math.Gain rho(k=1.288)
        annotation (Placement(transformation(extent={{-80,2},{-74,8}})));
      Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
        annotation (Placement(transformation(extent={{-80,-16},{-74,-10}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_m_flow_in=true,
        redeclare package Medium = Buildings.Media.Air,
        nPorts=1,
        use_T_in=true)
        annotation (Placement(transformation(extent={{-50,-10},{-38,2}})));
      Buildings.Fluid.Sources.Boundary_pT bou(redeclare package Medium =
            Buildings.Media.Air, nPorts=1) annotation (Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=180,
            origin={24,-4})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem2(
        redeclare package Medium = Buildings.Media.Air,
        tau=0.1,
        m_flow_nominal=0.1)
        annotation (Placement(transformation(extent={{2,-8},{8,0}})));
      Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
        annotation (Placement(transformation(
            extent={{-4,4},{4,-4}},
            rotation=-90,
            origin={-18,16})));
      Modelica.Blocks.Sources.RealExpression OP(y=pcm.heatCapacitor.T)
        annotation (Placement(transformation(extent={{54,-6},{68,6}})));
      Modelica.Blocks.Interfaces.RealOutput Tpcm
        annotation (Placement(transformation(extent={{100,-10},{120,10}})));
      Modelica.Blocks.Interfaces.RealInput Tair
        annotation (Placement(transformation(extent={{-186,-20},{-146,20}})));
      Modelica.Blocks.Math.Add add
        annotation (Placement(transformation(extent={{-106,32},{-96,42}})));
      Modelica.Blocks.Sources.Constant const(k=273.15)
        annotation (Placement(transformation(extent={{-134,36},{-126,44}})));
      Mixing_volume mixing_volume
        annotation (Placement(transformation(extent={{-16,-10},{-4,2}})));
      Modelica.Blocks.Sources.Constant Vair1(k=0.1388/50) "m3/s"
        annotation (Placement(transformation(extent={{-98,0},{-88,10}})));
      Modelica.Blocks.Sources.Constant Vair(k=0.1388) "m3/s"
        annotation (Placement(transformation(extent={{-98,18},{-88,28}})));
    equation
      connect(V_air1.y, pcm.Vair) annotation (Line(points={{-35.3,26},{-30,26},
              {-30,30.1},{-25.5833,30.1}},
                                        color={0,0,127}));
      connect(rho.y, boundary.m_flow_in) annotation (Line(points={{-73.7,5},{
              -64,5},{-64,0.8},{-51.2,0.8}}, color={0,0,127}));
      connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-73.7,
              -13},{-64,-13},{-64,-1.6},{-51.2,-1.6}}, color={0,0,127}));
      connect(boundary.ports[1], senTem1.port_a)
        annotation (Line(points={{-38,-4},{-28,-4}}, color={0,127,255}));
      connect(senTem2.port_b, bou.ports[1])
        annotation (Line(points={{8,-4},{20,-4}}, color={0,127,255}));
      connect(pcm.port_a, heatFlowSensor.port_a) annotation (Line(points={{-16.5,
              27.2},{-16.5,28},{-18,28},{-18,20}},       color={191,0,0}));
      connect(OP.y, Tpcm)
        annotation (Line(points={{68.7,0},{110,0}}, color={0,0,127}));
      connect(Tair, toKelvin.Celsius) annotation (Line(points={{-166,0},{-122,0},
              {-122,-13},{-80.6,-13}}, color={0,0,127}));
      connect(Tair, add.u2) annotation (Line(points={{-166,0},{-122,0},{-122,34},
              {-107,34}}, color={0,0,127}));
      connect(const.y, add.u1) annotation (Line(points={{-125.6,40},{-107,40}},
                                       color={0,0,127}));
      connect(add.y, pcm.Tair) annotation (Line(points={{-95.5,37},{-58,37},{
              -58,34},{-25.6667,34}},color={0,0,127}));
      connect(senTem2.port_a, mixing_volume.port_b)
        annotation (Line(points={{2,-4},{-4.6,-4}}, color={0,127,255}));
      connect(senTem1.port_b, mixing_volume.port_a)
        annotation (Line(points={{-22,-4},{-16,-4}},    color={0,127,255}));
      connect(heatFlowSensor.port_b, mixing_volume.Heatport_a) annotation (Line(
            points={{-18,12},{-18,10},{-10.24,10},{-10.24,2.36}}, color={191,0,
              0}));
      connect(Vair1.y, rho.u) annotation (Line(points={{-87.5,5},{-80.6,5}},
                                    color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},
                {100,100}})), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-140,-100},{100,100}})));
    end para_est;

    model Mixing_volume
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{80,-10},{100,10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a Heatport_a
        annotation (Placement(transformation(extent={{-14,96},{6,116}}),
            iconTransformation(extent={{-14,96},{6,116}})));
      Buildings.Fluid.Sensors.Temperature Tout(redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{58,0},{78,20}})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
        prescribedHeatFlow
        annotation (Placement(transformation(extent={{-22,4},{-4,22}})));
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort
        mixingVolumeHeatPort(
        m_flow_nominal=0.1,
        nPorts=2,
        redeclare package Medium = Air,
        V=1) annotation (Placement(transformation(extent={{18,4},{36,22}})));
      Modelica.Blocks.Sources.RealExpression Q_flow(y=Heatport_a.Q_flow)
        annotation (Placement(transformation(extent={{-52,8},{-40,18}})));
      Buildings.Fluid.Sensors.Temperature Tin(redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{-90,0},{-70,20}})));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
        prescribedTemperature
        annotation (Placement(transformation(extent={{-5,-5},{5,5}},
            rotation=90,
            origin={-5,61})));
      Modelica.Blocks.Sources.RealExpression T1(y=Tin.T)
                                                        "air Temperature"
        annotation (Placement(transformation(
            extent={{-6,-5},{6,5}},
            rotation=90,
            origin={-4,39})));
    equation
      connect(Tout.port, port_b)
        annotation (Line(points={{68,0},{90,0}}, color={0,127,255}));
      connect(mixingVolumeHeatPort.ports[1], Tout.port) annotation (Line(points={{25.2,4},
              {28,4},{28,0},{68,0}},          color={0,127,255}));
      connect(prescribedHeatFlow.port, mixingVolumeHeatPort.heatPort)
        annotation (Line(points={{-4,13},{18,13}}, color={191,0,0}));
      connect(Q_flow.y, prescribedHeatFlow.Q_flow)
        annotation (Line(points={{-39.4,13},{-22,13}}, color={0,0,127}));
      connect(port_a, Tin.port)
        annotation (Line(points={{-100,0},{-80,0}},color={0,127,255}));
      connect(Tin.port, mixingVolumeHeatPort.ports[2]) annotation (Line(points={{-80,0},
              {26,0},{26,4},{28.8,4}},          color={0,127,255}));
      connect(prescribedTemperature.port, Heatport_a)
        annotation (Line(points={{-5,66},{-4,66},{-4,106}}, color={191,0,0}));
      connect(prescribedTemperature.T, T1.y) annotation (Line(points={{-5,55},{
              -5,50.5},{-4,50.5},{-4,45.6}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                  Ellipse(
              extent={{-94,96},{88,-88}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Sphere,
              fillColor={170,213,255}),   Text(
              extent={{-154,112},{146,152}},
              textString="%name",
              lineColor={0,0,255})}), Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end Mixing_volume;

    annotation (
      uses(Buildings(version="6.0.0")));
  end PCM_HX;

  annotation (
    uses(Buildings(version="6.0.0"), Modelica(version="3.2.2"),
      OU44(version="2")),
              uses(Modelica(version="3.2.2"), Buildings(version="6.0.0"),
      IBPSA(version="3.0.0")));
end PCM;
