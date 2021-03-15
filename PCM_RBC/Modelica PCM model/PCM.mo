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
    parameter Modelica.SIunits.Distance H_c=0.0015;  // Height of air channel between CSM
    parameter Modelica.SIunits.Distance W_CSM=0.3; // Width of CSM plate
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
     //N_channels=m_PCM/2-1;
     N_channels=96;
     A_PCM=(m_PCM/2-1)*(L_CSM*W_CSM)*2;
     Vh=v_a*(N_channels/2*H_c*L_CSM);
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
        m_pcm=1,
        T(start=288.15))
        annotation (Placement(transformation(extent={{-38,52},{-24,66}})));
      convection_coefficient convection_coefficient1
        annotation (Placement(transformation(extent={{-86,-30},{-74,-14}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=15)
                annotation (Placement(transformation(
            extent={{-5,-5},{5,5}},
            rotation=90,
            origin={-33,21})));
      Modelica.Thermal.HeatTransfer.Components.Convection convection annotation (
          Placement(transformation(
            extent={{-7,7},{7,-7}},
            rotation=-90,
            origin={-33,-23})));
      Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor Tpcm
        annotation (Placement(transformation(extent={{-16,44},{-8,52}})));
      Modelica.Blocks.Math.Gain Area(k=0.135)
        annotation (Placement(transformation(extent={{-56,-28},{-46,-18}})));
      LiquidFraction_modify_2_2 Lf(Tmc=21.56 + 273.15, Tmh=22.32 + 273.15)
        annotation (Placement(transformation(extent={{6,36},{26,56}})));
      Modelica.Blocks.Interfaces.RealInput Vair annotation (Placement(
            transformation(extent={{-124,-34},{-100,-10}}),
                                                         iconTransformation(
              extent={{-122,-32},{-100,-10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (
          Placement(transformation(extent={{-10,-56},{8,-38}}),
            iconTransformation(extent={{-12,-58},{8,-38}})));
      Modelica.Blocks.Interfaces.RealInput Tair "K"
        annotation (Placement(transformation(extent={{-124,8},{-100,32}}),
            iconTransformation(extent={{-124,8},{-100,32}})));
    equation
      connect(convection.solid, thermalConductor.port_a)
        annotation (Line(points={{-33,-16},{-33,16}},color={191,0,0}));
      connect(thermalConductor.port_b, heatCapacitor.port) annotation (Line(
            points={{-33,26},{-33,52},{-31,52}}, color={191,0,0}));
      connect(heatCapacitor.port, Tpcm.port) annotation (Line(points={{-31,52},
              {-24,52},{-24,48},{-16,48}}, color={191,0,0}));
      connect(convection.Gc, Area.y)
        annotation (Line(points={{-40,-23},{-45.5,-23}}, color={0,0,127}));
      connect(Tpcm.T, Lf.Tpcm)
        annotation (Line(points={{-8,48},{-8,48.1},{4.5,48.1}}, color={0,0,127}));
      connect(Tair, Lf.Tair) annotation (Line(points={{-112,20},{-6,20},{-6,
              44.1},{4.5,44.1}},
                      color={0,0,127}));
      connect(Lf.H, heatCapacitor.H) annotation (Line(points={{27.5,48.9},{30,
              48.9},{30,62},{-48,62},{-48,61.24},{-39.4,61.24}}, color={0,0,127}));
      connect(Lf.dHdT, heatCapacitor.dHdT) annotation (Line(points={{27.5,43.9},
              {34,43.9},{34,66},{-50,66},{-50,57.6},{-39.4,57.6}}, color={0,0,
              127}));
      connect(convection_coefficient1.h, Area.u) annotation (Line(points={{-73.4,
              -21.04},{-60.7,-21.04},{-60.7,-23},{-57,-23}},
                                                     color={0,0,127}));
      connect(convection.fluid, port_a) annotation (Line(points={{-33,-30},{-32,
              -30},{-32,-47},{-1,-47}},
                                   color={191,0,0}));
      connect(Vair, convection_coefficient1.Vh) annotation (Line(points={{-112,
              -22},{-100,-22},{-100,-21.6},{-86.9,-21.6}},
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

    annotation (
      uses(Buildings(version="6.0.0")));
  end PCM_HX;

  package PCM_MPC

    model AHU_PCM_4stack
      PCM_HX.PCM_HX_partial_modify_2_media PCM1(Tinit=22.0)
        annotation (Placement(transformation(extent={{-64,12},{-44,32}})));
      PCM_HX.PCM_HX_partial_modify_2_media PCM2(Tinit=22.4)
        annotation (Placement(transformation(extent={{-26,10},{-6,30}})));
      PCM_HX.PCM_HX_partial_modify_2_media PCM3(Tinit=22.2)
        annotation (Placement(transformation(extent={{12,10},{32,30}})));
      PCM_HX.PCM_HX_partial_modify_2_media PCM4(Tinit=22.3)
        annotation (Placement(transformation(extent={{46,10},{66,30}})));
      Modelica.Blocks.Sources.RealExpression T_air(y=Te + 273.15)
        annotation (Placement(transformation(extent={{-86,18},{-72,30}})));
      Modelica.Blocks.Sources.RealExpression V_air1(y=gain1.y)
        annotation (Placement(transformation(extent={{-86,6},{-72,18}})));
      Modelica.Blocks.Sources.RealExpression V_air2(y=gain1.y)
        annotation (Placement(transformation(extent={{-44,6},{-30,18}})));
      Modelica.Blocks.Sources.RealExpression V_air3(y=gain1.y)
        annotation (Placement(transformation(extent={{-8,6},{6,18}})));
      Modelica.Blocks.Sources.RealExpression V_air4(y=gain1.y)
        annotation (Placement(transformation(extent={{28,6},{42,18}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem2(
                                 m_flow_nominal=0.1,
        redeclare package Medium = Air,
        T_start=295.15)
        annotation (Placement(transformation(extent={{-36,-20},{-30,-12}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
        m_flow_nominal=0.1,
        tau=1,
        redeclare package Medium = Air,
        T_start=294.05)
        annotation (Placement(transformation(extent={{-68,-20},{-62,-12}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem3(
                                 m_flow_nominal=0.1,
        redeclare package Medium = Air,
        T_start=295.25)
        annotation (Placement(transformation(extent={{2,-20},{8,-12}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem4(
        m_flow_nominal=0.1,
        tau=1,
        redeclare package Medium = Air,
        T_start=295.35)
        annotation (Placement(transformation(extent={{36,-20},{42,-12}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTem5(
                                 m_flow_nominal=0.1,
        redeclare package Medium = Air,
        T_start=295.65)
        annotation (Placement(transformation(extent={{70,-20},{76,-12}})));
      Modelica.Blocks.Sources.RealExpression T_air2(y=senTem2.T)
        annotation (Placement(transformation(extent={{-44,18},{-30,30}})));
      Modelica.Blocks.Sources.RealExpression T_air3(y=senTem3.T)
        annotation (Placement(transformation(extent={{-8,16},{6,28}})));
      Modelica.Blocks.Sources.RealExpression T_air4(y=senTem4.T)
        annotation (Placement(transformation(extent={{28,16},{42,28}})));
      Modelica.Blocks.Math.Gain rho(k=1.288)
        annotation (Placement(transformation(extent={{-116,-12},{-108,-4}})));
      Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
        annotation (Placement(transformation(extent={{-118,-36},{-108,-26}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_m_flow_in=true,
        nPorts=1,
        use_T_in=true,
        redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{-94,-22},{-82,-10}})));
      Buildings.Fluid.Sources.Boundary_pT bou1(
                                 nPorts=1, redeclare package Medium = Air)
                                           annotation (Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=180,
            origin={112,-16})));
      Modelica.Blocks.Sources.CombiTimeTable T_one_node_1(
        tableOnFile=true,
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/T_one_node_1.txt",
        timeScale=60,
        tableName="T_one_node")
        annotation (Placement(transformation(extent={{-54,-56},{-44,-46}})));
      Modelica.Blocks.Sources.CombiTimeTable T_pcm(
        tableOnFile=true,
        timeScale=60,
        tableName="Tpcm1",
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/Tpcm1.txt")
        annotation (Placement(transformation(extent={{16,-56},{6,-46}})));
      Modelica.Blocks.Sources.CombiTimeTable Tair_1node_4stacks_Matlab(
        tableOnFile=true,
        timeScale=60,
        tableName="Tair_1node_4stacks",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_one_node_simulation_results/Air_1node_4stacks.txt",
        columns={2,3,4,5}) "Air temperature from one node Matlab simulation"
        annotation (Placement(transformation(extent={{-54,-74},{-44,-64}})));

      Modelica.Blocks.Sources.CombiTimeTable Tpcm_1node_4stacks_Matlab(
        tableOnFile=true,
        timeScale=60,
        tableName="Tpcm_1node_4stacks",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_one_node_simulation_results/Tpcm_1node_4stacks.txt",
        columns={2,3,4,5})
        "pcm temperature from 1 node and 4 stacks Matlab simulation"
        annotation (Placement(transformation(extent={{16,-76},{6,-66}})));

      Modelica.Blocks.Sources.CombiTimeTable Tair_4stacks_meas(
        tableOnFile=true,
        timeScale=60,
        columns={2,3,4,5},
        tableName="Tair_4stacks_meas",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_measurements/Tair_4stacks_meas.txt")
        "Air temperature measurements 4 stacks"
        annotation (Placement(transformation(extent={{-54,-94},{-44,-84}})));
      Modelica.Blocks.Sources.CombiTimeTable Tpcm_4stacks_meas(
        tableOnFile=true,
        timeScale=60,
        columns={2,3,4,5},
        tableName="Tpcm_4stacks_meas",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_measurements/Tpcm_4tacks_meas.txt")
        "pcm temperature measurements 4stacks"
        annotation (Placement(transformation(extent={{16,-94},{6,-84}})));
      Modelica.Blocks.Interfaces.RealOutput Tsup "celsius degree"
        annotation (Placement(transformation(extent={{140,0},{160,20}})));
      Modelica.Blocks.Interfaces.RealOutput Vair_out "m3/s"
        annotation (Placement(transformation(extent={{138,-50},{158,-30}})));
      Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate(redeclare package
          Medium = PCM_HX.Air)
        annotation (Placement(transformation(extent={{80,-10},{92,-22}})));
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
        annotation (Placement(transformation(extent={{102,2},{118,18}})));
      Modelica.Blocks.Interfaces.RealInput VF "volume flowrate[m3/h]"
        annotation (Placement(transformation(extent={{-164,6},{-142,28}})));
      Modelica.Blocks.Interfaces.RealInput Te
        "ambient temperature[celsius degree]"
        annotation (Placement(transformation(extent={{-162,-42},{-140,-20}})));
      Modelica.Blocks.Math.Gain gain1(k=1/3600)
        annotation (Placement(transformation(extent={{-134,14},{-126,22}})));
      Modelica.Blocks.Math.Gain gain2(k=1/50)
        annotation (Placement(transformation(extent={{-132,-12},{-124,-4}})));
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort
        mixingVolumeHeatPort(
        redeclare package Medium = Buildings.Media.Air,
        m_flow_nominal=0.1,
        nPorts=2,
        V=0.135*0.0025)
             annotation (Placement(transformation(extent={{-56,-16},{-46,-6}})));
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort
        mixingVolumeHeatPort1(
        redeclare package Medium = Buildings.Media.Air,
        m_flow_nominal=0.1,
        nPorts=2,
        V=0.135*0.0025)
             annotation (Placement(transformation(extent={{-20,-16},{-10,-6}})));
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort
        mixingVolumeHeatPort2(
        redeclare package Medium = Buildings.Media.Air,
        m_flow_nominal=0.1,
        nPorts=2,
        V=0.135*0.0025)
             annotation (Placement(transformation(extent={{18,-16},{28,-6}})));
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort
        mixingVolumeHeatPort3(
        redeclare package Medium = Buildings.Media.Air,
        m_flow_nominal=0.1,
        nPorts=2,
        V=0.135*0.0025)
             annotation (Placement(transformation(extent={{52,-16},{62,-6}})));
    equation
      connect(T_air.y, PCM1.Tair) annotation (Line(points={{-71.3,24},{-61.6667,
              24}},                    color={0,0,127}));
      connect(V_air1.y, PCM1.Vair) annotation (Line(points={{-71.3,12},{-66,12},
              {-66,19.9},{-61.5833,19.9}}, color={0,0,127}));
      connect(V_air2.y, PCM2.Vair) annotation (Line(points={{-29.3,12},{-26,12},
              {-26,17.9},{-23.5833,17.9}}, color={0,0,127}));
      connect(V_air3.y, PCM3.Vair) annotation (Line(points={{6.7,12},{12,12},{
              12,17.9},{14.4167,17.9}}, color={0,0,127}));
      connect(V_air4.y, PCM4.Vair) annotation (Line(points={{42.7,12},{44,12},{
              44,17.9},{48.4167,17.9}}, color={0,0,127}));
      connect(T_air2.y, PCM2.Tair) annotation (Line(points={{-29.3,24},{-26,24},
              {-26,22},{-23.6667,22}}, color={0,0,127}));
      connect(T_air3.y, PCM3.Tair) annotation (Line(points={{6.7,22},{14.3333,
              22}},                 color={0,0,127}));
      connect(T_air4.y, PCM4.Tair) annotation (Line(points={{42.7,22},{44,22},{
              44,22},{48.3333,22}}, color={0,0,127}));
      connect(rho.y, boundary.m_flow_in) annotation (Line(points={{-107.6,-8},{
              -100,-8},{-100,-11.2},{-94,-11.2}},     color={0,0,127}));
      connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-107.5,
              -31},{-100,-31},{-100,-13.6},{-95.2,-13.6}}, color={0,0,127}));
      connect(boundary.ports[1], senTem1.port_a) annotation (Line(points={{-82,-16},
              {-68,-16}},                          color={0,127,255}));
      connect(senTem5.port_b, volumeFlowRate.port_a)
        annotation (Line(points={{76,-16},{80,-16}}, color={0,127,255}));
      connect(volumeFlowRate.port_b, bou1.ports[1])
        annotation (Line(points={{92,-16},{108,-16}}, color={0,127,255}));
      connect(volumeFlowRate.V_flow, Vair_out) annotation (Line(points={{86,
              -22.6},{86,-40},{148,-40}}, color={0,0,127}));
      connect(senTem5.T, fromKelvin.Kelvin) annotation (Line(points={{73,-11.6},
              {73,10},{100.4,10}}, color={0,0,127}));
      connect(fromKelvin.Celsius, Tsup)
        annotation (Line(points={{118.8,10},{150,10}}, color={0,0,127}));
      connect(Te, toKelvin.Celsius)
        annotation (Line(points={{-151,-31},{-119,-31}}, color={0,0,127}));
      connect(VF, gain1.u) annotation (Line(points={{-153,17},{-142.5,17},{
              -142.5,18},{-134.8,18}}, color={0,0,127}));
      connect(gain1.y, gain2.u) annotation (Line(points={{-125.6,18},{-120,18},
              {-120,4},{-136,4},{-136,-8},{-132.8,-8}}, color={0,0,127}));
      connect(gain2.y, rho.u)
        annotation (Line(points={{-123.6,-8},{-116.8,-8}}, color={0,0,127}));
      connect(senTem1.port_b, mixingVolumeHeatPort.ports[1])
        annotation (Line(points={{-62,-16},{-52,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort.ports[2], senTem2.port_a)
        annotation (Line(points={{-50,-16},{-36,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort.heatPort, PCM1.port_a) annotation (Line(
            points={{-56,-11},{-56,17.2},{-52.5,17.2}}, color={191,0,0}));
      connect(senTem2.port_b, mixingVolumeHeatPort1.ports[1])
        annotation (Line(points={{-30,-16},{-16,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort1.ports[2], senTem3.port_a)
        annotation (Line(points={{-14,-16},{2,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort1.heatPort, PCM2.port_a) annotation (Line(
            points={{-20,-11},{-18,-11},{-18,15.2},{-14.5,15.2}}, color={191,0,
              0}));
      connect(senTem3.port_b, mixingVolumeHeatPort2.ports[1])
        annotation (Line(points={{8,-16},{22,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort2.ports[2], senTem4.port_a)
        annotation (Line(points={{24,-16},{36,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort2.heatPort, PCM3.port_a) annotation (Line(
            points={{18,-11},{22,-11},{22,15.2},{23.5,15.2}}, color={191,0,0}));
      connect(senTem4.port_b, mixingVolumeHeatPort3.ports[1])
        annotation (Line(points={{42,-16},{56,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort3.ports[2], senTem5.port_a)
        annotation (Line(points={{58,-16},{70,-16}}, color={0,127,255}));
      connect(mixingVolumeHeatPort3.heatPort, PCM4.port_a) annotation (Line(
            points={{52,-11},{56,-11},{56,15.2},{57.5,15.2}}, color={191,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,
                -100},{140,100}})),
                              Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-140,-100},{140,100}})),
        experiment(StopTime=432000));
    end AHU_PCM_4stack;

    package Air =Buildings.Media.Air;

    model AHU_PCM_one_stack
      PCM_HX.PCM_HX_partial_modify_2_media PCM1(Tinit=22.0)
        annotation (Placement(transformation(extent={{-64,10},{-44,30}})));
      Modelica.Blocks.Sources.RealExpression T_air(y=Te + 273.15)
        annotation (Placement(transformation(extent={{-86,18},{-72,30}})));
      Modelica.Blocks.Sources.RealExpression V_air1(y=gain1.y)
        annotation (Placement(transformation(extent={{-86,6},{-72,18}})));
      Modelica.Blocks.Math.Gain rho(k=1.288)
        annotation (Placement(transformation(extent={{-116,-12},{-108,-4}})));
      Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
        annotation (Placement(transformation(extent={{-118,-36},{-108,-26}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1,
        redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{-92,-22},{-80,-10}})));
      Modelica.Blocks.Sources.CombiTimeTable T_one_node_1(
        tableOnFile=true,
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/T_one_node_1.txt",
        timeScale=60,
        tableName="T_one_node")
        annotation (Placement(transformation(extent={{-54,-56},{-44,-46}})));
      Modelica.Blocks.Sources.CombiTimeTable T_pcm(
        tableOnFile=true,
        timeScale=60,
        tableName="Tpcm1",
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/Tpcm1.txt")
        annotation (Placement(transformation(extent={{16,-56},{6,-46}})));
      Modelica.Blocks.Sources.CombiTimeTable Tair_1node_4stacks_Matlab(
        tableOnFile=true,
        timeScale=60,
        tableName="Tair_1node_4stacks",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_one_node_simulation_results/Air_1node_4stacks.txt",
        columns={2,3,4,5}) "Air temperature from one node Matlab simulation"
        annotation (Placement(transformation(extent={{-54,-74},{-44,-64}})));

      Modelica.Blocks.Sources.CombiTimeTable Tpcm_1node_4stacks_Matlab(
        tableOnFile=true,
        timeScale=60,
        tableName="Tpcm_1node_4stacks",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_one_node_simulation_results/Tpcm_1node_4stacks.txt",
        columns={2,3,4,5})
        "pcm temperature from 1 node and 4 stacks Matlab simulation"
        annotation (Placement(transformation(extent={{16,-76},{6,-66}})));

      Modelica.Blocks.Sources.CombiTimeTable Tair_4stacks_meas(
        tableOnFile=true,
        timeScale=60,
        columns={2,3,4,5},
        tableName="Tair_4stacks_meas",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_measurements/Tair_4stacks_meas.txt")
        "Air temperature measurements 4 stacks"
        annotation (Placement(transformation(extent={{-54,-94},{-44,-84}})));
      Modelica.Blocks.Sources.CombiTimeTable Tpcm_4stacks_meas(
        tableOnFile=true,
        timeScale=60,
        columns={2,3,4,5},
        tableName="Tpcm_4stacks_meas",
        fileName=
            "C:/Users/taoy/Desktop/Modelica PCM model/air_pcm_measurements/Tpcm_4tacks_meas.txt")
        "pcm temperature measurements 4stacks"
        annotation (Placement(transformation(extent={{16,-94},{6,-84}})));
      Modelica.Blocks.Interfaces.RealOutput Tsup "celsius degree"
        annotation (Placement(transformation(extent={{20,2},{40,22}})));
      Modelica.Blocks.Interfaces.RealOutput Vair_out "m3/s"
        annotation (Placement(transformation(extent={{20,-44},{40,-24}})));
      Modelica.Blocks.Interfaces.RealInput VF "volume flowrate[m3/h]"
        annotation (Placement(transformation(extent={{-164,8},{-142,30}})));
      Modelica.Blocks.Interfaces.RealInput Te
        "ambient temperature[celsius degree]"
        annotation (Placement(transformation(extent={{-162,-42},{-140,-20}})));
      Modelica.Blocks.Math.Gain gain1(k=1/3600)
        annotation (Placement(transformation(extent={{-134,14},{-126,22}})));
      Modelica.Blocks.Math.Gain gain2(k=1/50)
        annotation (Placement(transformation(extent={{-132,-12},{-124,-4}})));
      Modelica.Blocks.Interfaces.RealOutput H "liquid fraction"
        annotation (Placement(transformation(extent={{20,42},{40,62}})));
      Modelica.Blocks.Interfaces.RealOutput dHdT
        annotation (Placement(transformation(extent={{20,24},{40,44}})));
      Modelica.Blocks.Sources.RealExpression lf(y=PCM1.heatCapacitor.H)
        annotation (Placement(transformation(extent={{-6,44},{8,56}})));
      Modelica.Blocks.Sources.RealExpression lf_dT(y=PCM1.heatCapacitor.dHdT)
        annotation (Placement(transformation(extent={{-6,28},{8,40}})));
      Buildings.Fluid.Sensors.Temperature senTem2(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{-68,-16},{-58,-6}})));
      Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate(redeclare package
          Medium = Air)
        annotation (Placement(transformation(extent={{-30,-10},{-18,-22}})));
      Buildings.Fluid.Sources.Boundary_ph bou(
          nPorts=1, redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{12,-22},{0,-10}})));
      Buildings.Fluid.MixingVolumes.MixingVolume vol(
        V=0.45*0.3*0.015/2,
        m_flow_nominal=0.0036,
        nPorts=2,
        energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
        massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
        redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{-56,-16},{-48,-8}})));
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
        annotation (Placement(transformation(extent={{-8,6},{4,18}})));
      Modelica.Blocks.Sources.RealExpression Toutlet(y=Air.temperature_phX(
                p=101325,
                h=vol.hOut_internal,
                X={0.01,0.99}))
        annotation (Placement(transformation(extent={{-36,6},{-22,18}})));
      Buildings.Fluid.Sensors.Temperature senTem1(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{-46,-16},{-36,-6}})));
      Modelica.Blocks.Sources.RealExpression Tpcm(y=PCM1.heatCapacitor.T)
        annotation (Placement(transformation(extent={{-34,64},{-20,78}})));
      Modelica.Blocks.Interfaces.RealOutput Tpcm_1 "pcm temperature"
        annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
        annotation (Placement(transformation(extent={{-2,66},{8,76}})));
    equation
      connect(T_air.y, PCM1.Tair) annotation (Line(points={{-71.3,24},{-68,24},
              {-68,22},{-61.6667,22}}, color={0,0,127}));
      connect(V_air1.y, PCM1.Vair) annotation (Line(points={{-71.3,12},{-66,12},
              {-66,17.9},{-61.5833,17.9}}, color={0,0,127}));
      connect(rho.y, boundary.m_flow_in) annotation (Line(points={{-107.6,-8},{
              -100,-8},{-100,-11.2},{-92,-11.2}},     color={0,0,127}));
      connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-107.5,
              -31},{-100,-31},{-100,-13.6},{-93.2,-13.6}}, color={0,0,127}));
      connect(Te, toKelvin.Celsius)
        annotation (Line(points={{-151,-31},{-119,-31}}, color={0,0,127}));
      connect(VF, gain1.u) annotation (Line(points={{-153,19},{-142.5,19},{
              -142.5,18},{-134.8,18}}, color={0,0,127}));
      connect(gain1.y, gain2.u) annotation (Line(points={{-125.6,18},{-120,18},
              {-120,4},{-136,4},{-136,-8},{-132.8,-8}}, color={0,0,127}));
      connect(gain2.y, rho.u)
        annotation (Line(points={{-123.6,-8},{-116.8,-8}}, color={0,0,127}));
      connect(dHdT, dHdT)
        annotation (Line(points={{30,34},{30,34}}, color={0,0,127}));
      connect(lf.y, H)
        annotation (Line(points={{8.7,50},{20,50},{20,52},{30,52}},
                                                    color={0,0,127}));
      connect(lf_dT.y, dHdT)
        annotation (Line(points={{8.7,34},{30,34}}, color={0,0,127}));
      connect(boundary.ports[1], senTem2.port)
        annotation (Line(points={{-80,-16},{-63,-16}}, color={0,127,255}));
      connect(volumeFlowRate.V_flow, Vair_out) annotation (Line(points={{-24,
              -22.6},{-24,-34},{30,-34}}, color={0,0,127}));
      connect(volumeFlowRate.port_b, bou.ports[1])
        annotation (Line(points={{-18,-16},{0,-16}}, color={0,127,255}));
      connect(vol.heatPort, PCM1.port_a) annotation (Line(points={{-56,-12},{
              -56,15.2},{-52.5,15.2}}, color={191,0,0}));
      connect(fromKelvin.Celsius, Tsup)
        annotation (Line(points={{4.6,12},{30,12}}, color={0,0,127}));
      connect(Tsup, Tsup)
        annotation (Line(points={{30,12},{30,12}}, color={0,0,127}));
      connect(Toutlet.y, fromKelvin.Kelvin)
        annotation (Line(points={{-21.3,12},{-9.2,12}}, color={0,0,127}));
      connect(vol.ports[1], senTem1.port)
        annotation (Line(points={{-52.8,-16},{-41,-16}}, color={0,127,255}));
      connect(senTem1.port, volumeFlowRate.port_a)
        annotation (Line(points={{-41,-16},{-30,-16}}, color={0,127,255}));
      connect(senTem2.port, vol.ports[2])
        annotation (Line(points={{-63,-16},{-51.2,-16}}, color={0,127,255}));
      connect(Tpcm_1, fromKelvin1.Celsius) annotation (Line(points={{30,70},{26,
              70},{26,71},{8.5,71}}, color={0,0,127}));
      connect(fromKelvin1.Kelvin, Tpcm.y)
        annotation (Line(points={{-3,71},{-19.3,71}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,
                -100},{20,100}}), graphics={Rectangle(
              extent={{-106,42},{4,-42}},
              lineColor={28,108,200},
              lineThickness=1), Text(
              extent={{-76,2},{-26,-14}},
              lineColor={28,108,200},
              lineThickness=1,
              fontSize=10,
              textString="AHU
")}),                         Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-140,-100},{20,100}})),
        experiment(StopTime=432000));
    end AHU_PCM_one_stack;

    model AHU_PCM_one_stack_NoData
      PCM_HX.PCM_HX_partial_modify_2_media PCM1(Tinit=22.0)
        annotation (Placement(transformation(extent={{-66,10},{-46,30}})));
      Modelica.Blocks.Sources.RealExpression T_air(y=Te + 273.15)
        annotation (Placement(transformation(extent={{-86,18},{-72,30}})));
      Modelica.Blocks.Sources.RealExpression V_air1(y=gain1.y)
        annotation (Placement(transformation(extent={{-86,6},{-72,18}})));
      Modelica.Blocks.Math.Gain rho(k=1.288)
        annotation (Placement(transformation(extent={{-116,-12},{-108,-4}})));
      Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
        annotation (Placement(transformation(extent={{-118,-36},{-108,-26}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1,
        redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{-92,-22},{-80,-10}})));

      Modelica.Blocks.Interfaces.RealOutput Tsup "celsius degree"
        annotation (Placement(transformation(extent={{20,2},{40,22}})));
      Modelica.Blocks.Interfaces.RealOutput Vair_out "m3/s"
        annotation (Placement(transformation(extent={{20,-44},{40,-24}})));
      Modelica.Blocks.Interfaces.RealInput VF "volume flowrate[m3/h]"
        annotation (Placement(transformation(extent={{-164,8},{-142,30}})));
      Modelica.Blocks.Interfaces.RealInput Te
        "ambient temperature[celsius degree]"
        annotation (Placement(transformation(extent={{-162,-42},{-140,-20}})));
      Modelica.Blocks.Math.Gain gain1(k=1/3600)
        annotation (Placement(transformation(extent={{-134,14},{-126,22}})));
      Modelica.Blocks.Math.Gain gain2(k=1/96)
        annotation (Placement(transformation(extent={{-132,-12},{-124,-4}})));
      Modelica.Blocks.Interfaces.RealOutput H "liquid fraction"
        annotation (Placement(transformation(extent={{20,42},{40,62}})));
      Modelica.Blocks.Interfaces.RealOutput dHdT
        annotation (Placement(transformation(extent={{20,24},{40,44}})));
      Modelica.Blocks.Sources.RealExpression lf(y=PCM1.heatCapacitor.H)
        annotation (Placement(transformation(extent={{-6,44},{8,56}})));
      Modelica.Blocks.Sources.RealExpression lf_dT(y=PCM1.heatCapacitor.dHdT)
        annotation (Placement(transformation(extent={{-6,28},{8,40}})));
      Buildings.Fluid.Sensors.Temperature senTem2(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{-68,-16},{-58,-6}})));
      Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate(redeclare package
          Medium = Air)
        annotation (Placement(transformation(extent={{-30,-10},{-18,-22}})));
      Buildings.Fluid.Sources.Boundary_ph bou(
          nPorts=1, redeclare package Medium = Air)
        annotation (Placement(transformation(extent={{12,-22},{0,-10}})));
      Buildings.Fluid.MixingVolumes.MixingVolume vol(
        m_flow_nominal=0.0036,
        nPorts=2,
        energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
        massDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
        redeclare package Medium = Air,
        V=0.45*0.3*0.0015/2)
        annotation (Placement(transformation(extent={{-56,-16},{-48,-8}})));
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
        annotation (Placement(transformation(extent={{-8,6},{4,18}})));
      Modelica.Blocks.Sources.RealExpression Toutlet(y=Air.temperature_phX(
                p=101325,
                h=vol.hOut_internal,
                X={0.01,0.99}))
        annotation (Placement(transformation(extent={{-36,6},{-22,18}})));
      Buildings.Fluid.Sensors.Temperature senTem1(redeclare package Medium =
            Air)
        annotation (Placement(transformation(extent={{-46,-16},{-36,-6}})));
      Modelica.Blocks.Sources.RealExpression Tpcm(y=PCM1.heatCapacitor.T)
        annotation (Placement(transformation(extent={{-34,64},{-20,78}})));
      Modelica.Blocks.Interfaces.RealOutput Tpcm_1 "pcm temperature"
        annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
        annotation (Placement(transformation(extent={{-2,66},{8,76}})));
    equation
      connect(T_air.y, PCM1.Tair) annotation (Line(points={{-71.3,24},{-68,24},
              {-68,22},{-63.6667,22}}, color={0,0,127}));
      connect(V_air1.y, PCM1.Vair) annotation (Line(points={{-71.3,12},{-66,12},
              {-66,17.9},{-63.5833,17.9}}, color={0,0,127}));
      connect(rho.y, boundary.m_flow_in) annotation (Line(points={{-107.6,-8},{
              -100,-8},{-100,-11.2},{-92,-11.2}},     color={0,0,127}));
      connect(toKelvin.Kelvin, boundary.T_in) annotation (Line(points={{-107.5,
              -31},{-100,-31},{-100,-13.6},{-93.2,-13.6}}, color={0,0,127}));
      connect(Te, toKelvin.Celsius)
        annotation (Line(points={{-151,-31},{-119,-31}}, color={0,0,127}));
      connect(VF, gain1.u) annotation (Line(points={{-153,19},{-142.5,19},{
              -142.5,18},{-134.8,18}}, color={0,0,127}));
      connect(gain1.y, gain2.u) annotation (Line(points={{-125.6,18},{-120,18},
              {-120,4},{-136,4},{-136,-8},{-132.8,-8}}, color={0,0,127}));
      connect(gain2.y, rho.u)
        annotation (Line(points={{-123.6,-8},{-116.8,-8}}, color={0,0,127}));
      connect(dHdT, dHdT)
        annotation (Line(points={{30,34},{30,34}}, color={0,0,127}));
      connect(lf.y, H)
        annotation (Line(points={{8.7,50},{20,50},{20,52},{30,52}},
                                                    color={0,0,127}));
      connect(lf_dT.y, dHdT)
        annotation (Line(points={{8.7,34},{30,34}}, color={0,0,127}));
      connect(boundary.ports[1], senTem2.port)
        annotation (Line(points={{-80,-16},{-63,-16}}, color={0,127,255}));
      connect(volumeFlowRate.V_flow, Vair_out) annotation (Line(points={{-24,
              -22.6},{-24,-34},{30,-34}}, color={0,0,127}));
      connect(volumeFlowRate.port_b, bou.ports[1])
        annotation (Line(points={{-18,-16},{0,-16}}, color={0,127,255}));
      connect(vol.heatPort, PCM1.port_a) annotation (Line(points={{-56,-12},{
              -56,15.2},{-54.5,15.2}}, color={191,0,0}));
      connect(fromKelvin.Celsius, Tsup)
        annotation (Line(points={{4.6,12},{30,12}}, color={0,0,127}));
      connect(Tsup, Tsup)
        annotation (Line(points={{30,12},{30,12}}, color={0,0,127}));
      connect(Toutlet.y, fromKelvin.Kelvin)
        annotation (Line(points={{-21.3,12},{-9.2,12}}, color={0,0,127}));
      connect(vol.ports[1], senTem1.port)
        annotation (Line(points={{-52.8,-16},{-41,-16}}, color={0,127,255}));
      connect(senTem1.port, volumeFlowRate.port_a)
        annotation (Line(points={{-41,-16},{-30,-16}}, color={0,127,255}));
      connect(senTem2.port, vol.ports[2])
        annotation (Line(points={{-63,-16},{-51.2,-16}}, color={0,127,255}));
      connect(Tpcm_1, fromKelvin1.Celsius) annotation (Line(points={{30,70},{26,
              70},{26,71},{8.5,71}}, color={0,0,127}));
      connect(fromKelvin1.Kelvin, Tpcm.y)
        annotation (Line(points={{-3,71},{-19.3,71}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,
                -100},{20,100}}), graphics={Rectangle(
              extent={{-106,42},{4,-42}},
              lineColor={28,108,200},
              lineThickness=1), Text(
              extent={{-76,2},{-26,-14}},
              lineColor={28,108,200},
              lineThickness=1,
              fontSize=10,
              textString="AHU
")}),                         Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-140,-100},{20,100}})),
        experiment(StopTime=432000));
    end AHU_PCM_one_stack_NoData;

    model RBC_vair

      Modelica.Blocks.Interfaces.RealInput Tpcm
        annotation (Placement(transformation(extent={{-132,-8},{-100,24}})));
      Modelica.Blocks.Interfaces.RealInput Tzone
        annotation (Placement(transformation(extent={{-132,24},{-100,56}})));
      Modelica.Blocks.Interfaces.RealInput Tset
        annotation (Placement(transformation(extent={{-132,-68},{-100,-36}})));
      Modelica.Blocks.Interfaces.RealInput Tout
        annotation (Placement(transformation(extent={{-132,-40},{-100,-8}})));
      Modelica.Blocks.Interfaces.BooleanOutput y1
        "on off signal for pcm inlet vair to PCM"
        annotation (Placement(transformation(extent={{100,20},{120,40}})));
      Modelica.Blocks.Interfaces.BooleanOutput y2
        "on off signal for pcm outlet vair to zone "
        annotation (Placement(transformation(extent={{100,-8},{120,12}})));
      Modelica.Blocks.Interfaces.BooleanOutput y3
        "on off signal for Vair to zone directly"
        annotation (Placement(transformation(extent={{100,-36},{120,-16}})));
      Modelica.Blocks.Interfaces.RealOutput Vair
        annotation (Placement(transformation(extent={{100,-72},{120,-52}})));
    equation
      if Tzone>Tset then
        if Tpcm<Tout then
          y1=true;
          y2=true;
          y3=false;
          Vair=850;
        else //Tpcm>Tout
          y1=false;
          y2=false;
          y3=true;
          Vair=850;
        end if;
      else //Tzone<Tset
        if Tpcm<Tout then
          y1=false;
          y2=false;
          y3=false;
          Vair=0;
        else //Tpcm>Tout
          y1=true;
          y2=false;
          y3=false;
          Vair=850;
        end if;
      end if;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end RBC_vair;

    model MPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update
      extends Buildings.BaseClasses.BaseIcon;
      package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"});

      parameter Real TairInit=20;

      OU44.office_model.office_zone_MPC_est_para_KA_4MPC_controller_new
                 mPC_demo4_1(
        RInt=0.05,
        occ_gain=100/5,
        RExt=2,
        tmass=25,
        imass=17.89,
        shgc=7.94/4,
        Vi=486.5/4)
        annotation (Placement(transformation(extent={{88,-14},{120,18}})));
      Modelica.Blocks.Interfaces.RealInput solrad "solar radiation [w/m2]"
        annotation (Placement(transformation(extent={{-136,52},{-120,68}}),
            iconTransformation(extent={{-136,52},{-120,68}})));
      Modelica.Blocks.Interfaces.RealOutput Tin
        "indoor temperature [celsius degree]"
        annotation (Placement(transformation(extent={{160,0},{182,22}})));
      AHU_PCM_one_stack_NoData
                             aHU_PCM_one_stack_NoData
        annotation (Placement(transformation(extent={{-74,-2},{-46,26}})));
      Modelica.Blocks.Interfaces.RealOutput Tsup
        "supply ventilation temperature"
        annotation (Placement(transformation(extent={{160,24},{180,44}})));
      Modelica.Blocks.Interfaces.RealOutput H
        annotation (Placement(transformation(extent={{160,78},{180,98}})));
      Modelica.Blocks.Interfaces.RealOutput dHdT
        annotation (Placement(transformation(extent={{160,50},{180,70}})));
      Modelica.Blocks.Interfaces.RealOutput Qc
        "indoor temperature [celsius degree]"
        annotation (Placement(transformation(extent={{160,-90},{180,-70}}),
            iconTransformation(extent={{160,-90},{180,-70}})));
      Modelica.Blocks.Interfaces.RealOutput qc
        "indoor temperature [celsius degree]"
        annotation (Placement(transformation(extent={{160,-58},{180,-38}}),
            iconTransformation(extent={{160,-58},{180,-38}})));
      Modelica.Blocks.Math.Add add(k1=-1)
        annotation (Placement(transformation(extent={{110,-116},{120,-106}})));
      Modelica.Blocks.Math.Gain cp(k=1.030*1.29/3600) annotation (Placement(
            transformation(
            extent={{-4,-4},{4,4}},
            rotation=-90,
            origin={132,-82})));
      Modelica.Blocks.Math.MultiProduct multiProduct(nu=2)
        annotation (Placement(transformation(extent={{132,-112},{144,-100}})));
      Modelica.Blocks.Interfaces.RealOutput Qsup
        "supplied cooling or heating [kW]"
        annotation (Placement(transformation(extent={{160,-116},{180,-96}}),
            iconTransformation(extent={{160,-116},{180,-96}})));
      Modelica.Blocks.Interfaces.RealOutput Tpcm
        annotation (Placement(transformation(extent={{160,94},{180,114}})));
      Modelica.Blocks.Interfaces.RealInput Tout "solar radiation [w/m2]"
        annotation (Placement(transformation(extent={{-136,32},{-120,48}}),
            iconTransformation(extent={{-136,18},{-120,34}})));
      Modelica.Blocks.Interfaces.RealInput occ "solar radiation [w/m2]"
        annotation (Placement(transformation(extent={{-136,-22},{-120,-6}}),
            iconTransformation(extent={{-136,-22},{-120,-6}})));
      Modelica.Blocks.Interfaces.RealInput Vair "solar radiation [w/m2]"
        annotation (Placement(transformation(extent={{-136,-8},{-120,8}}),
            iconTransformation(extent={{-136,-44},{-120,-28}})));
      Modelica.Blocks.Interfaces.RealOutput Tinterior
        "indoor temperature [celsius degree]"
        annotation (Placement(transformation(extent={{160,-28},{182,-6}}),
            iconTransformation(extent={{160,-28},{182,-6}})));
      Modelica.Blocks.Sources.Constant const(k=0)
        annotation (Placement(transformation(extent={{-60,-36},{-50,-26}})));
      Modelica.Blocks.Logical.Switch switch1
        annotation (Placement(transformation(extent={{-38,-32},{-28,-22}})));
      Modelica.Blocks.Logical.Switch switch2
        annotation (Placement(transformation(extent={{-10,-16},{0,-6}})));
      Modelica.Blocks.Logical.Switch switch3
        annotation (Placement(transformation(extent={{-98,14},{-88,24}})));
      Modelica.Blocks.Sources.Constant const1(k=0)
        annotation (Placement(transformation(extent={{-88,-36},{-96,-28}})));
      Modelica.Blocks.Logical.Switch switch4
        annotation (Placement(transformation(extent={{18,2},{28,12}})));
      RBC_vair
          RBC_pcm
        annotation (Placement(transformation(extent={{-72,130},{-50,154}})));
      Modelica.Blocks.Sources.RealExpression realExpression(y=Tin)
        annotation (Placement(transformation(extent={{-116,140},{-96,160}})));
      Modelica.Blocks.Sources.RealExpression realExpression1(y=Tpcm)
        annotation (Placement(transformation(extent={{-116,110},{-96,130}})));
      Modelica.Blocks.Sources.RealExpression realExpression2(y=Tout)
        annotation (Placement(transformation(extent={{-116,92},{-96,112}})));
      Modelica.Blocks.Sources.Constant const2(k=24)
        annotation (Placement(transformation(extent={{-114,80},{-104,90}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression(y=RBC_pcm.y1)
        annotation (Placement(transformation(extent={{-118,38},{-98,58}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression1(y=RBC_pcm.y3)
        annotation (Placement(transformation(extent={{-60,38},{-40,58}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression2(y=RBC_pcm.y3)
        annotation (Placement(transformation(extent={{-66,-64},{-46,-44}})));
      Modelica.Blocks.Math.Gain gain(k=96)
        annotation (Placement(transformation(extent={{-34,-6},{-26,2}})));
      Modelica.Blocks.Math.Gain gain1(k=3600)
        annotation (Placement(transformation(extent={{-22,-6},{-14,2}})));
      Modelica.Blocks.Sources.BooleanExpression booleanExpression3(y=RBC_pcm.y2)
        annotation (Placement(transformation(extent={{-60,22},{-40,42}})));
      Modelica.Blocks.Interfaces.BooleanOutput y1
        annotation (Placement(transformation(extent={{160,142},{180,162}})));
      Modelica.Blocks.Interfaces.BooleanOutput y2
        annotation (Placement(transformation(extent={{160,132},{180,152}})));
      Modelica.Blocks.Interfaces.BooleanOutput y3
        annotation (Placement(transformation(extent={{160,122},{180,142}})));
      Modelica.Blocks.Interfaces.RealOutput Vair_out
        annotation (Placement(transformation(extent={{160,110},{180,130}})));
      Modelica.Blocks.Interfaces.RealOutput k
        "indoor temperature [celsius degree]" annotation (Placement(
            transformation(extent={{160,-136},{180,-116}}), iconTransformation(
              extent={{160,-90},{180,-70}})));
      Modelica.Blocks.Sources.RealExpression realExpression3(y=((Tout + 273.15)
             + (Tpcm + 273.15))/(Tsup + 273.15))
        annotation (Placement(transformation(extent={{12,-136},{32,-116}})));
      Modelica.Blocks.Discrete.ZeroOrderHold zeroOrderHold(samplePeriod=60)
        annotation (Placement(transformation(extent={{-92,148},{-86,154}})));
      Modelica.Blocks.Discrete.ZeroOrderHold zeroOrderHold1(samplePeriod=60)
        annotation (Placement(transformation(extent={{-90,138},{-84,144}})));
      Modelica.Blocks.Discrete.ZeroOrderHold zeroOrderHold2(samplePeriod=60)
        annotation (Placement(transformation(extent={{-90,100},{-84,106}})));
    equation
      connect(mPC_demo4_1.Solrad,solrad)  annotation (Line(points={{87.3905,
              14.2182},{80,14.2182},{80,60},{-128,60}},
                                            color={0,0,127}));
      connect(mPC_demo4_1.Tin, Tin) annotation (Line(points={{121.829,10.7273},
              {162,10.7273},{162,11},{171,11}},
                                color={0,0,127}));
      connect(aHU_PCM_one_stack_NoData.H, H) annotation (Line(points={{-44.25,
              19.28},{-16,19.28},{-16,88},{170,88}}, color={0,0,127}));
      connect(aHU_PCM_one_stack_NoData.dHdT, dHdT) annotation (Line(points={{-44.25,
              16.76},{0,16.76},{0,60},{170,60}},        color={0,0,127}));
      connect(mPC_demo4_1.qrad, qc) annotation (Line(points={{121.981,-8.03636},
              {121.981,-48},{170,-48}},color={0,0,127}));
      connect(mPC_demo4_1.Qrad, Qc) annotation (Line(points={{121.829,-5.27273},
              {154,-5.27273},{154,-80},{170,-80}}, color={0,0,127}));
      connect(cp.y, multiProduct.u[1])
        annotation (Line(points={{132,-86.4},{132,-103.9}}, color={0,0,127}));
      connect(add.y, multiProduct.u[2]) annotation (Line(points={{120.5,-111},{
              120.5,-108.1},{132,-108.1}}, color={0,0,127}));
      connect(multiProduct.y, Qsup)
        annotation (Line(points={{145.02,-106},{170,-106}}, color={0,0,127}));
      connect(aHU_PCM_one_stack_NoData.Tpcm_1, Tpcm) annotation (Line(points={{-44.25,
              21.8},{26,21.8},{26,104},{170,104}},        color={0,0,127}));
      connect(Tout, add.u2) annotation (Line(points={{-128,40},{-108,40},{-108,
              -114},{109,-114}}, color={0,0,127}));
      connect(Tout, aHU_PCM_one_stack_NoData.Te) annotation (Line(points={{-128,40},
              {-108,40},{-108,7.66},{-75.925,7.66}},     color={0,0,127}));
      connect(Tout, mPC_demo4_1.Tout) annotation (Line(points={{-128,40},{68,40},
              {68,9.56364},{87.5429,9.56364}},
                                            color={0,0,127}));
      connect(occ, mPC_demo4_1.Occ) annotation (Line(points={{-128,-14},{-82,
              -14},{-82,26},{56,26},{56,6.21818},{87.5429,6.21818}},
                                                             color={0,0,127}));
      connect(mPC_demo4_1.Tinterior, Tinterior) annotation (Line(points={{121.981,
              6.8},{156,6.8},{156,-17},{171,-17}},       color={0,0,127}));
      connect(Vair, switch1.u1) annotation (Line(points={{-128,0},{-60,0},{-60,
              -23},{-39,-23}}, color={0,0,127}));
      connect(const.y, switch1.u3)
        annotation (Line(points={{-49.5,-31},{-39,-31}}, color={0,0,127}));
      connect(switch1.y, switch2.u3) annotation (Line(points={{-27.5,-27},{
              -24.75,-27},{-24.75,-15},{-11,-15}}, color={0,0,127}));
      connect(switch2.y, mPC_demo4_1.Vair) annotation (Line(points={{0.5,-11},{
              52.25,-11},{52.25,-1.34545},{87.5429,-1.34545}},  color={0,0,127}));
      connect(switch2.y, cp.u) annotation (Line(points={{0.5,-11},{52,-11},{52,
              -52},{132,-52},{132,-77.2}}, color={0,0,127}));
      connect(switch3.y, aHU_PCM_one_stack_NoData.VF) annotation (Line(points={{-87.5,
              19},{-81.75,19},{-81.75,14.66},{-76.275,14.66}},        color={0,
              0,127}));
      connect(Vair, switch3.u1) annotation (Line(points={{-128,0},{-116,0},{
              -116,23},{-99,23}}, color={0,0,127}));
      connect(const1.y, switch3.u3) annotation (Line(points={{-96.4,-32},{-106,
              -32},{-106,15},{-99,15}}, color={0,0,127}));
      connect(switch4.y, Tsup) annotation (Line(points={{28.5,7},{54,7},{54,34},
              {170,34}}, color={0,0,127}));
      connect(switch4.y, mPC_demo4_1.Tsup) annotation (Line(points={{28.5,7},{
              54,7},{54,2.43636},{87.5429,2.43636}}, color={0,0,127}));
      connect(booleanExpression.y, switch3.u2) annotation (Line(points={{-97,48},
              {-94,48},{-94,30},{-104,30},{-104,19},{-99,19}}, color={255,0,255}));
      connect(booleanExpression1.y, switch4.u2) annotation (Line(points={{-39,48},
              {-14,48},{-14,7},{17,7}},   color={255,0,255}));
      connect(booleanExpression2.y, switch1.u2) annotation (Line(points={{-45,
              -54},{-42,-54},{-42,-27},{-39,-27}}, color={255,0,255}));
      connect(gain.y, gain1.u)
        annotation (Line(points={{-25.6,-2},{-22.8,-2}}, color={0,0,127}));
      connect(aHU_PCM_one_stack_NoData.Vair_out, gain.u) annotation (Line(
            points={{-44.25,7.24},{-44.25,1.62},{-34.8,1.62},{-34.8,-2}}, color=
             {0,0,127}));
      connect(gain1.y, switch2.u1) annotation (Line(points={{-13.6,-2},{-11,-2},
              {-11,-7}}, color={0,0,127}));
      connect(booleanExpression3.y, switch2.u2) annotation (Line(points={{-39,
              32},{-36,32},{-36,-10},{-11,-10},{-11,-11}}, color={255,0,255}));
      connect(aHU_PCM_one_stack_NoData.Tsup, switch4.u3) annotation (Line(
            points={{-44.25,13.68},{8,13.68},{8,3},{17,3}}, color={0,0,127}));
      connect(Tout, switch4.u1) annotation (Line(points={{-128,40},{14,40},{14,
              11},{17,11}}, color={0,0,127}));
      connect(RBC_pcm.y1, y1) annotation (Line(points={{-48.9,145.6},{146,145.6},
              {146,152},{170,152}},
                               color={255,0,255}));
      connect(RBC_pcm.y2, y2) annotation (Line(points={{-48.9,142.24},{150,
              142.24},{150,142},{170,142}},
                                    color={255,0,255}));
      connect(RBC_pcm.y3, y3) annotation (Line(points={{-48.9,138.88},{54.5,
              138.88},{54.5,132},{170,132}},
                                     color={255,0,255}));
      connect(RBC_pcm.Vair, Vair_out) annotation (Line(points={{-48.9,134.56},{
              126,134.56},{126,120},{170,120}},
                                           color={0,0,127}));
      connect(aHU_PCM_one_stack_NoData.Tsup, add.u1) annotation (Line(points={{-44.25,
              13.68},{8,13.68},{8,-108},{109,-108}},        color={0,0,127}));
      connect(realExpression3.y, k)
        annotation (Line(points={{33,-126},{170,-126}}, color={0,0,127}));
      connect(realExpression.y, zeroOrderHold.u) annotation (Line(points={{-95,
              150},{-94,150},{-94,151},{-92.6,151}}, color={0,0,127}));
      connect(zeroOrderHold.y, RBC_pcm.Tzone) annotation (Line(points={{-85.7,
              151},{-78,151},{-78,146.8},{-73.76,146.8}}, color={0,0,127}));
      connect(zeroOrderHold1.y, RBC_pcm.Tpcm) annotation (Line(points={{-83.7,
              141},{-80,141},{-80,142.96},{-73.76,142.96}}, color={0,0,127}));
      connect(zeroOrderHold1.u, realExpression1.y) annotation (Line(points={{
              -90.6,141},{-92,141},{-92,120},{-95,120}}, color={0,0,127}));
      connect(realExpression2.y, zeroOrderHold2.u) annotation (Line(points={{-95,102},
              {-92,102},{-92,103},{-90.6,103}},          color={0,0,127}));
      connect(zeroOrderHold2.y, RBC_pcm.Tout) annotation (Line(points={{-83.7,
              103},{-78,103},{-78,139.12},{-73.76,139.12}}, color={0,0,127}));
      connect(const2.y, RBC_pcm.Tset) annotation (Line(points={{-103.5,85},{-82,
              85},{-82,135.76},{-73.76,135.76}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,
                -140},{160,160}}), graphics={Rectangle(
              extent={{-120,162},{160,-138}},
              lineColor={0,140,72},
              fillColor={28,108,200},
              fillPattern=FillPattern.Forward), Text(
              extent={{-26,210},{66,124}},
              lineColor={28,108,200},
              fillColor={28,108,200},
              fillPattern=FillPattern.Forward,
              textStyle={TextStyle.Bold},
              textString="PCM_Zone
")}),                         Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-120,-140},{160,160}})),
        experiment(StopTime=432000, __Dymola_Algorithm="Dassl"));
    end MPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update;

    model test_RBC_PCMzone_update
      MPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update
        mPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update
        annotation (Placement(transformation(extent={{-36,-16},{-8,14}})));
      Modelica.Blocks.Sources.Pulse pulse4(
        period=86400,
        width=50,
        startTime=43200,
        offset=0,
        amplitude=1)
        annotation (Placement(transformation(extent={{-94,-8},{-88,-2}})));
      Modelica.Blocks.Sources.Pulse pulse1(
        period=86400,
        startTime=0,
        width=50,
        amplitude=0)
        annotation (Placement(transformation(extent={{-94,4},{-88,10}})));
      Modelica.Blocks.Math.Add signal
        annotation (Placement(transformation(extent={{-82,-2},{-78,2}})));
      Modelica.Blocks.Math.Gain solrad(k=400)
        annotation (Placement(transformation(extent={{-64,10},{-60,14}})));
      Modelica.Blocks.Math.Gain occ(k=7)
        annotation (Placement(transformation(extent={{-62,-6},{-58,-2}})));
      Modelica.Blocks.Math.Gain vair(k=850)
        annotation (Placement(transformation(extent={{-62,-14},{-58,-10}})));
      Modelica.Blocks.Sources.CombiTimeTable Tair_1(
        tableOnFile=true,
        timeScale=60,
        tableName="Tair1",
        fileName="C:/Users/taoy/Desktop/Modelica PCM model/Tair1.txt")
        annotation (Placement(transformation(extent={{-72,-28},{-68,-24}})));
      Modelica.Blocks.Math.Add signal1
        annotation (Placement(transformation(extent={{-54,2},{-50,6}})));
      Modelica.Blocks.Sources.Constant const2(k=0)
        annotation (Placement(transformation(extent={{-64,18},{-60,22}})));
      Modelica.Blocks.Sources.Pulse pulse2(
        period=86400,
        width=50,
        startTime=43200,
        offset=0,
        amplitude=25)
        annotation (Placement(transformation(extent={{-92,20},{-86,26}})));
      Modelica.Blocks.Sources.Pulse pulse3(
        period=86400,
        startTime=0,
        width=50,
        amplitude=12)
        annotation (Placement(transformation(extent={{-92,32},{-86,38}})));
      Modelica.Blocks.Math.Add signal2
        annotation (Placement(transformation(extent={{-80,26},{-76,30}})));
      Modelica.Blocks.Sources.Constant const1(k=0)
        annotation (Placement(transformation(extent={{-64,18},{-60,22}})));
      Modelica.Blocks.Sources.Constant const3(k=850)
        annotation (Placement(transformation(extent={{-64,-28},{-60,-24}})));
    equation
      connect(pulse1.y, signal.u1) annotation (Line(points={{-87.7,7},{-83.85,7},
              {-83.85,1.2},{-82.4,1.2}}, color={0,0,127}));
      connect(pulse4.y, signal.u2) annotation (Line(points={{-87.7,-5},{-83.85,
              -5},{-83.85,-1.2},{-82.4,-1.2}}, color={0,0,127}));
      connect(signal.y, solrad.u) annotation (Line(points={{-77.8,0},{-70,0},{
              -70,12},{-64.4,12}}, color={0,0,127}));
      connect(signal.y, occ.u) annotation (Line(points={{-77.8,0},{-70,0},{-70,
              -4},{-62.4,-4}}, color={0,0,127}));
      connect(signal.y, vair.u) annotation (Line(points={{-77.8,0},{-70,0},{-70,
              -12},{-62.4,-12}}, color={0,0,127}));
      connect(solrad.y,
        mPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update.solrad)
        annotation (Line(points={{-59.8,12},{-48,12},{-48,4},{-36.8,4}}, color=
              {0,0,127}));
      connect(occ.y, mPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update.occ)
        annotation (Line(points={{-57.8,-4},{-48,-4},{-48,-3.4},{-36.8,-3.4}},
            color={0,0,127}));
      connect(const2.y, signal1.u1) annotation (Line(points={{-59.8,20},{-56,20},
              {-56,5.2},{-54.4,5.2}}, color={0,0,127}));
      connect(signal1.y,
        mPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update.Tout)
        annotation (Line(points={{-49.8,4},{-48,4},{-48,0.6},{-36.8,0.6}},
            color={0,0,127}));
      connect(pulse3.y, signal2.u1) annotation (Line(points={{-85.7,35},{-81.85,
              35},{-81.85,29.2},{-80.4,29.2}}, color={0,0,127}));
      connect(pulse2.y, signal2.u2) annotation (Line(points={{-85.7,23},{-81.85,
              23},{-81.85,26.8},{-80.4,26.8}}, color={0,0,127}));
      connect(signal2.y, signal1.u2) annotation (Line(points={{-75.8,28},{-70,
              28},{-70,2.8},{-54.4,2.8}}, color={0,0,127}));
      connect(vair.y, mPC_demo_one_stack_KA_ZoneModel_RBC_PCM_zone_copy_update.Vair)
        annotation (Line(points={{-57.8,-12},{-48,-12},{-48,-5.6},{-36.8,-5.6}},
            color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=172800));
    end test_RBC_PCMzone_update;

  end PCM_MPC;

  annotation (
    uses(                            Modelica(version="3.2.2"),
      OU44(version="2"),
      Buildings(version="5.0.1")),
              uses(Modelica(version="3.2.2"), Buildings(version="6.0.0"),
      IBPSA(version="3.0.0")));
end PCM;
