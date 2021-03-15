clear all 
close all

% Data input file with the storage data needed for input and comparison:
load('29-06-17_04-07-17'); % 5 days - 120 hours - 7200 min. total.

% Number of hours for each timestep between each data point in input:
data_time_length=1/60; 

% Number of nodes:
n=10;

% Which data points to use in run:
start=11;                            % First datapoint from storage to use.
time=length(Speicher3_air)-2096-1;   % Number of datapoints to include.

% Input for storages Speicher2(SP25) and Speicher3(SP21):
% input1(1,1:time+1)=(Speicher2.S2F1(start:start+time));    
input1(1,1:time+1)=(Speicher3_air(start:start+time,1)); % Air input [C].
input2(1,1:4)=Speicher3_PCM(start,:);                   % Initial PCM temperature [C]
input2(1:n,:)=ones(n,4).*(input2)-(ones(4,n).*linspace(0,0.1,n))';                                          % OBS


% Number of panels per stack in the corresponding storage unit:
if Speicher2_PCM(start,1)==input2(1,1)
    N_CSM_per_stack=32;
else
    N_CSM_per_stack=25;
end

% Calling function for the model:
[l_f,T_PCM,Q_PCM,time__,T_out]=RT_model(input1(1,:),input2,data_time_length,n,N_CSM_per_stack);

% Plot of data from storage:
plot_RT_data3_better_resolution

% Plot of data from storage with model results:
ComparisonOfDataAndModel0_less_time


% Rubitherm model function:
function [l_f,T_PCM,Q_PCM,time__,T_out]=RT_model(data,T_PCM_i,data_time_length,n,N_CSM_per_stack)
% Script for SP21EK properties:
SP21EK_better_data

% Dimension and basic inputs:
m_CSM=2;                            % PCM mass in each CSM panel.
N_stacks=4;                         % Number of stacks in storage.
m_PCM=N_CSM_per_stack*m_CSM;        % Total mass of PCM in storage.    
time_step=3600*(data_time_length);  % Time step in seconds.
Vh=500;                             % Volume flow rate [m3/h]

% Properties for curve used for transistion:
Ac=0;
Kc=1;
Ah=0;
Kh=1;

% Specific heat capacity of air and PCM:
cp_PCM_sensible=2000;
cp_air=1004.9;


% Initial temperature used in ODE solver.
T1=T_PCM_i;
% Tolerance setting
options=odeset('RelTol',1e-3,'AbsTol',1e-6); % 1e-3, 1e-6 are default values.

% Values for transistion functions:
Actohs(1:n,1:4)=0;
Khtocs(1:n,1:4)=1;

% Predefinition of size to get MATLAB off my back:
lf_i=zeros(n,4);
mode=zeros(n,4);

% Defining initial temperature, liquid fraction etc for each node in each
% stack:
for ad=1:N_stacks    
    
    % Cooling of PCM initialy
    if T1(1:n,ad)>data(1)
        mode(1:n,ad)=2;
        lf_i(1:n,ad)=Ac+(Kc+Ac)./((1+exp(-Bc.*(T1(1:n,ad)-Mc))).^(1/vc));
        Khtocs(1:n,ad)=(lf_i(1:n,ad)-Ac).*(1+exp(-Bc.*(T1(1:n,ad)-Mc))).^(1/vc)+Ac;
    % Heating of PCM initialy:    
    else
        mode(1:n,ad)=6;
        lf_i(1:n,ad)=Ah+(Kh+Ah)./((1+exp(-Bh.*(T1(1:n,ad)-Mh))).^(1/vh));
        Actohs(1:n,ad)=(lf_i(1:n,ad)-Kh./(1+exp(-Bh.*(T1(1:n,ad)-Mh))).^(1/vh))./(1-1./(1+exp(-Bh.*(T1(1:n,ad)-Mh))).^(1/vh));
    end
end

% Loop for time period:
for t=1:length(data)
    % Minutes:
    t
    clear time Tout QPCM QPCM_accumulated TPCM lf t1 T_PCM1 H
    % Input for ODE solver:
    Tspan=[0 time_step];
    
    % Calculation of air density using linear function:
    rho_air=linear_air_density(data(t));
    
    % Mass flow rate per plate:
    ma=rho_air*Vh/((N_CSM_per_stack-1)*3600); % [kg/s]
    
    % Input for initial PCM temperature in ODE-solver for all nodes in all stacks:
    if t==1
        T11=T_PCM_i;
    else
        T11=reshape(T_PCM(end,:),n,N_stacks);
    end
    % Air temperature input:
    T_in=data(t);
    
    % PCM temperature, liquid fraction,the current curve etc. inputs for first stack:
    T1=T11(1:n,1);
    modes=mode(1:n,1);
    Actoh=Actohs(1:n,1);
    Khtoc=Khtocs(1:n,1);
    if t==1
    else
        lf_i=reshape(l_f(end,:),n,N_stacks); 
    end
    
    % Call of the ODE function:
    [t1,T_PCM1]=ode15s(@(t1,T_PCM1)f(t1,T_PCM1,T_in,Vh,m_PCM,lf_i(1:n,1),modes,Actoh,Khtoc,Mc,Mh,Bc,Bh,vc,vh,hf,rho_PCM_s,rho_PCM_l,k_PCM,n),Tspan,T1,options);
    % Script for the keeping track of the hysteresis curve in use:
    step_nodes2   
    
    mode(1:n,1)=modes(1:n,1);
    Actohs(1:n,1)=Actoh(1:n,1);
    Khtocs(1:n,1)=Khtoc(1:n,1);
    
    % PCM temperature and liquid fraction for next timestep:
    TPCM=[T_PCM1(1,1:n); T_PCM1(end,1:end)];
    lf=[H(1,1:n); H(end,1:n)];
    
    % Heat transfer in stack:
    QPCM_accumulated(:,1:n)=m_PCM/n.*cp_PCM_sensible.*(TPCM(:,1:n)-T1(1:n,1)')+m_PCM/n.*hf.*(lf(:,1)-lf_i(1:n,1)'); % Accumulated heat transfer
    
    QPCM=zeros(length(QPCM_accumulated(:,1))-1,n); % Pre-allocation
    % Heat transfer:
    for q=1:length(QPCM_accumulated(:,1))-1
        QPCM(q,1:n)=QPCM_accumulated(q+1,1:n)-QPCM_accumulated(q,1:n);
    end
    
    % Time passed for each data point:
    time(1:2,1)=[t1(1) t1(end)];
    
    % Air temperature after stack:
    Tout=[-QPCM((1:length(QPCM(:,1))),1)./((Vh*rho_air/3600).*(time(2)-time(1)).*cp_air)+ones(length(TPCM(:,1))-1,1).*T_in];
    
    % Stacks 2-4:
    for e=1:N_stacks-1
        % Initial liquid fraction and PCM temperature for ODE-solver:
        lfi(1:n,1)=lf_i(1:n,e+1);
        T1(1:n,1)=T11(1:n,e+1);
        if t==1
        else
            lf_i(1:n,e+1)=l_f(end,(n*e+1):(n*(e+1)))';
        end
        
        % Script for each time output for previous solver:
            clear T_PCM1 H
            % Initial inputs:
%             Tspan=[0 time(w+1)-time(w)];
            modes=mode(1:n,e+1);
            Actoh=Actohs(1:n,e+1);
            Khtoc=Khtocs(1:n,e+1);
            T_in=Tout(e);
            % Call of ODE-solver:
            [t1,T_PCM1]=ode15s(@(t1,T_PCM1)f(t1,T_PCM1,T_in,Vh,m_PCM,lfi(1:n,1),modes,Actoh,Khtoc,Mc,Mh,Bc,Bh,vc,vh,hf,rho_PCM_s,rho_PCM_l,k_PCM,n),Tspan,T1,options);
            
            % Output PCM temperature:
            TPCM(1,(n*e+1):(n*(e+1)))=T11(1:n,e+1);
            TPCM(2,(n*e+1):(n*(e+1)))=T_PCM1(end,1:n);
            
            % Script for the keeping track of the hysteresis curve in use:            
            step_nodes_stacks4
            % Results for next time step and saving results:
            mode(1:n,e+1)=modes(1:n,1);
            Actohs(1:n,e+1)=Actoh(1:n,1);
            Khtocs(1:n,e+1)=Khtoc(1:n,1);
            lf(1,(n*e+1):(n*(e+1)))=lf_i(1:n,e+1);
            lf(2,(n*e+1):(n*(e+1)))=H(end,1:n);
            QPCM_accumulated(1,(n*e+1):(n*(e+1)))=0;
            QPCM_accumulated(2,(n*e+1):(n*(e+1)))=m_PCM/n.*cp_PCM_sensible.*(T_PCM1(end,1:n)-T_PCM1(1,1:n))+m_PCM/n.*hf.*(lf(2,(n*e+1):(n*(e+1)))-lf_i(1:n,e+1)');
            
            QPCM(1,(n*e+1):(n*(e+1)))=QPCM_accumulated(2,(n*e+1):(n*(e+1)))-QPCM_accumulated(1,(n*e+1):(n*(e+1)));
            
            
            % Temperature after each of the stacks:
            Tout(e+1)=-QPCM(1,e+1)./((Vh*rho_air/3600).*(t1(end)-t1(1)).*cp_air)+T_in;
            % initial inputs for next time step:
            T1(1:n,1)=T_PCM1(end,1:n)';
            lfi(1:n,1)=H(end,1:n);
       
        
    end
    % Saving data:
    Tout(1,:)=Tout(length(time)-1,:);
    time_(1:length(time),t)=time_step*(t-1)+time;
    a=find(time_(:,t)==0);
    % Saving data:
    if t>1
        a=find(time_(:,t)==0);
        if isempty(a)==1
            a=length(time_(:,1))+1;
        end
        time__=[time__; time_(1:a-1,t)];
        T_PCM=[T_PCM;TPCM];
        T_out=[T_out;Tout];
        Q_PCM=[Q_PCM;QPCM];
        l_f=[l_f;lf];
        
    else
        time__=time_;
        T_PCM=TPCM;
        T_out=Tout;
        Q_PCM=QPCM;
        l_f=lf;
        
    end
    
    
end
% Script for plotting results from the nodes. Hysteresis curves etc.
RT_model_plot0_nodes
end

% Derivative function for ODE solver.
function dTdt=f(t1,T_PCM1,T_air,Vh,m_PCM,lf_i,modes,Actoh,Khtoc,Mc,Mh,Bc,Bh,vc,vh,hf,rho_PCM_s,rho_PCM_l,k_PCM,n)

% Constants:
cp_PCM_sensible=2000;                       % Specific heat capacity [J/kgK]
W_CSM=.3;                                   % Width of CSM plate [m]
L_CSM=0.45;                                 % Length of CSM plates [m]
D_CSM=0.015;                                % Thickness of CSM plates [m]
m_CSM=2;                                    % Mass of PCM in CSM plate [kg]
N_CSM=m_PCM/m_CSM;                          % Number of CSM plates
dx=D_CSM/(2*n);                                 % Step length between nodes [m]

% Values for hysteresis curves:
Ac=0;
Kc=1;
Ah=0;
Kh=1;

% Conductivity, dynamic viscosity, Prandtl number and density of air
% (assumed constant):
k=0.02364;
mu=0.00001729;
Pr=0.7344;
rho=1.288;
k_a=k(1);
my_a=mu(1);
rho_a=rho(1);
Pr_a=Pr(1);

% Calculating entrance region length and Reynolds number for heat transfer
% coefficient
H_c=5e-3;                           % Height of air channel between plates[m]
L_c=2*H_c;                          % Characteristic length [m]
N_channels=round(N_CSM)+1;          % Number of air channels
v_a=Vh/(N_channels*H_c*L_CSM*3600); % Air velocity in each channel [m/s]
Re=(rho_a*v_a*L_c)/my_a;            % Reynolds number
l=0.05*Re*Pr_a*L_c;                 % Entrance region length [m]

% Nusselt number in entrance region and the remaining:
Nu_entry=7.54+(0.03*(L_c/W_CSM)*Re*Pr_a)/(1+0.016*((L_c/W_CSM)*Re*Pr_a).^(2/3));
Nu_rest=7.54;

% Calculates average nusselt number in channel:
if l>W_CSM
    Nu=Nu_entry;
else
    Nu=Nu_entry*l/W_CSM+Nu_rest*(W_CSM-l)/W_CSM;
end

% Heat transfer coefficient:
U_PCM=Nu*k_a/L_c;       % Actually the h-value, same same but different.

% Prealocation:
dTdt=zeros(n,1);
H=zeros(n,1);

%% For loop for all of the nodes:
for q=1:n
    
    % Set of loops for determining whether which (hysteresis) curve each 
    % node is behaving as and finding the corresponding derivative of the
    % liquid fraction (for use in cp-term):
    % Boundary node connected to the air:
    if q==1
        
        % Cooling of the boundary node:
        if T_PCM1(q,1)>T_air
           
            % Continued cooling:
            if modes(q,1)==1
                H(q,1)=Ac+(Kc-Ac)/(exp(-Bc*(T_PCM1(q,1)-Mc))+1)^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
           
            % Change from heating to cooling:
            elseif modes(q,1)==5 || modes(q,1)==6
                Khtoc(q,1)=(lf_i(q,1)-Ac)*(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc)+Ac;  
                H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
            
            % Continued cooling after change:
            elseif modes(q,1)==2
                H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
            end
            
        % Heating of the boundary node:
        elseif T_PCM1(q,1)<=T_air
            
            % Continued heating:
            if modes(q,1)==5
                H(q,1)=Ah+(Kh-Ah)/(exp(-Bh*(T_PCM1(q,1)-Mh))+1)^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                
            % Change from cooling to heating:
            elseif modes(q,1)==1 || modes(q,1)==2
                Actoh(q,1)=(lf_i(q,1)-Kh/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh))/(1-1/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh));
                H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                
            % Continued heating after change:
            elseif modes(q,1)==6
                H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
            end
        end
        
        % Calculation of density, apparent heat capacity, thermal
        % diffusivity and dT/dt for ODE solver:
        rho_PCM=(rho_PCM_l+rho_PCM_s)/2;
        cp_PCM=(cp_PCM_sensible+dHdT*hf)*H(q,1)+(cp_PCM_sensible+dHdT*hf)*(1-H(q,1));
        alpha_PCM=k_PCM/(rho_PCM*cp_PCM);
        dTdt(q,1)=(2*U_PCM/(rho_PCM*cp_PCM*dx))*(T_air-T_PCM1(q,1))+((2*alpha_PCM)/(dx^2))*(T_PCM1(q+1,1)-T_PCM1(q,1));
        
        % Boundary node at the symmetry boundary (addiabatic):
    elseif q==n
        % Cooling of the boundary node:
        if T_PCM1(q,1)>T_PCM1(q-1,1)
            % Continued cooling:
            if modes(q,1)==1
                H(q,1)=Ac+(Kc-Ac)/(exp(-Bc*(T_PCM1(q,1)-Mc))+1)^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                % Change from heating to cooling:
            elseif modes(q,1)==5 || modes(q,1)==6
                Khtoc(q,1)=(lf_i(q,1)-Ac)*(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc)+Ac;
                H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                % Continued cooling after change:
            elseif modes(q,1)==2
                H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
            end
            % Heating of the boundary node:
        elseif T_PCM1(q,1)<=T_PCM1(q-1,1)
            % Continued heating:
            if modes(q,1)==5
                H(q,1)=Ah+(Kh-Ah)/(exp(-Bh*(T_PCM1(q,1)-Mh))+1)^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                % Change from cooling to heating:
            elseif modes(q,1)==1 || modes(q,1)==2
                Actoh(q,1)=(lf_i(q,1)-Kh/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh))/(1-1/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh));
                H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                % Continued heating after change:
            elseif modes(q,1)==6
                H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
            end
        end
        rho_PCM=(rho_PCM_l+rho_PCM_s)/2;
        cp_PCM=(cp_PCM_sensible+dHdT*hf)*H(q,1)+(cp_PCM_sensible+dHdT*hf)*(1-H(q,1));
        alpha_PCM=k_PCM/(rho_PCM*cp_PCM);
        dTdt(q,1)=((2*alpha_PCM)/(dx^2))*(T_PCM1(q-1,1)-T_PCM1(q,1));
    
        
    % All the enterior nodes:
    else
        % Maybe cooling of the interior node:
        if T_PCM1(q,1)>T_PCM1(q-1,1)
            
            % If true, the temperature gradient to the previous node is
            % higher than to the next node and it is therefore dominant for
            % the direction of the heat transfer, meaning that cooling is
            % happening:
            if T_PCM1(q,1)<T_PCM1(q+1,1) && abs(T_PCM1(q+1,1)-T_PCM1(q,1))<abs(T_PCM1(q-1,1)-T_PCM1(q,1))
                % Continued cooling:
                if modes(q,1)==1
                    H(q,1)=Ac+(Kc-Ac)/(exp(-Bc*(T_PCM1(q,1)-Mc))+1)^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                    % Change from heating to cooling:
                elseif modes(q,1)==5 || modes(q,1)==6
                    Khtoc(q,1)=(lf_i(q,1)-Ac)*(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc)+Ac;
                    H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                    % Continued cooling after change:
                elseif modes(q,1)==2
                    H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                end
                % The temperature gradient to the next node is higher than to
                % the previous node and it is therefore dominant for
                % the direction of the heat transfer, meaning that heating is
                % happening:
            else
                if modes(q,1)==5
                    H(q,1)=Ah+(Kh-Ah)/(exp(-Bh*(T_PCM1(q,1)-Mh))+1)^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                    % Change from cooling to heating:
                elseif modes(q,1)==1 || modes(q,1)==2
                    Actoh(q,1)=(lf_i(q,1)-Kh/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh))/(1-1/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh));
                    H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                    % Continued heating after change:
                elseif modes(q,1)==6
                    H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                end
            end
            rho_PCM=(rho_PCM_l+rho_PCM_s)/2;
            cp_PCM=(cp_PCM_sensible+dHdT*hf)*H(q,1)+(cp_PCM_sensible+dHdT*hf)*(1-H(q,1));
            alpha_PCM=k_PCM/(rho_PCM*cp_PCM);
            dTdt(q,1)=((2*alpha_PCM)/(dx^2))*(T_PCM1(q+1,1)+T_PCM1(q-1,1)-2*T_PCM1(q,1));
            % Maybe heating of the interior node:
        elseif T_PCM1(q,1)<=T_PCM1(q-1,1)
            % If true, the temperature gradient to the previous node is
            % higher than to the next node and it is therefore dominant for
            % the direction of the heat transfer, meaning that heating is
            % happening:
            if T_PCM1(q,1)>T_PCM1(q+1,1) && abs(T_PCM1(q+1,1)-T_PCM1(q,1))<abs(T_PCM1(q-1,1)-T_PCM1(q,1))
                % Continued heating:
                if modes(q,1)==5
                    H(q,1)=Ah+(Kh-Ah)/(exp(-Bh*(T_PCM1(q,1)-Mh))+1)^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                    % Change from cooling to heating:
                elseif modes(q,1)==1 || modes(q,1)==2
                    Actoh(q,1)=(lf_i(q,1)-Kh/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh))/(1-1/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh));
                    H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                    % Continued heating after change:
                elseif modes(q,1)==6
                    H(q,1)=Actoh(q,1)+(Kh-Actoh(q,1))/(1+exp(-Bh*(T_PCM1(q,1)-Mh)))^(1/vh);
                    dHdT=-(Bh*exp(Bh*(Mh - T_PCM1(q,1)))*(Ah - Kh))/(vh*(exp(Bh*(Mh - T_PCM1(q,1))) + 1)^(1/vh + 1));
                end
                % The temperature gradient to the next node is higher than to
                % previous node and it is therefore dominant for the direction
                % of the heat transfer, meaning that cooling is happening:
            else
                % Continued cooling:
                if modes(q,1)==1
                    H(q,1)=Ac+(Kc-Ac)/(exp(-Bc*(T_PCM1(q,1)-Mc))+1)^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                    % Change from heating to cooling:
                elseif modes(q,1)==5 || modes(q,1)==6
                    Khtoc(q,1)=(lf_i(q,1)-Ac)*(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc)+Ac;
                    H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                    % Continued cooling after change:
                elseif modes(q,1)==2
                    H(q,1)=Ac+(Khtoc(q,1)-Ac)/(1+exp(-Bc*(T_PCM1(q,1)-Mc)))^(1/vc);
                    dHdT=-(Bc*exp(Bc*(Mc - T_PCM1(q,1)))*(Ac - Kc))/(vc*(exp(Bc*(Mc - T_PCM1(q,1))) + 1)^(1/vc + 1));
                end
            end
            rho_PCM=(rho_PCM_l+rho_PCM_s)/2;
            cp_PCM=(cp_PCM_sensible+dHdT*hf)*H(q,1)+(cp_PCM_sensible+dHdT*hf)*(1-H(q,1));
            alpha_PCM=k_PCM/(rho_PCM*cp_PCM);
            dTdt(q,1)=((2*alpha_PCM)/(dx^2))*(T_PCM1(q+1,1)+T_PCM1(q-1,1)-2*T_PCM1(q,1));
        end
    end
end
end