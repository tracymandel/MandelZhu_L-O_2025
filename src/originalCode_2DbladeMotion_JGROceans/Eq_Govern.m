function [Y,Diff,Cdn,DragT,DragN,MaT,MaN,BT,BN,TT,TN,VBT,WT,AMN,VBN,WN]...
    =Eq_Govern(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m_unit,ma_unit,Dt,Ds,w0_unit,rhou_water,Cdt,Cdn,diameter,Area_cross,E,I,thick,T_wave,Umax)

%% ========================================================================
% copyright Longhuan Zhu 2020
% 
% For model details, see
% Zhu, L., Zou, Q.-P., Huguenard, K., and Fredriksson, D. W. (2020). Mechanisms for the Asymmetric
% Motion of Submerged Aquatic Vegetation in Waves: A Consistent-Mass Cable Model. Journal of
% Geophysical Research: Oceans, 125(2):e2019JC015517.
% 
% Primary code repository for Zhu et al. (2020) located at https://github.com/lzhu7/flexibleBladeDynamicsWavesCurrents
% ========================================================================


%%
% global Umax
KC=T_wave/diameter*Umax;
% KC=T_wave/diameter*0.206*1.234072045640771; %Luhar

%% Zhu et al. 2020 dataset used constant coefficients (see below)
Cdn=max(1.95,10*KC^(-1/3)); 
% Cdn = 1.95;

Re=Umax*diameter/1e-6;
Cdt = 0.074 * Re^(-1/5);
% Cdt = 0.1;

%%
if KC<20
    ma_unit_1=(1+(0.35*KC^(2/3)))*ma_unit;
else
    ma_unit_1=(1+(0.15*KC^(2/3)))*ma_unit;
end
ma_unit_2=(1+((KC-18)^2)/49)*ma_unit;
ma_unit=min(ma_unit_1,ma_unit_2);



%%
DragT=-1/2*rhou_water*Cdt*2*(diameter+thick)/4*( ...%drag tangent
        abs(X(9)-(Ux_wave_new(2)*cos(X(12))+Uz_wave_new(2)*sin(X(12))))*(X(9)-(Ux_wave_new(2)*cos(X(12))+Uz_wave_new(2)*sin(X(12)))) ...
        +abs(X(3)-(Ux_wave_new(1)*cos(X(6))+Uz_wave_new(1)*sin(X(6))))*(X(3)-(Ux_wave_new(1)*cos(X(6))+Uz_wave_new(1)*sin(X(6)))) ...
        +abs(X0(9)-(Ux_wave_old(2)*cos(X0(12))+Uz_wave_old(2)*sin(X0(12))))*(X0(9)-(Ux_wave_old(2)*cos(X0(12))+Uz_wave_old(2)*sin(X0(12)))) ...
        +abs(X0(3)-(Ux_wave_old(1)*cos(X0(6))+Uz_wave_old(1)*sin(X0(6))))*(X0(3)-(Ux_wave_old(1)*cos(X0(6))+Uz_wave_old(1)*sin(X0(6)))) ...
    );
MaT=m_unit/2/Dt*(X(9)+X(3)-X0(9)-X0(3))-m_unit/8/Dt*(X(10)+X(4)+X0(10)+X0(4))*(X(12)+X(6)-X0(12)-X0(6));%%%+ or -????
TT=-1/4*(X(11)*X(8)+X0(11)*X0(8)+X(5)*X(2)+X0(5)*X0(2));
BT=+1/2/Ds*(X(7)+X0(7)-X(1)-X0(1));
% VBT=rhou_water*diameter*thick*(...
%         1/2/Dt*(Ux_wave_new(2)*cos(X(12))+Ux_wave_new(1)*cos(X(6))-Ux_wave_old(2)*cos(X0(12))-Ux_wave_old(1)*cos(X0(6)))...
%         +1/2/Dt*(Uz_wave_new(2)*sin(X(12))+Uz_wave_new(1)*sin(X(6))-Uz_wave_old(2)*sin(X0(12))-Uz_wave_old(1)*sin(X0(6))) ...
%     ); %Virtual Buoyancy 0;%+
VBT=rhou_water*diameter*thick*(...
        1/2/Dt*(Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
        *1/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)))...
        +1/2/Dt*(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))...
        *1/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
    ); %Virtual Buoyancy 0;%+
VBTCable=0;
WT=-w0_unit/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)));

Y(1,1)=-MaT+BT+TT+DragT+VBT+WT;
% % Y(1,1)=-MaT+BT+TT+DragT+VBTCable+WT;

%%
DragN=-1/2*rhou_water*Cdn*diameter/4*( ...%drag normal
        abs(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12))))*(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12)))) ...
        +abs(X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6))))*(X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6)))) ...
        +abs(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12))))*(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12)))) ...
        +abs(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6))))*(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6)))) ...
    );

MaN=m_unit/2/Dt*(X(10)+X(4)-X0(10)-X0(4))+m_unit/8/Dt*(X(9)+X(3)+X0(9)+X0(3))*(X(12)+X(6)-X0(12)-X0(6));
BN=+1/2/Ds*(X(8)+X0(8)-X(2)-X0(2));
TN=+1/4*(X(11)*X(7)+X0(11)*X0(7)+X(5)*X(1)+X0(5)*X0(1));
AMN=-ma_unit/2/Dt*(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12)))...
    +X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6))) ...
    -(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12))))...
    -(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6))))); ...%added mass
    % added mass should be the derivative of the relative velocity, affect
    % a lot
AMNLuhar=-ma_unit/2/Dt*(X(10)+X(4)-X0(10)-X0(4))....
    -ma_unit/8/Dt*(X(9)+X(3)+X0(9)+X0(3))*(X(12)+X(6)-X0(12)-X0(6))....  %%DV/Dt=pV/pt+wxV
    -ma_unit/8/Dt*((Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
    *(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
    -(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))....
    *(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6))));...%%Luhar fomula
Diff=AMN-AMNLuhar;
VBNCable=-rhou_water*diameter*thick*( ...
        1/2/Dt*(Ux_wave_new(2)*sin(X(12))+Ux_wave_new(1)*sin(X(6))...
        -Ux_wave_old(2)*sin(X0(12))-Ux_wave_old(1)*sin(X0(6)))...
        -1/2/Dt*(Uz_wave_new(2)*cos(X(12))+Uz_wave_new(1)*cos(X(6))...
        -Uz_wave_old(2)*cos(X0(12))-Uz_wave_old(1)*cos(X0(6)))...
    );% ...%Virtual Buoyancy
VBN=-rhou_water*diameter*thick*( ...
    -1/2/Dt*(Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
    *1/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
    +1/2/Dt*(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))...
    *1/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)))...
    );% ...%Virtual Buoyancy
WN=-w0_unit/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)));
Y(2,1)=-MaN+BN+TN+DragN+AMN+VBN+WN;
% % Y(2,1)=-MaN+BN+TN+DragN+AMNLuhar+VBN+WN;
% % Y(2,1)=-MaN+BN+TN+DragN+AMN+VBNCable+WN;

%%
Y(3,1)=-1/2/Dt*(X(7)+X(1)-X0(7)-X0(1))+E*Area_cross/2/Ds*(X(9)+X0(9)-X(3)-X0(3))-E*Area_cross/4*(X(11)*X(10)+X0(11)*X0(10)+X(5)*X(4)+X0(5)*X0(4));
%%
Y(4,1)=-1/2/Dt*(X(12)+X(6)-X0(12)-X0(6))+1/2/Ds*(X(10)+X0(10)-X(4)-X0(4))+1/4*(X(11)*X(9)+X0(11)*X0(9)+X(5)*X(3)+X0(5)*X0(3));
%%
Y(5,1)=E*I/Ds*(X(11)-X(5))+1/2*(X(8)+X(2));
%%
Y(6,1)=1/Ds*(X(12)-X(6))-1/2*(X(11)+X(5));