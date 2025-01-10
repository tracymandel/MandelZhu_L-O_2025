function mainFunctionLinearWaveZhu(length,width,thick,rhou_structure,...
    E,rhou_water,Depth_water,Height_wave,T_wave,resultsFolder,validationCase,dataFileName,...
    num_T_period,n,Dt)
% Runs Longhuan's code to simulate blade motion
% ========================================================================
% copyright Longhuan Zhu 2020
% 
% For model details, see
% Zhu, L., Zou, Q.-P., Huguenard, K., and Fredriksson, D. W. (2020). Mechanisms for the Asymmetric
% Motion of Submerged Aquatic Vegetation in Waves: A Consistent-Mass Cable Model. Journal of
% Geophysical Research: Oceans, 125(2):e2019JC015517.
% 
% Primary code repository for Zhu et al. (2020) located at https://github.com/lzhu7/flexibleBladeDynamicsWavesCurrents
% ========================================================================

%% Inputs and model configurations


%constant parameters
g=9.81;             %[m/s^2] gravatational acceleration
tol=1e-8;           %Tolerance       

        
%hydrodynamic coefficients
Ca=1;%               %added mass coeff.
Cdn=nan;%            %drag coeff. Calculated in Eq_Govern.m with Luhar's formula
Cdt=nan;%            %friction coeff. Calculated in Eq_Govern.m with Luhar's formula
        
%control options
wave_option=1;      %0: measured water partical velocity
                    %1: linear wave theory
propagation=1;      %0: left (not available here for some reason)
                    %1: right
crest=0;            %0: crest (does not matter)
                    %1: trough
% settings to save data
data0 = 0;   % 1: save all variables
             % 0: save selected variables below
data1 = 'xPos'; % horizontal position of blade segments in one wave period [segment #, time step]
data2 = 'zPos'; % vertical position of blade segments [segment #, time step]
                % example: xPos(5, 2) is the horizontal position at the
                % fifth segment at time t=2*Dt.
data3 = 'tPos'; % time serious for one period for the position
data4 = 'Umax';
    
%% ======================== input is done here=============================
%% The following is the calculation
tic
%% initial calculation
%structure
Area_cross=width*thick;             %cross section area
I=width*thick^3/12;                 %second momentum for rectangular cross
m_unit=rhou_structure*Area_cross;   %mass for unit length
%         ma_unit=Ca*rhou_water*pi/4*width^2; %added mass for unit length

% this is one modification from the Zhu 2020 dataset
ma_unit=Ca*rhou_water*pi/4*width^2; %width*thick;%added mass for unit length
% ma_unit=Ca*rhou_water*width*thick; %added mass for unit length

w0_unit=(rhou_structure-rhou_water)*Area_cross*g;                   %wet weight for unit length

% keyboard 

%waves        
Omega_wave=2*pi/T_wave;             %angular frequency        
k_wave_number=wave_num(Omega_wave,Depth_water);%wave number (linear wave theory)
%check shallow water    
alpha_kh=k_wave_number*Depth_water;

%linear wave function
if wave_option
    Ux_wave=@(x,z,t) Height_wave/2*Omega_wave*cosh(k_wave_number.*z)/sinh(k_wave_number*Depth_water).*cos(k_wave_number.*x-Omega_wave.*t+crest*pi);
    Uz_wave=@(x,z,t) Height_wave/2*Omega_wave*sinh(k_wave_number.*z)/sinh(k_wave_number*Depth_water).*sin(k_wave_number.*x-Omega_wave.*t+crest*pi);
    Ux_max=@(x,z,t) Height_wave/2*Omega_wave*cosh(k_wave_number.*z)/sinh(k_wave_number*Depth_water);
end
%0.05;%
    %control
    T_total=num_T_period*T_wave;    %total time
    num_time=int32(T_total/Dt);     % make sure it's an integer...
    t=0:Dt:T_total;                 %time

    

%% Main Procedure

%initialize
    Ds=length/n;                    %segment length
    m(:,1)=zeros(n,1)+m_unit;       %segment mass
    ma_seg=ma_unit;                   %segment added mass
    w0(:,1)=zeros(n,1)+w0_unit;%*Ds;  %segment weight

    %keyboard
    %unknow vector: X=[T,Q,u,v,Omega,phi]
    T=zeros(n+1,num_time+1);        %tension
    Q=zeros(n+1,num_time+1);        %shear
    u=zeros(n+1,num_time+1);        %tangent velocity
    v=zeros(n+1,num_time+1);        %normal velocity
    Omega=zeros(n+1,num_time+1);    %curvature
    phi=zeros(n+1,num_time+1);      %angular to horizontal

%time=0;
    T(n+1,1)=0;                     %-rhou_water*g*(Depth_water-length)*Area_cross;            %0;
    for i=n+1:-1:2
        T(i-1,1)=T(i,1)+w0(i-1)*Ds;%+rhou_water*g*(Depth_water-(i-2)*Ds)*Area_cross; %pretension
    end
    phi(:,1)=zeros(n+1,1)+pi/2;     %vertical pi/2

       
    %solution definition X=[T,Q,u,v,Omega,phi]
    X=zeros(6*(n+1),1);
    for k=1:n+1                     %kth segment
        X(6*(k-1)+1)=T(k,1);
        X(6*(k-1)+2)=Q(k,1);
        X(6*(k-1)+3)=u(k,1);
        X(6*(k-1)+4)=v(k,1);
        X(6*(k-1)+5)=Omega(k,1);
        X(6*(k-1)+6)=phi(k,1);
    end
    
    %postion    
        DisX=zeros(n+1,num_time+1);     %x displacement for each point
        DisZ=zeros(n+1,num_time+1);     %z displacement for each point
        DisZ(1,1)=0;
        for k=2:n+1
                DisZ(k,1)=DisZ(k-1,1)+Ds;
        end
        
    %structure length--check
        Length_curvature=zeros(num_time+1,1);
        Length_curvature(1)=length;

    u_flow=zeros(n+1,num_time+1);       %water partical velocities: tangent direction
    v_flow=zeros(n+1,num_time+1);       %normal direction
%     global Umax
    Umax=zeros(n+1,num_time+1);
    %select wave theory
        if wave_option<1
            [u_flow(:,1),v_flow(:,1),Umax(:,1)]=wave_nearVege(DisX(:,1),DisZ(:,1),t(1),k_wave_number,Omega_wave);
%             [u_flow(:,1),~]=wave_nearVege(DisX(:,1),DisZ(:,1),t(1),k_wave_number,Omega_wave);
%             v_flow(:,1)=ones(n+1,1)*0.01;
        elseif wave_option<2
            u_flow(:,1)=Ux_wave(DisX(:,1),DisZ(:,1),t(1));
            v_flow(:,1)=Uz_wave(DisX(:,1),DisZ(:,1),t(1));
            Umax(:,1)=Ux_max(DisX(:,1),DisZ(:,1),t(1)); 
        end
        
CDN=zeros(n,num_time+1);
KC=zeros(n,num_time+1);
DragT=zeros(n,num_time+1);
DragN=zeros(n,num_time+1);
MaT=zeros(n,num_time+1);
MaN=zeros(n,num_time+1);
BT=zeros(n,num_time+1);
BN=zeros(n,num_time+1);
TT=zeros(n,num_time+1);
TN=zeros(n,num_time+1);
VBT=zeros(n,num_time+1);
WT=zeros(n,num_time+1);
AMN=zeros(n,num_time+1);
VBN=zeros(n,num_time+1);
WN=zeros(n,num_time+1);

%% time after 0

    j = 1;
    tolX = 0.01;
    tolZ = 0.01;
    xErr = 1;
    zErr = 1;
    t0=toc;

    while j < num_time+1 & (xErr>tolX || zErr>tolZ)

    % for j=2:num_time+1

        j = j+1;
        X0=X;

        %New position
            DisX(:,j)=DisX(:,j-1);
            DisZ(:,j)=DisZ(:,j-1);
            u_flow(:,j)=u_flow(:,j-1);
            v_flow(:,j)=v_flow(:,j-1);
            Umax(:,j)=Umax(:,j-1);

        %Substitude parameters into Governing Equations
            Eq_Govern_Done=@(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m,w0,Umax) Eq_Govern(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m,ma_seg,Dt,Ds,w0,rhou_water,Cdt,Cdn,width,Area_cross,E,I,thick,T_wave,Umax);

        [Y,KC(:,j),~,DragT(:,j),DragN(:,j),MaT(:,j),MaN(:,j),BT(:,j),BN(:,j),TT(:,j),TN(:,j),VBT(:,j),WT(:,j),AMN(:,j),VBN(:,j),WN(:,j)]=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));

        norm_Y=norm(Y);

        while norm_Y>tol
                
            %Form Jacobi
                J=J_value(@Jacob,Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,Umax(:,j));
            
            X=X-J\Y;

            %New position
                [DisX(:,j),DisZ(:,j),Length_curvature(j)]=Position(DisX(:,j-1),DisZ(:,j-1),X0,X,Dt,Ds,n);

            %select wave theory
                if wave_option<1
                    [u_flow(:,j),v_flow(:,j),Umax(:,j)]=wave_nearVege(DisX(:,j),DisZ(:,j),t(j),k_wave_number,Omega_wave);
                elseif wave_option<2  % LWT: wave_option = 1
                    u_flow(:,j)=Ux_wave(DisX(:,j),DisZ(:,j),t(j));
                    v_flow(:,j)=Uz_wave(DisX(:,j),DisZ(:,j),t(j));   
                    Umax(:,j)=Ux_max(DisX(:,j),DisZ(:,j),t(j));      
                end
                            
            [Y,KC(:,j),~,DragT(:,j),DragN(:,j),MaT(:,j),MaN(:,j),BT(:,j),BN(:,j),TT(:,j),TN(:,j),VBT(:,j),WT(:,j),AMN(:,j),VBN(:,j),WN(:,j)]=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));
               
            norm_Y=norm(Y);
        end

        %output results for jth time
            for k=1:n+1
                T(k,j)=X(6*(k-1)+1);
                Q(k,j)=X(6*(k-1)+2);
                u(k,j)=X(6*(k-1)+3);
                v(k,j)=X(6*(k-1)+4);
                Omega(k,j)=X(6*(k-1)+5);
                phi(k,j)=X(6*(k-1)+6);
            end
        %process report    
            % fprintf('T_wave= %4.1fs; Time is %8.4fs\n',T_wave,t(j)) 
        t1=toc;
        
        if rem(t(j),5)==0
            fprintf(['t=',num2str(t(j),'%5.2f'),'s; finish in ', num2str((t(end)-t(j))/t(j)*(t1-t0)), 's\n']) 
        end
    end
    
    
    %% check calculations
    %structural length check
        AveSum=mean(Length_curvature);
        AveSum2=sqrt(norm(Length_curvature)^2/(double(num_time)+1));
        if abs(AveSum2-length)/length<1e-3
        else
            warning('the blade length is heavily strectched due to small axial modulus')
        end

    % % figure for reconfiguration
    %     figure(1); clf; hold on
        T_num=T_wave/Dt;
    %     DT=round(T_num/64);
    % 
    %     order_T=1;        
    %     plot(DisX(n+1,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num),....
    %         DisZ(n+1,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num),'g-','LineWidth',2.5)
    % 
    %     for j=num_time+1-order_T*T_num:DT:num_time-(order_T-1)*T_num%1:DT:num_time+1;%1:DT:T_num+1;%
    %         plot(DisX(:,j),DisZ(:,j),'color',[39, 174, 96]/255,'Linewidth',0.1)
    %         hold on
    %     end
    % 
    %         xlabel('x(m)')
    %         ylabel('z(m)')
    %         axis equal;
    %         axis([-1 1 0 1]*length)
    %     grid on
    %     saveas(gcf,[resultsFolder,'\bladeMotion.png']);
        
    % % displacement of upper end (check steady state)
    %     figure(2); clf; hold on
    %     plot(t,DisX(n+1,:),t,DisZ(n+1,:),'LineWidth',2.5)
    %     legend('XDis-Upper end','ZDis-Upper end','location','best')
    %     Title=['Displacement:    T_w_a_v_e=',num2str(T_wave),'s  a_w_a_v_e=',num2str(Height_wave/2),'m; kh=',num2str(alpha_kh)];
    %     title(Title)
    %     xlabel('t(s)')
    %     ylabel('Displacement(m)')
    %     grid on
    %     saveas(gcf,[resultsFolder,'\tipDisplacement.png']);
        
%% save data
% save Data
mkdir(resultsFolder);

if data0
    save([resultsFolder,'/',validationCase,'/',dataFileName,'.mat']);
else
    order_T = 1;
    xPos = DisX(:,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num);
    zPos = DisZ(:,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num);
    tPos = (1:T_num)'*Dt;
    save([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'], data1, data2, data3, data4);
end


%%
toc          %ending time
% disp('===============END=================')2
