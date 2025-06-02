%% ========================================================================
%  Script Name:    EndoCtl_v1.m
%  Description:    Magnetic endoscope robot controller
%
%  Author:         Ruoyu Xu
%  Date:           2025-03-11
%  Version:        1.0
%
%  Example:        na
%      
%  Notes:          cal_inv_dyn_ipm.m(or cal_ipm_opt.m) is required,
%                  simulation with vrep akukamag_ipm.ttt
%
%  Dependencies:    Robotics System Toolbox
%
%  Copyright:      (c) Ruoyu Xu. All rights reserved.
% ========================================================================
clear;clc;close all;
%% vrep ini
[clientID, vrep ] = StartVrep(19999);
robot_names={'./LBRiiwa7R800/joint','./LBRiiwa7R800/link2_resp/joint','./LBRiiwa7R800/link3_resp/joint',...
    './LBRiiwa7R800/link4_resp/joint','./LBRiiwa7R800/link5_resp/joint','./LBRiiwa7R800/link6_resp/joint',...
    './LBRiiwa7R800/link7_resp/joint'};
robot_handles = get_handle_Joint(vrep,clientID,robot_names);
get_joint_positions_ini(clientID,vrep,robot_handles);
[~ ,baseframe] =vrep.simxGetObjectHandle (clientID,'./LBRiiwa7R800',vrep.simx_opmode_blocking);
[~,ipm_handles] = vrep.simxGetObjectHandle(clientID,'IPM',vrep.simx_opmode_blocking);
%% user para
robot = loadrobot('kukaIiwa7','DataFormat','column','Gravity',[0 0 -9.81]);
% set joint position limitation
q_max=[170;120;170;120;170;120;175]*pi/180;
dq_max=0.8*[85;85;100;75;130;130;135]*pi/180;
m_a=680;
m_c=10;
% Initial Pose
q_ini = [0; 20; 0; -90; 0; 70; 0] * pi / 180;
transform = getTransform(robot,q_ini,'iiwa_link_ee_kuka');
ik = inverseKinematics('RigidBodyTree',robot);
ik.SolverParameters.AllowRandomRestart = false;
weights = 0.5*[1 1 1 1 1 1];
initialGuess = q_ini;
taskFinal=transform;
taskFinal(1:3,4)=transform(1:3,4)+[0;0;0]; % adjust end position
q_ini = ik('iiwa_link_ee_kuka',taskFinal,weights,initialGuess); % ini joint position after adjustment
x_epm_ini=cal_xepm(robot,q_ini,m_a); % ini epm position
m_epm_hat=x_epm_ini(4:6)/norm(x_epm_ini(4:6));
x_ipm_ini = cal_ipm_opt(x_epm_ini,[x_epm_ini(1:2);x_epm_ini(3)-0.25;-m_c*m_epm_hat(1:2);m_c*m_epm_hat(3)]); % ini ipm position
x_ipm_c=x_ipm_ini; % the current ipm state
q_c=q_ini;
rotn0=taskFinal(1:3,3); % ini z vector of epm
%% vrep move to ini
set_robot_positions(vrep, clientID, robot_handles, q_ini)
vrep.simxSetObjectPosition(clientID,ipm_handles,baseframe,x_ipm_c(1:3)',vrep.simx_opmode_oneshot);
pause(1);

%% main loop
Ts = 0.005; % sim period
loop_rate = rateControl(1/Ts); % define timer loop
reset(loop_rate); % reset timer loop
cnt=0; % counter
tp=30; % sim time
t0=tic; % calculate initial time
dt=toc(t0);
while (dt<tp)
    dt=cnt*Ts;
    cnt=cnt+1;

    % desired position
    rr=0.2;
    ang=pi*dt/tp;
    dx=rr.*sin(ang);
    dy=rr.*(1-cos(ang));
    dz=0;
    drx=0;
    dry=30/57.3;
    drz=ang;
    p_ipm_d=x_ipm_ini(1:3)+[dx;dy;dz];
    h_ipm_d=eul2rotm([drz dry drx])*x_ipm_ini(4:6)/norm(x_ipm_ini(4:6));
    
    % feedback
    q_f=q_c; % robot joint position
    x_epm=cal_xepm(robot,q_f,m_a); % epm state
    m_epm_hat=x_epm(4:6)/norm(x_epm(4:6));
    x_ipm=x_ipm_c; % ipm state
    p_ipm=x_ipm(1:3); % ipm pos
    h_ipm=x_ipm(4:6)/norm(x_ipm(4:6)); % ipm heading
    p=x_ipm(1:3)-x_epm(1:3); % ipm epm relative pos
    m_epm=x_epm(4:6); % epm moment
    m_ipm=x_ipm(4:6); % ipm moment
    [PM_F,PM_tau] = cal_ipm_force(p,m_epm,m_ipm);
    ft=norm(PM_F(1:2));
    % Jacobian
    jacob_m = geometricJacobian(robot,q_f,'iiwa_link_ee_kuka'); 
    Jac_g=[jacob_m(4:6,:);jacob_m(1:3,:)]; % robot jacobian 6x7
    Jac_A=[eye(3) zeros(3); zeros(3) skew(m_epm_hat)']*Jac_g; % 6x7
    Jac_F = cal_ipm_jacob_ana(p,m_epm,m_ipm); % magnetic force jacobian 6x9
    J_FA=Jac_F(:,1:6)*[-eye(3) zeros(3); zeros(3) eye(3)]*Jac_A; % 6x9
%     J_FA=Jac_F*[-eye(3) zeros(3) zeros(3); zeros(3) eye(3) zeros(3); zeros(3,9)]*[Jac_A;zeros(3,7)]; % 6x9
    % controller
    err_p=p_ipm_d-p_ipm; % pos error
    err_h=cross(h_ipm,h_ipm_d); % heading error
    Km=diag([5 5 5 1 1 1]); % control parameter
    dq=J_FA'/(J_FA*J_FA'+0.00001*eye(6))*Km*[err_p;err_h]; % joint vel

    dq_rec1=LimitJointState(dq,dq_max); % limit joint vel

    dx_r=Jac_g*dq_rec1; % limit z
    dx_r(3)=0;
    dq_rec2=pinv(Jac_g)*dx_r;

    q_d=q_f+Ts*dq_rec2;

%     xend=getTransform(robot,q_d,'iiwa_link_ee_kuka');
%     an_s=acosd(rotn0'*xend(1:3,3));
%     if an_s>=30
%         dx_r2=Jac_g*dq_rec2; % limit z
%         dx_r2(3)=0;
%         dx_r2(5:6)=[0;0];
%         dq_rec2=pinv(Jac_g)*dx_r2;
%     end
%     q_d=q_f+Ts*dq_rec2;

    % sim
    q_c=q_d; % current joint state in simulation
    x_epm_c=cal_xepm(robot,q_c,m_a); % current epm state in simulation
%     x_ipm_c = cal_inv_dyn_ipm(x_epm_c,x_ipm); 

    x_ipm_cc = cal_ipm_opt(x_epm_c,x_ipm); % current ipm state in simulation

    % Three conditions
    flagc=1;
        % No block
    if flagc==1
        x_ipm_c=x_ipm_cc;
    elseif flagc==2
        % Enable block
        if dt>=5 && ft<=0.1
            x_ipm_c(4:6)=x_ipm_cc(4:6);       
        else
            x_ipm_c=x_ipm_cc;
        end
    elseif flagc==3
        % lose cotrol
        if dt>=5
            x_ipm_c(1:3)=x_ipm_c0(1:3)+[0.1;0.1;0];
            x_ipm_c(4:6)=x_ipm_c0(4:6);       
        else
            x_ipm_c=x_ipm_cc;
            x_ipm_c0=x_ipm_cc;
        end
    else
        x_ipm_c=x_ipm_cc;
    end


    T_robot=getTransform(robot,q_c,'iiwa_link_ee_kuka'); % end pose of manipulator


    m_ipm_c_hat=x_ipm_c(4:6)/norm(x_ipm_c(4:6));
    m_epm_c_hat=x_epm_c(4:6)/norm(x_epm_c(4:6));

    baseRipm_z=m_ipm_c_hat; % calculate ipm pose
    baseRipm_y=T_robot(1:3,2);
    baseRipm_x=cross(baseRipm_y,baseRipm_z);
    baseRipm=[baseRipm_x baseRipm_y baseRipm_z];

    set_robot_positions(vrep, clientID, robot_handles, q_d); % set robot in vrep
    vrep.simxSetObjectPosition(clientID,ipm_handles,baseframe,x_ipm_c(1:3)',vrep.simx_opmode_oneshot);
    vrep.simxSetObjectOrientation(clientID,ipm_handles,baseframe,rotm2eul(baseRipm,'xyz'),vrep.simx_opmode_oneshot);
    
    % record datau
    xipma(:,cnt)=x_ipm;
    xepma(:,cnt)=x_epm;
    qca(:,cnt)=q_f;
    xipmad(:,cnt)=[p_ipm_d;h_ipm_d];
    timea(:,cnt)=dt;
    f_maga(:,cnt)=[PM_F;PM_tau];
    Jaa(cnt)=det(J_FA*J_FA');

    %draw
    if mod(dt,0.01)==0
% %         hold off
% %         plot(xipma(1,:),xipma(2,:),'--r')
% %         hold on
% %         plot(xipmad(1,:),xipmad(2,:),'-k')
% %         drawnow
% 
%         hold off
%         plot(xipma(1,:),'--r')
%         hold on
%         plot(xipmad(1,:),'-k')
%         drawnow
    clc
    disp(['Time: ',num2str(dt), '; Error: ', num2str(1000*norm(err_p)), ' mm, ', num2str(acosd(h_ipm'*h_ipm_d)), ' deg.'])
    
    end
    

    waitfor(loop_rate); % loop period
end

% save('sim_EAST_block.mat','timea','xipmad','qca','xepma','xipma','f_maga')
%%
% x_epm_ini=cal_xepm(robot,q_ini,m_a)
%% subfun
function jacob_out = cal_ipm_jacob_ana(p,m_epm,m_ipm)
    u0= 4*pi*10^(-7);
    epm_norm = (m_epm'*m_epm)^0.5;
    ipm_norm = (m_ipm'*m_ipm)^0.5;
    p_norm = (p'*p)^0.5;
    m_epm_hat = m_epm/epm_norm;
    m_ipm_hat = m_ipm/ipm_norm;
    p_hat = p/p_norm;
    PM_D = 3*(p_hat*p_hat')-eye(3);
    PM_Z = eye(3)-5*(p_hat*p_hat');
    PM_G = eye(3)- p_hat*p_hat';
    F_p = 3*u0*epm_norm*ipm_norm/(4*pi*p_norm^5)*(m_epm_hat*m_ipm_hat'*PM_Z+m_ipm_hat*m_epm_hat'*PM_Z ...
        +m_ipm_hat'*m_epm_hat*PM_Z-5*(p_hat*p_hat')*m_epm_hat*m_ipm_hat'*PM_G ...
        -5*(p_hat*p_hat')*m_ipm_hat*m_epm_hat'*PM_G-5*m_ipm_hat'*(p_hat*p_hat')*m_epm_hat*PM_Z);
    F_ma = 3*u0*epm_norm*ipm_norm/(4*pi*p_norm^4)*(m_ipm_hat'*p_hat*eye(3)+m_ipm_hat*p_hat'+p_hat*m_ipm_hat'*PM_Z);
    F_mc = 3*u0*epm_norm*ipm_norm/(4*pi*p_norm^4)*(m_epm_hat'*p_hat*eye(3)+m_epm_hat*p_hat'+p_hat*m_epm_hat'*PM_Z);
    T_p = 3*u0*epm_norm*ipm_norm/(4*pi)*(skew(m_ipm_hat/p_norm^3)*(p_hat*m_epm_hat'*(PM_G/p_norm) ...
        +(PM_G/p_norm)*(p_hat'*m_epm_hat))+skew(PM_D*m_epm_hat)*(m_ipm_hat*p_hat'/p_norm^4));
    T_ma = u0*epm_norm*ipm_norm/(4*pi*p_norm^3)*skew(m_ipm_hat)*PM_D;
    T_mc = -u0*epm_norm*ipm_norm/(4*pi*p_norm^3)*skew(PM_D*m_epm_hat);
    jacob_out=[F_p F_ma F_mc;T_p T_ma T_mc];
end

function x_epm=cal_xepm(robot,q_f,m_a)
% cal state of epm
    transform = getTransform(robot,q_f,'iiwa_link_ee_kuka');
    pos_r=transform(1:3,4);
    rot_n=transform(1:3,3);
    p_epm=pos_r+0.09*rot_n; % tool position offset 9 cm
    m_epm=m_a*transform(1:3,1);
    x_epm=[p_epm;m_epm];
end

function  dq_out  = LimitJointState(dq_in,dqm)
    sca=abs(dq_in./dqm);
    sca_max=max(sca);
    if sca_max>1
        dq_out=dq_in/sca_max;
    else
        dq_out=dq_in;
    end
end

function S = skew(V)
	S = [
   	    0    -V(3)    V(2)
        V(3)      0    -V(1)
       -V(2)    V(1)      0
   ];
end
%% vrep subfunctions
% constructor
function [clientID, vrep ] = StartVrep(porta)
    addpath vrep_lib/;
    vrep = remApi('remoteApi');   % using the prototype file (remoteApiProto.m)
    vrep.simxFinish(-1);        % just in case, close all opened connections
    clientID = vrep.simxStart('127.0.0.1',porta,true,true,5000,5);% start the simulation
    if (clientID>-1)
        vrep.simxSetFloatingParameter(clientID,vrep.sim_floatparam_simulation_time_step,0.001,vrep.simx_opmode_oneshot_wait);
        vrep.simxSynchronous(clientID,true);
        [returnCode]=vrep.simxStartSimulation(clientID,vrep.simx_opmode_oneshot_wait);
        disp('remote API server connected successfully');
    else
        disp('failed connecting to remote API server');
        DeleteVrep(clientID, vrep); %call the destructor!
    end
end

function robot_handles = get_handle_Joint(vrep,clientID,robot_names)
% get robot handles
    len=length(robot_names);
    robot_handles=zeros(len,1);
    for i=1:len
        [~,robot_handles(i)] = vrep.simxGetObjectHandle(clientID,robot_names{i},vrep.simx_opmode_oneshot_wait);
    end  
end

function set_robot_positions(vrep, clientID, robot_handles, robot_positions)
    len=length(robot_handles);
    for i=1:len
        err = vrep.simxSetJointTargetPosition(clientID,robot_handles(i),robot_positions(i),vrep.simx_opmode_oneshot);
    end
end

function robot_positions = get_joint_positions_ini(clientID,vrep,robot_handles)
% run get_joint_position_ini for the first joint position feedback
    len=length(robot_handles);
    robot_positions = zeros(len,1);
    err = zeros(len,1);
    for j=1:len
        [err(j),robot_positions(j)]=vrep.simxGetJointPosition(clientID,robot_handles(j),vrep.simx_opmode_streaming);
        pause(0.1);
    end
end

function robot_positions = get_joint_positions(clientID,vrep,robot_handles)
    len=length(robot_handles);
    robot_positions = zeros(len,1);
    err = zeros(len,1);
    for j=1:len
        [err(j),robot_positions(j)]=vrep.simxGetJointPosition(clientID,robot_handles(j),vrep.simx_opmode_buffer);
    end
end

% destructor
function DeleteVrep(clientID, vrep)
    vrep.simxPauseSimulation(clientID,vrep.simx_opmode_oneshot_wait); % pause simulation
    vrep.simxStopSimulation(clientID,vrep.simx_opmode_oneshot_wait); % stop simulation
    vrep.simxGetPingTime(clientID);
    vrep.simxFinish(clientID);  % close the line if still open
    vrep.delete();              % call the destructor!
    disp('simulation ended');
end
