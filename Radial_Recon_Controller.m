clear;clc;close all;

[file,this_path] = uigetfile('.dat','Select Radial UTE File');
cd(this_path);
%% Read Raw Data
disp('Reading Data')
%This twix reading function is freely available online, but I made an edit
%to handle a "bug" in my sequence code. Briefly, the loop counting
%variables are stored as shorts, so the line counter gets screwed up if you
%collect too many projections. This function reads the data correctly in
%spite of this.
twix_obj = mapVBVD_pjn(file,'ignoreSeg');

if length(twix_obj)>1
    fid = double(twix_obj{2}.image());
    twix_obj = twix_obj{2};
else
    fid = double(twix_obj.image());
end

fid = permute(fid,[1,3,2]);

disp('Reading Data Complete')

%% Some Trajectory Parameters - These should match the default protocol
ImSize = 350; %Desired Recon size.
FOV = 350;
ADC_Dur = 960; %ADC in us
RUT = 140; %Ramp Up Time in us
Resolution = [FOV/ImSize FOV/ImSize FOV/ImSize]/1000;
gamma = 42.6;
MaxGrad = 13.980741193; %Max Gradient in mT/m from simulation ERR file
TR = 3.3;
%% Generate Trajectories - Kind of a Brute force method
NPro = size(fid,2);
Pts = size(fid,1);
Dur = ADC_Dur-RUT;

Grad = linspace(0,1,RUT);
Grad((RUT):ADC_Dur) = 1;

Grad = Grad * MaxGrad;

ADC_Dur = ADC_Dur*1e-6;
RUT = RUT*1e-6;
Dw = ADC_Dur/Pts;

Arm_untimed = cumtrapz(Grad);
Grad_Time = 0:1e-6:(ADC_Dur-1e-6);

Time = 0:Dw:(ADC_Dur-Dw);
Arm = interp1(Grad_Time,Arm_untimed,Time); %Now we are in mT*s/m

Arm = Arm*gamma/1000;

%% Rotate Trajectories - This assumes that "Golden Means" is selected

phi1 = 0.46557123;
phi2 = 0.6823278;
gs = 1;
gr = 1;
gp = 1;

r = zeros(1,NPro);
p = zeros(1,NPro);
s = zeros(1,NPro);
%Rotation code 
for i = 0:(NPro-1)
    kz = (i*phi1-floor(i*phi1))*2-1;
    ts = kz*gs;
    alpha = (i*phi2-floor(i*phi2))*2*pi;
    tr = sqrt(1-kz*kz)*cos(alpha)*gr;
    tp = sqrt(1-kz*kz)*sin(alpha)*gp;
    r(i+1) = tr;
    p(i+1) = tp;
    s(i+1) = -ts;
end

trajx_rot = zeros(length(Arm),NPro);
trajy_rot = trajx_rot;
trajz_rot = trajx_rot;
for i = 1:NPro
    trajx_rot(:,i) = r(i)*Arm';
    trajy_rot(:,i) = p(i)*Arm';
    trajz_rot(:,i) = s(i)*Arm';
end

traj = cat(3,trajx_rot,trajy_rot,trajz_rot);
traj = permute(traj,[3,1,2]);

kFOV_desired = 1./(Resolution);
kMax_desired = kFOV_desired/2;
max_k = max(kMax_desired); %Here, we are in 1/m
traj = traj/max_k/2;

%% Perform a trajectory gradient delay correction - Optional
% This parameter may need some manual editing - Gradient delay correction
% in units of dwell times (i.e. a value of 1 shifts the gradient shape in time by 1
% dwell time)
t_delay = 0;

traj = traj_delay_correction(traj,t_delay);

%% Clear some variables to make sure there is space in memory
clear trajx_rot;
clear trajy_rot;
clear trajz_rot;
clear Grad;
clear Grad_Time;
clear Gate_Raw;
clear Arm;
clear Arm_untimed;

%% Remove Descent to steady state - Optional
fid(:,1:1000,:) = [];
traj(:,:,1:1000) = [];

%% Check trajectory delay - Optional
% Use a small number of projections to reconstruct a quick image
nIter = 10;
n_sub = 5000;

tmpfid = fid(:,1:n_sub,:);
tmptraj = traj(:,:,1:n_sub);

[tmpfid,tmptraj] = ImTools.change_mat_size(tmpfid,tmptraj,ImSize/2,ImSize);

AllDCF = Recon.get_DCF(tmptraj,ImSize/2,nIter);
% Reconstruct
ncoils = size(fid,3);
Check_Image = zeros(ImSize/2,ImSize/2,ImSize/2,ncoils);
for i = 1:ncoils
    Check_Image(:,:,:,i) = Recon.mem_eff_recon(ImSize/2,tmpfid(:,:,i),tmptraj,AllDCF,i,ncoils);%base_noncart_recon(ImSize,fid_I,rtrajI);
end
%toc
Fin_Image = SOS_Coil_Combine(Check_Image);
imslice(Fin_Image/max(Fin_Image(:)));

clear tmpfid;
clear tmptraj;
%% Reconstruct all Projections
nIter = 10;
AllDCF = Recon.get_DCF(traj,ImSize,nIter);
% Reconstruct
% could also implement a coil combination step to make this faster
ncoils = size(fid,3);
Image = zeros(ImSize,ImSize,ImSize,ncoils);
for i = 1:ncoils
    Image(:,:,:,i) = Recon.mem_eff_recon(ImSize,fid(:,:,i),traj,AllDCF,i,ncoils);%base_noncart_recon(ImSize,fid_I,rtrajI);
end
%toc
Fin_Image = SOS_Coil_Combine(Image);
Fin_Image = Fin_Image/max(Fin_Image(:))*100;

% Can write out this image as nifti or DICOM, or display using your preferred method.
% This imslice tool is one I like:
imslice(Fin_Image);