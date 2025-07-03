clear all
close all
clc

%%
addpath('MNI\')

load vol2gm.mat
load('S07_texting.nirs','-mat');

nCh = size(SD.MeasList,1)/2;

%1) d: #Samples x 2*#Channels contains the raw intensity data. The first half 
% of the columns contains the signals for all channels at one wavelength,
% whereas the second half of the columns contains the signals for all channels 
% at the second wavelength.
%2) s: #Samples x #Conditions contains the information about timing of stimulus 
% presentation for the different conditions of the task
%3) SD: struct which contains all the information about the array, i.e., 
% source and detector positions, channels composition, wavelengths of the 
% system and unit of measure. The channel composition is contained in the 
% SD.MeasList matrix, where the rows correspond to channels and the columns 
% contain: source index, detector index, auxiliary column and wavelength type.
%4) aux: #Samples x #Auxiliary channels contains additional auxiliary signals
%5) t: #Samples x 1 is the time vector

%% 1) Plot array configuration
figure('name','Array configuration')
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    % src = starting point of the line
    % det = ending point of the line
end
xlabel('x')
ylabel('y')
zlabel('z')

%% 2) Compute SD distance and plot it as histogram
distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end
figure('name','SD distance Histogram')
histogram(distCh,20)
xlabel('SD distance [mm]')
ylabel('N of channels')

%% 3) Remove bad channels

dRange = [500 10^10];
SNRrange = 0;
SD.MeasListAct = removeNoisyChannels(d,dRange,SNRrange);

% Plot array configuration with bad channels highlighted
figure('name','3D configuration')
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if SD.MeasListAct(iCh)==1
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    else
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'r')
    end
end
xlabel('x')
ylabel('y')
zlabel('z')

%% 4) pre-Processing

% Intensity plot
% Visualize fNIRS data at the two wavelengths without removed channels
dgood = d;
dgood(:,SD.MeasListAct==0) = [];

%only good channels
figure('name','Original time series')
subplot(121)
plot(t,dgood(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
subplot(122)
plot(t,dgood(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')


% Convert to optical density changes
meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

% Visualize optical density changes of non-removed channels
dodDataRemoved = dodConv;
dodDataRemoved(:,SD.MeasListAct==0) = [];

figure('name','Original optical density')
subplot(121)
plot(t,dodDataRemoved(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')
subplot(122)
plot(t,dodDataRemoved(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

%%
addpath(genpath("C:\homer2"))

%% Wavelet

%the functions will know which channel is good and should be processed and which not
%Set iqr to 0.5
%iqr = 0.5; %threshold to detect outliers in wavelet details
%Run wavelet motion correction implemented in homer
%dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);
%save('dodWavelet.mat',"dodWavelet")
load("dodWavelet.mat")

waveletCorrection = meanValue.*exp(-dodWavelet);

% Visualize fNIRS data at the two wavelengths without removed channels
dgood_wavelet = waveletCorrection;
dgood_wavelet(:,SD.MeasListAct==0) = [];
%only good channel
figure('name','Wavelet process')
subplot(121)
plot(t,dgood_wavelet(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
subplot(122)
plot(t,dgood_wavelet(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')


% Plot wavelet corrected optical density data at first wavelength considering only good channels
dodWavGood = dodWavelet(:,SD.MeasListAct==1);
figure('name','Wavelet OD')
subplot(121)
plot(t,dodWavGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')
subplot(122)
plot(t,dodWavGood(:,end/2+1:end))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 2')

% Compare uncorrected vs. wavelet corrected data at each good channel
for iCh = 1:5 %nCh
    if SD.MeasListAct(iCh) == 1 % Plot only if it is a good channel
        figure('name','Comparison OD original-wavelet')
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodWavelet(:,iCh))
        title(num2str(iCh))
        legend('dod','dod Wavelet corrected')
        xlim([t(1) t(end)])
        pause
        close
    end
end

%% Spline motion correction
% Detect motion artifacts in signal
tMotion = 5; %window for motion artifacts detection
tMask = 10; %window of spline interpolation around the motion artifact
SDThresh = 5; %th on the std
AmpThresh = 0.15; %th on the apmlitude
fs = 1/(t(2)-t(1));

tIncMan = ones(length(t),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
% Motion detection technique. tIncCh is a matrix number of samples x twice n of
% channels which contains for each channel (column) the information about
% whether an artifact was present (0s) or not (1s). tInc is a vector which
% contains information on whether at that time sample in any of the channel
% was present an artifact (0s) or not (1s). tInc can therefore be obtained
% from tIncCh by setting to 0 every row that contains at least one 0. 
[tInc,tIncCh] = hmrMotionArtifactByChannel(dodConv, fs, SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);

% Spline interpolation
p = 0.99;
dodSpline = hmrMotionCorrectSpline(dodConv,t,SD,tIncCh,p);

splineCorrection = meanValue.*exp(-dodSpline);

% Visualize fNIRS data at the two wavelengths without removed channels
dgood_spline = splineCorrection;
dgood_spline(:,SD.MeasListAct==0) = [];
%only good channel
figure('name','Spline process')
subplot(121)
plot(t,dgood_spline(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
subplot(122)
plot(t,dgood_spline(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')


% Plot spline corrected optical density data in all good channels of first wavelength
dodSplineGood = dodSpline(SD.MeasListAct == 1);
figure('name','Spline OD')
subplot(121)
plot(t,dodSpline(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')
subplot(122)
plot(t,dodSpline(:,end/2+1:end))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 2')

% Compare uncorrected vs. wavelet corrected data at each good channel
for iCh = 1:7%nCh
    if SD.MeasListAct(iCh) == 1 % Plot only if it is a good channel
        figure('name','Comparison OD original-Spline')
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodSpline(:,iCh))
        title(num2str(iCh))
        legend('dod','dod Spline corrected')
        xlim([t(1) t(end)])
        pause
        close
    end
end


%% tPCA motion correction
varThresh = 0.80; % % of variance to remove
nIter = 5; % n of iterations

% Looking at the help of the function, we know that one of the outputs
% (tIncPCAbefore) contains information about detected motion artifacts
% before applying the PCA
[dodPCA,tIncPCAafter,svs,nSV,tIncPCAbefore] = hmrMotionCorrectPCArecurse(dodConv,fs,SD,tIncMan,tMotion,tMask,SDThresh,AmpThresh,varThresh,nIter);

PCAcorrection = meanValue.*exp(-dodPCA);

% Visualize fNIRS data at the two wavelengths without removed channels
dgood_pca = PCAcorrection;
dgood_pca(:,SD.MeasListAct==0) = [];
%only good channel
figure('name','tPCA process')
subplot(121)
plot(t,dgood_pca(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
subplot(122)
plot(t,dgood_pca(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

% Plot wavelet corrected optical density data at first wavelength considering only good channels
dodPCAGood = dodPCA(:,SD.MeasListAct==1);
figure('name','Comparison OD original-tPCA')
subplot(121)
plot(t,dodPCAGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')
subplot(122)
plot(t,dodPCAGood(:,end/2+1:end))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 2')

%compare optical density
for iCh = 1:5%nCh
    if SD.MeasListAct(iCh) == 1 % Plot only if it is a good channel
        figure(14);
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodPCA(:,iCh))
        title(num2str(iCh))
        legend('dod','dod PCA corrected')
        xlim([t(1) t(end)])
        pause
        close
    end
end

%% Band-pass filter
lowerCutOff = 0.01;
higherCutOff = 0.5; 
fs = 1/(t(2)-t(1));
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

Filt_signal = meanValue.*exp(-dodFilt);

% Visualize fNIRS data at the two wavelengths without removed channels
dgood_filt = Filt_signal;
dgood_filt(:,SD.MeasListAct==0) = [];
%only good channel
figure('name','Band-pass process')
subplot(121)
plot(t,dgood_filt(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
subplot(122)
plot(t,dgood_filt(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')


%comapare optical density
for iCh = 1:5%nCh
    if SD.MeasListAct(iCh) == 1 % Plot only if it is a good channel
        figure('name','Comparison OD wavelet-wavelet filtered')
        plot(t,dodWavGood(:,iCh))
        hold on;
        plot(t,dodFilt(:,iCh))
        title(num2str(iCh))
        legend('dod Wavelet','dod Wavelet filtered')
        xlim([t(1) t(end)])
        pause
        close
    end
end

% looking at the time series of the channel we identified a lot of
% MA spread all over the channels, so our intention was to apply PCA. PCA
% perform a great correction of the motion artifact, but it seems to erase
% also part of the signal. performing Wavelet is possible to correct the
% majority of MA and after the filtering step is possible to see how
% physiological noise, like Mayer Waves, is attenuated. Spline interpolation 
% seems to have poor effects since there aren't steps artifact or grate base drift.

%% Compute block-averaged hemodynamic response on the optical density data

tRange = [-2 40]; % range of timimg around stimulus to define a trial
sRange = fix(tRange*fs); % convert the time in seconds to samples
tHRF = tRange(1):1/fs:tRange(2); % time vector for the hemodynamic response (and trials)
dodAvg = zeros(length(tHRF),size(dodFilt,2),size(s,2)); % initialize the matrix that will contain our average hemodynamic response for each channel (for both wavelength) and condition
for iS = 1:size(s,2) % for each condition
    % Get the timing of stimulus presentation for that condition
    stimulusTiming = find(s(:,iS)==1); 
    % Initialize the matrix that will contain the single trial responses
    % for that condition for all channels at both wavelengths
    ytrial = zeros(length(tHRF),size(dodFilt,2),length(stimulusTiming));
    
    nTrial = 0;
    for iT = 1:length(stimulusTiming) % for each stimulus presented (for eacht trial)
        if (stimulusTiming(iT)+sRange(1))>=1 && (stimulusTiming(iT)+sRange(2))<=size(dodFilt,1) % Check that there are enough data pre and post stimulus (this is useful to check that the first stimulus is presented at least 2 seconds after the start of the acquisition and that the last stimulus has at least 18 seconds of data afterwards)
            nTrial = nTrial + 1;
            ytrial(:,:,nTrial) = dodFilt(stimulusTiming(iT)+[sRange(1):sRange(2)],:); % extract the trial from the dc data
        end
    end
    
    % Average trials (the fourth dimension of the ytrial matrix)
    dodAvg(:,:,iS) = mean(ytrial(:,:,1:nTrial),3);
    % Correct for the baseline
    for ii = 1:size(dodAvg,2) % for each channel (but you can do it directly without the for cycle)
        foom = mean(dodAvg(1:-sRange(1),ii,iS),1); % compute baseline as average of the signal in the -2:0 seconds time range
        dodAvg(:,ii,iS) = dodAvg(:,ii,iS) - foom; % subtract the baseline from the average hemodynamic responses
    end
end

chSel = [39 57 63 80];
figure('name', 'DOD block average')
for iCh = 1:length(chSel)
    subplot(2,2,iCh)
    plot(tHRF,squeeze(dodAvg(:,chSel(iCh),1)),'r','LineWidth',2)
    hold on;
    plot(tHRF,squeeze(dodAvg(:,chSel(iCh)+nCh,1)),'b','LineWidth',2)
    plot(tHRF,squeeze(dodAvg(:,chSel(iCh),2)),'g','LineWidth',2)
    plot(tHRF,squeeze(dodAvg(:,chSel(iCh)+nCh,2)),'m--','LineWidth',2)
    legend('OD first wl right hand','OD second wl rigth hand','OD first wl left hand','OD second wl left hand')
    title(num2str(chSel(iCh)))
    xlabel('Time [s]')
    ylabel('\DeltaHb [M]')
    xlim([tHRF(1) tHRF(end)])
    ylim([-0.5 0.5])
end

%%
load("S07_texting.jac", "-mat") %J = one sensitivity matrix fo each wl channels x nodes
%J.vol contains the Jacobian values for each node of the tetrahedral mesh.

addpath(genpath("C:\iso2mesh-master"))

%% 7) sensitivity matrix

load('MNI\HeadVolumeMesh.mat')

% Display whole array sensitivity on head model
HeadVolumeMesh.node(:,4) = (sum(J{1}.vol));
%.node = 3d posizion - value of segmentation (4th column)
%.elem  = indexes of the nodes composing a tetrahedron the fifth column can contain the value that can be associated to an element
% Display whole array sensitivity on GM volume mesh
figure('name', 'brain mesh all channels')
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
colorbar
title('all channels')
caxis([-3 0])

% Remove bad channels from Jacobian
for i = 1:length(SD.Lambda) % For each wavelength (we have two Js)
    tmp = J{i}.vol;
    JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end

HeadVolumeMesh.node(:,4) = (sum(JCropped{1}));

% Display whole array sensitivity on GM volume mesh
figure('name', 'brain mesh only good channels')
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
colorbar
title('only good channels')
caxis([-3 0])

%% 8) inverse problem 

load("MNI\GMSurfaceMesh.mat")

% Compute inverse of Jacobian (TIKONOV)
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our
% data
datarecon = -dodAvg;

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct (average dod)
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nCond = size(datarecon,3);
nWavs = length(SD.Lambda);

% Initialize final results matrices
hbo.vol = zeros(nFrames,nNodeVol,nCond);
hbr.vol = zeros(nFrames,nNodeVol,nCond);
hbo.gm = zeros(nFrames,nNodeGM,nCond);
hbr.gm = zeros(nFrames,nNodeGM,nCond);

% Obtain specific absorption coefficients for each wl
Eall = [];
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
end

% For each frame
for cond = 1:nCond
    for frame = 1:nFrames 
    
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
    
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);
    
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
    
        % Book-keeping and saving
        hbo.vol(frame,:,cond) = hbo_tmpVol;
        hbr.vol(frame,:,cond) = hbr_tmpVol;
        hbo.gm(frame,:,cond) = hbo_tmpGM;
        hbr.gm(frame,:,cond) = hbr_tmpGM;
    
    end
end

%% Plot reconstructed images
tRecon = [0 10 18]; % time point
baseline = abs(tRange(1)); % two seconds of baseline
sRecon = fix(tRecon*fs)+fix(baseline*fs); % Convert to samples
load greyJet % load colormap to make better images
hbO_image = cell(2,3); %first row first condition
hbR_image = cell(2,3); %second row second contidtion

for cond=1:nCond
    figure;
    sgtitle('condition '+string(cond))
    
    for instant = 1:length(tRecon)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbo.gm(sRecon(instant),:,cond);
        hbO_image{cond,instant} = GMSurfaceMesh.node;
        subplot(2,3,instant)
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        caxis([-0.03 0.03]) % Set the limit of the colorbar
        view([0 90]) % Set the view angle
        title(['HbO instant ' num2str(tRecon(instant)) 's'])
        colormap(greyJet) % set the loaded colormap
        hb = colorbar;
        hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
        axis off % remove axis
        
        GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(instant),:,cond);
        hbR_image{cond,instant} = GMSurfaceMesh.node;
        subplot(2,3,instant+3)
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        view([0 90])
        caxis([-0.03 0.03])
        title(['HbR instant ' num2str(tRecon(instant)) 's'])
        colormap(greyJet)
        hb = colorbar;
        hb.Label.String = {'\DeltaHbR [\muM]'};
        axis off
    end
end

%%
MeasListAct = SD.MeasListAct;
save('results_HW2_gruppo2','MeasListAct', 'hbO_image', 'hbR_image')


