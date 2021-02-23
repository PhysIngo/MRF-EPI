% T2 Map Reconstruction by Ingo Hermann, 2020-10-08
% Creale LUT for MRI scanner that MRF reconstruction can be performed
% directly through ICE on Siemens Scanner
% --------------------------------
% This scripts needs the user functions:
% MosaicOnOff.m
% openMaps.m
% makeMRFdictionaryAll.m

%% Config

% path where to find DICOM Files

mainp = 'D:\005 Matlab\MRF Scripts\VE11A\MRFPROSTATA_INVIVO';%DICOM4Group

% Create Maps with matlab? This is to check if the calcualted maps of the
% MRI are the same as the ones calcualted offline here in Matlab
calcMaps = true;
groupMatch = true;
generateLUT = false;
Mosaic = false;
Groups = 5000;

LutName = 'Prostate';
clear T2starList T1List

% generate the List of T1 and T2* times
T1List(1) = 1;
count = 2;
tmpV = 300;
while (tmpV < 3500)
    T1List(count) = tmpV;
    tmpV=tmpV.*1.05;
    count = count+1;
end
% tmpV = 5;

T2starList(1) = 1;
count = 2;
tmpV = 10;
while (tmpV < 500)
    T2starList(count) = tmpV;
    tmpV=tmpV.*1.05;
    count = count+1;
end

AngleCorrectList = 0.60:0.1:1.4;

% T1List=[500:200:3000];
% T2starList=[20:4:100];
% AngleCorrectList = 1;

OffResList=0;  %List of off-resonance values to be simulated [Hz]

% open the baseline images and the headers 
[tmpImg, tmpHead] = openMaps([mainp,'\*'],'number',1);

ImageSize = tmpHead.Width;
realTR = tmpHead.RepetitionTime;
realTE = tmpHead.EchoTime;

cut = 0;

% open the header and extract additional informations as lRepetitions...
fprintf('%s Reading in Dicom data\n', datestr(now));

% get list of DICOM files (should be either *.IMA or *.dcm)
filename = dir([mainp,'\*' ]);

% filename = filename(~ismember({filename.name},{'.','..'}));
filename = filename(~[filename.isdir]);

% read in dicoms information and extract relevant information for reco
FileNames = struct2cell(filename);

FileNames = strcat([mainp '\'], FileNames(1,:));
dcmInfos= cellfun(@(x) dicominfo(x), FileNames, 'UniformOutput', 0);
dcmInfos = cell2mat(dcmInfos);

dicomIma = fileread(FileNames{1,1});
allImas = strfind(dicomIma,'sWipMemBlock.alFree');
WipMemBlock = dicomIma(1,allImas(1):allImas(end)+30);

idx = strfind(dicomIma,'lRepetition');
idx(idx<8000) = [];
if dicomIma(idx(1)+17)<2
    meas = 1+str2double(dicomIma(idx(1)+17:idx(1)+19));
else
    meas = 1+str2double(dicomIma(idx(1)+17:idx(1)+18));
end

SSPat = strfind(dicomIma,'sWipMemBlock.alFree[3]');
if isempty(SSPat)
    SSPatString = 'Prep';
else
    SSPatString = '';    
end
Pauses = strfind(dicomIma,'sWipMemBlock.alFree[4]');
if isempty(Pauses)
    PauseDur = 0;
else
    PauseDur = 2500; 
end
Spoil = strfind(dicomIma,'sWipMemBlock.alFree[11]');
if isempty(Spoil)
    SpoilString = '';
else
    SpoilString = 'Spoil';    
end
Pulses = strfind(dicomIma,'sWipMemBlock.alFree[13]');
if isempty(Pulses)
    PulseNumber = 0;
else
    PulseNumber = str2double(dicomIma(Pulses+28)); 
    if PulseNumber<2 && PulseNumber~=0
        PulseNumber = str2double(dicomIma(Pulses+28:Pulses+29)); 
    end        
end
fprintf('The WIP parameters are: %s, %s, pause: %dms, IR: %d\n',...
    SSPatString,SpoilString,PauseDur,PulseNumber); 


%% read in Dicoms

% important to have this for either mosaic basline images or not...
% This might by replaced by the MosaicOnOff.m function...
if ~Mosaic
    [alltmpImg, head] = openMaps([mainp,'\*'],'all');
    SliceNum = round((size(dir(mainp),1)-2)/meas);
    SingleNum = ceil(sqrt(SliceNum));
    NewImageSize = ImageSize*ceil(sqrt(SliceNum));
%     MosaicImg = MosaicOnOff(alltmpImg,SliceNum);
    MosaicImg = zeros(NewImageSize,NewImageSize,meas);
    for measi=1:1:meas
        for x=1:1:SingleNum
            for y=1:1:SingleNum
                nn = y+SingleNum*(x-1)+(measi-1)*SliceNum;
                [row,col] = find(AllInstNumbers==nn);
                if y+SingleNum*(x-1)>SliceNum
                    continue;
                end
                MosaicImg(1+ImageSize*(x-1):ImageSize*(x),...
                    1+ImageSize*(y-1):ImageSize*(y),measi) = alltmpImg(:,:,col);
            end
        end
    end
else
    [alltmpImg, head] = openMaps([mainp,'\*'],'all');
end
  
Measurements = length(FileNames)-cut;
if ~Mosaic
    Measurements = meas;
end


% get slice location and so on here
% slice timing is important since it is a interleaved acquisition
numbers = zeros(length(dcmInfos)-cut,1);
SliceLoc = numbers;
for x1 = 1:length(dcmInfos)-cut
    numbers(x1) = dcmInfos(x1).InstanceNumber;
    SliceLoc(x1) = dcmInfos(x1).SliceLocation;  % slice location
end

SingleImageSize = dcmInfos(1).AcquisitionMatrix; % image dimensions
if SingleImageSize(1) == 0
    SingleImageSize(1) = SingleImageSize(2);
    SingleImageSize(2) = 0;
    SingleImageSize(4) = SingleImageSize(3);
    SingleImageSize(3) = 0;
end
% timing of each slice
AC = [char(dcmInfos.AcquisitionDate) char(dcmInfos.AcquisitionTime)];
AcquisitionTime = datenum(AC, 'yyyymmddHHMMSS.FFF');

% read in dicom data
[~,sortInds] = sort(numbers);
DCMIMAGES = cellfun(@(x) double(dicomread(x)), FileNames, 'UniformOutput', 0);
DCMIMAGES = cat(3,DCMIMAGES{sortInds});
if ~Mosaic
    DCMIMAGES = MosaicImg;
end
% DCMIMAGES=squeeze(DCMIMAGES(:,:,cut+1:end));

SliceLoc = SliceLoc(sortInds);
AcquisitionTime = AcquisitionTime(sortInds);

% get number of slices
slices = unique(SliceLoc);
% get the timing for SMS and Mulstilice readout, as they order the images weirdly
% try
%     SliceTiming = dcmInfos(1).Private_0019_1029;
% catch
%     SliceTiming = 1;
% endNumSlices = length(slices);
if ~Mosaic
    NumSlices = 1;
    nSliceLoc = SliceLoc(1)+SliceLoc.*0;
    nslices = slices(1)+slices.*0;
    clear SliceLoc;
    SliceLoc = nSliceLoc(1:meas);
%         SliceTiming = [1:SliceNum];
    SliceTiming = reshape([SliceNum/2+1:SliceNum; 1:SliceNum/2],1,SliceNum);
%     SliceTiming = [31 1 32 2 33 3 34 4 35 5 36 6 37 7 38 8 39 9 40 10 41 11 42 12 43 13 44 14 45 15 46 16 47 17 48 18 49 19 50 20 51 21 52 22 53 23 54 24 55 25 56 26 57 27 58 28 59 29 60 30]';
else
    SliceTiming = reshape([SliceNum/2+1:SliceNum; 1:SliceNum/2],1,SliceNum);
%     SliceTiming = [31 1 32 2 33 3 34 4 35 5 36 6 37 7 38 8 39 9 40 10 41 11 42 12 43 13 44 14 45 15 46 16 47 17 48 18 49 19 50 20 51 21 52 22 53 23 54 24 55 25 56 26 57 27 58 28 59 29 60 30]';
end
SliceTimingUinque = unique(SliceTiming);
[atmp,btmp] = sort(SliceTimingUinque);
for x1 = 1:length(atmp)
    SliceTiming(SliceTiming == atmp(x1)) = btmp(x1);
end
% slice timing cause it is a interleaved acquisition
NumSlices = length(slices);
if ~Mosaic
    NumSlices = 1;
    nSliceLoc = SliceLoc(1);
    nslices = slices(1);
    clear SliceLoc;
    SliceLoc = nSliceLoc(1:meas);
    SliceTiming = reshape([SliceNum/2+1:SliceNum; 1:SliceNum/2],1,SliceNum);
end

MRFBilder = zeros(size(DCMIMAGES,1),size(DCMIMAGES,2),size(DCMIMAGES,3)/NumSlices, NumSlices);
[d1,d2,d3] = size(MRFBilder);
% for 3D reshape
for iloop = 1:NumSlices
    MRFBilder(:,:,:,iloop) = DCMIMAGES(:,:,SliceLoc == slices(iloop));
    TRReadout3D(:, iloop) = AcquisitionTime(SliceLoc == slices(iloop));
end

% get TR timing of each slice
TRAfterExcitation = str2double(cellstr(datestr(TRReadout3D(1,:) - min(TRReadout3D(1,:)), 'ss.FFF')))*1000;
for x1 = 1:size(TRReadout3D,2)
    TRReadout3D(1:end-1,x1) = str2double(cellstr(datestr(TRReadout3D(2:end,x1) - TRReadout3D(1:end-1,x1), 'ss.FFF')))*1000;
end


SliceTimingUinque = unique(SliceTiming);
[atmp,btmp] = sort(SliceTimingUinque);
for x1 = 1:length(atmp)
    SliceTiming(SliceTiming == atmp(x1)) = btmp(x1);
end

% if you use SMS you only have a fraction of the dicitonaries
SMSPatINFO = dcmInfos(1).Private_0051_1011;
try
    iPat = str2double(SMSPatINFO(2:3));
    SMSFactor = str2double(SMSPatINFO(5:end));
catch
    iPat = str2double(SMSPatINFO(2:end));
    SMSFactor = 1;
end

% get the acqusistion order of the slices
NumSlicesSMS = length(SliceTiming);
if ~Mosaic
    NumSlicesSMS = SliceNum;
end
NumOfSlices = NumSlicesSMS/SMSFactor;
numX = size(DCMIMAGES,1) / SingleImageSize(1);
numY = size(DCMIMAGES,2) / SingleImageSize(4);
SliceAcquistionOrder = zeros(size(DCMIMAGES,1), size(DCMIMAGES,2));

% calculate the SliceAcquistionOrder
ct = 0;
for x1 = 1:numX
    for x2 = 1:numY
        ct = ct+1;
        if ct <= NumSlicesSMS
            xi = ((x1-1)*SingleImageSize(1)+1):x1*SingleImageSize(1);
            yi = ((x2-1)*SingleImageSize(4)+1):x2*SingleImageSize(4);
            SliceAcquistionOrder(xi,yi) = SliceTiming(ct);
        end
    end
end
TE_min = dcmInfos(1,1).EchoTime;

MRFBilderNew = MRFBilder;
%% calculate dictionaries / LUTs and save them for a scanenr
if generateLUT
    clear allDic;
    fprintf('%s Starting to generate %d LUTs\n', datestr(now), NumOfSlices);

    % LUTs will be saved in the "..\DICOM\LUTxx'
    svpath = strcat(mainp,'\',LutName, num2str(NumOfSlices));
    [~,~ ] = mkdir(svpath);

    NTR = meas; %length(TRReadout3D);

    % TEs from cpp
    TEs_3_new = [0.0, 3.0, 5.1, 6.2, 5.2, 2.8, 0.4, 0.3, 3.7, 9.7, 16.0, 19.4, 17.9, 4.7, 9.8, 21.1, 30.4, 33.1, 27.4, 4.8, ...
        0.0, 5.2, 18.8, 34.8, 45.6, 34.2, 17.3, 3.5, 7.2, 46.0, 60.5, 30.6, 0.5, 0.9, 3.5 ...
        0.0, 3.0, 5.1, 6.2, 5.2, 2.8, 0.4, 0.3, 3.7, 9.7, 16.0, 19.4, 17.9, 4.7, 9.8, 21.1, 30.4, 33.1, 27.4, 4.8, ...
        0.0, 5.2, 18.8, 34.8, 45.6, 34.2, 17.3, 3.5, 7.2, 46.0, 60.5, 30.6, 0.5, 0.9, 3.5 ...
        0.0, 3.0, 5.1, 6.2, 5.2, 2.8, 0.4, 0.3, 3.7, 9.7, 16.0, 19.4, 17.9, 4.7, 9.8, 21.1, 30.4, 33.1, 27.4, 4.8, ...
        0.0, 5.2, 18.8, 34.8, 45.6, 34.2, 17.3, 3.5, 7.2, 46.0, 60.5, 30.6, 0.5, 0.9, 3.5 ...
        0.0, 3.0, 5.1, 6.2, 5.2, 2.8, 0.4, 0.3, 3.7, 9.7, 16.0, 19.4, 17.9, 4.7, 9.8, 21.1, 30.4, 33.1, 27.4, 4.8, ...
        0.0, 5.2, 18.8, 34.8, 45.6, 34.2, 17.3, 3.5, 7.2, 46.0, 60.5, 30.6, 0.5, 0.9, 3.5];
    TE_new = TEs_3_new + TE_min;


    % from CPP - in Matlab format
    FlipAngles_List_new= [ 0, 2, 4, 6, 8, 10, 12, 14, 16, 18,...
        20, 22, 24, 26, 28, 30, 32, 34, 36, 38, ...
        40, 42, 44, 46, 48, 50, 52, 54, 56, 58,...
        60, 62, 64, 66, 68, 70, 72, 74, 76, 78,...
        80, 82, 84, 86, 88, ...
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18,...
        20, 22, 24, 26, 28, 30, 32, 34, 36, 38, ...
        40, 42, 44, 46, 48, 50, 52, 54, 56, 58,...
        60, 62, 64, 66, 68, 70, 72, 74, 76, 78,...
        80, 82, 84, 86, 88, ... 
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18,...
        20, 22, 24, 26, 28, 30, 32, 34, 36, 38, ...
        40, 42, 44, 46, 48, 50, 52, 54, 56, 58,...
        60, 62, 64, 66, 68, 70, 72, 74, 76, 78,...
        80, 82, 84, 86, 88, ... 
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18,...
        20, 22, 24, 26, 28, 30, 32, 34, 36, 38, ...
        40, 42, 44, 46, 48, 50, 52, 54, 56, 58,...
        60, 62, 64, 66, 68, 70, 72, 74, 76, 78,...
        80, 82, 84, 86, 88 ];
    NFR = length(FlipAngles_List_new);
    SelectFlipAngles_new = [
        23, 21, 25, 28, 29, 29, 27, 24, 19, 18, 23, 28, 32, 32, 19, 20, 30, 38, 43, 41, ...
        36, 27, 17, 20, 25, 32, 33, 34, 33, 30, 27, 23, 19, 18, 21,...
        23, 21, 25, 28, 29, 29, 27, 24, 19, 18, 23, 28, 32, 32, 19, 20, 30, 38, 43, 41, ...
        36, 27, 17, 20, 25, 32, 33, 34, 33, 30, 27, 23, 19, 18, 21,...
        23, 21, 25, 28, 29, 29, 27, 24, 19, 18, 23, 28, 32, 32, 19, 20, 30, 38, 43, 41, ...
        36, 27, 17, 20, 25, 32, 33, 34, 33, 30, 27, 23, 19, 18, 21,...
        23, 21, 25, 28, 29, 29, 27, 24, 19, 18, 23, 28, 32, 32, 19, 20, 30, 38, 43, 41, ...
        36, 27, 17, 20, 25, 32, 33, 34, 33, 30, 27, 23, 19, 18, 21];

    SelectFlipAngles(1:NTR) = SelectFlipAngles_new(1:NTR);
    FlipAngles_List(1:NFR) = FlipAngles_List_new(1:NFR);
    TE(1:NTR) = TE_new(1:NTR);
    TEs_3(1:NTR) = TEs_3_new(1:NTR);
    
    % each slice get it's own dictionary
    figure;hold on;
    
    clear newPulses newTimes;
    for SL = 1:NumOfSlices

        fprintf('%s Generate LUT: slice %d of %d\n', datestr(now), SL, NumOfSlices);
        TR = TRReadout3D(:,1);
        if NumOfSlices == 1
            TR = TEs_3 + realTR;
        else
            TR = TEs_3 + (realTR-20)/NumOfSlices;
        end
        RFpulses = FlipAngles_List(SelectFlipAngles+1)./180*pi;

        % calc RefDic
        RefDic = combvec(T1List,T2starList,OffResList, AngleCorrectList);
        RefDic(:,RefDic(1,:) < RefDic(2,:)) =  [];
        
        % always show the dictionary
        ShowDictionary = true;
        if ShowDictionary
            [D, TRsim, alltmr, allpuls]= makeMRFdictionaryAll(RFpulses, ...
                RefDic ,TR , TE, NumOfSlices, SL, iPat, SMSFactor,SSPatString,SpoilString,...
                'Pulses', PulseNumber, 'Pause', PauseDur,'Change',0);
        end
        
        % show the simulated rf pulses 
        precisionFactor = 10;
        if NumOfSlices>20
            precisionFactor = 1;
        end
        if NumOfSlices<200
            alltimes = 1/precisionFactor:1/precisionFactor:round(alltmr(end),0)+1;
            allpulses = zeros(length(alltimes),1);
            firstCount = 1;
            for tmrCount = 1:1:round(alltmr(end)+1,0)*precisionFactor
                if length(allpuls)>=firstCount
                    if tmrCount == round(alltmr(firstCount)*precisionFactor)
                        allpulses(tmrCount) = allpuls(firstCount);
                        firstCount = firstCount+1;
                    end
                end
            end
            cMap = parula(NumOfSlices);
                    hold on;
            if NumOfSlices<20
                bar(alltimes(allpulses~=0),allpulses(allpulses~=0),0.025,'FaceColor',cMap(SL,:),'EdgeColor',cMap(SL,:));
                bar(alltimes(allpulses>1.65),allpulses(allpulses>1.65),0.025,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
                bar(alltimes(allpulses==0.85),allpulses(allpulses==0.85),0.025,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
            else
                bar(alltmr,allpuls,0.025,'FaceColor',cMap(SL,:),'EdgeColor',cMap(SL,:));
                bar(alltmr(allpuls>1.65),allpuls(allpuls>1.65),0.025,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
            end
        end
        
        newTimes(SL,:) = alltmr(1,:);
        newPulses(SL,:) = allpuls(1,:);
        
        
        Dictionary = abs(D);

        % only matchi if the relevant slice is there
        Mtmp = reshape(abs(MRFBilder), d1*d2,d3);
        Mtmp = Mtmp(SliceAcquistionOrder == SL,:);

        % create the LUTs for the scanner
        % the LUTs need to be saved in an anatomical order, but the dictionaries
        % are calcualted in a temporal order
        % anatomica order: Slice 1 = slices that is "up the top of the brain",
        % slice 2 = 1 slice "futher down"
        % temporal order: Slice 1 = first slices that is measured, Slice 2 = 2nd
        % slice that is measured
        % therefore we need to map the temporal onto the anatomical slice order
        idxAnatomical = find(SliceTiming == SL);

        D = single(abs(D));
        % LUT values
        ttt = sqrt(sum(D.^2, 2));
        
        % save all LUT to txt files in the folder defined in the beginning
        sname2 = strcat(svpath, '\LUTFactor_', num2str(idxAnatomical(1)-1), '.txt');
        fileID = fopen(sname2,'w');
        for i = 1:size(ttt,1)
            fprintf(fileID, '%f ', ttt(i,1));
            fprintf(fileID, '\r\n');
        end
        fclose(fileID);
        
        % when you have SMS factor 3
        if SMSFactor == 3
            sname3 = strcat(svpath, '\LUTFactor_', num2str(20+idxAnatomical(1)-1), '.txt');
            fileID3 = fopen(sname3,'w');
            for i = 1:size(ttt,1)
                fprintf(fileID3, '%f ', ttt(i,1));
                fprintf(fileID3, '\r\n');
            end
            fclose(fileID3);

            sname3 = strcat(svpath, '\LUTFactor_', num2str(40+idxAnatomical(1)-1), '.txt');
            fileID3 = fopen(sname3,'w');
            for i = 1:size(ttt,1)
                fprintf(fileID3, '%f ', ttt(i,1));
                fprintf(fileID3, '\r\n');
            end
            fclose(fileID3);
        end
        
        D = D./repmat(sqrt(sum(D.^2, 2)), 1, Measurements); % normalize dictionarz
        c = horzcat(RefDic',D);
        allDic(:,:,SL) = c;

        % write the LUTs
        sname = strcat(svpath, '\LUT_', num2str(idxAnatomical(1)-1), '.txt');
        fileID = fopen(sname,'w');
        for i = 1:size(c,1)
            for ii = 1:size(c,2)
                fprintf(fileID, '%f ', c(i,ii));
            end
            fprintf(fileID, '\r\n');
        end
        fclose(fileID);
        
        % when you have SMS factor 3
        if SMSFactor == 3
            sname = strcat(svpath, '\LUT_', num2str(20+idxAnatomical(1)-1), '.txt');
            fileID = fopen(sname,'w');
            for i = 1:size(c,1)
                for ii = 1:size(c,2)
                    fprintf(fileID, '%f ', c(i,ii));
                end
                fprintf(fileID, '\r\n');
            end

            sname = strcat(svpath, '\LUT_', num2str(40+idxAnatomical(1)-1), '.txt');
            fileID = fopen(sname,'w');
            for i = 1:size(c,1)
                for ii = 1:size(c,2)
                    fprintf(fileID, '%f ', c(i,ii));
                end
                fprintf(fileID, '\r\n');
            end
        end
            
            
            
            
    end

     %% _______________________ Generate GroupMatch Files ______________________%
     % this is only for fast group matching 
     clear allnc2;
    if (groupMatch)
        SLStart = 1;
        for SL = SLStart:NumOfSlices
            dGroups = Groups;  
                            
            numberOfMatchingGroups = dGroups/10;
            [dx, dy] = size(D);
            aGroupNumber = zeros(dx,1); % Per D signal, which groupt it's in
            aGroupSize = ceil(dx/dGroups); % # lenght D / # groups
            missCol = abs(dx-dGroups*aGroupSize);
        %     dGroups = dGroups - ceil( missCol/aGroupSize)+1;
            if dx-dGroups*aGroupSize<0
                dGroups = ceil(dx/aGroupSize)-1;
            end
            if missCol == 0
                dGroups = dGroups-1;
            end
            missCol = aGroupSize - abs(dx-dGroups*aGroupSize);
            NumOfGroups = dGroups;
            aMeanSignals = zeros(dGroups,dy); % Average Signal of each group
            cGroupDict = cell(aGroupSize,3); % data saved each group: 1 = Dk, 2 = VK, 3 = numDic
            numDic = 1:dx;
            idxAnatomical = find(SliceTiming == SL);

            RefDic(:,:) =  squeeze(allDic(:,1:4,SL)');
            D(:,:) = squeeze(allDic(:,5:end,SL));
            fprintf('%s Make the template signals of %s\n', datestr(now), num2str(idxAnatomical(1)-1))
            DicLoop = D;
            for x1 = 1:dGroups
                tmpVal = D(aGroupNumber==0,:);
                aMeanSignals(x1,:) = tmpVal(ceil(rand*(size(tmpVal,1)-1))+1,:);
                aMeanSignals(x1,:) = aMeanSignals(x1,:)./repmat(sqrt(sum(aMeanSignals(x1,:).^2,2)),1,size(aMeanSignals(x1,:),2));
                dotpro=abs(reshape(DicLoop(DicLoop>=0),[],Measurements)*aMeanSignals(x1,:)');
                [~,index] = sort(dotpro); 
                if (size(aGroupNumber(aGroupNumber == 0),1) < aGroupSize)
                    aGroupNumber(index(end-size(aGroupNumber(aGroupNumber == 0),1)+1:end)) = x1;
                else
                    aGroupNumber(index(end-aGroupSize+1:end)) = x1;
                end
                DicLoop(index(end-aGroupSize+1:end),:) = 0;
            end
            aMeanSignals = aMeanSignals./repmat(sqrt(sum(aMeanSignals.^2,2)),1,size(aMeanSignals,2));


            newD = D;
            nRefDic = RefDic;

            fprintf('%s Write GroupMatch to file of %s\n', datestr(now), num2str(idxAnatomical(1)-1))
            
            sname = strcat(svpath,'\','LUT_', num2str(idxAnatomical(1)-1),'.txt');
            fileID2 = fopen(sname,'w');
            old = 1;
            for x1 = 1:dGroups

                newD = single(abs(D(aGroupNumber==x1,:)));
                newD = newD./repmat(sqrt(sum(newD.^2, 2)), 1, Measurements); % normalize dictionarz
                c1 = horzcat(nRefDic(:,aGroupNumber==x1)',newD);
                
                savesz = size(aGroupNumber(aGroupNumber == x1),1);
                for ii = 1:savesz
                    for jj = 1:size(c1,2)
                        fprintf(fileID2, '%f ', c1(ii,jj));
                    end
                    fprintf(fileID2, '\r\n');
                end

                allc1(old:old+savesz-1,:,SL) = c1(:,:);
                old = x1*savesz+1;

            end
            fclose(fileID2);

            nc = aMeanSignals*aMeanSignals';
            [~, nc2] = sort(1-nc,2);
            sname = strcat(svpath,'\','GroupCorr_', num2str(idxAnatomical(1)-1), '_', num2str(NumOfSlices), '_noGroup_',num2str(NumOfGroups),'.txt');
            fileID3 = fopen(sname,'w');        
            for ii = 1:size(aMeanSignals,1)
                % amount of matching groups
                for jj = 1:numberOfMatchingGroups
                    fprintf(fileID3, '%0.0f ', nc2(ii,jj)-1);
                end
                fprintf(fileID3, '\r\n');
            end
            fclose(fileID3);

            sname = strcat(svpath,'\','GroupMean_', num2str(idxAnatomical(1)-1), '_', num2str(NumOfSlices), '_noGroup_',num2str(NumOfGroups),'.txt');
            fileID4 = fopen(sname,'w');        
            for ii = 1:size(aMeanSignals,1)
                for jj = 1:size(aMeanSignals,2)
                    fprintf(fileID4, '%f ', aMeanSignals(ii,jj));
                end
                fprintf(fileID4, '\r\n');
            end
            fclose(fileID4);
    %         allMeanSignals(:,:,SL) = aMeanSignals(:,:);
            allMeanSignals(1:numberOfMatchingGroups,:,SL) = aMeanSignals(1:numberOfMatchingGroups,:);
            allnc2(:,1:numberOfMatchingGroups,SL) = nc2(:,1:numberOfMatchingGroups);
        end
    end
    fprintf('%s Finished calculating LUTs \n', datestr(now));
else% otherwise load the LUT if it already exists!!!
    % loading LUT
    dGroups = Groups;
    fprintf('%s Starting loading LUTs \n', datestr(now));
    numberOfMatchingGroups = dGroups/10;% make 5 groups -- this needs to be adjusted
    clear allnc2 allDic;
    idxAnatomical = 1;
    svpath = strcat(mainp,'\',LutName, num2str(NumOfSlices));
    name = strcat(svpath,'\LUT_',num2str(idxAnatomical-1),'.txt');
    c=dlmread(name);
    dx = size(c,1);

    aGroupSize = ceil(dx/dGroups); % # lenght D / # groups
    dGroups = dx/aGroupSize;
    NumOfGroups = dGroups;
    for SL = 1:NumOfSlices
        idxAnatomical = find(SliceTiming == SL);
        fprintf('%s loading LUT_%d \n', datestr(now),idxAnatomical-1);
        svpath = strcat(mainp,'\',LutName, num2str(NumOfSlices));
        name = strcat(svpath,'\LUT_',num2str(idxAnatomical-1),'.txt');
        c=dlmread(name); 
        allDic(:,:,SL) = c;
        if groupMatch == 1
            name = strcat(svpath,'\GroupMean_',num2str(idxAnatomical-1),'_',num2str(NumOfSlices),'_noGroup_',num2str(NumOfGroups),'.txt');
            aMeanSignals=dlmread(name); 
            allMeanSignals(:,:,SL) = aMeanSignals(1:numberOfMatchingGroups,:);
            name = strcat(svpath,'\GroupCorr_',num2str(idxAnatomical-1),'_',num2str(NumOfSlices),'_noGroup_',num2str(NumOfGroups),'.txt');
            nc2=dlmread(name); 
            allnc2(:,1:numberOfMatchingGroups,SL) = nc2(:,1:numberOfMatchingGroups);  
            allc1 = allDic;
        end
    end

    fprintf('%s Finished loading LUTs \n', datestr(now));

end
MRFBilderNew = MRFBilder;
%% Match the D with the data if maps should be generated
% idea: match the D with the dicionary (so no real measurement) to
% analyze a wide specrum of T1 and T2* values

% perform some denoising here if you want
MRFBilderNew = MRFBilder;

% apply some denoising if you want:
% but additionally needs Denoising_Ingo.m, MPdenoising.m
GaussF = false;
DeepL = true;
MedianF = false;
MPPCA = true;

if MedianF
    for i=1:1:meas
        MRFBilderNew_Median(:,:,i) = medfilt2(squeeze(MRFBilder(:,:,i)), [3 3]);
    end
    fprintf('%s Finished with Median Filter way \n', datestr(now))
end
    
if GaussF
    MRFBilderNew_Gauss = Denoising_Ingo(MRFBilder, 0.5, 'Gauss');
    fprintf('%s Finished with Gaussian Filter \n', datestr(now))
end

if DeepL
    net = denoisingNetwork('DnCNN');
    for netCounter = 1:1:meas
        noisyI = uint8(MRFBilder(:,:,netCounter));
        denoiseI = denoiseImage(noisyI,net);
        MRFBilderNew_Deep(:,:,netCounter) = double(denoiseI);
    %     figure
    %     imshowpair(noisyI, MRFBilderNew_Gauss(:,:,netCounter),'montage');
    end
    fprintf('%s Finished with Deep Learning Filter \n', datestr(now))
end

if MPPCA
    for i=1:1:35
        tmp = MosaicOnOff(MRFBilderNew(:,:,i),NumOfSlices);
        MRFBilderNew_All(:,:,:,i) = tmp;
    end
    [denoised] = MPdenoising(MRFBilderNew_All);
    for i=1:1:35
        tmp =  denoised(:,:,:,i);
        MRFBilderNew(:,:,i) = MosaicOnOff(tmp,NumOfSlices);
    end
end


tic;

% make the matching in the traditional way since group matching is not much
% fast in Matlab and SVD also not since we only have 35 measurements

offSet = 0;
saveMaps = true;
fprintf('%s Matching the dicitionary on the traditional way\n', datestr(now))
tic
[d1,d2,d3] = size(MRFBilderNew(:,:,1+offSet:end));
T1Map = zeros(d1*d2,1);
T2Map = T1Map;
FitMap = T1Map;
B1Map = T1Map;
NTR = 35;
d3 = NTR-offSet;
for loop = 1:NumOfSlices
    fprintf('%s Matching slice %d of %d\n', datestr(now), loop, NumOfSlices);
    % Use magnitude of dictioanry
    Dictionary = abs(allDic(:,5+offSet:NTR+4,loop));
    RefDic = allDic(:,1:4,loop)';
    D = Dictionary;
    
    % only match if the relevant slice is there
    Mtmp = reshape(abs(MRFBilderNew(:,:,1+offSet:NTR)), d1*d2,d3);
    Mtmp = Mtmp(SliceAcquistionOrder == loop,:);
    M = Mtmp;
    
    
    [dd1,dd2] = size(M);
    % [dd1, dd2] = size(D);

    D = single(abs(D));
    D = D./repmat(sqrt(sum(D.^2,2)), 1,dd2);

    M = single(abs(M));
    M = M./repmat(sqrt(sum(M.^2,2)), 1, dd2);

    
    Mtmp = M;
    % calculate max memory (a bit of a gess, but better than nothing)
    CorrelationCoeff = zeros(dd1,1);
    MaxCorrVal = ones(dd1,1);

    [MaxCorrVal, CorrelationCoeff] = MatchMRF(D(:,:), M(:,:));
    T1Map(SliceAcquistionOrder == loop) = RefDic(1,MaxCorrVal);
    T2Map(SliceAcquistionOrder == loop) = RefDic(2,MaxCorrVal);
    FitMap(SliceAcquistionOrder == loop) = CorrelationCoeff;
    B1Map(SliceAcquistionOrder == loop) = RefDic(4,MaxCorrVal);
    
end
timetra = toc;
fprintf('%s Traditional way took: %0.2f\n', datestr(now), timetra)
T1Map = reshape(T1Map, d1,d2);
T2Map = reshape(T2Map, d1,d2);
FitMap = reshape(FitMap, d1,d2);
B1Map = reshape(B1Map, d1,d2);
%     FitMap(FitMap<0.98) = 0;
    
figure;imagesc(T1Map,[0 2000]);colormap;colorbar;
figure;imagesc(T2Map,[0 120]);colormap;colorbar;
figure;imagesc(B1Map,[min(B1Map(:)) max(B1Map(:))]);colormap;colorbar;
figure;imagesc(FitMap,[0.98 1]);colormap;colorbar;

%     T1Map = T1Map.*FitMap;
%     T2Map = T2Map.*FitMap;
%     B1Map = B1Map.*FitMap;

eval([sprintf('maps%0.0d',NTR),'.T1Map = T1Map;']);
eval([sprintf('maps%0.0d',NTR),'.T2Map = T2Map;']);
eval([sprintf('maps%0.0d',NTR),'.B1Map = B1Map']);
eval([sprintf('maps%0.0d',NTR),'.FitMap = FitMap']);
eval([sprintf('maps%0.0d',NTR),'.MRFBilderNew = MRFBilderNew']);

toc


T1 = MosaicOnOff(T1Map,NumOfSlices);
T2 = MosaicOnOff(T2Map,NumOfSlices);
B1 = MosaicOnOff(B1Map,NumOfSlices);
for i=1:1:NumOfSlices
%     idxAnatomical = find(SliceTiming == i);
    metadata = head(i,1);
    metadataAcquistionNumber = i;
    metadata.SeriesDescription = ['T1'];
    metadata.WindowCenter = 1000;
    metadata.WindowWidth = 1000;
    metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'1'];
    mkdir(metadata.SeriesDescription);
    metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
    dicomwrite(squeeze(T1(:,:,i))./256^2, [metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);

    
    metadata = head(i,1);
    metadata.SeriesDescription = ['T2'];
    metadata.WindowCenter = 100;
    metadata.WindowWidth = 100;
    metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'2'];
    mkdir(metadata.SeriesDescription);
    metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
    dicomwrite(squeeze(T2(:,:,i))./256^2, [metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);

    
    metadata = head(i,1);
    metadata.SeriesDescription = ['B1'];
    metadata.WindowCenter = 1000;
    metadata.WindowWidth = 500;
    metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'3'];
    mkdir(metadata.SeriesDescription);
    metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
    dicomwrite(squeeze(B1(:,:,i).*1000)./256^2, [metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);
end



