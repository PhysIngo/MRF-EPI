% T2 Map Reconstruction by Ingo Hermann, 2021-01-22
% Creale LUT for MRI scanner that MRF reconstruction can be performed
% directly through ICE on Siemens Scanner
% --------------------------------
% This script needs the user functions:
% MosaicOnOff.m
% openMaps.m
% makeMRFdictionaryAll.m

%% Config

% path where to find DICOM Files



% Create Maps with matlab? This is to check if the calcualted maps of the
% MRI are the same as the ones calcualted offline here in Matlab
calcMaps = true;
generateLUT = true;


isScanner = 4;   %1 for Mannheim, 2 for Barcelona, 3 for Trio, 4 for Maastricht
if isScanner == 1
    iPat = 3;SMSFactor = 1;NumOfSlices = 60;
    LutName = 'Mannheim';
    mainp = 'D:\005 Matlab\MRF Scripts\DicomMannheim\';
elseif isScanner == 2
	iPat = 2;SMSFactor = 3;NumOfSlices = 20;
    LutName = 'Barcelona';
    mainp = 'D:\005 Matlab\MRF Scripts\DicomBarcelona\';%DICOM4Group
elseif isScanner == 3
    NumOfSlices = 20;
    LutName = 'Trio';
elseif isScanner == 4
    iPat = 2;SMSFactor = 1;NumOfSlices = 31;
    LutName = 'Mannheim'; %Maastricht
    mainp = 'D:\004 Messungen\MRF_Knee_20211211_Scannexus\DICOM\10211_TEST_KNEE_21_02_11-10_07_06-STD-1_3_12_2_1107_5_2_43_67003\SEBASTIAN_WEINGAERTNER_PRISMA_DEVELOPMENT_20210211_100749_242000\MRF_KNEE_POST_FIRST_LOAD_1_0014\';
end
    
clear T2starList T1List
prec = 0.05;

% generate the List of T1 and T2* times
T1List(1) = 1;
count = 2;
tmpV = 300;
while (tmpV < 3500)
    T1List(count) = tmpV;
    tmpV=tmpV.*(1.00+prec);
    count = count+1;
end
% tmpV = 5;

T2starList(1) = 1;
count = 2;
tmpV = 10;
while (tmpV < 500)
    T2starList(count) = tmpV;
    tmpV=tmpV.*(1.00+prec);
    count = count+1;
end

AngleCorrectList = 0.60:0.1:1.4;

% overwrite the old list to save time
% this is a small list for test reasons
% T1List=[500:50:1500 1600:200:3000];
% T2starList=[20:2:100 120:20:300];
% AngleCorrectList = 0.8:0.2:2.4;

% T1List=[500:100:2000];
% T2starList=[20:5:60 70:10:150];
% AngleCorrectList = 0.8:0.2:1.2;

OffResList=0;  %List of off-resonance values to be simulated [Hz]

% open the baseline images and the headers 
[baseline,head,dcm,meta] = openMaps([mainp,'*'],'all');
baselineMRF = reshape(baseline,size(baseline,1),size(baseline,2),[],size(baseline,3));

ImageSize = head(1,1).Width;
realTR = head(1,1).RepetitionTime;
realTE = head(1,1).EchoTime;

cut = 0;

% open the header and extract additional informations as lRepetitions...
fprintf('%s Reading in Dicom data\n', datestr(now));

clear MRFBilder;
MRFBilderNew = squeeze(baselineMRF);
if length(size(MRFBilderNew))==3
    for i=1:1:35
        MRFBilder(:,:,:,i) = MosaicOnOff(MRFBilderNew(:,:,i),NumOfSlices*SMSFactor);
    end
end

MRFBilderNew = MRFBilder;

meas = size(MRFBilderNew,4);
if isScanner == 3
    NumOfSlices = size(MRFBilderNew,3);

    for i=1:1:size(baseline,3)
        a(i)= head(i,1).SliceLocation;
    end
    NumOfSlices = length(unique(a));
    meas = size(baseline,3)/NumOfSlices;
    MRFBilder = reshape(baseline,size(baseline,1),size(baseline,2),NumOfSlices,meas);
end
%% calculate dictionaries / LUTs and save them for a scanenr
if generateLUT
    clear allDic;
    fprintf('%s Starting to generate %d LUTs\n', datestr(now), NumOfSlices);

    % LUTs will be saved in the "..\DICOM\LUTxx'
    svpath = strcat(mainp,LutName, num2str(NumOfSlices));
    [~,~ ] = mkdir(svpath);

    NTR = meas; %length(TRReadout3D);

    if isScanner == 1 || isScanner == 2 || isScanner == 4 %Mannheim & Barcelona
        TEs_new = [0.0, 3.0, 5.1, 6.2, 5.2, 2.8, 0.4, 0.3, 3.7, 9.7, 16.0, 19.4, 17.9, 4.7, 9.8, 21.1, 30.4, 33.1, 27.4, 4.8, ...
            0.0, 5.2, 18.8, 34.8, 45.6, 34.2, 17.3, 3.5, 7.2, 46.0, 60.5, 30.6, 0.5, 0.9, 3.5];
        
        FlipAngles_List_new = [
            23, 21, 25, 28, 29, 29, 27, 24, 19, 18, 23, 28, 32, 32, 19, 20, 30, 38, 43, 41, ...
            36, 27, 17, 20, 25, 32, 33, 34, 33, 30, 27, 23, 19, 18, 21 ];
    elseif isScanner == 3 %Trio
        TEs_new = [0.0, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.6, 2.6, 2.6, 2.6, 2.5, 2.3, 2.2, 2, 1.8, 1.5, 1.3, 1, ...
            0.7, 0.5, 0.3, 0.2, 0, 0, 0, 0.1, 0.3, 0.6, 0.9, 1.3, 1.8, 2.4, 3, 3.7, 4.3, 5, 5.7, 6.3, ...
            6.9, 7.4, 7.8, 8.1, 8.2, 8.3, 8.2, 7.9, 7.5, 7, 6.4, 5.7, 4.9, 4, 3.2, 2.4, 1.6, 1, 0.5, 0.1, ...
            0, 0.1, 0.5, 1.1, 2, 3.1, 4.6, 6.2, 8.1, 10.1, 12.3, 14.5, 16.7, 18.9, 20.9, 22.7, 24.2, 25.5, 26.3, 26.8, ...
            26.8, 26.3, 25.4, 24.1, 22.3, 20.2, 17.8, 15.3, 12.6, 9.9, 7.3, 4.9, 2.9, 1.3, 0.3, 0, 0.4, 1.5, 3.5, 6.3, ...
            10, 14.3, 19.4, 25, 31.1, 37.4, 43.9, 50.4, 56.5, 62.2, 67.3, 71.5, 74.8, 76.9, 77.8, 77.4, 75.8, 72.8, 68.6, 63.2, ...
            57, 50, 42.4, 34.7, 27, 19.8, 13.2, 7.6, 3.4, 0.8, 0, 1.3, 4.6, 10.2, 18.1, 28, 39.9, 53.5, 68.6, 84.7, ...
            101.4, 118.4, 135, 150.7, 165.1, 177.7, 188, 195.7, 200.3,0.0, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.6];
        
        FlipAngles_List_new= [ 8, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 14, 14, 14, 15, 14, 14, 14, ...
            14, 13, 13, 12, 11, 11, 10, 9, 8, 7, 6, 4, 3, 2, 0, 3, 4, 6, 7, 8, ...
            10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 18, 18, 18, 18, 17, 17, 16, 16, 15, 14, ...
            13, 12, 11, 9, 8, 6, 5, 4, 2, 0, 3, 5, 8, 10, 13, 15, 17, 19, 21, 23, ...
            24, 26, 27, 28, 28, 29, 29, 29, 29, 28, 27, 26, 25, 24, 22, 21, 19, 16, 14, 12, ...
            10, 7, 5, 2, 0, 3, 4, 5, 6, 7, 9, 10, 10, 11, 13, 14, 14, 15, 16, 17, ...
            17, 17, 18, 18, 18, 19, 19, 19, 18, 18, 18, 18, 17, 17, 16, 15, 15, 14, 13, 12, ...
            11, 11, 10, 8, 7, 6, 5, 4, 2, 0, 2, 3, 4, 5, 6, 6, 7, 8, 9];
    end
    TE_new = TEs_new + realTE;
   
    FlipAngles_List(1:meas) = 2.*FlipAngles_List_new(1:meas);
    TE(1:meas) = TE_new(1:meas);
    
    % each slice get it's own dictionary
    figure;hold on;
    
    clear newPulses newTimes;
    for SL = 1:NumOfSlices

        fprintf('%s Generate LUT: slice %d of %d\n', datestr(now), SL, NumOfSlices);
        TR = TE+(realTR-20.05)/(NumOfSlices);
        if isScanner == 4
            TR = TR - realTE;
        end
        RFpulses = FlipAngles_List./180*pi;

        % calc RefDic
        RefDic = combvec(T1List,T2starList,OffResList, AngleCorrectList);
        RefDic(:,RefDic(1,:) < RefDic(2,:)) =  [];
        
        if isScanner == 1  %Mannheim
            [D, TRsim, alltmr, allpuls] = makeMRFdictionary_EPI_T2Star_Mannheim_s(RFpulses, ...
                RefDic ,TR , TE, NumOfSlices, SL, iPat, SMSFactor);
        elseif isScanner == 2   %Barcelona
            [D, TRsim] = makeMRFdictionary_EPI_T2Star_Barcelona_s(RFpulses, ...
                RefDic ,TR , TE, NumOfSlices, SL, iPat, SMSFactor);
            alltmr = [];allpuls=[];
        elseif isScanner == 3   %Trio
            [D, TRsim, alltmr, allpuls] = makeMRFdictionaryAll_Trio_s(RFpulses, ...
                RefDic ,TR , TE, NumOfSlices, SL);
        elseif isScanner == 4   %Maastricht
            [D, TRsim, alltmr, allpuls] = makeMRFdictionary_EPI_T2Star_Barcelona_s(RFpulses, ...
                RefDic ,TR , TE, NumOfSlices, SL, iPat, SMSFactor);
        end
        
        % show the simulated rf pulses, only available for Mannheim data...
        if isScanner == 1 || isScanner == 4
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
        end
        
        Dictionary = abs(D);

        D = single(abs(D));
        % LUT values
        ttt = sqrt(sum(D.^2, 2));
        %take care here from the weird slice ordering with iPat
        if isScanner == 2
            SliceOrder =[2 4 6 8 10 20 12 14 16 18 3 5 7 9 11 1 13 15 17 19];
        elseif isScanner == 4
            SliceOrder = [2:2:NumOfSlices 1:2:NumOfSlices];
        else
            SliceOrder = [2:2:NumOfSlices 1:2:NumOfSlices];
        end
        Anatomical = SliceOrder(SL);
        
        % save all LUT to txt files in the folder defined in the beginning
        sname2 = strcat(svpath, '\LUTFactor_', num2str(Anatomical), '.txt');
        fileID = fopen(sname2,'w');
        for i = 1:size(ttt,1)
            fprintf(fileID, '%f ', ttt(i,1));
            fprintf(fileID, '\r\n');
        end
        fclose(fileID);
        
        
        D = D./repmat(sqrt(sum(D.^2, 2)), 1, meas); % normalize dictionarz
        c = horzcat(RefDic',D);
        allDic(:,:,Anatomical) = c;

        % write the LUTs
        sname = strcat(svpath, '\LUT_', num2str(Anatomical), '.txt');
        fileID = fopen(sname,'w');
        for i = 1:size(c,1)
            for ii = 1:size(c,2)
                fprintf(fileID, '%f ', c(i,ii));
            end
            fprintf(fileID, '\r\n');
        end
        fclose(fileID);
    end


else% otherwise load the LUT if it already exists!!!
    % loading LUT
    
    fprintf('%s Starting loading LUTs \n', datestr(now));
    clear allnc2 allDic;
   
    for SL = 1:NumOfSlices
        fprintf('%s loading LUT_%d \n', datestr(now),SL);
%         svpath = strcat(mainp,LutName, num2str(NumOfSlices));
        name = strcat(svpath,'\LUT_',num2str(SL),'.txt');
        c=dlmread(name); 
        allDic(:,:,SL) = c;
    end

    fprintf('%s Finished loading LUTs \n', datestr(now));

end
MRFBilderNew = MRFBilder;






%% Match the D with the data if maps should be generated
% idea: match the D with the dicionary (so no real measurement) to
% analyze a wide specrum of T1 and T2* values

% perform some denoising here if you want
MRFBilderNew = MRFBilder;
% figure;imagesc(MosaicOnOff(max(MRFBilder,[],4)))
% figure;plot(squeeze(MRFBilderNew(151,82,1:31,:))');colororder(jet(40));

figure;hold on;
theSlc = 1;
for theSlc = 1:1:31
    clf;
    Dictionary = abs(allDic(:,5:35+4,theSlc));
    D = Dictionary;    
    D = single(abs(D));
    D = D./repmat(sqrt(sum(D.^2,2)), 1, 35);
    M = reshape(squeeze(MRFBilderNew(:,:,theSlc,1:35)),[],35);  
    M = single(abs(M));
    M = M./repmat(sqrt(sum(M.^2,2)), 1, 35);
    M = reshape(M, d1, d2, 35);
    
    fill([1:35 35:-1:1],[min(abs(D(:,1:35)),[],1) max(abs(D(:,35:-1:1)),[],1)],'b');
    hold on;plot(squeeze(M(151,82,:))','k-','LineWidth',2);ylim([0 1]);
    pause(0.5);
end
% apply some denoising if you want:

MPPCA = false;

if MPPCA
    [MRFBilderNew] = MPdenoising(MRFBilderNew);
end

isDicom = false;         %save as Dicom
tic;

% make the matching in the traditional way since group matching is not much
% fast in Matlab and SVD also not since we only have 35 measurements

saveMaps = true;
fprintf('%s Matching the dicitionary on the traditional way\n', datestr(now))
tic
[d1,d2,d3,d4] = size(MRFBilderNew);
T1 = zeros(d1*d2,d3);
T2 = T1;
FitR = T1;
B1 = T1;

NTR = meas;
% for NTR=1:1:160
for loop = 1:NumOfSlices*SMSFactor
    fprintf('%s Matching slice %d of %d\n', datestr(now), loop, NumOfSlices);
    % only match if the relevant slice is there
    M = reshape(squeeze(MRFBilderNew(:,:,loop,1:NTR)),[],NTR);   
    
    % Use magnitude of dictioanry
    theSLC = 1+mod(loop-1,NumOfSlices);
%     [~,theSLC] = find(SliceOrder==loop);
    Dictionary = abs(allDic(:,5:NTR+4,theSLC));
    RefDic = allDic(:,1:4,theSLC)';
    D = Dictionary;    
    
    [dd1,dd2] = size(M);
    
    D = single(abs(D));
    D = D./repmat(sqrt(sum(D.^2,2)), 1,dd2);

    M = single(abs(M));
    M = M./repmat(sqrt(sum(M.^2,2)), 1, dd2);
    
    % calculate max memory (a bit of a gess, but better than nothing)
    CorrelationCoeff = zeros(dd1,1);
    MaxCorrVal = ones(dd1,1);

    [MaxCorrVal, CorrelationCoeff] = MatchMRF(D(:,:), M(:,:));
    T1(:,loop) = RefDic(1,MaxCorrVal);
    T2(:,loop) = RefDic(2,MaxCorrVal);
    B1(:,loop) = RefDic(4,MaxCorrVal);
    FitR(:,loop) = CorrelationCoeff;
end
timetra = toc;
fprintf('%s Traditional way took: %0.2f\n', datestr(now), timetra)
T1Map = reshape(T1, d1,d2,d3);
T2Map = reshape(T2, d1,d2,d3);
B1Map = reshape(B1, d1,d2,d3);
FitMap = reshape(FitR, d1,d2,d3);


tresh = 0;
mask = max(MRFBilderNew(:,:,:,:),[],4);
mask(mask<tresh) = 0;mask(mask>0) = 1;
T1Map = T1Map.*mask;
T2Map = T2Map.*mask;
B1Map = B1Map.*mask;
FitMap = FitMap.*mask;
%     FitMap(FitMap<0.98) = 0;
% myMaps(:,:,:,NTR) = T2Map;
% myMaps2(:,:,:,NTR) = T1Map;
% myMaps3(:,:,:,NTR) = B1Map;
% myMaps4(:,:,:,NTR) = FitMap;
% end
    
figure;imagesc(MosaicOnOff(T1Map),[0 2000]);colormap;colorbar;
figure;imagesc(MosaicOnOff(T2Map),[0 120]);colormap;colorbar;
figure;imagesc(MosaicOnOff(B1Map),[0.5 1.5]);colormap;colorbar;
figure;imagesc(MosaicOnOff(FitMap),[0.98 1]);colormap;colorbar;

%     T1Map = T1Map.*FitMap;
%     T2Map = T2Map.*FitMap;
%     B1Map = B1Map.*FitMap;

toc
if isDicom
    for i=1:1:NumOfSlices
    %     idxAnatomical = find(SliceTiming == i);
        metadata = head(i,1);
        metadataAcquistionNumber = i;
        metadata.SeriesDescription = ['T1'];
        metadata.WindowCenter = 1000;
        metadata.WindowWidth = 1000;
        metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'1'];
        mkdir([mainp,metadata.SeriesDescription]);
        metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
        dicomwrite(squeeze(T1Map(:,:,i))./256^2, [mainp,metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);


        metadata = head(i,1);
        metadata.SeriesDescription = ['T2'];
        metadata.WindowCenter = 100;
        metadata.WindowWidth = 100;
        metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'2'];
        mkdir([mainp,metadata.SeriesDescription]);
        metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
        dicomwrite(squeeze(T2Map(:,:,i))./256^2, [mainp,metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);


        metadata = head(i,1);
        metadata.SeriesDescription = ['B1'];
        metadata.WindowCenter = 1000;
        metadata.WindowWidth = 500;
        metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'3'];
        mkdir([mainp,metadata.SeriesDescription]);
        metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
        dicomwrite(squeeze(B1Map(:,:,i).*1000)./256^2, [mainp,metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);

        metadata = head(i,1);
        metadata.SeriesDescription = ['Fit'];
        metadata.WindowCenter = 1000;
        metadata.WindowWidth = 100;
        metadata.SeriesInstanceUID = [metadata.SeriesInstanceUID,'4'];
        mkdir([mainp,metadata.SeriesDescription]);
        metadata.Filename = [metadata.SeriesDescription,'slc',num2str(round(40+metadata.SliceLocation))];
        dicomwrite(squeeze(FitMap(:,:,i).*1000)./256^2, [mainp,metadata.SeriesDescription,'\',metadata.Filename,'.dcm'], metadata);
    end
end


%% make some plots to understand the data
figure;

dataVec = squeeze(D(1:29,:))';
refVec = squeeze(RefDic(2,1:29));
subplot(1,2,1);plot(dataVec);colororder(parula(64));
c = colorbar;c.Label.String = 'Increasing T_1';

dataVec = squeeze(D(1:29:1470,:))';
refVec = squeeze(RefDic(2,1:29:1470));
subplot(1,2,2);plot(dataVec);colororder(parula(64));
c = colorbar;c.Label.String = 'Increasing T_2*';
% c.Limits = [refVec(1) refVec(end)];

