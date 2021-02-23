% T2 Map Reconstruction by Ingo Hermann, 2020-10-08
% Reconstruction for T2 Map for a 60 slice and three weigthing T2 weigthed
% sequence if all images are in the following folder
% This script was made for the MRF-EPI MS study to reconstruct from the 3
% PD-TSE T2 weighted images a T2Map
% --------------------------------
% This scripts needs the user functions:
% MosaicOnOff.m
% openMaps.m
% T2Fit.m

Outfolder = 'G:\Studien\MRF\Messdaten\';
allSub = ls(Outfolder);
delNames = zeros(size(allSub,1),1);
for count = 1:1:size(allSub,1)
    tmpName = strtrim(allSub(count,:));
    folder = strtrim(allSub(count,:));
    if exist([Outfolder,folder,'\T2Map_Base\'] ,'dir')~=7 || ~contains(folder,'Healthy')  || ~contains(folder,'notFinished') %||...
%         exist([Outfolder,folder,'\T2Map_quant\'] ,'dir')==7% ||...
%         ~contains(tmpName,'MRF') 
        delNames(count) = count;
    end
end
delNames(delNames==0) = [];
allSub(delNames,:) = [];


for bigLoop = 1:1:size(allSub,1)
    folder = strtrim(allSub(bigLoop,:));

    fprintf('%s, start loading images from %s\n',datestr(now),folder);
    slcs = 60;
    for i=1:1:slcs*3
        [img, tmp] = openMaps([Outfolder,folder,'\T2Map_Base\*'],'num',i);
        TE(i) = tmp.EchoTime;
        SN(i) = tmp.SeriesNumber;
        SL(i) = tmp.SliceLocation;
        AN(i) = tmp.AcquisitionNumber;
        IN(i) = tmp.InstanceNumber;
        heads(i,1) = tmp;
        allImages(:,:,i) = img(9:end-8,9:end-8);
    end
    [~, order] = sort(SL);
    allImagesOrder = allImages(:,:,order);
    sz = size(allImagesOrder);
    allImage = reshape(allImagesOrder,sz(1),sz(2),slcs,3);
%     data = zeros(1920,1920,3);
    data(:,:,1) = MosaicOnOff(squeeze(allImagesOrder(:,:,1:3:end)),slcs);
    data(:,:,2) = MosaicOnOff(squeeze(allImagesOrder(:,:,2:3:end)),slcs);
    data(:,:,3) = MosaicOnOff(squeeze(allImagesOrder(:,:,3:3:end)),slcs);

    fprintf('%s, finished loading images from %s\n',datestr(now),folder);
    % get all the TEs and the image size and the max and min of every slice 
    weight1 = MosaicOnOff(squeeze(data(:,:,1)),slcs);
    weight2 = MosaicOnOff(squeeze(data(:,:,2)),slcs);
    weight3 = MosaicOnOff(squeeze(data(:,:,3)),slcs);
    
    TEs = unique(TE);
    ps = size(data,1);
    maxData = max(data,[],3);
    minData = min(data,[],3);
    diff = maxData-minData;

    % define difference Vector to speed up reconstruction time 
    res = zeros(ps*ps,1);   
    T2 = zeros(ps*ps,3);
    dataVec = reshape(data,ps*ps,3);
    diffVec = reshape(diff,ps*ps,1);
    weight1Vec = dataVec(:,1);
    weight2Vec = dataVec(:,2);
    weight3Vec = dataVec(:,3);

    % The fit whic can take some time ~ 1500 seconds
    % parameters for the lsqcurve fit
    opts = optimset('lsqcurvefit');
    opts = optimset(opts,'Display','off'); 
    opts = optimset(opts,'Algorithm','levenberg-marquardt');
    opts = optimset(opts,'MaxIter',100);

    pixelTresh = 40;

    tic;
    parfor i=1:1:ps*ps
        if mod(i,ps*ps/100)==0
            fprintf('%s, %0.0f percent fitted\n',datestr(now),round(i*100/(ps*ps),0));
        end
        if diffVec(i)<pixelTresh
            continue;
        else
            [T2(i,:), ~, res(i)] = T2_Fit(TEs', dataVec(i,:)', opts, 'Parameter',2);
        end
    end
    toc
    T2Map = reshape(T2(:,3),ps,ps);
    resMap = reshape(res,ps,ps);
    T2Map(T2Map>4000) = 4000;
    T2Map(T2Map<0) = 0;
    mkdir([Outfolder,folder,'\T2Map_quant\']);
    T2s = MosaicOnOff(T2Map,slcs);
    
    % this save every single images as a single dicom !!!
    % I might use this more often :)
    mkdir([Outfolder,folder,'\T2Map_quant_single\']);
    for i=1:1:slcs
        tmp = heads(i,1);
        tmp.SeriesDescription = 'T2Map';
        tmp.WindowCenter = 100;
        tmp.WindowWidth = 200;
        tmp.SeriesInstanceUID = [tmp.SeriesInstanceUID,'3'];
        dicomwrite(uint16(squeeze(T2s(:,:,i))),[Outfolder,folder,'\T2Map_quant_single\T2Map_slc',num2str(SL(i)),'.dcm'],tmp);
    end
    dicomwrite(uint16(T2Map),[Outfolder,folder,'\T2Map_quant\T2Map.dcm']);
    figure;imagesc(T2Map,[0 200]);
    save([Outfolder,folder,'\T2Map.mat'],'T2Map');
    
    
    theHead = tmp(1,1);
    tmpT2 = MosaicOnOff(allImages(:,:,1:60));
    mkdir ([Outfolder,folder,'\T2']);
    Tname = [Outfolder,folder,'\T2\',folder,'_T2.dcm'];
    dicomwrite(tmpT2./256^2, Tname, theHead);
    dicm2nii(Tname, [Outfolder,folder,'\'], '.nii.gz',[folder,'_T2'])

end % big loop