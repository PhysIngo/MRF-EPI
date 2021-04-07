function [dict, TRsim, varargout]= makeMRFdictionary_EPI_T2Star_Mannheim_s(RFpulses, RefDic ,TR , TE, NumOfSlices, SliceNum, iPat, SMSFactor, varargin)



T1 = RefDic(1,:);
T2star = RefDic(2,:);
AngleCorrection = RefDic(4,:);


rf = abs(RFpulses'*AngleCorrection);  % acount for angle correction
% rph=angle(RFpulses'*AngleCorrection);

d1 = size(RefDic,2);    % number of fingerprints
d2 = length(RFpulses);  % length of fingerprints
dict = zeros(d1,d2);    % output dictionary
M0 = zeros(3,d1);       % temporarz magnetization vector, [mx,my,mz]
M0(3,:) = 1;


% m_lInitialPrepScans-----m_lSliceAccelPrepScans----m_lPATPrepScans----------m_lSliceAccelDummyScans-------m_lAdjPrepScans--------Imaging----------
% m_lInitialDummyScans ---m_lSliceAccelRefScans								 m_lSliceAccelDummyScans																				   
% ----------n----------------------1---------------------1/0---------------------------m-------------------------1/0-------------------------------
% 																																		   
% --------dummy--------------------DAQ-------------------DAQ--------------------------dummy---------------------DAQ-----------------DAQ------------
% 																																		   
% --------all slices------------all slices------------all slices---------------reduced slices--------------reduced slices------reduced slices----
% 
% in case of SMS:
% m_lInitialDummyScans      = 0
% m_lSliceAccelRefScans     = 1 per slice (so if 4 slices and SMS 4 = 16)
% m_lPATPrepScans           = 1 per slice * (Grappa=3->3 Grappa=2->1 )
% m_lSliceAccelDummyScans   = 3sec/TR + 1 
% m_lAdjPrepScans           = 0



% if SMSFactor > 1
%     for i = 1:SMSFactor*NumOfSlices
%         if (mod(i-1, NumOfSlices)+1) == SliceNum
%             M0 = rx_flip3D(M0, rf(1,:), rph(1,:));
%         end
%         M0 = relx3D(M0, TR(1), T1, T2); % on average (NumOfSlices*SMSFactor) waiting time
%     end
%     M0(1:2,:) = 0;
% end

% account for m_lSliceAccelRefScans
if iPat == 2
    k = 0;
elseif iPat == 3
    k = 2;
else
    k = 0;
end
% for x1 = (1 + k)
% %     at the beginning there is a 20ms pause fore some reason...
%    M0 = relx3D(M0, 20, T1, T2);
%    tmr = tmr + 20;
%     for i = 1:NumOfSlices
%         if (mod(i-1, NumOfSlices)+1) == SliceNum
%             M0 = rx_flip3D(M0, rf(1,:), rph(1,:));
%             fprintf('iPAT %d, Time %2.1f, FA  %2.0f \n', x1, tmr, rf(1,1)*180/pi)
%         end
%         M0 = relx3D(M0, TR(1), T1, T2); % on average (NumOfSlices*SMSFactor) waiting time
%         tmr = tmr + TR(1);
%         M0(1:2,:) = 0;
%     end
% end

% PINGO 
counter = 1;
tmr = 29.3;
InvOffset = 20.64;
SROffset = 1;
beginPause = 20; 
if max(strcmp(varargin,'Pulses'))
    idx = 1 + find(strcmp(varargin,'Pulses'));
    AmountOfPulses = varargin{1,idx};
else
    AmountOfPulses = 0;
end
    
if NumOfSlices == 1
    repetitions = 0;
    tmr = 29.17;
    InvPulses = 1;
elseif NumOfSlices == 4
    tmr = 30.1;
    repetitions = 2;
    InvPulses = 2;
elseif NumOfSlices == 8
    tmr = 30.1;
    repetitions = 2;
    InvPulses = 2;
elseif NumOfSlices == 16
    tmr = 30.1;
    repetitions = 3;
    InvPulses = 4;
elseif NumOfSlices == 24
    tmr = 30.1;
    repetitions = 3;
    InvPulses = 6;
elseif NumOfSlices == 60
    repetitions = 2;
    InvPulses = 15;
else
    repetitions = 5;
    InvPulses = 15;
end

if max(strcmp(varargin,'Change'))
    idx = 1 + find(strcmp(varargin,'Change'));
    change = varargin{1,idx};
else
    change = 0;
end

if max(strcmp(varargin,'Repetitions'))
    idx = 1 + find(strcmp(varargin,'Repetitions'));
    repetitions = varargin{1,idx};
    k = 0;
end
% try to put 0.3 here but it's not working ...

if ~max(strcmp(varargin,'Prep'))
    repetition = 0;
end
tmr=tmr+0;
for loop = 1:1:(1 + k)+repetitions
    for puls = 1:1:NumOfSlices
%         modify to see all the pulses
        if puls == SliceNum
            M0 = rx_flip3D(M0, rf(1,:));
%             save time points
            alltmr (counter) = tmr;
            allpuls (counter) = RFpulses(1);
            counter = counter + 1;
        end
        M0 = relx3D(M0, TR(1), T1, T2star); % on average (NumOfSlices*SMSFactor) waiting time
        tmr = tmr + TR(1);
        M0(1:2,:) = 0;
        if (puls==NumOfSlices && puls>1)
            M0 = relx3D(M0, beginPause, T1, T2star); % on average (NumOfSlices*SMSFactor) waiting time
            tmr = tmr + beginPause;
        end
        M0(1:2,:) = 0;
    end
end

% if we include a 10 sec pause after the ref scans, set magnetiazation to 1
if SMSFactor > 1
    M0(1:2,:) = 0;
    M0(3,:) = 1;
end

if max(strcmp(varargin,'Kill'))
    M0(1:2,:) = 0;
    M0(3,:) = 0;
end
% Offset after inversion puls
TRsim = zeros(d2,1);
ct = 0;

% I include this to one loop below so we can skip all of this
if (1==0)
    if max(strcmp(varargin,'Sat'))
        M0 = rx_flip3D(M0, pi/2);
    %     M0(1:2,:) = 0;
        M0(3,:) = 0;
        % save time points
        alltmr (counter) = tmr+1/2;
        allpuls (counter) = 0.85;
        counter = counter + 1;
        % relax after saturation
        M0 = relx3D(M0, SROffset, T1, T2star);
        ct = ct+InvOffset;

        tmr = tmr + SROffset;
    else
        % +++ PINGO +++ INVERSION RECOVERY
        %180 global inversion pulse just at the beginning
        if max(strcmp(varargin,'Pause')) && change==0
            idx = 1 + find(strcmp(varargin,'Pause'));
            pauseDuration = varargin{1,idx};
            M0 = relx3D(M0, pauseDuration, T1, T2star);
            tmr = tmr + pauseDuration;
        end

        M0 = relx3D(M0, InvOffset/4, T1, T2star);
        if (change == 0)
            M0 = rx_flip3D(M0, pi);

            % save time points
            alltmr (counter) = tmr+InvOffset/4;
            allpuls (counter) = 1.7;
            counter = counter + 1;

        end

        % relax after inverstion, it takes 22ms
        M0 = relx3D(M0, InvOffset/4*3, T1, T2star);
        ct = ct+InvOffset;

        tmr = tmr + InvOffset;
        % +++ PINGO +++ INVERSION RECOVERY
    end
end


invPulsCC = 1;
for i=1:d2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% only for not 60 slices MRF so only KIDNEY and so on ... %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if  i > 1 && AmountOfPulses > 0
        if mod(i-1,AmountOfPulses)==0 && i<21
            if max(strcmp(varargin,'SR'))
                M0 = rx_flip3D(M0, pi/2);
    %             M0(1:2,:) = 0;
                M0(3,:) = 0;
                % save time points
                alltmr (counter) = tmr+1/2;
                allpuls (counter) = 0.85;
                counter = counter + 1;
                % relax after saturation
                M0 = relx3D(M0, SROffset, T1, T2star);
                ct = ct+InvOffset;
                tmr = tmr + SROffset;
            else
                if max(strcmp(varargin,'Spoil')) &&  i>16 
                    tmr = tmr + 1.0;
                    M0 = rx_flip3D(M0, pi/2);
                    M0(3,:) = 0;
                    % save time points
                    alltmr (counter) = tmr;
                    allpuls (counter) = 0.85;
                    counter = counter + 1;           
                    M0 = relx3D(M0, 6.3, T1, T2star);
                    M0(1,:) = 0;
                    M0(2,:) = 0;
                    tmr = tmr + 6.3;
                else
                    if max(strcmp(varargin,'IRP'))
                        if max(strcmp(varargin,'Pause'))
                            idx = 1 + find(strcmp(varargin,'Pause'));
                            pauseDuration = varargin{1,idx};
                            M0 = relx3D(M0, pauseDuration, T1, T2star);
                            tmr = tmr + pauseDuration;
                        end
                    end
                    % +++ PINGO +++ INVERSION RECOVERY
                    %180 global inversion pulse just at the beginning            
                    M0 = relx3D(M0, InvOffset/4, T1, T2star);
                    M0 = rx_flip3D(M0, pi);

                    % save time points
                    alltmr (counter) = tmr+InvOffset/4;
                    allpuls (counter) = 1.7;
                    counter = counter + 1;

                    % relax after inverstion, it takes 22ms
                    M0 = relx3D(M0, InvOffset/4*3, T1, T2star);
                    ct = ct+InvOffset;

                    tmr = tmr + InvOffset;
                    % +++ PINGO +++ INVERSION RECOVERY
                end
            end
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% only for the 60 slices MRF so only brain ms study ... %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if i==36
        M0 = zeros(3,d1);       % temporarz magnetization vector, [mx,my,mz]
        M0(3,:) = 1;
    end
    for x1 =1:NumOfSlices  
        % modify to see all the pulses
%         if (mod((i-1)*NumOfSlices+x1,4+2*NumOfSlices)==0) && InvPulses > 1 && ...
%                 invPulsCC <InvPulses && NumOfSlices > 15
        if i>35
            if (mod(i-1-35,2))==0
                if ( ((x1-1) == (round((i-1-35)/2)*4)  ))
                    M0 = rx_flip3D(M0, pi);
                    alltmr (counter) = tmr;
                    allpuls (counter) = 1.7;
                    counter = counter + 1;   
                    tmr = tmr + 15;
                    invPulsCC = invPulsCC+1;
                end
            end
        else
            
            if (mod(i-1,2))==0
                if (change + (round((i-1)/2)*4) > NumOfSlices)
                    nx1 = x1+NumOfSlices;
                else
                    nx1 = x1;
                end %|| (nx1-1) == (round((i-1)/2)*4-35)
                if ( ((nx1-1-change) == (round((i-1)/2)*4)  ) && invPulsCC<=NumOfSlices/4)
                    M0 = rx_flip3D(M0, pi);
                    alltmr (counter) = tmr;
                    allpuls (counter) = 1.7;
                    counter = counter + 1;   
                    tmr = tmr + 15;
                    invPulsCC = invPulsCC+1;
                end
            end
        end
        if(x1 == SliceNum)
            % save time points  
            alltmr (counter) = tmr;
            allpuls (counter) = RFpulses(i);
            counter = counter + 1;   
        end
        if(x1 == SliceNum)
            %fprintf('Excitation %2.0f:, Time %2.1f, FA %2.1f, TE %2.1f, TR %2.1f\n', i, tmr, rf(i,1)*180/pi, TE(i), TR(i));
            tmr = tmr + TR(i);
            TRsim(i) = ct;
            ct=0;
            % if it's the slice we're interested in
            % excitation
            M0 = rx_flip3D(M0, rf(i,:));
            
        
            % relax by TE
            M0 = relx3D(M0, TE(i), T1, T2star);
            % readout
            dict(:,i)=(M0(1,:)+1i.*M0(2,:))';
            % time between TE and TR-End
            % relax
            M0 = relx3D(M0, TR(i) - TE(i), T1, T2star);
            % spoil transversal magnatization
            M0(1:2,:) = 0;
            ct = ct + TR(i);

        else % (x1 == SliceNum)
            % if it's not the slice we're interested in, just relax
            M0 = relx3D(M0, TR(i), T1, T2star);
            ct = ct + TR(i);
            tmr = tmr + TR(i);
        end
        if (x1>0 && x1 == SliceNum)
            M0 = relx3D(M0, beginPause, T1, T2star); % on average (NumOfSlices*SMSFactor) waiting time
            tmr = tmr + beginPause;
        end
    end % (x1 == SliceNum)
end % x1 =1:NumOfSlices

varargout{1} = alltmr;
varargout{2} = allpuls;
varargout{3} = M0;
end
% --- PINGO ---

% for i=1:d2
%     % at the beginning there is a 20ms pause fore some reason...
%     M0 = relx3D(M0, 20, T1, T2);
%     ct = ct + 20;
%     tmr = tmr + 20;
%     for x1 =1:NumOfSlices 
%         % global IR pulse is run before first slice in first measurement
%         %if (x1 == 1 && i == 1)
%         if(mod(i-1,2) == 0)
%             if (floor((i-1) / 2)*4 == (x1-1))
%                 % 180 global inversion pulse
%                 M0 = relx3D(M0, InvOffset/4, T1, T2);
%                 M0 = rx_flip3D(M0, pi, rph(1,:));
%                 % relax after inverstion, it takes 22ms
%                 M0 = relx3D(M0, InvOffset/4*3, T1, T2);
%                 ct = ct+InvOffset;
%                 %fprintf('IRPuls, Time %2.1f, Repetition %d, Slice %d \n', tmr + 5, i, x1);
%                 tmr = tmr + InvOffset;
%             end
%         end
%                
%         
%         if(x1 == SliceNum)
%             %fprintf('Excitation %2.0f:, Time %2.1f, FA %2.1f, TE %2.1f, TR %2.1f\n', i, tmr, rf(i,1)*180/pi, TE(i), TR(i));
%             tmr = tmr + TR(i);
%             TRsim(i) = ct;
%             ct=0;
%             % if it's the slice we're interested in
%             % excitation
%             M0 = rx_flip3D(M0, rf(i,:), rph(i,:));
%             % relax by TE
%             M0 = relx3D(M0, TE(i), T1, T2);
%             % readout
%             dict(:,i)=(M0(1,:)+1i.*M0(2,:))';
%             % time between TE and TR-End
%             % relax
%             M0 = relx3D(M0, TR(i) - TE(i), T1, T2);
%             % spoil transversal magnatization
%             M0(1:2,:) = 0;
%             ct = ct + TR(i);
%         
%         else % (x1 == SliceNum)
%             % if it's not the slice we're interested in, just relax
%             M0 = relx3D(M0, TR(i), T1, T2);
%             ct = ct + TR(i);
%              tmr = tmr + TR(i);
%         end % (x1 == SliceNum)
%     end % x1 =1:NumOfSlices
% end