function [dict, TRsim, varargout]= makeMRFdictionary_EPI_T2Star_Barcelona_s(RFpulses, RefDic ,TR , TE, NumOfSlices, SliceNum, iPat, SMSFactor)

tmr = 25;

T1 = RefDic(1,:);
T2 = RefDic(2,:);
AngleCorrection = RefDic(4,:);


rf = abs(RFpulses'*AngleCorrection);  % acount for angle correction
% rph=angle(RFpulses'*AngleCorrection);

d1 = size(RefDic,2);    % number of fingerprints
d2 = length(RFpulses);  % length of fingerprints
dict = zeros(d1,d2);    % output dictionary
M0 = zeros(3,d1);       % temporarz magnetization vector, [mx,my,mz]
M0(3,:) = 1;

InvOffset = 20.64;
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



counter = 1;
if SMSFactor > 1
    for i = 1:SMSFactor*NumOfSlices
        if (mod(i-1, NumOfSlices)+1) == SliceNum
            M0 = rx_flip3D(M0, rf(1,:));
        end
        M0 = relx3D(M0, TR(1), T1, T2); % on average (NumOfSlices*SMSFactor) waiting time
    end
    M0(1:2,:) = 0;
end

% account for m_lSliceAccelRefScans
if iPat == 2
    k = 2;
elseif iPat == 3
    k = 2;
else
    k = 0;
end
for x1 = 1:(1 + k)
    % at the beginning there is a 20ms pause fore some reason...
   % M0 = relx3D(M0, 20, T1, T2);
   % tmr = tmr +20;
    for i = 1:NumOfSlices
        if (mod(i-1, NumOfSlices)+1) == SliceNum
            M0 = rx_flip3D(M0, rf(1,:));
            
            alltmr (counter) = tmr;
            allpuls (counter) = RFpulses(1);
            counter = counter + 1;   
            %fprintf('iPAT %d, Time %2.1f, FA  %2.0f \n', x1, tmr, rf(1,1)*180/pi)
        end
        M0 = relx3D(M0, TR(1), T1, T2); % on average (NumOfSlices*SMSFactor) waiting time
        tmr = tmr + TR(1);
        M0(1:2,:) = 0;
    end
end



% if we include a 10 sec pause after the ref scans, set magnetiazation to 1
% if SMSFactor > 1 && NumOfSlices > 10
%     M0(1:2,:) = 0;
%     M0(3,:) = 1;
% end

TRsim = zeros(d2,1);
ct = 0;
for i=1:d2
    % at the beginning there is a 20ms pause fore some reason...
    M0 = relx3D(M0, 20, T1, T2);
    ct = ct + 20;
    tmr = tmr + 20;
    for x1 =1:NumOfSlices 
        % global IR pulse is run before first slice in first measurement
        %if (x1 == 1 && i == 1)
        if(mod(i-1,2) == 0)
            if (floor((i-1) / 2)*4 == (x1-1))
                % 180 global inversion pulse
                M0 = relx3D(M0, 20.64/4, T1, T2);
                M0 = rx_flip3D(M0, pi);
                
                alltmr (counter) = tmr+InvOffset/4;
                allpuls (counter) = 1.7;
                counter = counter + 1;

                % relax after inverstion, it takes 22ms
                M0 = relx3D(M0, 20.64/4*3, T1, T2);
                ct = ct+20.64;
                %fprintf('IRPuls, Time %2.1f, Repetition %d, Slice %d \n', tmr + 5, i, x1);
                tmr = tmr + 20.64;
            end
        end
               
        
        if(x1 == SliceNum)
            %fprintf('Excitation %2.0f:, Time %2.1f, FA %2.1f, TE %2.1f, TR %2.1f\n', i, tmr, rf(i,1)*180/pi, TE(i), TR(i));
            tmr = tmr + TR(i);
            TRsim(i) = ct;
            ct=0;
            % if it's the slice we're interested in
            % excitation
            M0 = rx_flip3D(M0, rf(i,:));
            
            alltmr (counter) = tmr;
            allpuls (counter) = RFpulses(i);
            counter = counter + 1;
                
            % relax by TE
            M0 = relx3D(M0, TE(i), T1, T2);
            % readout
            dict(:,i)=(M0(1,:)+1i.*M0(2,:))';
            % time between TE and TR-End
            % relax
            M0 = relx3D(M0, TR(i) - TE(i), T1, T2);
            % spoil transversal magnatization
            M0(1:2,:) = 0;
            ct = ct + TR(i);
        
        else % (x1 == SliceNum)
            % if it's not the slice we're interested in, just relax
            M0 = relx3D(M0, TR(i), T1, T2);
            ct = ct + TR(i);
             tmr = tmr + TR(i);
        end % (x1 == SliceNum)
    end % x1 =1:NumOfSlices
varargout{1} = alltmr;
varargout{2} = allpuls;
varargout{3} = M0;
end