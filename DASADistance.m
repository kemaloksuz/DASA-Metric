function [ DASA,LocErr,Type1Err,Type2Err] = DASADistance(GroundTruth, EstimatedStates, c, p)

%DEFINITION: DASA Distance function computes the DASA error along with its components 
%            at each time interval between the Ground Truth and Estimated States.

%INPUT: GroundTruth                 : MXKXN true realization of target states.
%                                     M is the dimension of state space
%                                     K is the total time interval
%                                     N is maximum number of targets at some
%                                     time interval
%                              
%       Domain of Ground Truth      : If there exists a target it is the state of the target; else -1
%
%       Estimated States            : MXKXN' estimated realization of target states.
%                                     N' is maximum number of estimated targets
%
%       Domain of Estimated States  : If there is an estimation assign the estimation; else -1
%       
%
%       c                           : Cut-off length
%       Domain and Size of c        : If c is scalar, then penalty is same
%       for all types of errors. If c is a 4X1 vector, then corresponding
%       penalty for each component is used. Then, the structure is as follows:
%                               c= [c Loc, c alpha, c beta, upperbound]'
%
%       p                           : DASA Metric Order, it is a scalar.
%                                  
%
%OUTPUT:DASA            : 1XK vector consisting of Total DASA errors.
%       LocErr          : 1XK vector consisting of DASA Localization errors.
%       Type1Err        : 1XK vector consisting of DASA Type 1 errors.
%       Type2Err        : 1XK vector consisting of DASA Type 2 errors.

%EXAMPLE:
%       
% GroundTruth(:,:,1) =
% 
%   362.0623  362.6301  362.9970  363.5144  364.0466  364.3454  364.9183  365.5057  365.9434  366.4977
%  -455.8079 -452.1482 -448.5619 -445.0821 -441.5779 -437.9930 -434.5992 -430.9577 -427.3920 -423.6721
% 
% 
% GroundTruth(:,:,2) =
% 
%   485.9775  480.8408  475.8439  470.8610  465.7639  460.8014  456.1248  451.7636  447.1491  442.6162
%   456.5848  458.8090  461.3339  463.6244  466.0539  468.4297  470.6152  472.7711  474.9831  477.3244
% 
% 
% GroundTruth(:,:,3) =
% 
%  -336.2130 -336.7359 -337.7273 -338.7342 -339.7409 -340.6016 -341.3556 -342.3669 -343.0680 -343.8559
%   -64.3663  -60.0341  -55.6847  -51.2141  -46.7820  -42.5332  -38.5079  -34.6482  -30.4335  -26.4386
% 
% 
% GroundTruth(:,:,4) =
% 
%    -1.0000   -1.0000 -278.9907 -280.9238 -282.8936 -285.0170 -287.0806 -288.8406 -290.6403 -292.6236
%    -1.0000   -1.0000 -244.1233 -245.8744 -247.6065 -249.5122 -251.4246 -253.5385 -255.4720 -257.3146
%
% EstimatedStates(:,:,1) =
% 
%   363.0000  478.9444 -338.5653 -339.0986 -339.6671 -341.9257 -341.9198 -345.1394 -344.9768 -345.7858
%  -456.0000  459.8737  -55.7887  -50.5572  -45.8694  -41.9101  -39.0135  -33.5995  -29.4367  -25.5599
% 
% 
% EstimatedStates(:,:,2) =
% 
%   484.0000 -336.9853  364.5326  365.8345  366.4009 -634.5137  467.7418  462.3439  367.5562  367.6631
%   457.0000  -59.6003 -445.6498 -442.2585 -436.0311  463.2302  465.8146  466.1866 -426.3090 -423.2922
% 
% 
% EstimatedStates(:,:,3) =
% 
%  -336.0000  363.8509   -1.0000  504.7148   -1.0000  366.7831   -1.0000  366.2407  454.5981 -288.7736
%   -64.0000 -450.6612   -1.0000  270.7727   -1.0000 -434.1454   -1.0000 -427.9305  468.0099 -258.5361
% 
% 
% EstimatedStates(:,:,4) =
% 
%    -1.0000   -1.0000   -1.0000   -1.0000   -1.0000  463.5634   -1.0000 -285.2318 -702.4216  446.9392
%    -1.0000   -1.0000   -1.0000   -1.0000   -1.0000  470.2829   -1.0000 -254.8227 -388.1714  470.3628
% 
% 
% EstimatedStates(:,:,5) =
% 
%    -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000 -287.2036   -1.0000
%    -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000 -257.7218   -1.0000
% p=2
% c=10
%
% Total DASA Error  :[1.3139,1.7012,7.2724,7.9236,7.6998,6.8403,8.9507,6.8288,7.3630,4.7293]
% Localization Error:[1.3139,1.7012,2.4033,2.6382,4.3098,3.3641,0.7576,3.3248,2.9031,4.7293]
% Type 1 Error      :[0,0,0,5.0000,0,5.0000,5.0000,5.0000,6.3246,0]
% Type 2 Error      :[0,0,7.0711,7.0711,7.0711,5.0000,8.6603,5.0000,4.4721,0]
%
% If c=[10,10,20,10]'
%
% Total DASA Error  :[1.3139,1.7012,9.0086,9.1150,9.1496,8.1697,9.6384,8.1637,8.3352,4.7293]
% Localization Error:[1.3139,1.7012,2.4033,2.6382,4.3098,3.3641,0.7576,3.3248,2.9031,4.7293]
% Type 1 Error      :[0,0,0,5.0000,0,5.0000,5.0000,5.0000,6.3246,0]
% Type 2 Error      :[0,0,14.1421,14.1421,14.1421,10.0000,17.3205,10.0000,8.9443,0]
%
% by Kemal OKSUZ at Bogazici University on 28th May 2016
%                           revised on 31th October 2016
%
% For more questions, mail to kemal.oksz@gmail.com
%
% REFERENCE: 
% For Hungarian Algorithm, we use Yi Cao implementation that may be found in:
% http://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3-
% This is the implementation of "Munkres Assignment Algorithm, Modified for Rectangular Matrices", 
% http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html

    %Check the structure of c parameter and decide where all types of
    %errors are penalized equally or not.
    [flag, ~]=size(c);
    if flag==1
        c=ones(4,1)*c;
    elseif flag~=4
        error('c must be either a scalar of a 4x1 vector')
    end
    
    %Find the total length of the process
    [~,K,~]=size(GroundTruth(1,:,1));
    
    %Initialize Data Structures
    DASA=zeros(1,K);
    LocErr=zeros(1,K);
    Type1Err=zeros(1,K);
    Type2Err=zeros(1,K);
    
    %Initialize delta that is a small value
    delta=10e-4;
    
    %For each time interval calculate the DASA error in this loop
    for i=1:K
        %Find the number of truth states and corresponding  
        %indexes for ith period
        [~,~,m]=size(GroundTruth(1,i,:));
        Tracks=[];
        counter=0;
        for j=1:m
            if (GroundTruth(1,i,j)~=-1)
                Tracks=[Tracks, j];
                counter=counter+1;
            end
        end  
        m=counter; %m is the number of Ground truths
        
        %Find the number of estimated states and corresponding  
        %indexes for ith period        
        [~,~,n]=size(EstimatedStates(1,i,:));
        EstimatedTracks=[];
        counter=0;        
        for j=1:n
            if (EstimatedStates(1,i,j)~=-1)
                EstimatedTracks=[EstimatedTracks, j];
                counter=counter+1;
            end
        end
        n=counter; %n is the number of Estimated targets
        
        %Initialize mhat and nhat. This corresponds to 1st step of the
        %computation algorithm in the paper.
        if m>n
            mhat=m-n;
            nhat=0;
        else
            nhat=n-m;
            mhat=0;
        end
        
        %2nd step of the computation algorithm:Initialize l to be the maximum of m and n.
        l=max(m,n);
        
        %3rd step of the computation algorithm: Initialize cost matrix and
        %assign the appropriate costs. Notice that for d(x,y) Euclidean
        %distance is used.
        C=zeros(n,m);
        for j=1:n
            for k=1:m
                C(j,k)=(norm(EstimatedStates(:,i,EstimatedTracks(1,j))-GroundTruth(:,i,Tracks(1,k))));
                if C(j,k)>c(1,1)
                    C(j,k)=c(1,1)+delta;
                end

            end
        end       
        
        %4th step of the computation algorithm: Run Hungarian Algorithm in
        %order to find out the best assignment between ground truth and
        %estimated states sets
        [assignment,~] = Hungarian(C);
        
        %5th step of the computation algorithm: We determine the number of
        %association that is more distant than c of localization. Then
        %increment nhat and mhat and find out omega.        
        [~,assignnum]=size(assignment);
        loc=0;
        omega=0;
        for j=1:assignnum
            if assignment(1,j) ~= 0
                if C(j,assignment(1,j))<=c(1,1)
                    omega=omega+1;
                    loc=loc+C(j,assignment(j))^p;
                else
                    mhat=mhat+1;
                    nhat=nhat+1;
                end
            end
        end
        
        %Calculate the Localization Component if it is defined. 
        if omega==0
            LocErr(1,i)=NaN;
        else
            LocErr(1,i)=(loc/omega)^(1/p);          
        end
        
        %In order to catch the error in which the localization error is
        %undefined we also calculate scaled localization component of
        %Equation 9. We directly exploit this value while calculating the
        %total DASA value below.
        MetricDASA=(loc/l)^(1/p);
        
        % Calculate Type 1 and Type 2 Error Components
        Type1Err(1,i)=c(2,1)*((nhat/l)^(1/p));
        Type2Err(1,i)=c(3,1)*((mhat/l)^(1/p));  
        
        %Calculate Unnormalized Total DASA Error
        DASA(1,i)=(MetricDASA^p+Type1Err(1,i)^p+Type2Err(1,i)^p).^(1/p);
        
        %Check whether the flexible adaptation of the metric is used and
        %calculate the normalization constant.
        if flag==1
            Z=(((m+nhat)/l)^(1/p));
        else
            Z=(((omega*(c(1,1)^p)+nhat*(c(2,1)^p)+mhat*(c(3,1)^p))/l)^(1/p))/c(4,1);
        end
        
        %Calculate the total DASA value by normalizing it.
        DASA(1,i)=DASA(1,i)/Z;
    end
end

