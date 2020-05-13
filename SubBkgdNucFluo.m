%% Do background subtraction for nuclear channel 

% Want to go thru each time frame and take the avg fluo intensity of the
% image where nuclei are not present, and then subtract this value from the
% nuclear fluo intensity. These corrected values will be saved in the
% CompiledNuclei structure for later use. 

% load file location info
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%%
% Go thru each time frame 
for tt = 1:length(ElapsedTime) 
    %At each time frame go thru and get the pixel x,y coordinates of each
    %nucleus
    AllCents=[];
    AllRadii=[];
    AllNucMasks=0;
    for nn=1:length(CompiledNuclei) 
        FramesExist = [CompiledNuclei(nn).Frames];
        if ismember(tt,FramesExist) %if nucleus exists at this time point
            Idx = find(FramesExist == tt);
            NucCent = [CompiledNuclei(nn).xPos(Idx),CompiledNuclei(nn).yPos(Idx)];
            AllCents = [AllCents;NucCent];
            AllRadii = [AllRadii; CompiledNuclei(nn).Radius(Idx)];
        end
    end
    %open the nuclear channel image at this time point and subtract out the
    %nuclear regions. Need to make circles centered on the coordinates.
    %getmidpointcircle is downloaded from Matlab filexchange. It uses the
    %mid-point circle algorithm. 
    OrigIm = imread([PreProcPath filesep Prefix filesep Prefix '-His_' iIndex(tt,3) '.tif']);
    %get size of image - this should be the same for each time point but 
    [m,n] = size(OrigIm);
    % go thru each nucleus present at this time frame, subtract from the
    % original image using imsubtract
    for nn=1:size(AllCents,1)
        %Get nuclear circumference
        [x,y] = getmidpointcircle(AllCents(nn,1),AllCents(nn,2),AllRadii(nn));
        %create mask of filled in circle  
        NucMask = poly2mask(x,y,m,n);
        
        AllNucMasks = AllNucMasks + NucMask; %image of all nuclei
    end
    %Logical index to remove where nuclei are 
    BkgdIm = OrigIm; 
    %If nucleus is present, is not background
    BkgdIm(AllNucMasks==1) = 0; 
    %Just testing this rn 
    BkgdFluo(tt) = mean(mean(BkgdIm));
end

%For now take avg background across the whole movie and subtract this at
%each time point from each nucleus
AvgBkgd = mean(BkgdFluo);
for nn = 1:length(CompiledNuclei)
    CompiledNuclei(nn).FluoBkgdSub = [CompiledNuclei(nn).FluoTimeTrace - AvgBkgd];
    CompiledNuclei(nn).BkgdSubVal = AvgBkgd;
end

%Save the updated CompiledNuclei structure 
% save('CompiledNuclei','CompiledNuclei');