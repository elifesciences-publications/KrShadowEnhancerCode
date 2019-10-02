%% Create SpotDiff structure to record fluorescence data grouped by the corresponding nucleus 
%% nc14 find spots associated with each nucleus
% Get spots specifically from nc14
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
nc14time=(max([CompiledParticles_14.NucEnd]))-nc14;
SpotCount=[];


% Use average AP position of spots to put into a discrete AP bin
for tt=1:length(CompiledParticles_14)
    APEstm(tt)=round(CompiledParticles_14(tt).MeanAP,3);
    for jj=1:length(APbinID)
        if APEstm(tt) < APbinID(jj)
            CompiledParticles_14(tt).APBin=APbinID(jj);
            break;
        end
    end
end 

%Find the spots associated with each nucleus 
    for ii=1:length(nuclei)
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    CompiledParticles_14(ii).FullFrames=zeros(1,nc14time);  %make field as long as number of frames in nc14
    if length(temp)==1
        SpotCount(ii,1)=1;
    end
    
    
    %need to only look at frames both exist in 
    
    if length(temp)==2
        SpotCount(ii,2)=1;  
    end
    
        
   end

%% Compile information about the two alleles per nucleus across the embryo during nc14 into SpotDiff structure  

 LastFrame=max([CompiledParticles_14.NucEnd]);
 nc14frames=[nc14:LastFrame];
 
for ii=1:length(nuclei)      
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));  
    CompiledParticles_14(temp(1)).nc14Frame=zeros(1,length(ElapsedTime));       %make frame value sit at that index position in the .nc14Frames structure
    %Entering the frames where the 1st spot of current nucleus has
    %transcription
    for qq=1:length(CompiledParticles_14(temp(1)).Frame)
        CompiledParticles_14(temp(1)).nc14Frame(CompiledParticles_14(temp(1)).Frame(qq))=CompiledParticles_14(temp(1)).Frame(qq);
    end
    counter=0;      %Also make nc14Fluo field same length as nc14Frames so can have 0 values for frames where schintz (nucleus) exists but not fluor 
    CompiledParticles_14(temp(1)).nc14Fluo=nan(1,length(ElapsedTime));
    
    %Go thru the frames where there was fluorescence and enter those values
    %in the new .nc14Fluo field at the column # of that frame
    for jj=1:length(CompiledParticles_14(temp(1)).nc14Frame)
        if any(CompiledParticles_14(temp(1)).nc14Frame(jj))
            counter=counter+1;
            CompiledParticles_14(temp(1)).nc14Fluo(jj)=CompiledParticles_14(temp(1)).Fluo(counter);
        end
    end
       CompiledParticles_14(temp(1)).nc14Fluo(CompiledParticles_14(temp(1)).nc14Fluo==0)=nan; %Fill in current 0's as nan's so can later fill in frames where nucleus exists but no spot as 0s

    %Repeat for 2nd transcription spot of the nucleus (if existant)
    if length(temp)==2
        CompiledParticles_14(temp(2)).nc14Frame=zeros(1,length(ElapsedTime));       %make frame value sit at that index position in the .nc14Frames structure
    for zz=1:length(CompiledParticles_14(temp(2)).Frame)
        CompiledParticles_14(temp(2)).nc14Frame(CompiledParticles_14(temp(2)).Frame(zz))=CompiledParticles_14(temp(2)).Frame(zz);
    end
    counter2=0;
        
    CompiledParticles_14(temp(2)).nc14Fluo=nan(1,length(ElapsedTime));
             
    for jj=1:length(CompiledParticles_14(temp(2)).nc14Frame)
        
        if any(CompiledParticles_14(temp(2)).nc14Frame(jj)) %(~isempty(CompiledParticles_14(temp(2)).nc14Frame(jj))) & (CompiledParticles_14(temp(2)).nc14Frame(jj)~=0)
            counter2=counter2+1;
            CompiledParticles_14(temp(2)).nc14Fluo(jj)=CompiledParticles_14(temp(2)).Fluo(counter2);
        end
    end
    
    CompiledParticles_14(temp(2)).nc14Fluo(CompiledParticles_14(temp(2)).nc14Fluo==0)=nan;
               
    end
  
    
        for jj=1:length(nc14frames)
            tempa=[];tempb=[];tempa2=[];tempb2=[];
            %Find what frames in nc14 this nucleus exists in
            tempa=find([schnitzcells(CompiledParticles_14(temp(1)).Nucleus).frames]==nc14frames(jj));
            %Find what frames the 1st transcriptional spot of this nucleus
            %has fluo values in
            tempb=find([CompiledParticles_14(temp(1)).nc14Frame]==nc14frames(jj));

        %schnitz exists but no spot at that frame
            if (~isempty(tempa)) & (isempty(tempb))
                CompiledParticles_14(temp(1)).nc14Frame(nc14frames(jj))=1;
            end
        end
        %Indicate frames where the nucleus exists but no transcription spot
        for jj=1:length(CompiledParticles_14(temp(1)).nc14Frame)   %For nuclei where found exist in schnitz but no particle, have fluo value of 0 (other indices where no frame in compiledpar and not id-ed in schnitz should be nan's rn
            if CompiledParticles_14(temp(1)).nc14Frame(jj)==1
                CompiledParticles_14(temp(1)).nc14Fluo(jj)=0;
            end
        end
        %Repeat for 2nd transcriptional spot if 2 alleles
        if length(temp)==2
            for jj=1:length(nc14frames)
                tempa2=find([schnitzcells(CompiledParticles_14(temp(2)).Nucleus).frames]==nc14frames(jj));
       
                tempb2=find([CompiledParticles_14(temp(2)).nc14Frame]==nc14frames(jj));
                    if (~isempty(tempa2)) & (isempty(tempb2))
                        CompiledParticles_14(temp(2)).nc14Frame(nc14frames(jj))=1;
                    end 
            end
       

            for jj=1:length(CompiledParticles_14(temp(2)).nc14Frame)
                if CompiledParticles_14(temp(2)).nc14Frame(jj)==1
                    CompiledParticles_14(temp(2)).nc14Fluo(jj)=0;
                end
            end
        end 
        
       
        %Creating SpotDiff structure to record the nc14 fluor values by
        %each nucleus, along with other info about the nucleus
        DiffArray=[];
        for jj=1:length(CompiledParticles_14(temp(1)).nc14Frame)
%              
            DiffArray(1,jj)=CompiledParticles_14(temp(1)).nc14Fluo(jj); 
            SpotDiff(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;  
            SpotDiff(ii).OriginalParticle=CompiledParticles_14(temp(1)).OriginalParticle;
            %Note that fluorescence error is a single value for each spot
            SpotDiff(ii).Err1=CompiledParticles_14(temp(1)).FluoError*ones(size(CompiledParticles_14(temp(1)).nc14Frame)); %calculate error for each spot for later reference/plotting
            SpotDiff(ii).Frame1=[CompiledParticles_14(temp(1)).Frame];
            SpotDiff(ii).nc14Frame1=[CompiledParticles_14(temp(1)).nc14Frame];
            SpotDiff(ii).APBin=CompiledParticles_14(temp(1)).APBin;
            SpotDiff(ii).Nucleus=CompiledParticles_14(temp(1)).Nucleus;
            SpotDiff(ii).SpotOne=[DiffArray(1,:)];
            SpotDiff(ii).TotalmRNAOne=[CompiledParticles_14(temp(1)).TotalmRNA];
            % Added 8/28/19 to track AP position of particles
            SpotDiff(ii).xPos1 = [CompiledParticles_14(temp(1)).xPos];
            SpotDiff(ii).yPos1 = [CompiledParticles_14(temp(1)).yPos];
            
        % Add field of smoothed trace for noise calcs 1/22/19
            if length(ElapsedTime)==length(SpotDiff(ii).SpotOne)
                SpotDiff(ii).SmoothSpotOne=smooth(ElapsedTime,SpotDiff(ii).SpotOne,0.1,'lowess');
                for bb=1:length(SpotDiff(ii).SpotOne)
                    if isnan(SpotDiff(ii).SpotOne(bb)) %If nucleus doesn't exist, again turn that frame into a nan value
                        SpotDiff(ii).SmoothSpotOne(bb)=nan;
                    end
                    if SpotDiff(ii).SmoothSpotOne(bb) <0  %1/29/19 get rid of negs from smoothing
                         SpotDiff(ii).SmoothSpotOne(bb)=0;
                    end
                end
            else
                SpotDiff(ii).SmoothSpotOne=nan;
            end
        
        end
        %Record same values for 2nd allele if exists
        if length(temp)==2
            for qq=1:length(CompiledParticles_14(temp(2)).nc14Fluo)
                DiffArray(2,qq)=CompiledParticles_14(temp(2)).nc14Fluo(qq);
            end
        
            SpotDiff(ii).Original2ndParticle=CompiledParticles_14(temp(2)).OriginalParticle;

            SpotDiff(ii).Err2=CompiledParticles_14(temp(2)).FluoError*ones(size(CompiledParticles_14(temp(2)).nc14Frame));
            SpotDiff(ii).Frame2=[CompiledParticles_14(temp(2)).Frame];
            SpotDiff(ii).nc14Frame2=[CompiledParticles_14(temp(2)).nc14Frame];
            SpotDiff(ii).TotalmRNATwo=[CompiledParticles_14(temp(2)).TotalmRNA];
            SpotDiff(ii).SpotTwo=[DiffArray(2,:)];
            
            %Added 8/28/19 to keep track of X,Y location of particles
            SpotDiff(ii).xPos2 = [CompiledParticles_14(temp(2)).xPos];
            SpotDiff(ii).yPos2 = [CompiledParticles_14(temp(2)).yPos];
            
         %Add smoothed trace field for noise calcs 1/22/19
            if length(ElapsedTime)==length(SpotDiff(ii).SpotTwo)
                SpotDiff(ii).SmoothSpotTwo=[smooth(ElapsedTime,SpotDiff(ii).SpotTwo,0.1,'lowess')];
                for bb=1:length(SpotDiff(ii).SpotTwo)
                    if isnan(SpotDiff(ii).SpotTwo(bb))
                       SpotDiff(ii).SmoothSpotTwo(bb)=nan;
                    end
                    if SpotDiff(ii).SmoothSpotTwo <0
                       SpotDiff(ii).SmoothSpotTwo=0;
                    end
                end
           else
              SpotDiff(ii).SmoothSpotTwo=nan;
           end
         
           SpotDiff(ii).AbsDiffBtwn=[abs(DiffArray(1,:)-DiffArray(2,:))];
        
           SpotDiff(ii).SpotCorr=[corrcoef(DiffArray','Rows','complete')];
    end
end

%Get rid of all empty field values 
for ss=1:length(SpotDiff)
    if isempty(SpotDiff(ss).SpotOne)
        SpotDiff(ss).SpotOne=nan;
    end
    if isfield(SpotDiff, 'SpotTwo')
        if isempty(SpotDiff(ss).SpotTwo)
            SpotDiff(ss).SpotTwo=nan;
        end
        if isempty(SpotDiff(ss).Original2ndParticle)
           SpotDiff(ss).Original2ndParticle=nan;
        end
        if isempty(SpotDiff(ss).Err2)
           SpotDiff(ss).Err2=nan;
        end
        if isempty(SpotDiff(ss).Frame2)
           SpotDiff(ss).Frame2=nan;
        end
        if isempty(SpotDiff(ss).nc14Frame2)
           SpotDiff(ss).nc14Frame2=nan;
        end
        if isempty(SpotDiff(ss).AbsDiffBtwn)
           SpotDiff(ss).AbsDiffBtwn=nan;
        end
        if isempty(SpotDiff(ss).SpotCorr)
           SpotDiff(ss).SpotCorr=nan;
        end
    end
end

%% Spot Correlation
% this section finds the correlation of activity between the two alleles in
% a single nucleus. Not really needed as later scripts (ie Total noise via
% covariance values) look into this but here if you want it

for j=1:length(SpotDiff)
    if ~isempty(SpotDiff(j).MeanAP)
        
    APEstm(j)=round(SpotDiff(j).MeanAP,3);
    
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            SpotDiff(j).APBin=APbinID(jj);
            break;
        end
    end
    end
    
end

 for ii=1:length(SpotDiff)
      if isempty(SpotDiff(ii).APBin)
          SpotDiff(ii).APBin=0;
      end
  end

Bins=unique([SpotDiff.APBin]);%want to find mean correlation of spots at each APbin for fewer points/easier to see graph when combining
APCorr=zeros(1,length(APbinID));
if isfield(SpotDiff, 'SpotTwo')
for zz=1:length(Bins)
    temp3=find([APbinID==Bins(zz)]);
    temp2=find([SpotDiff.APBin]==Bins(zz)); 
    for qq=1:length(temp2)
        if ~isempty(SpotDiff(temp2(qq)).SpotCorr)
            if length(SpotDiff(temp2(qq)).SpotCorr)==1
                APCorr(qq,temp3)=[SpotDiff(temp2(qq)).SpotCorr(1)];
            else
    APCorr(qq,temp3)=[SpotDiff(temp2(qq)).SpotCorr(1,2)];
            end
        end
    end 
end
APCorr(APCorr==0)=nan;
MeanAPCorr=nanmean(APCorr,1);


end
% 
 %save('SpotCorrelationAdj','SpotDiff')
% save('MeanAPCorrelationAdj','MeanAPCorr')
% save('CompiledParticles_nc14', 'CompiledParticles_14')


