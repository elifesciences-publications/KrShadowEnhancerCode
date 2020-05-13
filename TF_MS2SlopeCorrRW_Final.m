%% Track mRNA production and TF-FP levels
%Code to load data sets once have more 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList = {'KrBoth_BcdGFP', 'KrDist_BcdGFP'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).Particles(1).APbinID];
    for ee=1:NEmbryos
        %load needed data
        PrefixName = Data(ee).Particles.Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        load(filename);
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        %Group data by AP bins
        for aa=1:length(APbinID)
            NucAP=[];
            RawNucAP=[];
            SpotAP=[];
            RawSpotAP=[];
            CorrAP=[];
            d_NucAP =[];
            d_SpotAP=[];
            
            ProductionAP=find([BurstProperties.APBin]==APbinID(aa));
            %If no nuclei at AP bin, fill nan's length of nc14 movie - all
            %will be 100 frames
            if isempty(ProductionAP)
                NucAP=[NucAP;nan(1,101)];
                RawNucAP=[RawNucAP;nan(1,101)];
                CorrAP = [CorrAP; nan];
                SpotAP=[SpotAP;nan(1,101)];
                RawSpotAP=[RawSpotAP;nan(1,101)];
                d_NucAP=[d_NucAP; nan(1,100)];
                d_SpotAP =[d_SpotAP; nan(1,100)];
                
            else
                for zz=1:length(ProductionAP)
                    if isempty(BurstProperties(ProductionAP(zz)).NucFluo)
                        NucAP = [NucAP; nan(1,101)];
                        RawNucAP=[RawNucAP; nan(1,101)];
                        SpotAP = [SpotAP; nan(1,101)];
                        RawSpotAP =[RawSpotAP; nan(1,101)];
                        CorrAP = [CorrAP; nan];
                        d_NucAP = [d_NucAP; nan(1,100)]; %Value needs to be 1 less bc when do diff, end up with vector w 1 fewer entries
                        d_SpotAP =[d_SpotAP; nan(1,100)];
                        
                    
                    elseif length(BurstProperties(ProductionAP(zz)).NucFluo) < 101
                    %Deal w some embryos being slightly less than 50
                    %minutes of tracking
                    
                        %BurstProperties(ProductionAP(zz)).SmoothNucFluo = [BurstProperties(ProductionAP(zz)).SmoothNucFluo, nan(1,(101-length(BurstProperties(ProductionAP(zz)).SmoothNucFluo)))]';
                        %Use photobleach corrected data 3/27/20
                        %Add in statement to make all oriented same way (ie
                        %row vs column vector) 
                        [m,n] = size(BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr)
                        %Make into row vector
                        if n==1
                            BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr = [BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr]';
                        end
                        BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr = [BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr, nan(1,(101-length(BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr)))]';
                        %pre photobleach corrected nuclear data 
                        [m3,n3] = size(BurstProperties(ProductionAP(zz)).SmoothNucFluo)
                        %Make into row vector
                        if n3==1
                            BurstProperties(ProductionAP(zz)).SmoothNucFluo = [BurstProperties(ProductionAP(zz)).SmoothNucFluo]';
                        end
                        BurstProperties(ProductionAP(zz)).SmoothNucFluo = [BurstProperties(ProductionAP(zz)).SmoothNucFluo, nan(1,(101-length(BurstProperties(ProductionAP(zz)).SmoothNucFluo)))]';
                        %do same for spot data
                        [m2,n2] = size(BurstProperties(ProductionAP(zz)).SmoothTrace)
                        if n2==1
                            BurstProperties(ProductionAP(zz)).SmoothTrace = [BurstProperties(ProductionAP(zz)).SmoothTrace]';
                        end
                        BurstProperties(ProductionAP(zz)).SmoothTrace = [BurstProperties(ProductionAP(zz)).SmoothTrace, nan(1,(101-length(BurstProperties(ProductionAP(zz)).SmoothTrace)))]';
                        
                        %Also do for raw spot data 
                        [m4,n4] = size(BurstProperties(ProductionAP(zz)).RawTrace)
                        if n4 > 1 %???
                            BurstProperties(ProductionAP(zz)).RawTrace = [BurstProperties(ProductionAP(zz)).RawTrace];
                        end
                        BurstProperties(ProductionAP(zz)).RawTrace = [BurstProperties(ProductionAP(zz)).RawTrace, nan(1,(101-length(BurstProperties(ProductionAP(zz)).RawTrace)))];
                    end
                %concatenate all nuclei in that AP bin to one long matrix
                    %Use the smoothed nuclear fluorescence
                    if ~isempty(BurstProperties(ProductionAP(zz)).NucFluo)
                        NucAP = [NucAP;[BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr]'];
                        RawNucAP = [RawNucAP;[BurstProperties(ProductionAP(zz)).SmoothNucFluo]'];
                    
                        SpotAP = [SpotAP;[BurstProperties(ProductionAP(zz)).SmoothTrace]'];
                        RawSpotAP = [RawSpotAP; [BurstProperties(ProductionAP(zz)).RawTrace]];
                        %Save d(fluo)/dt info
                        %Need to check length of Elapsed Time 
                        if (length(ElapsedTime) - nc14) < 100
                            ElapsedTime = [ElapsedTime,nan(1,(100-(length(ElapsedTime)-nc14)))];
                        end
                        d_NucAP = [d_NucAP;[diff([BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr])./([diff(ElapsedTime(nc14:(nc14+100)))]')]'];
                        d_SpotAP = [d_SpotAP;[diff(BurstProperties(ProductionAP(zz)).SmoothTrace)./([diff(ElapsedTime(nc14:(nc14+100)))]')]'];
                        
                    
                    %Go ahead and do the correlations here too so doesn't
                    %have to be done in plotting section 
                        CorrInfo = corrcoef([BurstProperties(ProductionAP(zz)).SmoothNucFluoBleachCorr],[BurstProperties(ProductionAP(zz)).SmoothTrace],'Rows','pairwise');
                        CorrInfo = CorrInfo(2);
                        CorrAP = [CorrAP; CorrInfo];
                    end
                end
                
            end
            NucFluo(cc).Embryo(ee).APBin(aa).AllNucTraces = NucAP;
            NucFluo(cc).Embryo(ee).APBin(aa).AllSpotTraces=SpotAP;
            NucFluo(cc).Embryo(ee).APBin(aa).AllRawSpotTraces = RawSpotAP;
            NucFluo(cc).Embryo(ee).APBin(aa).NucSpotCorrs=CorrAP;
            NucFluo(cc).Embryo(ee).APBin(aa).d_NucTraces = d_NucAP;
            NucFluo(cc).Embryo(ee).APBin(aa).d_SpotTraces = d_SpotAP;
            NucFluo(cc).Embryo(ee).APBin(aa).AllRawNucTraces = RawNucAP;  %pre-photobleach corrected nuclear traces
        end
        %clear CompiledParticles 
        %clear BurstProperties
    end
    %Combine embryo data into one big matrix for whole construct 
    for aa=1:length(APbinID)
        AllEmbNucAP=[];
        AllEmbSpotAP=[];
        AllEmbRawSpotAP=[];
        AllEmbCorrAP=[];
        Alld_NucAP =[];
        Alld_SpotAP=[];
        AllEmbRawNucAP=[];
        for ee=1:NEmbryos
            AllEmbNucAP = [AllEmbNucAP; NucFluo(cc).Embryo(ee).APBin(aa).AllNucTraces];
            AllEmbSpotAP = [AllEmbSpotAP;NucFluo(cc).Embryo(ee).APBin(aa).AllSpotTraces];
            AllEmbCorrAP = [AllEmbCorrAP;NucFluo(cc).Embryo(ee).APBin(aa).NucSpotCorrs];
            Alld_NucAP = [Alld_NucAP; NucFluo(cc).Embryo(ee).APBin(aa).d_NucTraces];
            Alld_SpotAP = [Alld_SpotAP; NucFluo(cc).Embryo(ee).APBin(aa).d_SpotTraces];
            AllEmbRawNucAP = [AllEmbRawNucAP; NucFluo(cc).Embryo(ee).APBin(aa).AllRawNucTraces];
            AllEmbRawSpotAP = [AllEmbRawSpotAP; NucFluo(cc).Embryo(ee).APBin(aa).AllRawSpotTraces];
        end
        NucFluo(cc).APBin(aa).AllNucTraces = AllEmbNucAP;
        NucFluo(cc).APBin(aa).AllSpotTraces = AllEmbSpotAP;
        NucFluo(cc).APBin(aa).AllCorrs = AllEmbCorrAP;
        NucFluo(cc).APBin(aa).Alld_NucTraces = Alld_NucAP;
        NucFluo(cc).APBin(aa).Alld_SpotTraces = Alld_SpotAP;
        NucFluo(cc).APBin(aa).AllRawNucTraces = AllEmbRawNucAP;
        NucFluo(cc).APBin(aa).AllRawSpotTraces = AllEmbRawSpotAP;
    end
end


EggLength = APbinID.*100;
FigDirect =[DropboxFolder filesep 'Figures' filesep 'KrResubmission'];
DistalColor=[1 64 172]./255;
BothColor=[52 119 71]./255;
Colorz =[BothColor; DistalColor];
colormap(Colorz)

% Save the data 
save([DropboxFolder filesep 'Constructs' filesep 'TF_MS2Corr'],'NucFluo');
%% TF vs MS2 slope all AP together, compare Distal vs SE distribution, raw TF trace 

% need to initalize the TrueCorrelation matrix so don't add a bunch of 0's
% to make same length at the Distal column since had more embryos for that
% construct 
TrueCorrelation = nan(700,2);

for cc = 1:length(ConstructList)
    AllNucVals=[];
    AllSpotVals=[];
    for aa = [16:22]%1:length(APbinID)
        AllNucVals = [AllNucVals;[NucFluo(cc).APBin(aa).AllRawNucTraces(:,1:100)]]; %limit so same # values as slopes
        AllSpotVals = [AllSpotVals;[NucFluo(cc).APBin(aa).Alld_SpotTraces]];
    end
    %Calc correlation btwn TF and MS2
    for ss = 1:size(AllSpotVals,1) 
            RealCorr = corrcoef([AllSpotVals(ss,:)],[AllNucVals(ss,:)],'Rows','pairwise');
            TrueCorrelation(ss,cc) = RealCorr(2);
            clear RealCorr 
    end
    %Graph distribution
    figure
    histogram([TrueCorrelation(:,cc)]);
    xlabel('Bcd-MS2 slope correlation');
    ylabel('# nuclei');
    title({ConstructList{cc},'avg corr no correction',num2str(nanmean([TrueCorrelation(:,cc)]))});
    MedValues(cc) = nanmedian(TrueCorrelation(:,cc));
end

%Compare distribution btwn constructs
[h,p] = kstest2([TrueCorrelation(:,1)],[TrueCorrelation(:,2)]);

% plot the 2 distributions together - pretty histograms 
fontsize=10;
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;

set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);

figure
Y.ShadowPair = [TrueCorrelation(:,1)];
Y.Distal = [TrueCorrelation(:,2)];
colormap(Colorz);
%nolines makes so doesn't plot the points corresponding to bin edges.
%samebins make bins same width.
nhist(Y,'noerror','smooth','numbers','color','colormap','nolines','samebins','binfactor',0.1);

xlabel('Bcd-MS2 slope correlation');
ylabel('number of nuclei');
title({'no correction',num2str(p)});
%legend

%Save figure
saveas(gcf, [FigDirect filesep 'DistalvSECorrDistrib_NoCorrection','.svg'],'svg');


%% Generate example Bcd - MS2 slope traces 

% need to initalize the TrueCorrelation matrix so don't add a bunch of 0's
% to make same length at the Distal column since had more embryos for that
% construct 
TrueCorrelation = nan(700,2);
ElapsedTimeUse = ([1:101] .* 0.5);
MS2Color =[192 15 0]./255;

for cc = 1:length(ConstructList)
    
    for aa = [16:22]%1:length(APbinID)
        y = datasample([1:length(NucFluo(cc).APBin(aa).AllCorrs)],3)
        for ss = 1:length(y)
            %Only plot if the correlation is >= the median 
            if NucFluo(cc).APBin(aa).AllCorrs(y(ss)) >= ((nanmedian(NucFluo(cc).APBin(aa).AllCorrs))-(nanmedian(NucFluo(cc).APBin(aa).AllCorrs)*0.5))
                figure
                yyaxis left 
                plot(ElapsedTimeUse,[NucFluo(cc).APBin(aa).AllRawNucTraces(y(ss),:)],'Color',BothColor,'LineWidth',2.5,'DisplayName','BcdGFP');
                %ylim([(nanmax([NucFluo(cc).APBin(aa).AllRawNucTraces(ss,:)])*0.5), (nanmax([NucFluo(cc).APBin(aa).AllRawNucTraces(ss,:)])*1.5)]);
                ylabel('BcdGFP fluorescence (AU)');
                
                set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
                
                yyaxis right 
                plot(ElapsedTimeUse, [NucFluo(cc).APBin(aa).AllSpotTraces(y(ss),:)],'Color',MS2Color,'LineWidth',2.5,'DisplayName','MS2');
                ylabel('MS2 fluorescence (AU)');
                set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
                xlabel('time into nc14 (min)');
                xlim([0,50]);
                %Make axis colors black 
                ax = gca;
                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = 'k';
                
                legend
                
                title([ConstructList{cc},num2str(aa),'%EL',num2str(ss),'spot']);
            end
        end
        
    end
end

%Something to run when you have a ex you want to use and want to add points
%for slope calculations 
hold on 
tt = [36 37]; %Time points to use
yyaxis right
plot(ElapsedTimeUse(tt),NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt),'o','LineWidth',1.5,'Color','k','HandleVisibility','off');
%Add vertical dashed line
plot([ElapsedTimeUse(tt(1)) ElapsedTimeUse(tt(1))],[ 0 NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1))],'LineWidth',1.5,'LineStyle','--','Color','k','HandleVisibility','off');
plot([ElapsedTimeUse(tt(2)) ElapsedTimeUse(tt(2))],[0 NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2))],'LineWidth',1.5,'LineStyle','--','Color','k','HandleVisibility','off');
% Add horizontal dashed line
yyaxis right 
plot([ElapsedTimeUse(tt(1)) ElapsedTimeUse(end)], [NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1)) NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1))],'Color','k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off');
plot([ElapsedTimeUse(tt(2)) ElapsedTimeUse(end)], [NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2)) NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2))],'Color','k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off');
%Also indicate the Bcd point 
yyaxis left 
plot(ElapsedTimeUse(tt(1)), NucFluo(cc).APBin(aa).AllRawNucTraces(ss,tt(1)),'*','LineWidth',1.5,'Color','k','HandleVisibility','off');


%% Example MS2 vs MS2 slope traces
 

% need to initalize the TrueCorrelation matrix so don't add a bunch of 0's
% to make same length at the Distal column since had more embryos for that
% construct 
TrueCorrelation = nan(700,2);
ElapsedTimeUse = ([1:101] .* 0.5);
MS2Color =[192 15 0]./255;
SlopeColor =[152 10 0]./255;

for cc = 1:length(ConstructList)
    
    for aa = [16:22]%1:length(APbinID)
        y = datasample([1:length(NucFluo(cc).APBin(aa).AllCorrs)],3)
        for ss = 1:length(y)
            %Only plot if the correlation is >= the median 
            if NucFluo(cc).APBin(aa).AllCorrs(y(ss)) >= ((nanmedian(NucFluo(cc).APBin(aa).AllCorrs))-(nanmedian(NucFluo(cc).APBin(aa).AllCorrs)*0.5))
                figure
                yyaxis left 
                plot(ElapsedTimeUse(1:100),[NucFluo(cc).APBin(aa).AllSpotTraces(y(ss),1:100)],'Color',MS2Color,'LineWidth',2.5,'DisplayName','MS2 trace');
                %ylim([(nanmax([NucFluo(cc).APBin(aa).AllRawNucTraces(ss,:)])*0.5), (nanmax([NucFluo(cc).APBin(aa).AllRawNucTraces(ss,:)])*1.5)]);
                ylabel('MS2 fluorescence (AU)');
                
                set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
                
                yyaxis right 
                plot(ElapsedTimeUse(1:100), [NucFluo(cc).APBin(aa).Alld_SpotTraces(y(ss),:)],'Color',SlopeColor,'LineWidth',2.5,'LineStyle',':','DisplayName','Slope');
                ylabel('slope of MS2 fluorescence (AU)');
                set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
                xlabel('time into nc14 (min)');
                xlim([0,50]);
                %Make axis colors black 
                ax = gca;
                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = 'k';
                
                legend
                
                title([ConstructList{cc},num2str(aa),'%EL',num2str(ss),'spot']);
            end
        end
        
    end
end

%Something to run when you have a ex you want to use and want to add points
% %for slope calculations 
hold on 
tt = [36 37]; %Time points to use
yyaxis left
plot(ElapsedTimeUse(tt),NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt),'o','LineWidth',1.5,'Color','k','HandleVisibility','off');
%Add vertical dashed line
plot([ElapsedTimeUse(tt(1)) ElapsedTimeUse(tt(1))],[ 0 NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1))],'LineWidth',1.5,'LineStyle','--','Color','k','HandleVisibility','off');
plot([ElapsedTimeUse(tt(2)) ElapsedTimeUse(tt(2))],[0 NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2))],'LineWidth',1.5,'LineStyle','--','Color','k','HandleVisibility','off');
% Add horizontal dashed line
yyaxis left 
plot([ElapsedTimeUse(1) ElapsedTimeUse(tt(1))], [NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1)) NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(1))],'Color','k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off');
plot([ElapsedTimeUse(1) ElapsedTimeUse(tt(2))], [NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2)) NucFluo(cc).APBin(aa).AllSpotTraces(ss,tt(2))],'Color','k','LineWidth',1.5,'LineStyle','--','HandleVisibility','off');
%Also indicate the Bcd point 
yyaxis right 
plot(ElapsedTimeUse(tt(1)), NucFluo(cc).APBin(aa).Alld_SpotTraces(ss,tt(1)),'*','LineWidth',1.5,'Color','k','HandleVisibility','off');
% Add a horizontal line for the slope value too 
plot([ElapsedTimeUse(tt(1)) ElapsedTimeUse(end)], [NucFluo(cc).APBin(aa).Alld_SpotTraces(ss,tt(1)) NucFluo(cc).APBin(aa).Alld_SpotTraces(ss,tt(1))],'Color','k','LineWidth',1.5,'LineStyle',':','HandleVisibility','off');
