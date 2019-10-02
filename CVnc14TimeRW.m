%% calculate coefficient of variance across time of nc14 for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDistEmpty';'KrProxEmpty';'KrBothEmpty';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C';'Kr2xDistdBcd';'Kr2xProxdHb';'KrProx17C';'Kr2xProx17C';'Kr2xDist32C';'Kr2xDist17C';'Kr2xDistEmpty';'Kr2xDistEmpty32C';'KrInvSE';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

NumberNuclei=nan(20,41,[length(ConstructList)]);
% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    APTable=[];
    ConTimeAvg=[];
    AllNucTimeTable=[];
    AllTotmRNA=[];
    AllNucStrengthTable=[];
    ConFluoMeanz=[];
    AllVariance=[];
    AllSD=[];
    AllMedStrength=[];
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
     
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(Filename);

        APstuff=[SpotDiff.APBin];

        for aa=1:length(APbinID)
            APsubset=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            %Spotstuff=[APsubset.SpotOne];
            % Use smoothed data 20190123 RW
            Spotstuff=[APsubset.SmoothSpotOne];
            
            APSpotsubset=APsubset(~isempty(Spotstuff));
            APNucs(aa)=length(APsubset);
        end
        MostAPNucs=max(APNucs);
        MostAPNucs=2*MostAPNucs;
        TimeTable=nan(90,41); %changed to make all consistent size across embryos 2/6/19
        TotmRNATable=nan(90,41);
        TotmRNATable=nan(90,41);
        StrengthTable=nan(90, 41);
        MeanzTable=nan(90,41);
        VarianceTable=nan(90,41);
        SDTable=nan(90,41);
        MedStrengthTable=nan(90,41);
        for aa=1:length(APbinID)
            APsubset=[];
%           
            APsubset=SpotDiff(APstuff==APbinID(aa));
            if ~isempty(APsubset)
               
            for bb=1:length(APsubset)
                NumberNuclei(ee,aa,cc)=length(APsubset);
                TimeLength=[];
                TimeLength=sum(~isnan(APsubset(bb).SmoothSpotOne));   %Only want to look at frames where nucleus exisits (i.e 0 or number values)
                if TimeLength==1 %Ignore spots with only 1 frame of existance since we can't calculate the noise across time of that spot
                    AvgFluo=nan;
                    MedianFluo=nan;
                else
                AvgFluo=nanmean([APsubset(bb).SmoothSpotOne]);
                MedianFluo=nanmedian([APsubset(bb).SmoothSpotOne]);
                
                end
                VarFluo=nanstd([APsubset(bb).SmoothSpotOne]);
                VarianceFluo=nanvar([APsubset(bb).SmoothSpotOne]);
                
            TimeTable(bb,aa)=(VarFluo/AvgFluo); %Table of CV values for this embryo
            if isempty(APsubset(bb).TotalmRNAOne) %Record total mRNA produced by each spot in same format so can compare CV to expression lvl
                TotmRNATable(bb,aa)=nan;
            else
            TotmRNATable(bb,aa)=APsubset(bb).TotalmRNAOne;
            end
            MeanzTable(bb,aa)=AvgFluo;
            VarianceTable(bb,aa)=VarianceFluo;
            StrengthTable(bb,aa)=(VarianceFluo/AvgFluo);
            SDTable(bb,aa)=VarFluo;
            MedStrengthTable(bb,aa)=(VarianceFluo/MedianFluo);
            end
            % Do the same calculations for the 2nd alleles 
            for bb=1:length(APsubset)
                if isfield(APsubset,'SpotTwo') & (~isempty(APsubset(bb).SpotTwo))
                    TimeLength2=[];
                TimeLength2=(sum(~isnan(APsubset(bb).SmoothSpotTwo)));
                if TimeLength2==1
                    AvgFluo2=nan;
                    MedianFluo2=nan;
                else
                AvgFluo2=nanmean([APsubset(bb).SmoothSpotTwo]);
                MedianFluo2=nanmedian([APsubset(bb).SmoothSpotTwo]);
                end
                VarFluo2=nanstd([APsubset(bb).SmoothSpotTwo]);
                VarianceFluo2=nanvar([APsubset(bb).SmoothSpotTwo]);
                TimeTable(bb+(length(APsubset)),aa)=(VarFluo2/AvgFluo2);
                if isempty(APsubset(bb).TotalmRNATwo)
                    TotmRNATable(bb+(length(APsubset)),aa)=nan;
                else
                TotmRNATable(bb+(length(APsubset)),aa)=APsubset(bb).TotalmRNATwo;
                end
                MeanzTable(bb+(length(APsubset)),aa)=AvgFluo2;
                VarianceTable(bb+(length(APsubset)),aa)=VarianceFluo2;
                StrengthTable(bb+(length(APsubset)),aa)=(VarianceFluo2/AvgFluo2);
                SDTable(bb+(length(APsubset)),aa)=VarFluo2;
                MedStrengthTable(bb+(length(APsubset)),aa)=(VarianceFluo2/MedianFluo2);
                end
            end
            
            end
        end
        %Organize into a unified structure, including breaking down by
        %embryo
        WholeTimeTable(cc).Embryo(ee).TimeTable=TimeTable;
        WholeTimeTable(cc).Embryo(ee).TotalmRNATable=TotmRNATable;
        WholeTimeTable(cc).Embryo(ee).MeanFluo=MeanzTable;
        WholeTimeTable(cc).Embryo(ee).Variance=VarianceTable;
        WholeTimeTable(cc).Embryo(ee).FanoFactors=StrengthTable;
        WholeTimeTable(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        WholeTimeTable(cc).Embryo(ee).SD=SDTable;
        WholeTimeTable(cc).Embryo(ee).MedFano=MedStrengthTable;
        %Create arrays to hold all of the data from all embryos of a given
        %construct
        AllNucTimeTable=[AllNucTimeTable;TimeTable];
        AllTotmRNA=[AllTotmRNA;TotmRNATable];
        AllNucStrengthTable=[AllNucStrengthTable; StrengthTable];
        ConTimeAvg=[ConTimeAvg;(nanmean(TimeTable))];
        ConFluoMeanz=[ConFluoMeanz; MeanzTable];
        AllVariance=[AllVariance;VarianceTable];
        AllSD=[AllSD; SDTable];
        AllMedStrength=[AllMedStrength; MedStrengthTable];
    end
    %Organize the construct-wide data
    WholeTimeTable(cc).ConstructAvg=nanmean(ConTimeAvg);
    WholeTimeTable(cc).AllAvgFluo=ConFluoMeanz;
    WholeTimeTable(cc).AllVariance=AllVariance;
    WholeTimeTable(cc).AllNucs=AllNucTimeTable;
    WholeTimeTable(cc).TotalmRNA=AllTotmRNA;
    WholeTimeTable(cc).AllFanoFactor=AllNucStrengthTable;
    WholeTimeTable(cc).AvgAllNucs=nanmean([WholeTimeTable(cc).AllNucs]);
    WholeTimeTable(cc).MedAllNucs=nanmedian([WholeTimeTable(cc).AllNucs]);
    WholeTimeTable(cc).SDAllNucs=nanstd([WholeTimeTable(cc).AllNucs]); %SD of all of the spots of a construct 
    WholeTimeTable(cc).AllSD=AllSD; %SD of each spot across nc14 
    WholeTimeTable(cc).AllMedStrength=AllMedStrength;
    
     for aa=1:length(APbinID)
         WholeTimeTable(cc).SEAllNucs(aa)=(WholeTimeTable(cc).SDAllNucs(aa))/(sqrt(sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa)))));
     end
     WholeTimeTable(cc).Conf95All=(WholeTimeTable(cc).SEAllNucs).*1.95; 
     WholeTimeTable(cc).ConstructName=ConstructList{cc};
end

% Get rid of places where only have one data point for a whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa))) ==1
            WholeTimeTable(cc).AllNucs([find(~isnan(WholeTimeTable(cc).AllNucs(:,aa)))],aa)=nan;
            WholeTimeTable(cc).AvgAllNucs(aa)=nan;
            WholeTimeTable(cc).MedAllNucs(aa)=nan;
        end
        NumberNuclei(cc,aa)=sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa)));
    end
end
%% Plotting 
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
DoubProxColor=[215 183 58] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;
Colors(7).Color=DistalColor;
Colors(8).Color=ProxColor;
Colors(9).Color=BothSepColor;
Colors(10).Color=BothColor;
Colors(11).Color=DoubProxColor;
Colors(12).Color=DistalColor;
Colors(13).Color=ProxColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;
Colors(18).Color=DoubDistColor;
Colors(19).Color=DoubProxColor;
Colors(20).Color=ProxColor;
Colors(21).Color=DoubProxColor;
Colors(22).Color=DoubDistColor;
Colors(23).Color=DoubDistColor;
Colors(24).Color=DoubDistColor;
Colors(25).Color=DoubDistColor;
Colors(26).Color=BothColor;
Colors(27).Color=DistalColor;
Colors(28).Color=DistalColor;
Colors(29).Color=DistalColor;

fontsize=15;
fontname='Helvetica';
x_width=3; y_width=2.25;
xSize = 6; ySize = 4.5; xLeft = 0.1; yTop = 0.1;

EggLength=APbinID .* 100;
% Set folder to save figures
FigDirect=[DropboxFolder filesep 'Figures'];
% Save structure containing data to create figures
save([DropboxFolder filesep 'Constructs' filesep 'TemporalCVData'],'WholeTimeTable');

%% Compare temporal CV as f(x) of egg length duplicated vs shadow enhancer pair
figure 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs, WholeTimeTable(4).Conf95All, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).Conf95All, 'Color',Colors(5).Color,'LineWidth',2.5);
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
xlabel('% egg length');
xlim([0 100]);
ylim([0 4]);
ylabel('coefficient of variation');
%legend('2x Proximal', 'Both');
print('-painters',[FigDirect filesep '2EnhancervBothTempCVSmoothed'],'-dsvg');
%% Compare Distal @ endogenous spacing
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle',':');
hold on
errorbar(EggLength, WholeTimeTable(27).AvgAllNucs, WholeTimeTable(27).Conf95All, 'Color',Colors(27).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
xlabel('% egg length');
xlim([0 100]);
ylim([0 4]);
ylabel('coefficient of variation');
print('-painters',[FigDirect filesep 'EndogvsOrigDistTempCV'],'-dsvg');
%Zoom to center 20%
ZoomXSize=3; ZoomYSize=2.25;
set(gcf, 'PaperUnits','inches','PaperPosition',[xLeft yTop ZoomXSize ZoomYSize]);
xlabel('% egg length');
xlim([40 60]);
ylim([0 3]);
ylabel('coefficient of variation');
%legend('2x Proximal', 'Both');
print('-painters',[FigDirect filesep 'EndogDistvsOrigTempCVZoom'],'-dsvg');

%% Temporal CV as a function of mean expression (fluorescence) 
for cc=[1:6,27]%length(ConstructList)
    figure
    h=scatter(WholeTimeTable(cc).AllAvgFluo(:),WholeTimeTable(cc).AllNucs(:),[],Colors(cc).Color);
    title(ConstructList{cc});
    ylabel('temporal CV');
    xlabel('mean fluorescence');
    ylim([0 5]);
    xlim([0 45000]);
    set(gca, 'FontSize', fontsize, 'FontName', fontname);
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    print('-painters',[FigDirect filesep ConstructList{cc} 'CVvsMeanFluo'],'-dsvg');
end