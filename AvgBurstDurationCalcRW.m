%% Calculate average duration of transcriptional bursts for each construct by AP position 

%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'KrBoth17C';'Kr2xDist32C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'Kr2xDistEmpty';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?, y/n','s');
SlopeUse=input('Want to use slope calculations?, y/n','s');

%Count for each construct
AvgDuration=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgDurationCon=[];
    firsttime=1;
    ConDurAllAP=[];
    ConDurSE=[];
    ConDurSD=[];
    for ee=1:NEmbryos 
        PrefixName=Data(ee).Prefix;
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        %seperate out by AP bin
        DurAllAP=[];
        for aa=1:length(APbinID)
            DurAP=[];
            DurationAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(DurationAP)
                DurAP=[DurAP; nan];
            else
            for bb=1:length(DurationAP)
                if ~isempty(BurstProperties(DurationAP(bb)).Duration)
                DurAP=[DurAP;[BurstProperties(DurationAP(bb)).Duration]'];  %put all durations at a given AP value in a column going down
                else
                    DurAP=[DurAP; nan];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(DurAP)
                DurAllAP(bb,aa,ee)=DurAP(bb);
            end
            DurAllAP(DurAllAP==0)=nan;
            
            DurationSD(ee,aa,cc)=nanstd(DurAllAP(:,aa,ee));
        %7.16.18
        DurationSE(ee,aa,cc)=DurationSD(ee,aa,cc)/sqrt(sum(~isnan(DurAP))); % n= # spots with bursts
            clear DurationAP
        end
        AvgDurAllAP(cc).EmbryosDur(ee).MeanDuration=nanmean(DurAllAP(:,:,ee));
        AvgDurAllAP(cc).EmbryosDur(ee).Durs=(DurAllAP(:,:,ee));
        AvgDurAllAP(cc).EmbryosDur(ee).SD=DurationSD(ee,:,cc);
        AvgDurAllAP(cc).EmbryosDur(ee).SE=DurationSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(DurAllAP,3)
            ConDurAllAP=[ConDurAllAP; DurAllAP(:,:,bb)];
        end
        for bb=1:size(DurationSD,3)
            ConDurSD=[ConDurSD; DurationSD(:,:,bb)];
        end
        for bb=1:size(DurationSE,3)
            ConDurSE=[ConDurSE;DurationSE(:,:,bb)];
        end
        
    end
        AvgDurAllAP(cc).AvgDur=nanmean(ConDurAllAP);  %Avg duration of all embryos of a construct
        AvgDurAllAP(cc).ConSD=nanmean(ConDurSD);
        AvgDurAllAP(cc).ConSE=nanmean(ConDurSE);
        AvgDurAllAP(cc).AllDurs=[ConDurAllAP];
        AvgDurAllAP(cc).AllSD=nanstd(ConDurAllAP);
        for aa=1:length(APbinID)
            %Calc SE per AP bin since there are different number of spots
            %in different AP bins 
            AvgDurAllAP(cc).AllSE(aa)=([AvgDurAllAP(cc).AllSD(aa)])/(sqrt(sum(~isnan(AvgDurAllAP(cc).AllDurs(:,aa)))));
        end
        AvgDurAllAP(cc).All95Conf=([AvgDurAllAP(cc).AllSE]).*1.95;

end
%Get rid of single data points of a whole AP bin
for cc=1:length(ConstructList)
    AvgDurAllAP(cc).Construct=ConstructList{cc}
    for aa=1:length(APbinID)
        if sum(~isnan(AvgDurAllAP(cc).AllDurs(:,aa)))==1
            AvgDurAllAP(cc).AllSD(aa)=nan;
            AvgDurAllAP(cc).AvgDur(aa)=nan;
            AvgDurAllAP(cc).AllSE(aa)=nan;
            AvgDurAllAP(cc).All95Conf(aa)=nan;
        end
    end
end

%% Plot mean of each embryo at specific AP position for each construct

Egglength=APbinID.*100; %Nicer to look at % egg length vs fractions

%Define the colors of your constructs
DistalColor=[1 64 172]./255;
DistalEmptyColor=[8 210 238] ./ 255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
BothEmptyColor=[12 250 100] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;


Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=DoubProxColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor;
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;
Colors(18).Color=DoubDistColor;
Colors(19).Color=ProxColor;
Colors(20).Color=DoubDistColor;
Colors(21).Color=DoubProxColor;
Colors(22).Color=DoubProxColor;
Colors(23).Color=DoubDistColor;
Colors(24).Color=BothColor;
Colors(25).Color=BothColor;
Colors(26).Color=BothColor;
Colors(27).Color=DistalColor;
Colors(28).Color=DistalColor;
Colors(29).Color=DistalColor;
Colors(30).Color=DistalColor;

%Display/font settings
fontsize=12;
fontname='Helvetica';
x_width=3; y_width=2.25;
xSize = 6; ySize = 4.5; xLeft = 0.1; yTop = 0.1;
%Folder to save to 
FigDirect=[DropboxFolder filesep 'Figures'];
% Save structure needed to create all graphs
save([DropboxFolder filesep 'Constructs' filesep 'BurstDurationSlope'],'AvgDurAllAP');

%% singles/duplicated enhancers vs SE pair
%Either single vs SE pair
figure
errorbar(Egglength,AvgDurAllAP(27).AvgDur,AvgDurAllAP(27).All95Conf,'Color',Colors(27).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'SinglesvBothDuration_ED'],'-dsvg');

%doubles vs SE
figure
errorbar(Egglength,AvgDurAllAP(6).AvgDur,AvgDurAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(7).AvgDur,AvgDurAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min');
xlim([0 100]);

%Singles only
figure
errorbar(Egglength,AvgDurAllAP(27).AvgDur,AvgDurAllAP(27).All95Conf,'Color',Colors(27).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'SinglesDuration_ED'],'-dsvg');
%% Distal @ Endogenous 
%Distal at endogneous spacing vs proximal to promoter 
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle',':');
hold on 
errorbar(Egglength,AvgDurAllAP(27).AvgDur,AvgDurAllAP(27).All95Conf,'Color',Colors(27).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
ylim([0 5]);
print( [FigDirect filesep 'EndogvsOrigDistDuration'],'-dpdf');
