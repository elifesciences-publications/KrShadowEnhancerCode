
%% Determine and graph total mRNA produced during nc14 by construct

% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info - respond 'y' 
ncUse=input('Want to only use nc14? y/n','s');
SlopeUse=input('Want to use Slope calculations? y/n','s');

%Organize all integrated fluorescence recordings for each construct
AvgmRNAProd=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgmRNAProdCon=[];
    firsttime=1;
    Timez=[];
    ConmRNAProdAllAP=[];
    ConmRNAProdSE=[];
    ConProdSD=[];
    EmbsArray=[];
     mRNAProdAllAP=[];
     mRNAProdErrorAllAP=[];
% Load the data for each embryo
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
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        
        %seperate out by AP bin
        for aa=1:length(APbinID)
            ProdAP=[];
            ProdErrorAP=[];
            ProductionAP=find([BurstProperties.APBin]==APbinID(aa)); %find the nuclei with bursts 
            if isempty(ProductionAP)
                ProdAP=[ProdAP; nan];
            else
            for bb=1:length(ProductionAP)
                ProdAP=[ProdAP;BurstProperties(ProductionAP(bb)).TotalmRNA];  %put all mRNA outputs at a given AP value in a column going down
                ProdErrorAP=[ProdErrorAP;BurstProperties(ProductionAP(bb)).TotalmRNAError];           
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(ProdAP)
                mRNAProdAllAP(bb,aa,ee)=ProdAP(bb);
            end
            for bb=1:length(ProdErrorAP)
                mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
            end
            mRNAProdAllAP(mRNAProdAllAP==0)=nan;  %Remove 0's recorded in BurstProperties to indicate nucleus exist without that spot of transcription active
            mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            ProductionSD(ee,aa,cc)=nanstd(mRNAProdAllAP(:,aa,ee));
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(sum(~isnan(ProdAP))); %n=# of nuclei that produce mRNA
            clear ProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        AvgProdAllAP(cc).nc14Time(ee)=ElapsedTime(end)-ElapsedTime(nc14);
        
         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).AllProds=[mRNAProdAllAP(:,:,ee)];
         Timez=[Timez,(ElapsedTime(end)-ElapsedTime(nc14))];
    end
    %Combine the data from all embryos of a construct
    for bb=1:size(mRNAProdAllAP,3)
            ConmRNAProdAllAP=[ConmRNAProdAllAP; mRNAProdAllAP(:,:,bb)];
        end
        for bb=1:size(ProductionSD,3)
            ConProdSD=[ConProdSD; ProductionSD(:,:,bb)];
        end
        for bb=1:size(ProductionSE,3)
            ConmRNAProdSE=[ConmRNAProdSE;ProductionSE(:,:,bb)];
        end

        AvgProdAllAP(cc).nc14Time=mean(Timez);
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        AvgProdAllAP(cc).AllSD=nanstd([AvgProdAllAP(cc).AllProds]);
        for aa=1:length(APbinID)
            %Perform SE calc by AP bin as there are different # of data
            %points in different AP bins
            AvgProdAllAP(cc).AllSE(aa)=(([AvgProdAllAP(cc).AllSD(aa)])./sqrt(sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))));
        end
        AvgProdAllAP(cc).All95Conf=(AvgProdAllAP(cc).AllSE).*1.95;
        AvgProdAllAP(cc).EmbryoMeans=[EmbsArray];
        clear mRNAProdAllAP mRNAProdErrorAllAP;
end

% Get rid of single values for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))==1
            AvgProdAllAP(cc).AllSD(aa)=nan;
            AvgProdAllAP(cc).AvgProd(aa)=nan;
            AvgProdAllAP(cc).AllSE(aa)=nan;
            AvgProdAllAP(cc).All95Conf(aa)=nan;
        end
    end
end


    

%% Plotting 
EggLength=APbinID.*100;
FractInfo=[DropboxFolder filesep 'Constructs' filesep 'FractON.mat'];

% 3/3/2019 updated/corrected Frnap value to be 379 instead of 377
% MS2Conversion =F1  which for Wunderlich Nikon scope is = 1338

MS2Conversion=(379*((1.275+4.021)/1.5)); %Frnap * Telongation ((length of MS2 +rest of transcript) / elongation rate)
HbMS2Conversion=(379*((1.275+4.021)/1.5));
for cc=1:length(ConstructList)
    AvgProdAllAP(cc).AvgmRNACounts=[[AvgProdAllAP(cc).AvgProd]./MS2Conversion];
end

% Set construct colors
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
Colors(16).Color='k';
Colors(17).Color=DoubDistColor;
Colors(18).Color=DistalColor;
Colors(19).Color=DoubDistColor;
Colors(20).Color=BothColor;
Colors(21).Color=DoubDistColor;
Colors(22).Color=ProxColor;
Colors(23).Color=DoubDistColor;
Colors(24).Color=DoubProxColor;
Colors(25).Color=DoubProxColor;
Colors(26).Color=DoubDistColor;
Colors(27).Color=DoubDistColor;
Colors(28).Color=BothColor;
Colors(29).Color=BothColor;
Colors(30).Color=DoubDistColor;
Colors(31).Color=BothColor;
Colors(32).Color=DistalColor;
Colors(33).Color=DistalColor;
Colors(34).Color=DistalColor;
% Set font/display parameters
FontUsed=input('Want larger font?','s');
if FontUsed=='y'
    fontsize=15;
else
fontsize=10;
end
fontname='Helvetica';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 6; ySize = 4.5; xLeft = 0.1; yTop = 0.1;

FigDirect=[DropboxFolder filesep 'Figures'];
%% AP bins of max expression each construct
% Record the AP bin of max expression for each construct to have for other
% calculations (ie noise)
ExpressionBins=cell(2,(length(ConstructList)));
for cc=1:length(ConstructList)
    AvgProdAP=[nanmean(AvgProdAllAP(cc).AllProds)];
    APBintoUse=find((nanmean(AvgProdAllAP(cc).AllProds)==max(nanmean(AvgProdAllAP(cc).AllProds))));
    ExpressionBins{1,cc}=APBintoUse;
    ExpressionBins{2,cc}=ConstructList{cc};
end
save([DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins'],'ExpressionBins');

%% Create figures for core constructs
for cc=[32,2,6,7,9];
    %figure
 errorbar(EggLength,AvgProdAllAP(cc).AvgProd,AvgProdAllAP(cc).All95Conf,'Color',Colors(cc).Color,'LineWidth',2.5);
 hold on 
 yyaxis right
 plot(EggLength, (AvgProdAllAP(cc).AvgProd./MS2Conversion), 'Color', Colors(cc).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 950000]);
 print( [FigDirect filesep ConstructList{cc} 'TotalmRNA'],'-dpdf');
end

%% Singles/duplicated vs shadow enhancer pair
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd),(AvgProdAllAP(7).All95Conf),'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 print( [FigDirect filesep 'DuplicatedvsBothTotalmRNA'],'-dpdf');

 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd),(AvgProdAllAP(2).All95Conf),'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
 print( [FigDirect filesep 'SinglesvsBothTotalmRNA_ED'],'-dpdf');
 
 %Just singles
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd),(AvgProdAllAP(2).All95Conf),'Color',Colors(2).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
 print( [FigDirect filesep 'SinglesTotalmRNA_ED'],'-dpdf');
 %% Compare Endogenous Distal spacing
figure
  %yyaxis left
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineStyle',':','LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color,'LineStyle',':','LineWidth',2.5);
 plot(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion), 'Color', Colors(32).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 900000]);
print( [FigDirect filesep 'EndogvsOrigDistTotalmRNA'],'-dpdf');