% Load chosen enhancer constructs to compare 
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'Kr2xDist32C';'Kr2xProx32C';'KrBoth32C'; 'KrDist17C';'KrBoth17C';'Kr2xDist17C';'KrProx17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

ConCorrMean=[];
MeanAPCorr=[];
MeanAllCorr=[];
APCorr=[];

for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgCorrCon=[];
    firsttime=1;
    ConCorrAllAP=[];
    ConCorrSE=[];
    ConCorrSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        
        filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(filename);
        
        %seperate out by AP bin
        CorrAllAP=[];
        for aa=2:length(APbinID) %skip APbin 0
            CorrAP=[];
            CorrelationAP=find([SpotDiff.APBin]==APbinID(aa));
            if isempty(CorrelationAP)
                CorrAP=[CorrAP; nan];
            else
            for bb=1:length(CorrelationAP)  
                if ~isfield(SpotDiff, 'SpotCorr')
                    
                    CorrAP=[CorrAP; nan];
                elseif isempty(SpotDiff(CorrelationAP(bb)).SpotCorr)
                    CorrAP=[CorrAP;nan];
                elseif length([SpotDiff(CorrelationAP(bb)).SpotCorr])==1
                    CorrAP=[CorrAP;SpotDiff(CorrelationAP(bb)).SpotCorr(1)];
                else
                CorrAP=[CorrAP;SpotDiff(CorrelationAP(bb)).SpotCorr(1,2)];  
                end
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(CorrAP)
                CorrAllAP(bb,aa,ee)=CorrAP(bb);
            end
            CorrAllAP(CorrAllAP==0)=nan;
            
            CorrelationSD(ee,aa,cc)=nanstd(CorrAllAP(:,aa,ee));
        
        CorrelationSE(ee,aa,cc)=CorrelationSD(ee,aa,cc)/sqrt(sum(~isnan(CorrAP)));
            clear CorrelationAP
        end
        %Compile all correlation data for a construct in one long column
        %to put in structure
        for bb=1:size(CorrAllAP,3)
            ConCorrAllAP=[ConCorrAllAP; CorrAllAP(:,:,bb)];
        end
        for bb=1:size(CorrelationSD,3)
            ConCorrSD=[ConCorrSD; CorrelationSD(:,:,bb)];
        end
        for bb=1:size(CorrelationSE,3)
            ConCorrSE=[ConCorrSE;CorrelationSE(:,:,bb)];
        end
        AvgCorrAllAP(cc).EmbryoCorr(ee).AvgCorr=nanmean(CorrAllAP(:,:,ee));
        AvgCorrAllAP(cc).EmbryoCorr(ee).SE=CorrelationSE(ee,:,cc);
        AvgCorrAllAP(cc).EmbryoCorr(ee).SD=CorrelationSD(ee,:,cc);
    end
        AvgCorrAllAP(cc).AvgCorr=nanmean(ConCorrAllAP);  
        AvgCorrAllAP(cc).CorrSD=nanmean(ConCorrSD);
        AvgCorrAllAP(cc).CorrSE=nanmean(ConCorrSE);
        AvgCorrAllAP(cc).AllCorrs=[ConCorrAllAP];
        AvgCorrAllAP(cc).AllSD=nanstd(AvgCorrAllAP(cc).AllCorrs);
        %calculate SE by AP bin since different number of tracked spots in
        %different AP bins 
        for aa=1:length(APbinID)
        AvgCorrAllAP(cc).AllSE(aa)=(AvgCorrAllAP(cc).AllSD(aa))/sqrt(sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa))));
        AvgCorrAllAP(cc).NumNucUsed(aa)=sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa)));
        end
        AvgCorrAllAP(cc).All95Conf=(AvgCorrAllAP(cc).AllSE) .* 1.95;

end

%Remove data relying on only one data point for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
    if sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa)))==1
        AvgCorrAllAP(cc).AllSD(aa)=nan;
        AvgCorrAllAP(cc).AllSE(aa)=nan;
        AvgCorrAllAP(cc).All95Conf(aa)=nan;
        AvgCorrAllAP(cc).AvgCorr(aa)=nan;
    end
    end
end

%% Plotting 
%Plot mean of each embryo for each construct
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[149 188 114] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[149 188 114] ./ 255;
BothEmptyColor=[12 250 100] ./ 255;
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
Colors(10).Color=DoubDistColor;
Colors(11).Color=DoubProxColor;
Colors(12).Color=BothColor;
Colors(13).Color=DistalColor;
Colors(14).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures'];

EggLength=APbinID .* 100;
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;


%% Graph allele correlation as fx of egg length for homozygous and heterozygous embryos

figure 
errorbar(EggLength, AvgCorrAllAP(1).AvgCorr, AvgCorrAllAP(1).All95Conf, 'Color',Colors(1).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(2).AvgCorr, AvgCorrAllAP(2).All95Conf, 'Color',Colors(2).Color, 'LineWidth',2.5); 
errorbar(EggLength, AvgCorrAllAP(3).AvgCorr, AvgCorrAllAP(3).All95Conf, 'Color',Colors(3).Color, 'LineWidth',2.5); 
%legend('Distal', 'Proximal', 'Separated');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('allele correlation');
xlim([0 100]);
%print( [FigDirect filesep 'SinglesvSepCorrelation'],'-dsvg');

%interpolate missing points in proximal (bc only had 1 value point) 
F=fillmissing(AvgCorrAllAP(2).AvgCorr,'linear');
PointsKept=[~isnan(AvgCorrAllAP(2).AvgCorr)];
F(1:(find(PointsKept,1,'first')-1))=nan;
F((find(PointsKept,1,'last')+1):end)=nan;
figure 
errorbar(EggLength, AvgCorrAllAP(1).AvgCorr, AvgCorrAllAP(1).All95Conf, 'Color',Colors(1).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, F, AvgCorrAllAP(2).All95Conf, 'Color',Colors(2).Color, 'LineWidth',2.5); 
errorbar(EggLength, AvgCorrAllAP(3).AvgCorr, AvgCorrAllAP(3).All95Conf, 'Color',Colors(3).Color, 'LineWidth',2.5); 
%legend('Distal', 'Proximal', 'Separated');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('allele correlation');
xlim([0 100]);
print( [FigDirect filesep 'SinglesvSepCorrelation'],'-dsvg');


%% Save avg correlation vectors
DistalCorr=[AvgCorrAllAP(1).AvgCorr];
ProximalCorr=[AvgCorrAllAP(2).AvgCorr];
SepCorr=[AvgCorrAllAP(3).AvgCorr];
save([DropboxFolder filesep 'Constructs' filesep 'DistalAlleleCorr'],'DistalCorr');
save([DropboxFolder filesep 'Constructs' filesep 'ProximalAlleleCorr'],'ProximalCorr');
save([DropboxFolder filesep 'Constructs' filesep 'SepAlleleCorr'],'SepCorr');

