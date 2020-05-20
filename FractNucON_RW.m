%% Fraction of nuclei active as f(x) of egg length, time 
%load constructs
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth'}%,'KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','KrSEChrom3','Kr2xDistChrom3','Kr2xProxChrom3','Kr2xDistChrom3_Empty','Kr4_1xBcd','KrSEChr3_Empty','Kr3_1xBcd','Kr4_1xBcd'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};

[SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2, OutputFolder...
 ]= readMovieDatabase('2017-08-03-mKr1_E1')  %just any random dataset to give us the dropbox folder location

% Go thru each construct
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    FluoThreshold = 379; %Fluorescence value of 1 mRNA added 12/2/19
    
    % Go through each embryo
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        nc14 = Data(ee).nc14;
        % Load SpotDiff data and originial schnitzcells 
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat']
        FileName2=[DropboxFolder filesep PrefixName filesep 'APDetection.mat']
        FileName3=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat']
        load(Filename);
        load(CompPars);
        load(FileName2);
        load(FileName3);
        APstuff=[SpotDiff.APBin];
        
        % Put nucleus field in Schnitz 
        for ss=1:length(schnitzcells)
            schnitzcells(ss).Nucleus=ss;
        end
        %Limit to schnitz to nuclei that exist in nc14 
        schnitzcells_n14=[schnitzcells];
        for ss=1:length(schnitzcells_n14)
            if (max(schnitzcells_n14(ss).frames)) < nc14
                schnitzcells_n14(ss).Nucleus=nan;
            end
        end
        %Need to add AP bin info to schnitzcells
        
        %first translate x,y coordinates to AP position 
        %Angle between the x-axis and the AP-axis
        if exist('coordPZoom', 'var')
            APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
        else
            error('coordPZoom not defined. Was AddParticlePosition.m run?')
        end
        APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
        NucPosition=[];


        for i=1:length(schnitzcells_n14)     %find the AP position of each nucleus across time 
            for j=1:length(schnitzcells_n14(i).frames)
            
                %Angle between the x-axis and the particle using the A position as a
                %zero
                Angles=atan((schnitzcells_n14(i).ceny(j)-coordAZoom(2))./(schnitzcells_n14(i).cenx(j)-coordAZoom(1)));
            
                %Distance between the points and the A point
                Distances=sqrt((coordAZoom(2)-schnitzcells_n14(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells_n14(i).cenx(j)).^2);
                APPositions=Distances.*cos(Angles-APAngle);
                NucTimeFrame=schnitzcells_n14(i).frames(j);  %making the columns the time frame, each row is a nucleus 
                NucPosition(i,NucTimeFrame)=APPositions/APLength;
            
            end
        end
        
        for ii=1:size(NucPosition,1)
            schnitzcells_n14(ii).APPos=[NucPosition(ii,:)];
            schnitzcells_n14(ii).MeanAP=nanmean(schnitzcells_n14(ii).APPos);
        end
        
        %estimate AP bin
        for j = 1:length(schnitzcells_n14)
            APEstm(j) = round(schnitzcells_n14(j).MeanAP,2);
            for jj = 1:length(APbinID)
                if APEstm(j) < APbinID(jj) 
                    schnitzcells_n14(j).APbin = APbinID(jj)
                    break;
                end
            end
        end 
        NucAPStuff = [schnitzcells_n14.APbin];
        
        for aa=1:length(APbinID)
            %Find Nuclei in each AP bin 
            APsubset = SpotDiff(APstuff == APbinID(aa));
            NucSubset = schnitzcells_n14(NucAPStuff == APbinID(aa));
            APmRNA=[];
            APFrames=[];
            if isempty(NucSubset)
                NucFrames = nan(1,100);
                SpotmRNA = nan(1,100);
                APmRNA = [APmRNA; SpotmRNA];
                APFrames = [APFrames; NucFrames];
                %continue
            elseif ~isempty(APsubset) 
                for ss=1:length(APsubset)
                    SpotmRNA =[APsubset(ss).SmoothSpotOne]';
                    %find the frames w fluorescence
                    %ONFrames = [find(APsubset(ss).SmoothSpotOne >1)];
                    % set minimal threshold as FRNAP 12/2/19
                    ONFrames = [find(APsubset(ss).SmoothSpotOne > FluoThreshold)];
                    OFFFrames = [find(APsubset(ss).SmoothSpotOne <= FluoThreshold)]; %find frames where no fluo
                    %NucFrames = zeros(1,(length(ElapsedTime)-nc14));
                    %Change to make length of elapsedtime bc that's what
                    %the spotdiff field of smoothtrace is 
                    NucFrames = nan(1,(length(ElapsedTime))); %changed to nan 12/2/19
                    %Mark which frames had fluorescence
                    NucFrames(ONFrames) = 1;
                    %Mark which frames existed but OFF
                    NucFrames(OFFFrames) = 0;
                    %Restrict to nc14
                    NucFrames = NucFrames(nc14:end);
                    SpotmRNA = SpotmRNA(nc14:end);
                    % Make all of them 100 points long (50 minutes) 
                    if (length(NucFrames) < 100) | (length(SpotmRNA) <100)
                        NucFrames= [NucFrames, nan(1,(100 -(length(NucFrames))))];
                        SpotmRNA = [SpotmRNA, nan(1,(100 -(length(SpotmRNA))))];
                    end
                    if (length(NucFrames) > 100) | (length(SpotmRNA) >100)
                        NucFrames = NucFrames(1:100);
                        SpotmRNA = SpotmRNA(1:100);
                    end
                     
                    APmRNA = [APmRNA; SpotmRNA];
                    APFrames = [APFrames; NucFrames];
                    
                end
            elseif (~isempty(NucSubset) & isempty(APsubset))
                %Fill with zeros if no active nuclei
                for nn = 1:length(NucSubset) 
                    NucFrames = zeros(1,100);
                    SpotmRNA = zeros(1,100);
                    APmRNA = [APmRNA; SpotmRNA];
                    APFrames = [APFrames; NucFrames];
                    
                end
            end
            
            FractAct(cc).Embryo(ee).APbin(aa).TimeFrames = APFrames;
            FractAct(cc).Embryo(ee).APbin(aa).TimeFracts = (nansum(FractAct(cc).Embryo(ee).APbin(aa).TimeFrames)./(sum(~isnan(FractAct(cc).Embryo(ee).APbin(aa).TimeFrames))));
            if length(FractAct(cc).Embryo(ee).APbin(aa).TimeFracts) == 1
                FractAct(cc).Embryo(ee).APbin(aa).TimeFracts = [nan(1,100)];
            end
            FractAct(cc).Embryo(ee).APbin(aa).TimeProd = APmRNA;
            
        end
        
    end
    for aa=1:length(APbinID)
        AllFract =[];
        AllProd=[]
        EmbFract =[];
        for ee=1:NEmbryos
            
            AllFract = [AllFract; [FractAct(cc).Embryo(ee).APbin(aa).TimeFrames]];
            AllProd = [AllProd; [FractAct(cc).Embryo(ee).APbin(aa).TimeProd]];
            EmbFract =[EmbFract; FractAct(cc).Embryo(ee).APbin(aa).TimeFracts];
        end
        FractAct(cc).APbin(aa).EmbFracts =EmbFract;
        %Calculate things based on each embryo being a replicate
        FractAct(cc).APbin(aa).AvgEmbFracts = nanmean(FractAct(cc).APbin(aa).EmbFracts);
        FractAct(cc).APbin(aa).SDEmbFracts = nanstd(FractAct(cc).APbin(aa).EmbFracts);
        FractAct(cc).APbin(aa).SEEmbFracts = ((FractAct(cc).APbin(aa).SDEmbFracts)./NEmbryos);
        FractAct(cc).APbin(aa).EmbFracts95CI = FractAct(cc).APbin(aa).SEEmbFracts .* 1.95;
        FractAct(cc).APbin(aa).AllTimeFrames = AllFract;
        % Actually divide number nuclei on by total number of nuclei
        FractAct(cc).APbin(aa).AllFracts = nansum(FractAct(cc).APbin(aa).AllTimeFrames)./(sum(~isnan(FractAct(cc).APbin(aa).AllTimeFrames)))
        FractAct(cc).APbin(aa).AllTimeProd = AllProd;
        FractAct(cc).APbin(aa).AvgFractON = nanmean(FractAct(cc).APbin(aa).AllFracts);
        FractAct(cc).APbin(aa).AvgProd = nanmean(FractAct(cc).APbin(aa).AllTimeProd);
        FractAct(cc).APbin(aa).SDFractON = nanstd(FractAct(cc).APbin(aa).AllFracts);
        FractAct(cc).APbin(aa).SDProd = nanstd(FractAct(cc).APbin(aa).AllTimeProd);
        %Number of nuclei per bin per time 
        NumbNuc = sum(~isnan(FractAct(cc).APbin(aa).AllTimeFrames));
        FractAct(cc).APbin(aa).Fract95CI = ((FractAct(cc).APbin(aa).SDFractON)./sqrt(NumbNuc));
        FractAct(cc).APbin(aa).Prod95CI = ((FractAct(cc).APbin(aa).SDProd)./sqrt(NumbNuc));
        
    end
    for tt = 1:100
        TimeCenterFrame = [];
        TimeCI =[];
        for aa =1:length(APbinID)
            TimeCenterFrame = [TimeCenterFrame,[FractAct(cc).APbin(aa).AvgEmbFracts(tt)]];
            TimeCI = [TimeCI, FractAct(cc).APbin(aa).EmbFracts95CI(tt)]; %FIX
            
        end
        FractAct(cc).TimeFract(tt).TimeAPFract = [TimeCenterFrame];
        FractAct(cc).TimeFract(tt).TimeAPCI = [TimeCI];
    end
            
            
end
%% Plotting
EggLength = APbinID.*100;
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
PDistalColor = [129 161 214]./255;


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
Colors(35).Color=BothColor;
Colors(36).Color = DoubDistColor;
Colors(37).Color = DoubProxColor;
%Colors(38).Color = DoubDistColor;
Colors(38).Color = DoubDistColor;
Colors(39).Color = DoubDistColor;
Colors(40).Color = BothColor;
Colors(41).Color=BothColor;

% Set font/display parameters
FontUsed=input('Want larger font?','s');
if FontUsed=='y'
    fontsize=15;
else
fontsize=10;
end
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;

FigDirect=[DropboxFolder filesep 'Figures'];

%% Compare singles 
figure 
errorbar([0.5:0.5:50],FractAct(1).APbin(21).AllFracts, FractAct(1).APbin(21).Fract95CI, 'Color',Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(2).APbin(21).AllFracts, FractAct(2).APbin(21).Fract95CI, 'Color',Colors(2).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
title([num2str(EggLength(21)), '%', ' ','egg length']);

% 35% Egglength
figure 
errorbar([0.5:0.5:50],FractAct(1).APbin(15).AllFracts, FractAct(1).APbin(15).Fract95CI, 'Color',Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(2).APbin(15).AllFracts, FractAct(2).APbin(15).Fract95CI, 'Color',Colors(2).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
title([num2str(EggLength(15)), '%', ' ','egg length']);

% 60% egg length
figure 
errorbar([0.5:0.5:50],FractAct(1).APbin(25).AllFracts, FractAct(1).APbin(25).Fract95CI, 'Color',Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(2).APbin(25).AllFracts, FractAct(2).APbin(25).Fract95CI, 'Color',Colors(2).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
title([num2str(EggLength(25)), '%', ' ','egg length']);

% Across AP at 30min into nc14 
figure
errorbar(EggLength,FractAct(1).TimeFract(60).TimeAPFract,FractAct(1).TimeFract(60).TimeAPCI,'Color',Colors(1).Color, 'LineWidth',3.5);
hold on 
errorbar(EggLength,FractAct(2).TimeFract(60).TimeAPFract,FractAct(2).TimeFract(60).TimeAPCI,'Color',Colors(2).Color, 'LineWidth',3.5);
errorbar(EggLength,FractAct(9).TimeFract(60).TimeAPFract,FractAct(9).TimeFract(60).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',3.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(30), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvSingsFractON_30min','.pdf'],'pdf');

% 10 min into nc14
figure
errorbar(EggLength,FractAct(1).TimeFract(20).TimeAPFract,FractAct(1).TimeFract(20).TimeAPCI,'Color',Colors(1).Color, 'LineWidth',3.5);
hold on 
errorbar(EggLength,FractAct(2).TimeFract(20).TimeAPFract,FractAct(2).TimeFract(20).TimeAPCI,'Color',Colors(2).Color, 'LineWidth',3.5);
errorbar(EggLength,FractAct(9).TimeFract(20).TimeAPFract,FractAct(9).TimeFract(20).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',3.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(10), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvSingsFractON_10min','.pdf'],'pdf');

% 20 min into nc14
figure
errorbar(EggLength,FractAct(1).TimeFract(40).TimeAPFract,FractAct(1).TimeFract(40).TimeAPCI,'Color',Colors(1).Color, 'LineWidth',3.5);
hold on 
errorbar(EggLength,FractAct(2).TimeFract(40).TimeAPFract,FractAct(2).TimeFract(40).TimeAPCI,'Color',Colors(2).Color, 'LineWidth',3.5);
errorbar(EggLength,FractAct(9).TimeFract(40).TimeAPFract,FractAct(9).TimeFract(40).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',3.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(20), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvSingsFractON_20min','.pdf'],'pdf');

% 40 min into nc14 
figure
errorbar(EggLength,FractAct(1).TimeFract(80).TimeAPFract,FractAct(1).TimeFract(80).TimeAPCI,'Color',Colors(1).Color, 'LineWidth',3.5);
hold on 
errorbar(EggLength,FractAct(2).TimeFract(80).TimeAPFract,FractAct(2).TimeFract(80).TimeAPCI,'Color',Colors(2).Color, 'LineWidth',3.5);
errorbar(EggLength,FractAct(9).TimeFract(80).TimeAPFract,FractAct(9).TimeFract(80).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',3.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(40), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvSingsFractON_40min','.pdf'],'pdf');
%% Compare doubles vs SE
figure 
errorbar([0.5:0.5:50],FractAct(6).APbin(21).AllFracts, FractAct(6).APbin(21).Fract95CI, 'Color',Colors(6).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(7).APbin(21).AllFracts, FractAct(7).APbin(21).Fract95CI, 'Color',Colors(7).Color, 'LineWidth', 2.5);
errorbar([0.5:0.5:50],FractAct(9).APbin(21).AllFracts, FractAct(9).APbin(21).Fract95CI, 'Color',Colors(9).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
xlim([0 40]);
title([num2str(EggLength(21)), '%', ' ','egg length']);

% 35% egglength 
figure 
errorbar([0.5:0.5:50],FractAct(6).APbin(15).AllFracts, FractAct(6).APbin(15).Fract95CI, 'Color',Colors(6).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(7).APbin(15).AllFracts, FractAct(7).APbin(15).Fract95CI, 'Color',Colors(7).Color, 'LineWidth', 2.5);
errorbar([0.5:0.5:50],FractAct(9).APbin(15).AllFracts, FractAct(9).APbin(15).Fract95CI, 'Color',Colors(9).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
xlim([0 40]);
title([num2str(EggLength(15)), '%', ' ','egg length']);

% 60% egglength
figure 
errorbar([0.5:0.5:50],FractAct(6).APbin(25).AllFracts, FractAct(6).APbin(25).Fract95CI, 'Color',Colors(6).Color, 'LineWidth', 2.5);
hold on 
errorbar([0.5:0.5:50],FractAct(7).APbin(25).AllFracts, FractAct(7).APbin(25).Fract95CI, 'Color',Colors(7).Color, 'LineWidth', 2.5);
errorbar([0.5:0.5:50],FractAct(9).APbin(25).AllFracts, FractAct(9).APbin(25).Fract95CI, 'Color',Colors(9).Color, 'LineWidth', 2.5);
ylabel('fraction of nuclei active');
xlabel('time into nc14 (min)');
xlim([0 40]);
title([num2str(EggLength(25)), '%', ' ','egg length']);

% Across AP at 30min into nc14 
figure
errorbar(EggLength,FractAct(6).TimeFract(60).TimeAPFract,FractAct(6).TimeFract(60).TimeAPCI,'Color',Colors(6).Color, 'LineWidth',4.5);
hold on 
errorbar(EggLength,FractAct(7).TimeFract(60).TimeAPFract,FractAct(7).TimeFract(60).TimeAPCI,'Color',Colors(7).Color, 'LineWidth',4.5);
errorbar(EggLength,FractAct(9).TimeFract(60).TimeAPFract,FractAct(9).TimeFract(60).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',4.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(30), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvDoubFractON_30min','.pdf'],'pdf');

% 20 min into nc14
figure
errorbar(EggLength,FractAct(6).TimeFract(40).TimeAPFract,FractAct(6).TimeFract(40).TimeAPCI,'Color',Colors(6).Color, 'LineWidth',4.5);
hold on 
errorbar(EggLength,FractAct(7).TimeFract(40).TimeAPFract,FractAct(7).TimeFract(40).TimeAPCI,'Color',Colors(7).Color, 'LineWidth',4.5);
errorbar(EggLength,FractAct(9).TimeFract(40).TimeAPFract,FractAct(9).TimeFract(40).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',4.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80]);
title([num2str(20), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvDoubFractON_20min','.pdf'],'pdf');

% 40 min into nc14 
figure
errorbar(EggLength,FractAct(6).TimeFract(80).TimeAPFract,FractAct(6).TimeFract(80).TimeAPCI,'Color',Colors(6).Color, 'LineWidth',4.5);
hold on 
errorbar(EggLength,FractAct(7).TimeFract(80).TimeAPFract,FractAct(7).TimeFract(80).TimeAPCI,'Color',Colors(7).Color, 'LineWidth',4.5);
errorbar(EggLength,FractAct(9).TimeFract(80).TimeAPFract,FractAct(9).TimeFract(80).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',4.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80])
title([num2str(40), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvDoubFractON_40min','.pdf'],'pdf');

% 10 min into nc14 
figure
errorbar(EggLength,FractAct(6).TimeFract(20).TimeAPFract,FractAct(6).TimeFract(20).TimeAPCI,'Color',Colors(6).Color, 'LineWidth',4.5);
hold on 
errorbar(EggLength,FractAct(7).TimeFract(20).TimeAPFract,FractAct(7).TimeFract(20).TimeAPCI,'Color',Colors(7).Color, 'LineWidth',4.5);
errorbar(EggLength,FractAct(9).TimeFract(20).TimeAPFract,FractAct(9).TimeFract(20).TimeAPCI,'Color',Colors(9).Color, 'LineWidth',4.5);
ylabel('fraction of nuclei active');
xlabel('% egg length');
xlim([20 80])
title([num2str(10), 'minutes into nc14']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SEvDoubFractON_10min','.pdf'],'pdf');