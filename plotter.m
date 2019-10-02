% Simulates model under different 
% assumptions and creates figures,
% -- Inputs:
% model: proximal, distal, shadow, etc (see switch 
% statement below)
% experimental: 1 if you want to plot experimental plots,
% 0 otherwise.
% figCounter: keeps track of figure number 
% bestVals1: parameter set for proximal
% bestVals2: parameter set for distal
% -- Outputs:
% figCounter: returns a track of figure number
% to the top level function. 
% Alvaro Fletcher
% Last updated: 9/01/19
function figCounter = plotter(model, ...
    experimental, figCounter, bestVals1, bestVals2)


    % number of nuclei per bin
    extra.numSpots = 80;
    % number of consecutive simulations to compute
    % confidence interval
    maxIter = 5;

    % defines size of figures to be plotted
    xLength = 900;
    yLength = 700;

    degIn1 = bestVals1(1);
    koffIn1 = bestVals1(2);
    konIn1 = bestVals1(3);
    nIn1 = bestVals1(4);
    scalingIn1 = bestVals1(5);
    sdIn1 = bestVals1(6);

    degIn2 = bestVals2(1);
    koffIn2 = bestVals2(2);
    konIn2 = bestVals2(3);
    nIn2 = bestVals2(4);
    scalingIn2 = bestVals2(5);
    sdIn2 = bestVals2(6);


    % defaults
    fixedTF = false;
    withoutBursts = false;
    duplicatedSim = false;
    temperature = false;
    hot = false;
    shadow = false;  
    onlyCorr = false;
    tempDiff = false;
    distAtProx = true;
    % we use the duplicated functions to get 
    % the correlation 'for free'. However,
    % we only want to calculate it when we 
    % are not simulating enhancer competition
    corrBool = true;                            
    
    switch model
        case 'proxfixed'
            fixedTF = true;
            model = 'proximal';
            colorOfChoice =  [238 123 23]./255;
        case 'distfixed'
            fixedTF = true;
            model = 'distal';
            colorOfChoice =   [1 64 172]./255;
        case 'shadowfixed'
            model = 'shadow';
            fixedTF = true;
            shadow = true;
            onlyCorr = true;
            colorOfChoice =   [52 119 71]./255;    
        case 'proxnequal1'
            withoutBursts = true;
            model = 'proximal';
            colorOfChoice =  [238 123 23]./255;
        case 'distnequal1'
            withoutBursts = true;
            model = 'distal';
            colorOfChoice =   [1 64 172]./255;
        case 'shadownequal1'
            model = 'shadow';
            withoutBursts = true;
            shadow = true;
            onlyCorr = true;
            colorOfChoice =   [52 119 71]./255;  
        case 'prox2x'
             duplicatedSim = true;
             model = 'proximal';
             corrBool = false;
             colorOfChoice =  [215, 183, 58]./255;
        case 'dist2x'
            duplicatedSim = true;
            model = 'distal';
            corrBool = false;
            colorOfChoice =  [73 184 253] ./ 255;
        case 'shadow'
             shadow = true;
             corrBool = false;
             colorOfChoice =   [52 119 71]./255;
        case 'proxhot'
            hot = true;
            model = 'proximal';
            temperature = true;
            corrBool = false;
            colorOfChoice =  [255 50 23]./255;
        case 'disthot'
            hot = true;
            model = 'distal';
            temperature = true;
            corrBool = false;
            distAtProx = false;
            colorOfChoice =   [116 64 172]./255;
        case 'proxcold'
            hot = false;
            model = 'proximal';
            temperature = true;
            corrBool = false;
            colorOfChoice=  [238 123 185]./255; % pink
        case 'distcold'
            hot = false;
            model = 'distal';
            temperature = true;
            corrBool = false;
            distAtProx = false;
            colorOfChoice =   [0 0 255]./255; 
        case 'proxhotdiff'
            hot = true;
            model = 'proximal';
            temperature = true;
            corrBool = false;
            tempDiff = true;
            colorOfChoice =  [255 50 23]./255;
        case 'disthotdiff'
            hot = true;
            model = 'distal';
            temperature = true;
            corrBool = false;
            tempDiff = true;
            distAtProx = false;
            colorOfChoice =   [116 64 172]./255;
        case 'proxcolddiff'
            hot = false;
            model = 'proximal';
            temperature = true;
            corrBool = false;
            tempDiff = true;
            colorOfChoice=  [238 123 185]./255; % pink
        case 'distcolddiff'
            hot = false;
            model = 'distal';
            temperature = true;
            corrBool = false;
            tempDiff = true;
            distAtProx = false;
            colorOfChoice =   [0 0 255]./255; 
        case 'shadowcorr'
            shadow = true;
            model = 'shadow';
            onlyCorr = true;
            colorOfChoice =   [52 119 71]./255;
        case 'edist'
            model = 'distal';
            distAtProx = false;
            colorOfChoice =  [208 76 240] ./255;
            corrBool = false;
        case 'proximal'          
            colorOfChoice =  [238 123 23]./255;
        case 'distal'
            colorOfChoice =   [1 64 172]./255; % blue distal

    end
    
    % parameters for simulation
    upperEgg = 80; 
    lowerEgg = 20; 

    if corrBool
        % we start at 2 because of .,.., 
        % loading correlation data
        files = dir(strcat(pwd,'/allCorr'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allCorr/',files(myFiles).name));
        end

        DistalCorr = myStruct(1).files.DistalCorr;
        ProximalCorr = myStruct(2).files.ProximalCorr;
        ShadowCorr = myStruct(3).files.SepCorr;

    end

    if temperature

        if hot

            files = dir(strcat(pwd,'/allHotData'));
            for myFiles = 3:length(files)
                myStruct(myFiles - 2).files = load(strcat(pwd,'/allHotData/',files(myFiles).name));
            end


            DistBurstSize = myStruct(1).files.Dist32CBurstSize;
            DistDuration = myStruct(2).files.Dist32CDuration;
            DistFreq = myStruct(3).files.Dist32CFreq;
            DistIF = myStruct(5).files.Dist32CTotalProd;


            ProxBurstSize = myStruct(6).files.Prox32CBurstSize;
            ProximalDuration = myStruct(7).files.Prox32CDuration;
            ProxFreq = myStruct(8).files.Prox32CFreq;
            ProxIF = myStruct(10).files.Prox32CTotalProd;


        else

            files = dir(strcat(pwd,'/allColdData'));
            for myFiles = 3:length(files)
                myStruct(myFiles - 2).files = load(strcat(pwd,'/allColdData/',files(myFiles).name));
            end

            DistBurstSize = myStruct(1).files.Dist17CBurstSize;
            DistDuration = myStruct(2).files.Dist17CDuration;
            DistFreq = myStruct(3).files.Dist17CFreq;
            DistIF = myStruct(5).files.Dist17CTotalProd;


            ProxBurstSize = myStruct(6).files.Prox17CBurstSize;
            ProximalDuration = myStruct(7).files.Prox17CDuration;
            ProxFreq = myStruct(8).files.Prox17CFreq;
            ProxIF = myStruct(10).files.Prox17CTotalProd;

        end


    elseif duplicatedSim

        files = dir(strcat(pwd,'/allDoubData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allDoubData/',files(myFiles).name));
        end

        % the addition of endog here is meaningless but it 
        % helps keep the code simple
        DistBurstSize = myStruct(1).files.DoubDistBurstSize;
        DistDuration = myStruct(2).files.DoubDistDuration;
        DistFreq = myStruct(3).files.DoubDistFreq;
        DistIF = myStruct(4).files.DoubDistTotalProd;
        DistCV =  myStruct(9).files.DoubDistTempCV;

        ProxBurstSize = myStruct(5).files.DoubProxBurstSize;
        ProximalDuration = myStruct(6).files.DoubProxDuration;
        ProxFreq = myStruct(7).files.DoubProxFreq;
        ProxIF = myStruct(8).files.DoubProxTotalProd;
        ProxCV =  myStruct(10).files.DoubProxTempCV;


    elseif shadow

        files = dir(strcat(pwd,'/allShadowData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allShadowData/',files(myFiles).name));
        end

        ShadowBurstSize = myStruct(1).files.BothBurstSize;
        ShadowDuration = myStruct(2).files.BothDuration;
        ShadowFreq = myStruct(3).files.BothFreq;
        ShadowIF = myStruct(4).files.BothTotalProd;
        ShadowCV = myStruct(5).files.BothTempCV;

    % loading single enhancer data
    else     

        files = dir(strcat(pwd,'/allData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allData/',files(myFiles).name));
        end

        ProxBurstSize = myStruct(9).files.ProxBurstSize;
        ProximalDuration = myStruct(10).files.ProximalDuration;
        ProxFreq = myStruct(11).files.ProxFreq;
        ProxIF = myStruct(13).files.ProxTotalProd;
        ProxCV = myStruct(16).files.ProxTempCV;
        % fixing the duration vectors that are not the right dims 
        ProximalDuration = ProximalDuration';

        if distAtProx

            files = dir(strcat(pwd,'/allDistalData'));
            for myFiles = 3:length(files)
                myStruct(myFiles - 2).files = load(strcat(pwd,'/allDistalData/',files(myFiles).name));
            end

            DistBurstSize = myStruct(1).files.DistBurstSize;
            DistDuration = myStruct(2).files.DistalDuration;
            DistFreq = myStruct(3).files.DistFreq;
            DistIF = myStruct(5).files.DistTotalProd; 
            DistCV =  myStruct(6).files.DistTempCV; 

            DistDuration = DistDuration';

        else

            DistBurstSize = myStruct(3).files.EndogDistBurstSize;
            DistDuration = myStruct(4).files.EndogDistDuration;
            DistFreq = myStruct(5).files.EndogDistFreq;
            DistIF = myStruct(7).files.EndogDistTotalProd;
            DistCV = myStruct(15).files.EndogDistTempCV;

        end

    end

    % we are dividing raw fluorescence data by F_rnap or F_1 since our model
    % operates in transcript count.
    F1 = 1338;
    switch model
        case 'proximal'  
            prox = true;
            expFreq = ProxFreq;
            expDur = ProximalDuration;
            if corrBool
                expCorr = ProximalCorr;
            end
            expSize = ProxBurstSize/F1;
            expIF = ProxIF/F1;
            if ~temperature
                expCV = ProxCV;
            end
        case 'distal'  
            prox = false; 
            expFreq = DistFreq;
            expDur = DistDuration; 
            if corrBool
                expCorr = DistalCorr;
            end
            expSize = DistBurstSize/F1;
            expIF = DistIF/F1;
            if ~temperature
                expCV = DistCV;
            end
        case 'shadow'
            prox = false;
            expFreq =  ShadowFreq;
            expDur = ShadowDuration;
            if corrBool
                expCorr = ShadowCorr;
            end
            expSize = ShadowBurstSize/F1;
            expIF = ShadowIF/F1;
            expCV = ShadowCV;
        otherwise  
            disp('You gave a nonexistent option.')  
            return
    end

    extra.alpha = 1.9540; 
    extra.endSpan = 50;
    enhancer.mean = 50;
    extra.numBinsRWdata = 26;
    extra.eggLength = linspace(lowerEgg, upperEgg, extra.numBinsRWdata);
    extra.ElapsedTime = 0:0.5:extra.endSpan;

    if prox
       if tempDiff && temperature
            if hot 
               enhancer.rFull =  [0  0  0  0  0  0  0  0 ...     
               77.4100   67.1404   74.3840   90.3404  103.1877 ...
               98.2185   87.2952   76.7283   69.6812   65.7044 ...
               56.2718  0  0  0  0  0  0  0];
            else
                enhancer.rFull = [0  0  0  0  0  0  0  0  0  0 ... 
                   73.8782  104.0331  105.6240  101.5717   93.0417 ...
                   83.3756   79.1392   68.0422   70.5296   71.0012 ...
                   0  0  0  0  0  0];
            end
        % regular value of r for single prox
        else
           enhancer.rFull =  [0   0   67.0257   60.0482   57.4959   51.5140   53.9837   53.1280 ...
               60.4217   90.3269  104.0844  114.0616  119.2144  101.4468   87.6735   79.4335 ...
               68.7479   63.0808   86.9607   0   0   0   0   0 ...
               0   0];
        end
    % proximal at distal case 
    elseif distAtProx
        enhancer.rFull = [0   0   0   0   0   0   0   0   58.8609 ...
           74.4993   86.2019  111.3550  141.3504  131.6084  127.4489  115.2295   98.5590 ...
           95.3566   81.3611   71.7599   67.8933   59.8787   60.7270   59.1140   62.3468 ... 
           64.2655 ];
    % the endogenous distal case
    else

        if tempDiff && temperature
            if hot
                enhancer.rFull = [0  0  0  0  0  0  0 48.6836   62.6932   67.8274   73.2909  101.7584 ...
                  135.9031  115.3749  123.7736  113.7628   94.9080 ...
                   82.5860   67.0451   63.8155   58.2935   56.8102 ...
                   54.3845   53.3738   38.7997 0];
            else
                enhancer.rFull = [0  0  0  0  0  0  0  0  0  0 ... 
                  109.4928  143.5294  158.0884  158.3929  151.3082 ...
                  141.1055  128.9081  117.1859   98.6541   83.0483 ...
                   79.4257   71.1059   69.2333   59.6465   65.1885 ...
                   67.8755];
            end
        % regular value of distal for single distal
        else            
            enhancer.rFull = [0    0    0    0    0     ...
               40.4270   47.7268   56.6218   57.3002  72.3034  95.2149  123.3258  134.0092  ... 
               132.5551  121.6160 102.2820   87.8931   71.4522   65.2000   58.7315 ...
               56.1441   50.0337   0   0   0  0];
        end
    end

    % setting the r of the proximal and distal at their endogenous locations if
    % annealing to shadow
    if shadow 
        % when doing correlation we want proximal and distal at proximal location
        if corrBool
           enhancer.rFull1 =  [0   0   67.0257   60.0482   57.4959   51.5140   53.9837   53.1280 ...
               60.4217   90.3269  104.0844  114.0616  119.2144  101.4468   87.6735   79.4335 ...
               68.7479   63.0808   86.9607   0   0   0   0   0 ...
               0   0];
           enhancer.rFull2 = [0   0   0   0   0   0   0   0   58.8609 ...
               74.4993   86.2019  111.3550  141.3504  131.6084  127.4489  115.2295   98.5590 ...
               95.3566   81.3611   71.7599   67.8933   59.8787   60.7270   59.1140   62.3468 ... 
               64.2655 ];
        % when doing all else we want r for edistal
        else
           enhancer.rFull1 =  [0   0   67.0257   60.0482   57.4959   51.5140   53.9837   53.1280 ...
               60.4217   90.3269  104.0844  114.0616  119.2144  101.4468   87.6735   79.4335 ...
               68.7479   63.0808   86.9607   0   0   0   0   0 ...
               0   0];
           enhancer.rFull2 = [0    0    0    0    0     ...
               40.4270   47.7268   56.6218   57.3002  72.3034  95.2149  123.3258  134.0092  ... 
               132.5551  121.6160 102.2820   87.8931   71.4522   65.2000   58.7315 ...
               56.1441   50.0337   0   0   0  0];
        end
    % for a duplicated distal case we want the r of the distal at prox followed by the 
    % r of the endogenous distal
    elseif duplicatedSim && ~prox
        enhancer.rFull1 = [0   0   0   0   0   0   0   0   58.8609 ...
           74.4993   86.2019  111.3550  141.3504  131.6084  127.4489  115.2295   98.5590 ...
           95.3566   81.3611   71.7599   67.8933   59.8787   60.7270   59.1140   62.3468 ... 
           64.2655 ];
        enhancer.rFull2 = [0    0    0    0    0     ...
           40.4270   47.7268   56.6218   57.3002  72.3034  95.2149  123.3258  134.0092  ... 
           132.5551  121.6160 102.2820   87.8931   71.4522   65.2000   58.7315 ...
           56.1441   50.0337   0   0   0  0];
    end

    if duplicatedSim

        enhancer.k_on1 = konIn1;
        enhancer.k_off1 = koffIn1; 
        enhancer.k_on2 = konIn2;
        enhancer.k_off2 = koffIn2; 

        enhancer.degTF = degIn1;
        enhancer.tf = nIn1;
        enhancer.scaling = scalingIn1;
        enhancer.sd =  sdIn1;

    % using the lowest energy sets from single annealing
    elseif shadow 
            % 1 = prox in sifdc
            % 2 = dist in sifdc    
            enhancer.k_on1 = konIn1;
            enhancer.k_off1 = koffIn1; 
            enhancer.k_on2 = konIn2;
            enhancer.k_off2 = koffIn2; 

            enhancer.degTF1 = degIn1;
            enhancer.tf1 = nIn1;
            enhancer.scaling1 =  scalingIn1;
            enhancer.sd1 = sdIn1;

            enhancer.degTF2 = degIn2;
            enhancer.tf2 = nIn2;
            enhancer.scaling2 = scalingIn2;
            enhancer.sd2 =  sdIn2;

    % conditions for single case given by analysis    
    else

        enhancer.degTF = degIn1;
        enhancer.tf = nIn1;
        enhancer.scaling = scalingIn1;
        enhancer.sd =  sdIn1;
        enhancer.k_on = konIn1;
        enhancer.k_off = koffIn1; 

    end

    % determines the width of the lines in the simulated plots
    width = 2; 
    % determines the width of the lines in the experimental plots
    widthExp = 1;
    % determines the size of the font in the plots
    fontSize = 10;  
    zVal = 1.96;
    eggLength = extra.eggLength;
    numBins = extra.numBinsRWdata;
    props = ["amps","freq","dur","corrs","IF","CV","sizes"];

    simMat = zeros(length(props),numBins,maxIter);  
    % dimension with the number of iterations
    itersDim = 3;

    % indexes of each property on the list above
    freqi = 2;
    duri = 3;
    corri = 4;
    ifi = 5;
    cvi = 6;
    sizei = 7; 
     
    for plotIter = 1:maxIter
        simMat(:,:,plotIter) = mainfunction(enhancer, extra);
    end

    % getting the means across all iterations
    meanGiven = nanmean(simMat(:,:,:),itersDim);
    % 95% confidence interval
    confGiven = zVal * nanstd(simMat(:,:,:),0,itersDim)/sqrt(maxIter);

    if corrBool 
        figure(100)

        errorbar(eggLength,meanGiven(corri,:),confGiven(corri,:),'Color',colorOfChoice, 'LineWidth',width)
        hold on
        xlim([0 100])
        yticks([-0.2 0 0.25 0.5 0.75 1])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Allele Correlation','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

    end

    if ~onlyCorr

        figure(101)

        errorbar(eggLength,meanGiven(cvi,:),confGiven(cvi,:),'Color',colorOfChoice, 'LineWidth',width)
        hold on
        axis([0 100 0 7])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Coefficient of Variation','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label


        figure(102)

        subplot(2,2,1)

        errorbar(eggLength, meanGiven(ifi,:),confGiven(ifi,:),'Color',colorOfChoice, 'LineWidth',width) 
        hold on
        xlim([0 100])
        ylim([0 inf])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Integrated Fluorescence','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

        subplot(2,2,2)

        errorbar(eggLength, meanGiven(freqi,:),confGiven(freqi,:),'Color',colorOfChoice, 'LineWidth',width) 
        hold on
        xlim([0 100])
        ylim([0 inf])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Bursts per Minute','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

        subplot(2,2,3)

        errorbar(eggLength, meanGiven(sizei,:),confGiven(sizei,:),'Color',colorOfChoice, 'LineWidth',width) 
        hold on
        xlim([0 100])
        ylim([0 inf])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Fluorescence Intensity','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

        subplot(2,2,4)

        errorbar(eggLength, meanGiven(duri,:),confGiven(duri,:),'Color',colorOfChoice, 'LineWidth',width) 
        hold on
        xlim([0 100])
        ylim([1 4])
        xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
        ylabel('Burst duration (min)','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

        set(gcf, 'Position',  [100, 100, xLength, yLength])

    end

    if experimental 

        if onlyCorr

            corrColor =   [52 119 71]./255;
            figure (51)
            errorbar(eggLength,meanGiven(corri,:),confGiven(corri,:),'Color',colorOfChoice, 'LineWidth',width)
            hold on
            xlim([0 100])
            yticks([-0.2 0 0.25 0.5 0.75 1])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
            ylabel('Allele Correlation','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label
    
            hold on

            figure(51)
            % makes 41 evenly spaced points from 0 to 100
            eggLength = linspace(0,100,41);
            plot(eggLength, expCorr, ':','Color', corrColor, 'LineWidth',widthExp)
            xlim([0 100])

        else
            
            numRows = 2;
            numCols = 2;

            figure(figCounter)

            plottingCount = 1;

            % integrated fluorescence
            subplot(numRows, numCols, plottingCount)
            
            hold on
            
            % extra.eggLength describes unit length steps from the 20 percentile length of the egg
            % to the 80th percentile. 
            errorbar(eggLength, meanGiven(ifi,:),confGiven(ifi,:),'Color',colorOfChoice, 'LineWidth',width) 
            hold on
            xlim([0 100])
            ylim([0 inf])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
            ylabel('Integrated Fluorescence','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label


            plottingCount =  plottingCount + 1;


            % plotting frequencies
            subplot(numRows, numCols, plottingCount)

            hold on

            errorbar(eggLength, meanGiven(freqi,:),confGiven(freqi,:),'Color',colorOfChoice, 'LineWidth',width) 
            hold on
            xlim([0 100])
            ylim([0 inf])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
            ylabel('Bursts per Minute','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label



            plottingCount =  plottingCount + 1;
             

            % plotting durations
            subplot(numRows, numCols, plottingCount)

            hold on


            errorbar(eggLength, meanGiven(duri,:),confGiven(duri,:),'Color',colorOfChoice, 'LineWidth',width) 
            hold on
            xlim([0 100])
            ylim([1 4])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
            ylabel('Burst duration (min)','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

            plottingCount =  plottingCount + 1;

            % plotting burst size
            subplot(numRows, numCols, plottingCount)

            hold on

            errorbar(eggLength, meanGiven(sizei,:),confGiven(sizei,:),'Color',colorOfChoice, 'LineWidth',width) 
            hold on
            xlim([0 100])
            ylim([0 inf])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
            ylabel('Fluorescence Intensity','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label

            set(gcf, 'Position',  [100, 100, xLength, yLength])
            hold on

            figure(50)

            errorbar(eggLength,meanGiven(cvi,:),confGiven(cvi,:),'Color',colorOfChoice, 'LineWidth',width)
            hold on
            axis([0 100 0 7])
            xlabel('% egg length','fontweight','bold','fontsize',fontSize,'FontName','Arial') % x-axis label
            ylabel('Coefficient of Variation','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label
           
            hold on

            % correlation comes for free when both can bind at once
            if corrBool
                    if prox
                        corrColor = [238 123 23]./255;
                    else
                        corrColor =   [1 64 172]./255;                 
                    end     
                    figure (51)
                    errorbar(eggLength,meanGiven(corri,:),confGiven(corri,:),'Color',colorOfChoice, 'LineWidth',width)
                    hold on
                    xlim([0 100])
                    yticks([-0.2 0 0.25 0.5 0.75 1])
                    xlabel('% egg length','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % x-axis label
                    ylabel('Allele Correlation','fontweight','bold','fontsize',fontSize, 'FontName', 'Arial') % y-axis label
                    hold on

            else 
                    % placeholder
                    corrColor = [255 255 255]./255;

            end

            % color of experimental data 
            experimentalColor =  colorOfChoice;

            plottingCount = 1;
            % makes 41 evenly spaced points from 0 to 100
            eggLength = linspace(0,100,41);
            % IF
            figure(figCounter)
            subplot(numRows, numCols, plottingCount)
            plot(eggLength, expIF,':','Color', experimentalColor, 'LineWidth',widthExp)
            xlim([0 100])
            plottingCount =  plottingCount + 1;
            % Freq
            subplot(numRows, numCols, plottingCount)
            plot(eggLength, expFreq, ':','Color', experimentalColor, 'LineWidth',widthExp)
            xlim([0 100])
            plottingCount =  plottingCount + 1;

            % Dur
            subplot(numRows, numCols, plottingCount)
            plot(eggLength, expDur, ':','Color', experimentalColor, 'LineWidth',widthExp)
            xlim([0 100])
            plottingCount =  plottingCount + 1;
            % sizes
            subplot(numRows, numCols, plottingCount)
            plot(eggLength, expSize, ':','Color', experimentalColor, 'LineWidth', widthExp)
            xlim([0 100])

            set(gcf, 'Position',  [100, 100, xLength, yLength])
            % CV


            if ~temperature 
                figure(50)
                plot(eggLength, expCV,':','Color', experimentalColor, 'LineWidth',widthExp)
                xlim([0 100])
            end

            if corrBool
                figure(51)
                plot(eggLength, expCorr, ':','Color', corrColor, 'LineWidth',widthExp)
                xlim([0 100])
            end

            figCounter = figCounter + 1;

        end

    end 


    % Generates a vector of values that correspond to 
    % the values of mRNA production rates based on the assumption
    % of using a bell curve as production rate across the egg.
    % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % prodTF: the production rate of the TF
    % Alvaro Fletcher
    function prodTF = TFproduction(enhancer, extra)

        % range of egg length for beta and gamma
        tfProdLength = extra.eggLength;

        prodTF = enhancer.scaling * mygauss(tfProdLength, enhancer.mean, enhancer.sd);

        function fGauss = mygauss(gaussVar, muGauss, sGauss)

            p1Gauss = -.5 * (( (gaussVar)  - muGauss)/sGauss) .^ 2;
            p2Gauss = (sGauss * sqrt(2 * pi));
            fGauss = exp(p1Gauss) ./ p2Gauss; 

        end 

    end

    % An equivalent version of TFproduction but for two TFs.
    % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % prodTF1: the production rate of the TF 1
    % prodTF2: the production rate of the TF 2
    % Alvaro Fletcher
    function [prodTF1, prodTF2] = TFproductionShadow(enhancer, extra)

        % range of egg length for beta and gamma
        tfProdLength = extra.eggLength;

        prodTF1 = enhancer.scaling1 * mygauss(tfProdLength, enhancer.mean, enhancer.sd1);
        prodTF2 = enhancer.scaling2 * mygauss(tfProdLength, enhancer.mean, enhancer.sd2);

        function fGauss = mygauss(gaussVar, muGauss, sGauss)

            p1Gauss = -.5 * (( (gaussVar)  - muGauss)/sGauss) .^ 2;
            p2Gauss = (sGauss * sqrt(2 * pi));
            fGauss = exp(p1Gauss) ./ p2Gauss; 

        end 

    end  


    % Runs simulations of the model.
    % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % simMatTemp: struct containing burst properties and/or allele
    % correlation of simulated data
    % Alvaro Fletcher
    function simMatTemp = mainfunction(enhancer, extra)

        holderEggLength = length(extra.eggLength);
        holderElapsedTime = extra.ElapsedTime;
        holderNumSpots = extra.numSpots;
        
        % we ignore bins with less than one data point
        ignoreThreshold = 1;  

        % function handles for parfor loop
        fcnSingleEnhancer = @singleEnhancer;
        fcnDupEnhancer = @duplicatedEnhancer;
        fcnDiscrete = @crndiscrete;
        fcnShadow = @shadowPair;
        fcnFixedTF = @singleEnhancerFixed;
        fcnFixedTFDoub = @duplicatedEnhancerFixed;
        fcnFixedTFShadow = @shadowPairFixed;

        if shadow
            [prodTF1, prodTF2] = TFproductionShadow(enhancer, extra);
            prodTF = zeros(1,extra.numBinsRWdata);
        else
            prodTF1 = zeros(1,extra.numBinsRWdata);
            prodTF2 = zeros(1,extra.numBinsRWdata);
            prodTF = TFproduction(enhancer, extra);
        end

        % shadow case
        if shadow || (duplicatedSim && ~prox)
            rateRNA = zeros(1,extra.numBinsRWdata);
            rateRNA1 = enhancer.rFull1;
            rateRNA2 = enhancer.rFull2;
        else
            rateRNA1 = zeros(1,extra.numBinsRWdata);
            rateRNA2 = zeros(1,extra.numBinsRWdata);   
            rateRNA = enhancer.rFull;
        end

        numSpotsPar = extra.numSpots;

        % holds the mRNA values for SpotDiff
        tempSimR = zeros(holderEggLength, length(extra.ElapsedTime), holderNumSpots);
        tempCorr = zeros(1,holderEggLength);
        tempIF = zeros(1,holderEggLength);
        tempCV = zeros(1,holderEggLength);

        % going over different bins in the egg
        parfor i = 1:holderEggLength

            production = prodTF(i);
            myParR = rateRNA(i);
      
            % limits the number of interactions for the Gillespie algorithm.
            % However, we want all simulations to be limited by the time span of
            % nc14 and just choose a very large number for maxIterations
            maxIterations = 150000;

            % this is only for single enhancer cases that fit to all burst properties
            if ~(duplicatedSim || shadow)

                % the locations of the mRNA and transcription factor for our network
                RNAloc = 3;

                parTempSimR = zeros(length(holderElapsedTime), holderNumSpots);
                parTempIF = zeros(1, holderNumSpots);
                parTempCV = zeros(1, holderNumSpots);

                % going over all the spots in the bin
                for spot = 1:numSpotsPar

                    if fixedTF
                        % defining the network for single enhancer with
                        % fixed TF
                        CRN = feval(fcnFixedTF, enhancer, extra, production, myParR);       
                    else 
                        % defining the network for single enhancer with 
                        % fluctuating TFs
                        CRN = feval(fcnSingleEnhancer, enhancer, extra, production, ... 
                            myParR);
                    end

                    % runs the Gillespie algorithm 
                    expN = feval(fcnDiscrete, CRN, maxIterations);
                    
                    % saving the mRNA concentration over time for this simulation
                    R = expN(:,RNAloc); 

                    parTempSimR(:,spot) = R;
                    parTempIF(spot) = trapz(holderElapsedTime, R);
                    parTempCV(spot) = std(R)/mean(R);

                end

                tempSimR(i,:,:) = parTempSimR; 
                tempIF(i) = nanmean(parTempIF);
                tempCV(i) =  nanmean(parTempCV);

            end

            if duplicatedSim || shadow || corrBool

                % picking the production rate of TF for this bin
                production1 = prodTF1(i);
                myParR1 = rateRNA1(i);
                production2 = prodTF2(i);
                myParR2 = rateRNA2(i);

                % note that specifying "prox" and "dist" does not necessarily mean 
                % the proximal enhancer or distal enhancer but rather whichever 
                % enhancers are found in the proximal and distal locations
                RNAlocProx = 5;
                RNAlocDist = 8;

                if ~corrBool
                    parTempSimR = zeros(length(holderElapsedTime), holderNumSpots);
                    parTempIF = zeros(1, holderNumSpots);
                    parTempCV = zeros(1, holderNumSpots);
                end

                parTempCorr = zeros(1, holderNumSpots);
                for spot = 1:numSpotsPar

                    if shadow
                        if fixedTF
                            CRN = feval(fcnFixedTFShadow,enhancer, extra, production1, ...
                                production2, myParR1, myParR2);                          
                        else
                            CRN = feval(fcnShadow, enhancer, extra, production1, ...
                                production2, myParR1, myParR2);
                        end  
                    % duplicated case or correlation
                    else
                        if fixedTF
                            CRN = feval(fcnFixedTFDoub, enhancer, extra, production, myParR);        
                        else
                            if prox
                                CRN = feval(fcnDupEnhancer, enhancer, extra, production, myParR, myParR);
                            % we have r values for distal and edistal
                            else
                                % but for correlation we are calculating it between two distals at the proximal location
                                if corrBool
                                     CRN = feval(fcnDupEnhancer, enhancer, extra, production, myParR, myParR);
                                else                             
                                     CRN = feval(fcnDupEnhancer, enhancer, extra, production, myParR1, myParR2);
                                end
                            end
                        end
                    
                    end

                    % runs the Gillespie algorithm 
                    expN = feval(fcnDiscrete,CRN,maxIterations);      

                    % the RNA locations for the correlation models
                    % that include R1 and R2 (RNA produced by enhancer 1
                    % and enhancer 2 respectively) are selected per 
                    % choice of model at the beginning
                    proxR = expN(:,RNAlocProx);
                    distR = expN(:,RNAlocDist);

                    if ~corrBool

                        R = proxR + distR;

                        parTempSimR(:,spot) = R;
                        parTempIF(spot) = trapz(holderElapsedTime, R);
                        parTempCV(spot) = std(R)/mean(R);

                    else 
                        % storing the correlation between R1 and R2 
                        corrholder = corrcoef(proxR,distR);
                        parTempCorr(spot) = corrholder(1,2);

                    end
                   
                    
                end

                % storing all mRNA data and IF averages into our main 
                % structure
                if ~corrBool
                     tempSimR(i,:,:) = parTempSimR; 
                     tempIF(i) = nanmean(parTempIF);
                     tempCV(i) =  nanmean(parTempCV);
                else
                     tempCorr(i) = nanmean(parTempCorr);
                end

            end    

        end

        simData.corr = tempCorr;
        simData.CV = tempCV;
        simData.IF = tempIF;

        % building SpotDiff and simData.corr after the parallelizations
        SpotDiff = struct('SpotOne',cell(1, holderEggLength * extra.numSpots), ...
                    'APBin',cell(1, holderEggLength * extra.numSpots));

        postCounter = 1;
        for i = 1:length(extra.eggLength)
            for j = 1:extra.numSpots
                SpotDiff(postCounter).SpotOne = tempSimR(i,:,j);
                SpotDiff(postCounter).APBin = extra.eggLength(i);
                postCounter = postCounter + 1;
            end
        end  

        % obtaining data from burst calling algorithm
        BurstProperties = modSlopeBurstCallingRW(SpotDiff,extra.ElapsedTime);

        % Initializing...
        simData.amps = [];
        simData.sizes = [];
        simData.freqs = [];
        simData.durs = [];
        simData.uniqueBins = [];

        % storing all unique values in the AP bin field of BurstProperties
        uniqueAPBins = unique([BurstProperties(1:end).APBin]);

        % removing unnecessary NaNs
        uniqueAPBins = uniqueAPBins(~isnan(uniqueAPBins));

        % putting all the AP bin data into vector form
        APBinData = [BurstProperties(1:end).APBin];

        % sometimes not all AP bins are represented or they can be out of order
        for i = 1:length(uniqueAPBins)

            % finding all indexes in BurstProperties for the current AP bin
            currentAPBinVals =  find(APBinData == uniqueAPBins(i));

            % checking number of empty rows, it suffices to check amplitude
            % since the existence of any burst property implies the 
            % existence of all the others
       
            % tracks the number of empty rows
            nonEmptyCounter = 0;

            % in case nothing happened for every single bin then no 
            % burst property field in the struct is created (unlikely)
            if isfield(BurstProperties, 'BurstSize')
            
                for j = currentAPBinVals
                    
                    if ~isempty(BurstProperties(j).BurstSize) 
                         nonEmptyCounter = nonEmptyCounter + 1;       
                    end
                    
                end

            end

            % we ignore AP bins with not enough activity
            if  nonEmptyCounter > ignoreThreshold

                 simData.amps(end + 1) = ... 
                  nanmean([BurstProperties(currentAPBinVals).BurstAmplitude]); 
                simData.sizes(end + 1) = ... 
                      nanmean([BurstProperties(currentAPBinVals).BurstSize]); 
                simData.durs(end + 1) = ...
                      nanmean([BurstProperties(currentAPBinVals).Duration]);   
                simData.freqs(end + 1) = ...
                       nanmean([BurstProperties(currentAPBinVals).Frequency]); 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);

            else 
                simData.amps(end + 1) = nan;  
                simData.sizes(end + 1) = nan;       
                simData.durs(end + 1) = nan;
                simData.freqs(end + 1) = nan; 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);
           
            end 

        end

        amps = simData.amps;
        freqs = simData.freqs;
        durs = simData.durs;
        IF = simData.IF;
        sizes = simData.sizes;
        CV = simData.CV;
        corr = simData.corr;

        simMatTemp = [amps;freqs;durs;corr;IF;CV;sizes];

    end

    % Defines a single enhancer network to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR: production rate of mRNA
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = singleEnhancer(params, extra, prodTF, simR)

        if withoutBursts
            n1 = 1; 
        else
            % number of TFs per burst
            n1 = params.tf;
        end

        % E = enhancer 1, T = TF, C = complex, R=RNA

        %  E + T  <-> C -> C + R
        %  R -> 0
        %  n*T <- 0 <- T
        
        reaction=[...
            %   E  C  R  T    E  C  R  T   
                1  0  0  1    0  1  0  0  
                0  1  0  0    1  0  0  1  
                0  1  0  0    0  1  1  0     
                0  0  1  0    0  0  0  0  
                0  0  0  0    0  0  0  n1  
                0  0  0  1    0  0  0  0       
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with m entries
        CRN.params = [params.k_on params.k_off simR extra.alpha prodTF params.degTF];  

        %	        E   C   R   T   
        CRN.init = [1   0   0   0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end


    % Defines a single enhancer network with fixed TF numbers
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR: production rate of mRNA
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = singleEnhancerFixed(params, extra, prodTF, simR)

        % we need to convert the gaussian value to an integer
        numTFfixed = round(prodTF);

        % E = enhancer 1, T = TF, C = complex, R=RNA

        %  E + T  <-> C -> C + R
        %  R -> 0
        reaction=[...
            %   E  C  R  T    E  C  R  T   
                1  0  0  1    0  1  0  0  
                0  1  0  0    1  0  0  1  
                0  1  0  0    0  1  1  0     
                0  0  1  0    0  0  0  0      
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with m entries
        CRN.params = [params.k_on params.k_off simR extra.alpha];  

        %           E   C   R   T   
        CRN.init = [1   0   0   numTFfixed];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end


    % Defines a double enhancer network 
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = duplicatedEnhancer(params, extra, prodTF, simR1, simR2)

        if withoutBursts
            n1 = 1; 
        else
            % number of TFs per burst
            n1 = params.tf;
        end

        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  -> C1 -> C1 + R
        %  E2 + T1  -> C2 -> C2 + R
        %  R -> 0
        %  n*T1 <- 0 <- T1

        % it is convenient to leave T2 so that locations of R1
        % and R2 in the  and shadow case are the same
        
        reaction=[...
          % E1 E2 C1 C2 R1 T1 T2 R2   E1 E2 C1 C2 R1 T1  T2 R2 
            1  0  0  0  0  1  0  0    0  0  1  0  0  0   0  0 
            0  0  1  0  0  0  0  0    1  0  0  0  0  1   0  0 
            0  0  1  0  0  0  0  0    0  0  1  0  1  0   0  0 
            0  1  0  0  0  1  0  0    0  0  0  1  0  0   0  0 
            0  0  0  1  0  0  0  0    0  1  0  0  0  1   0  0        
            0  0  0  1  0  0  0  0    0  0  0  1  0  0   0  1
            0  0  0  0  1  0  0  0    0  0  0  0  0  0   0  0 
            0  0  0  0  0  0  0  0    0  0  0  0  0  n1  0  0 
            0  0  0  0  0  1  0  0    0  0  0  0  0  0   0  0 
            0  0  0  0  0  0  0  1    0  0  0  0  0  0   0  0       
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with reactionsNum entries

        % the r's will be the same for correlation purposes and for 2x proximal.
        % however, they are different in 2x distal due to having the data to calculate
        % r of edistal vs. distal 
        if duplicatedSim
            CRN.params=[params.k_on1 params.k_off1 simR1 params.k_on2 params.k_off2 ...
              simR2 extra.alpha prodTF params.degTF extra.alpha];   
        else 
            CRN.params=[params.k_on params.k_off simR1 params.k_on params.k_off ...
              simR1 extra.alpha prodTF params.degTF extra.alpha];  
        end

        %         E1 E2 C1 C2 R1 T1 T2 R2 
        CRN.init=[1  1  0  0  0  0  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

        
    end


    % Defines a double enhancer network with fixed TF number
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % NOTE: this network is only used to obtain allele correlation
    % between homozygotes, as such, we only require one value of r.
    % Alvaro Fletcher
    function CRN = duplicatedEnhancerFixed(params, extra, prodTF, simR)

        % we need to convert the gaussian value to an integer
        numTFfixed = round(prodTF);

        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  -> C1 -> C1 + R
        %  E2 + T1  -> C2 -> C2 + R
        %  R -> 0

        % it is convenient to leave T2 so that locations of R1
        % and R2 in the duplicated and shadow case are the same
        
        reaction=[...
          % E1 E2 C1 C2 R1 T1 T2 R2   E1 E2 C1 C2 R1 T1  T2 R2 
            1  0  0  0  0  1  0  0    0  0  1  0  0  0   0  0 
            0  0  1  0  0  0  0  0    1  0  0  0  0  1   0  0 
            0  0  1  0  0  0  0  0    0  0  1  0  1  0   0  0 
            0  1  0  0  0  1  0  0    0  0  0  1  0  0   0  0 
            0  0  0  1  0  0  0  0    0  1  0  0  0  1   0  0        
            0  0  0  1  0  0  0  0    0  0  0  1  0  0   0  1
            0  0  0  0  1  0  0  0    0  0  0  0  0  0   0  0 
            0  0  0  0  0  0  0  1    0  0  0  0  0  0   0  0       
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with reactionsNum entries

        CRN.params=[params.k_on params.k_off simR params.k_on params.k_off ...
          simR extra.alpha extra.alpha];  

        %         E1 E2 C1 C2 R1 T1         T2 R2 
        CRN.init=[1  1  0  0  0  numTFfixed  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

        
    end



    % Defines a double enhancer network 
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF1: production rate of TF 1
    % prodTF2: production rate of TF 2
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN  = shadowPair(params, extra, prodTF1, prodTF2, simR1, simR2)


        if withoutBursts
            n1 = 1;
            n2 = 1; 
        else
            % number of TFs per burst
            n1 = params.tf1;
            n2 = params.tf2;
        end


        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  <-> C1 -> C1 + R
        %  E2 + T2 <-> C2 -> C2 + R
        %  R -> 0
        %  n1*T1 <- 0 <- T1
        %  n2*T2 <- 0 <- T2 
        
    
        reaction=[...
       %  E1 E2 C1 C2 R1 T1 T2 R2   E1 E2 C1 C2 R1 T1  T2  R2  
          1  0  0  0  0  1  0  0    0  0  1  0  0  0   0   0   
          0  0  1  0  0  0  0  0    1  0  0  0  0  1   0   0   
          0  0  1  0  0  0  0  0    0  0  1  0  1  0   0   0   
          0  1  0  0  0  0  1  0    0  0  0  1  0  0   0   0   
          0  0  0  1  0  0  0  0    0  1  0  0  0  0   1   0         
          0  0  0  1  0  0  0  0    0  0  0  1  0  0   0   1   
          0  0  0  0  1  0  0  0    0  0  0  0  0  0   0   0   
          0  0  0  0  0  0  0  0    0  0  0  0  0  n1  0   0   
          0  0  0  0  0  1  0  0    0  0  0  0  0  0   0   0    
          0  0  0  0  0  0  0  0    0  0  0  0  0  0   n2  0   
          0  0  0  0  0  0  1  0    0  0  0  0  0  0   0   0   
          0  0  0  0  0  0  0  1    0  0  0  0  0  0   0   0   
     
        ];


        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum, 1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum, speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with m entries
        CRN.params = [params.k_on1 params.k_off1 simR1 params.k_on2 ...
        params.k_off2 simR2 extra.alpha ...
        prodTF1 params.degTF1 prodTF2 params.degTF2 extra.alpha];  

        %         E1 E2 C1 C2 R1 T1 T2 R2
        CRN.init=[1  1  0  0  0  0  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end

    % Defines a double enhancer network with fixed TF numbers
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF1: production rate of TF 1
    % prodTF2: production rate of TF 2
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN  = shadowPairFixed(params, extra, prodTF1, prodTF2, simR1, simR2)

        % we need to convert the gaussian value to an integer
        numTFfixed1 = round(prodTF1);
        numTFfixed2 = round(prodTF2);
        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  <-> C1 -> C1 + R
        %  E2 + T2 <-> C2 -> C2 + R
        %  R -> 0
    
        reaction=[...
       %  E1 E2 C1 C2 R1 T1 T2 R2   E1 E2 C1 C2 R1 T1  T2  R2  
          1  0  0  0  0  1  0  0    0  0  1  0  0  0   0   0   
          0  0  1  0  0  0  0  0    1  0  0  0  0  1   0   0   
          0  0  1  0  0  0  0  0    0  0  1  0  1  0   0   0   
          0  1  0  0  0  0  1  0    0  0  0  1  0  0   0   0   
          0  0  0  1  0  0  0  0    0  1  0  0  0  0   1   0         
          0  0  0  1  0  0  0  0    0  0  0  1  0  0   0   1   
          0  0  0  0  1  0  0  0    0  0  0  0  0  0   0   0   
          0  0  0  0  0  0  0  1    0  0  0  0  0  0   0   0   
     
        ];


        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum, 1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum, speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with m entries
        CRN.params = [params.k_on1 params.k_off1 simR1 params.k_on2 ...
        params.k_off2 simR2 extra.alpha extra.alpha];  

        %         E1 E2 C1 C2 R1 T1           T2           R2
        CRN.init=[1  1  0  0  0  numTFfixed1  numTFfixed2  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end
   

    % This function is given aGE CRN data structure (see samplereaction files), 
    % and optionally an initial state.  
    % It carries out the Gillespie algorithm for a fixed number of steps and 
    % produces a plot.  
    % This file has been corrected to have the right reaction propensities when
    % more than one identical molecule is involved.  For instance, if 2A->B,
    % and A=10, then the propensity is A*(A-1)=90, not A^2=100.  
    % --Inputs:
    % CRN: the network to be simulated
    % maxiterations: the maximum number of iterations
    % -- Outputs: 
    % expN: a vector containing simulated mRNA counts at 30 second intervals
    %  German Enciso
    % Last updated: 2009-11-14
    %
    % With some modifications by Alvaro Fletcher on 07/23/2018:
    % 1. Fixed a bug with the limiting time span of the simulation. 
    % 2. This version of crndiscrete also generates data as is recorded
    % experimentally, in this case, every 30 second intervals. 
    function expN = crndiscrete(CRN,maxiterations)

        nGE=CRN.speciesNum;   % number of species
        mGE=CRN.reactionsNum;   % number of reactions
        source=CRN.source;   % source is for source
        target=CRN.target;   % target is for target
        paramsGE=CRN.params;   %  all the parameters in the system
        initGE=CRN.init;
        tspan=CRN.tspan;

        %default maximum number of iterations
        if nargin==1    
            maxiterations=500*nGE*mGE; 
        end

        % make sure initGE is aGE column vector
        temp=size(initGE);

        if temp(1)==1 

            initGE=initGE'; 

        end

        N(1,:) = initGE;   % the initial state;
        tGE(1) = tspan(1);  %  the initial time;


        % rachel records concentrations of different species
        % at specific times so we have uniform recordings for all
        % simulations

        expN(1,:) = N(1,:);
        myTime = 2;

        for iGE = 1:maxiterations-1

            xGE = N(iGE,:);  % vector of concentrations
      
            % going over every reaction
            for k = 1:mGE

                pow=source(k,:);   % vector of powers
                reactionrate(k) = paramsGE(k) * prod( xGE.^pow );

                % correct reaction rate in case that any of the powers is larger than one    
                % going over every species
                for el = 1:nGE

                    if xGE(el) ~= 0

                        for alpha = 2:source(k,el)

                              reactionrate = reactionrate * (xGE(el) - alpha + 1)/xGE(el);

                        end

                    end

                end   

            end

            rTGE = sum(reactionrate); % Total rate
            u_time = rand;         % Uniformly distributed number

            % time(iGE) is the lapse of time before the next reaction
            time(iGE) = -log(u_time)/rTGE;

            timeHolderGE = tGE(iGE) + time(iGE);
                
            % preventing time from going above tspan(2)
            if timeHolderGE > tspan(2)  
                

                 for mycount = floor(tGE(iGE)) + 0.5:0.5:tspan(2) 
                     
                    % setting the remainder of expN to the last value recorded 
                    expN(myTime,:) = N(iGE,:);      
                    myTime = myTime + 1;


                 end   

                 break

            else 
                    
                  tGE(iGE+1) = timeHolderGE;

                  % we choose 0.5 to resemble the experimental data from Rachel

                  for mycount = floor(tGE(iGE)) + 0.5:0.5:floor(tGE(iGE+1))
                      
                        expN(myTime,:) = N(iGE,:);      
                        myTime = myTime + 1;

                  end

            end



            % choose a number numGE defining which reaction took place
            uGE=rTGE*rand;
            aGE=0;

            % mGE is the number of reactions
            for kk=1:mGE

                aGE=aGE+reactionrate(kk);

                if uGE <= aGE 
                        
                      % we choose the kk^th reaction to take place
                      numGE=kk; 
                      break  

                end

            end
            
            N(iGE+1,:) = N(iGE,:) + target(numGE,:) - source(numGE,:);

        end  

    end

    % Calulates burst properties (among other things) 
    % from experimental or simulated data.
    % --Inputs:
    % SpotDiff: a struct containing simulated data for all bins
    % ElapsedTime: a vector containing the time at which each point  
    % was recorded in the traces. 
    % --Outputs:
    % BurstProperties: a struct containing burst properties
    % Rachel Waymack. 
    %
    % Modifications done by Alvaro Fletcher on 07/23/2018:
    % 1. Made the script into a function.
    % 2. Deleted calculation for SpotTwo (allele correlation 
    % is calculated directly from the simulated data)
    % 3. RNAScale and nc14 parameters set to 1 (model operates in transcript count). 
    function [BurstProperties] = modSlopeBurstCallingRW(SpotDiff,ElapsedTime)

        % MOD: added RNA scale to ON threshold and OFF threshold
        RNAscale = 1;
        % nc14 length just set it to 1
        nc14 = 1;
        
        %Slope thresholds
        ONThreshold = RNAscale;   %6/20/18 thinking ~377AU=1 polymerase so 1500 ~ 5ish polymerases
        OFFThreshold = -RNAscale; 

        for n=1:length(SpotDiff)
            if length(SpotDiff(n).SpotOne)==length(ElapsedTime)
               %try limiting frames to nc14 to reduce extreme smoothing
               % MOD: got rid of the nc14 variable here
               % did I put the 100 here?
               if length(ElapsedTime) >= 100
               SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:(nc14+100)),SpotDiff(n).SpotOne(nc14:(nc14+100)),0.1,'lowess');
               else
                   ShortTrace='y';
                   SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotOne(nc14:end),0.1,'lowess');
               end
               %SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotOne(nc14:end),0.1,'lowess');
            %fluorescence should never be negative but sometimes with lots of 0s
            %looks like the smoothing f(x) makes it so
               SmoothParticles(n).Smoothed(SmoothParticles(n).Smoothed<0)=0;
               if exist('ShortTrace') & ShortTrace=='y'
               SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:end))]');
               else
               SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:(nc14+100)))]');
               end
        %Slope thresholds
               %OFFThresholdSlope=dx<=Q;
        % AboveLine=[SmoothParticles(ss).Smoothed];
        % AboveLine(~ONThresholdSlope)=nan;
        % AboveLine(OFFThresholdSlope)=nan;
        % BelowLine=[SmoothParticles(ss).Smoothed];
        % BelowLine(ONThresholdSlope)=nan;
        % BelowLine(~OFFThresholdSlope)=nan;
        StartFramePool=[];
        StartSlopePool=[];
        EndSlopePool=[];
        Starts=[find(SmoothParticles(n).RawSlope>=ONThreshold)];
        Starts=min(Starts);
        StartFramePool=[StartFramePool,Starts];
        EndPts=[find(SmoothParticles(n).RawSlope<=OFFThreshold)];
        EndPts=min(EndPts);
        EndFramePool=[1000,EndPts];
        if ~isempty(EndPts) & (~isempty(Starts))
        while EndFramePool(end) ~= EndFramePool(end-1)
            for ss=EndFramePool(end):length(SmoothParticles(n).RawSlope)
                if SmoothParticles(n).RawSlope(ss) >=ONThreshold
                    %if StartFramePool(end)~=ss
                    StartFramePool=[StartFramePool,ss];
                    StartSlopePool=[StartSlopePool, SmoothParticles(n).RawSlope(ss)];
                    break
                    %end
                end
            end
            for ss=StartFramePool(end):(length(SmoothParticles(n).RawSlope))
                if SmoothParticles(n).RawSlope(ss) <=OFFThreshold
                    EndFramePool=[EndFramePool, ss];
                    EndSlopePool=[EndSlopePool,SmoothParticles(n).RawSlope(ss)];
                    break
                end
            end
            if sum(SmoothParticles(n).RawSlope(EndFramePool(end):end)>=ONThreshold)==0
        %         EndFramePool(end+1)=EndFramePool(end);
        %         StartFramePool(end+1)=StartFramePool(end);
                PotentialPeaks(n).Starts=[StartFramePool];
                PotentialPeaks(n).Ends=[EndFramePool(2:end)];
                break
            elseif sum(SmoothParticles(n).RawSlope(StartFramePool(end):end)<=OFFThreshold)==0
        %         EndFramePool(end+1)=EndFramePool(end);
        %         StartFramePool(end+1)=StartFramePool(end);
                PotentialPeaks(n).Starts=[StartFramePool];

                PotentialPeaks(n).Ends=[EndFramePool(2:end)];
                break
            end
        end
        if length(PotentialPeaks(n).Ends)>1 & (PotentialPeaks(n).Ends(end)==PotentialPeaks(n).Ends(end-1))
            PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(1:end-1)];
        end
        if length(PotentialPeaks(n).Starts)>1 & (PotentialPeaks(n).Starts(end)==PotentialPeaks(n).Starts(end-1))
            PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(1:end-1)];
        end
        BurstProperties(n).SmoothTrace=SmoothParticles(n).Smoothed;
        %Get rid of repeated 1st and 2nd ON frames in case an OFF came before ON in
        %original trace - added 7/27/18
        if (length(PotentialPeaks(n).Starts)>=2)& (PotentialPeaks(n).Starts(1)==PotentialPeaks(n).Starts(2))
            PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(2:end)];
        end
        % If a trace starts with an OFF frame before an ON frame, need to ignore it
        % 7/27/18
        if PotentialPeaks(n).Ends(1) < PotentialPeaks(n).Starts(1)
            if length(PotentialPeaks(n).Ends)==1
                PotentialPeaks(n).Ends=nan;
                PotentialPeaks(n).Starts=nan;
            else
            PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(2:end)];
            end
        end
        if (isnan(PotentialPeaks(n).Ends)) & (isnan(PotentialPeaks(n).Starts))
            continue
        end
        %If start frames > end frames, ignore the last start since we can't say how
        %long the burst lasts
        if length(PotentialPeaks(n).Starts) > (length(PotentialPeaks(n).Ends))
            BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts(1:end-1));
            BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
            BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end-1))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
            BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts(1:end-1)];
        else
        BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts);
        BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
        BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
        BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts];
        end

        BurstProperties(n).OFFFrames=[PotentialPeaks(n).Ends];

        if PotentialPeaks(n).Ends(1)==1
            BurstProperties(n).BurstAmplitude(1)=0;
        else
            BurstProperties(n).BurstAmplitude(1)=max(SmoothParticles(n).Smoothed([1:PotentialPeaks(n).Ends(1)]));
        end

            for pp=2:length(PotentialPeaks(n).Ends)
            BurstProperties(n).BurstAmplitude(pp)=max(SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(pp-1):PotentialPeaks(n).Ends(pp)]));%[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
            BurstProperties(n).BurstAmplitude(pp)=(BurstProperties(n).BurstAmplitude(pp))-(SmoothParticles(n).Smoothed(PotentialPeaks(n).Starts(pp)));
            end

        BurstProperties(n).BurstAmplitude=BurstProperties(n).BurstAmplitude(BurstProperties(n).BurstAmplitude~=0); %get rid of 1st 0 bc above started at pp=2
        
        for pp=1:(length(BurstProperties(n).ONFrames)-1)
        BurstProperties(n).BurstSize(pp)=(trapz([BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1)],BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1))));
        end
        if isfield(BurstProperties, 'BurstSize')
        BurstProperties(n).BurstSize(end+1)=(trapz([BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end)], BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end))));
        end
        % if PotentialPeaks(n).Ends(1)==1
        %     BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(2:end)]-1)]';
        % else
        % BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
        % end
       % temp=find([CompiledParticles.Nucleus]==SpotDiff(n).Nucleus);
       % BurstProperties(n).TotalmRNAError=CompiledParticles(temp(1)).TotalmRNAError;

        BurstProperties(n).FirstTimeOn=ElapsedTime(PotentialPeaks(n).Starts(1));
        BurstProperties(n).NBursts=length(BurstProperties(n).BurstAmplitude);
        %Frequency set btwn time of first ON and 50min into nc14 
        if ((length(ElapsedTime))-nc14)<100
            BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(end)-ElapsedTime(PotentialPeaks(n).Starts(1)));
        else
        BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(nc14+100)-ElapsedTime(PotentialPeaks(n).Starts(1)));
        end
        BurstProperties(n).FrequencyAct=(BurstProperties(n).NBursts)/(ElapsedTime(PotentialPeaks(n).Ends(end))-ElapsedTime(PotentialPeaks(n).Starts(1)));
        BurstProperties(n).APBin=SpotDiff(n).APBin;
       % BurstProperties(n).Nucleus=SpotDiff(n).Nucleus;
        %BurstProperties(n).Interburst=
      %  BurstProperties(n).TotalmRNA=SpotDiff(n).TotalmRNAOne;
        BurstProperties(n).FractON=sum(BurstProperties(n).Duration)/(length(ElapsedTime(nc14:end))); 
        BurstProperties(n).TotalElapsedTime=ElapsedTime(end)-ElapsedTime(nc14);

        else 
            
            BurstProperties(n).APBin = SpotDiff(n).APBin;
            
        end
            end
        end
        
        for bb=1:length(BurstProperties)
            if isempty(BurstProperties(bb).APBin)
                BurstProperties(bb).APBin=nan;
                % BurstProperties(bb).Nucleus=nan;
                % MOD: added uniqueID
                BurstProperties(bb).uniqueID = nan;
            end

        end
        clear ShortTrace
        %save('BurstPropertiesSlope','BurstProperties');
    end

end