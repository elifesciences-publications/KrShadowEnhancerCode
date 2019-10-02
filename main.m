% Runs simulated annealing under different assumptions. 
% --Inputs: 
% strat: Either sifd (fitting to size, integrated fluorescence,
% frequency, and duration) or sifdc (fitting to size, 
% integrated fluorescence, frequency, duration, and correlation)
% model: proximal,distal, shadow, proxfixed, distfixed,
% prox2x, dist2x, proxnequal1, distnequal1, prox2xboot,
% dist2xboot, shadowboot, disthot, distcold, proxhot, proxcold
% proxhotdiff, proxcolddiff, disthotdiff, distcolddiff
% init: initial temperature for annealing
% random: 1 for random initial conditions, 0 for 
% analysis-derived initial conditions
% degIn1 - ... - sdIn1: parameters/initial conditions of enhancer at 
% proximal location (only for double enhancer fittings).
% degIn2 - ... - sdIn2: parameters/initial conditions of enhancer at 
% proximal location2 (only for double enhancer fittings).
% NOTE - input all parameters as strings 
% --Outputs:
% storage: the fitted parameter set
% out: the error associated with the fitted parameter set
% energies: the SSE over time
% iters: the total number of iterations 
% times: the amount of time spent per iteration
% NOTE - all the above quantities will be written to 
% text files in the same directory.
% Last updated: 9/01/19
function [storage, out, energies, iters, times] = ...
 main(model,strat,init,random, ...
    degIn1,koffIn1,konIn1,nIn1,scalingIn1,sdIn1, ...
    degIn2,koffIn2,konIn2,nIn2,scalingIn2,sdIn2)

    tic

    if nargin == 4
        degIn1 = '0';
        koffIn1 = '0';
        konIn1 = '0';
        nIn1 = '0';
        scalingIn1 = '0';
        sdIn1 = '0';
        degIn2 = '0';
        koffIn2 = '0';
        konIn2 = '0';
        nIn2 = '0';
        scalingIn2 = '0';
        sdIn2 = '0';

    % 2x proximal enhancer runs
    elseif nargin == 10

        degIn2 = '0';
        koffIn2 = '0';
        konIn2 = '0';
        nIn2 = '0';
        scalingIn2 = '0';
        sdIn2 = '0';

    end

    degIn1 = str2double(degIn1);
    koffIn1 = str2double(koffIn1);
    konIn1 = str2double(konIn1);
    nIn1 = str2double(nIn1);
    scalingIn1 = str2double(scalingIn1);
    sdIn1 = str2double(sdIn1);

    degIn2 = str2double(degIn2);
    koffIn2 = str2double(koffIn2);
    konIn2 = str2double(konIn2);
    nIn2 = str2double(nIn2);
    scalingIn2 = str2double(scalingIn2);
    sdIn2 = str2double(sdIn2); 

    % so that randomness is generated from the current time
    rng('shuffle')
    
    % number of nuclei per bin
    extra.numSpots = 80;
    % number of tries per temperature
    tries = 30;
    % time limit of 24 hours = 86400
    timeLim = 86400;
    % maximum number off tries in total in the anneal
    maxRejectAnneal = Inf;
    % cooling and stopping temperature for the annealing
    cool = 0.8;
    % fraction of the original temperature at which annealing 
    stoppingFraction  = 10;  

    % defaults
    fixedTF = false;
    withoutBursts = false;
    duplicatedAnneal = false;
    bootstrapping = false;
    temperature = false;
    hot = false;
    shadow = false; 
    tempDiff = false;
    % determines the rate of production to be used: 
    % endogenous distal or distal                                
    if (strcmp(strat,'sifd') && strcmp(model,'distal')) 
        distAtProx = false;           
    else 
        distAtProx = true;
    end

    switch model
        % correlation experiments (run with sifdc)
        case 'proxfixed'
            fixedTF = true;
            model = 'proximal';
        case 'distfixed'
            fixedTF = true;
            model = 'distal';
        case 'proxnequal1'
            withoutBursts = true;
            model = 'proximal';
        case 'distnequal1'
            withoutBursts = true;
            model = 'distal';
        % double experiments (run with sifd)
        case 'prox2x'
             duplicatedAnneal = true;
             model = 'proximal';
        case 'dist2x'
             duplicatedAnneal = true;
             model = 'distal';
        case 'prox2xboot'
            bootstrapping = true;
            duplicatedAnneal = true;
            model = 'proximal';
        case 'dist2xboot'
            bootstrapping = true;
            duplicatedAnneal = true;
            model = 'distal';
        case 'shadowboot'
            bootstrapping = true;
            shadow = true;
            model = 'shadow';
        case 'shadow'
            shadow = true;
        % temperature experiments (run with sifd)
        case 'proxhot'
            hot = true;
            model = 'proximal';
            temperature = true;
        case 'disthot'
            hot = true;
            model = 'distal';
            temperature = true;
            distAtProx = false;
        case 'proxcold'
            hot = false;
            model = 'proximal';
            temperature = true;
        case 'distcold'
            hot = false;
            model = 'distal';
            temperature = true;
            distAtProx = false;
        case 'proxhotdiff'
            hot = true;
            model = 'proximal';
            temperature = true;
            tempDiff = true;
        case 'disthotdiff'
            hot = true;
            model = 'distal';
            temperature = true;
            tempDiff = true;
            distAtProx = false;
        case 'proxcolddiff'
            hot = false;
            model = 'proximal';
            temperature = true;
            tempDiff = true;
        case 'distcolddiff'
            hot = false;
            model = 'distal';
            temperature = true;
            tempDiff = true;
            distAtProx = false;
        % annealing sifd to distal enhancer data
        % (distal at proximal, default with sifd
        % is endogenous distal)
        case 'distspecial'
            model = 'distal';
            distAtProx = true;      
    end

    % storing args in their original string form
    randomOrig = random;
    initOrig = init; 
    % arguments to compiled matlab code are passed as
    % strings
    random = str2double(random);
    init = str2double(init);
    
    % range for simulation
    upperEgg = 80; 
    lowerEgg = 20; 

    % we start at 2 because of pointers '.' and '..', 
    % loading correlation data
    if ~strcmp(strat,'sifd')
        files = dir(strcat(pwd,'/allData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allData/',files(myFiles).name));
        end

        DistalCorr = myStruct(1).files.DistalCorr;
        ProximalCorr = myStruct(14).files.ProximalCorr;
        corrBool = true;
    else 
        corrBool = false;
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
            DistInterBurst = myStruct(4).files.Dist32CInterBurst;
            DistIF = myStruct(5).files.Dist32CTotalProd;


            ProxBurstSize = myStruct(6).files.Prox32CBurstSize;
            ProximalDuration = myStruct(7).files.Prox32CDuration;
            ProxFreq = myStruct(8).files.Prox32CFreq;
            ProxInterBurstTime = myStruct(9).files.Prox32CInterBurst;
            ProxIF = myStruct(10).files.Prox32CTotalProd;


        else

            files = dir(strcat(pwd,'/allColdData'));
            for myFiles = 3:length(files)
                myStruct(myFiles - 2).files = load(strcat(pwd,'/allColdData/',files(myFiles).name));
            end

            DistBurstSize = myStruct(1).files.Dist17CBurstSize;
            DistDuration = myStruct(2).files.Dist17CDuration;
            DistFreq = myStruct(3).files.Dist17CFreq;
            DistInterBurst = myStruct(4).files.Dist17CInterBurst;
            DistIF = myStruct(5).files.Dist17CTotalProd;

            ProxBurstSize = myStruct(6).files.Prox17CBurstSize;
            ProximalDuration = myStruct(7).files.Prox17CDuration;
            ProxFreq = myStruct(8).files.Prox17CFreq;
            ProxInterBurstTime = myStruct(9).files.Prox17CInterBurst;
            ProxIF = myStruct(10).files.Prox17CTotalProd;

        end


    elseif duplicatedAnneal

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

        ProxBurstSize = myStruct(5).files.DoubProxBurstSize;
        ProximalDuration = myStruct(6).files.DoubProxDuration;
        ProxFreq = myStruct(7).files.DoubProxFreq;
        ProxIF = myStruct(8).files.DoubProxTotalProd;

    elseif shadow

        files = dir(strcat(pwd,'/allShadowData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allShadowData/',files(myFiles).name));
        end

        ShadowBurstSize = myStruct(1).files.BothBurstSize;
        ShadowDuration = myStruct(2).files.BothDuration;
        ShadowFreq = myStruct(3).files.BothFreq;
        ShadowIF = myStruct(4).files.BothTotalProd;

    % loading single enhancer data
    else     

        files = dir(strcat(pwd,'/allData'));
        for myFiles = 3:length(files)
            myStruct(myFiles - 2).files = load(strcat(pwd,'/allData/',files(myFiles).name));
        end

        ProxBurstSize = myStruct(9).files.ProxBurstSize;
        ProximalDuration = myStruct(10).files.ProximalDuration;
        ProxFreq = myStruct(11).files.ProxFreq;
        ProxInterBurstTime = myStruct(12).files.ProxInterBurstTime;
        ProxIF = myStruct(13).files.ProxTotalProd;
        % fixing the duration vectors that are not the right dimensions
        ProximalDuration = ProximalDuration';

        % loading data of distal at proximal or of endogenous distal
        if distAtProx

            files = dir(strcat(pwd,'/allDistalData'));
            for myFiles = 3:length(files)
                myStruct(myFiles - 2).files = load(strcat(pwd,'/allDistalData/',files(myFiles).name));
            end

            DistBurstSize = myStruct(1).files.DistBurstSize;
            DistDuration = myStruct(2).files.DistalDuration;
            DistFreq = myStruct(3).files.DistFreq;
            DistInterBurst = myStruct(4).files.DistInterBurstTime;
            DistIF = myStruct(5).files.DistTotalProd; 

            DistDuration = DistDuration';

        else

            DistBurstSize = myStruct(3).files.EndogDistBurstSize;
            DistDuration = myStruct(4).files.EndogDistDuration;
            DistFreq = myStruct(5).files.EndogDistFreq;
            DistInterBurst = myStruct(6).files.EndogDistInterBurst;
            DistIF = myStruct(7).files.EndogDistTotalProd;

        end

    end


    % getting the means of burst properties for initial conditions
    if ~shadow && ~duplicatedAnneal
        proxDurIdeal = nanmean(ProximalDuration);
        proxBetweenBurstsIdeal = nanmean(ProxInterBurstTime);
        distDurIdeal = nanmean(DistDuration);
        distBetweenBurstsIdeal = nanmean(DistInterBurst);
    end

    % we are dividing raw fluorescence data by F_rnap or F_1 since our model
    % operates in transcript count.
    F1 = 1338;
    % analysis for ICs
    switch model
        case 'proximal'  
            prox = true;
            if ~duplicatedAnneal
                burstLength = proxDurIdeal;
                betweenBursts = proxBetweenBurstsIdeal;
            end
            expFreq = ProxFreq;
            expDur = ProximalDuration;
            if corrBool
                 expCorr = ProximalCorr;
            end
            expSize = ProxBurstSize/F1;
            expIF = ProxIF/F1;
        case 'distal'  
            prox = false; 
            if ~duplicatedAnneal
                burstLength = distDurIdeal;
                betweenBursts = distBetweenBurstsIdeal;
            end
            expFreq = DistFreq;
            expDur = DistDuration; 
            if corrBool
               expCorr = DistalCorr;
            end
            expSize = DistBurstSize/F1;
            expIF = DistIF/F1;
        case 'shadow'
            prox = false;
            expFreq =  ShadowFreq;
            expDur = ShadowDuration;
            expSize = ShadowBurstSize/F1;
            expIF = ShadowIF/F1;
        otherwise  
            disp('You gave a nonexistent option.')  
            return
    end


    % Note that each of the data vectors have 41 entries in
    % them therefore each bin corresponds to around 2.43% in egg length. 
    % This implies that the first and last 20% of the egg length
    % are described by the first and last 8 bins of those 41 bins.
    % bins. Therefore, we only focus on bins 8 to 33.
    lowestBin = 8;
    highestBin = 33; 

    yFreq = expFreq(lowestBin:highestBin);
    yDur = expDur(lowestBin:highestBin);
    ySize = expSize(lowestBin:highestBin);
    yIF = expIF(lowestBin:highestBin);

    if corrBool
        yCorr = expCorr(lowestBin:highestBin);
    end
    % holds the maximum values for each burst property for 
    % subsequent normalizing
    yMaxGlobal =  [max(ySize) max(yIF) max(yFreq) max(yDur)];
    % indexes of different params in the vector
    nIndex = 1;
    scalingIndex = 2;
    sdIndex = 3;
    degIndex = 4; 
    konIndex = 5;
    koffIndex = 6;
    cIndex = 7;

    % initial conditions for the annealing. For all the code
    % [name of parameter] followed by a 1 corresponds to the 
    % enhancer in the proximal location while if it is followed
    % by a 2 it corresponds to the location of the distal. 
    if random 
        if duplicatedAnneal || shadow 
            k_onIC1 = 10.^((2--3).*rand + -3);
            k_offIC1 = 10.^((2--3).*rand + -3); 
            k_onIC2 = 10.^((2--3).*rand + -3);
            k_offIC2 = 10.^((2--3).*rand + -3); 

            if shadow
                nIC1 = floor(10^(2*rand + 1));
                scalingIC1 = 10.^((2--3).*rand + -3); 
                sdIC1 = randi(20);
                degIC1 = 10.^((2--3).*rand + -3);

                nIC2 = floor(10^(2*rand + 1));
                scalingIC2 = 10.^((2--3).*rand + -3); 
                sdIC2 = randi(20);
                degIC2 = 10.^((2--3).*rand + -3);
            else
                nIC = floor(10^(2*rand + 1));
                scalingIC = 10.^((2--3).*rand + -3); 
                sdIC = randi(20);
                degIC = 10.^((2--3).*rand + -3);
            end
        else 
            nIC = floor(10^(2*rand + 1));
            scalingIC = 10.^((2--3).*rand + -3); 
            sdIC = randi(20);
            degIC = 10.^((2--3).*rand + -3);
            k_onIC = 10.^((2--3).*rand + -3);
            k_offIC = 10.^((2--3).*rand + -3);
        end 
    % feed ICs from annealing in the single case
    elseif duplicatedAnneal 
        if prox
            % 1 = prox in sifdc
            k_onIC1 = konIn1;
            k_offIC1 = koffIn1; 
            k_onIC2 = k_onIC1;
            k_offIC2 = k_offIC1; 

        else
            % 1 = dist in sifdc
            % 2 = edist in sifd
            k_onIC1 = konIn1;
            k_offIC1 = koffIn1; 
            k_onIC2 = konIn2;
            k_offIC2 = koffIn2; 

        end

        % automatically these to be those of distal 
        % in sifdc
        degIC = degIn1;
        nIC = nIn1;
        scalingIC = scalingIn1;
        sdIC =  sdIn1;
    % using the lowest energy sets from single annealing
    elseif shadow 
            % 1 = prox in sifdc
            % 2 = dist in sifdc    
            k_onIC1 = konIn1;
            k_offIC1 = koffIn1; 
            k_onIC2 = konIn2;
            k_offIC2 = koffIn2; 

            degIC1 = degIn1;
            nIC1 = nIn1;
            scalingIC1 =  scalingIn1;
            sdIC1 = sdIn1;

            degIC2 = degIn2;
            nIC2 = nIn2;
            scalingIC2 = scalingIn2;
            sdIC2 =  sdIn2;

    % conditions for single case given by analysis    
    else

        nIC = 10;
        scalingIC = 40; 
        sdIC = 6;
        degIC =  1;  
        k_onIC = 1/betweenBursts;
        k_offIC = 1/burstLength; 

    end


    % annealing will stop after reaching a fraction of the original 
    % temperature
    stop = init/stoppingFraction;

    % setting the conditions for annealing
    def = struct(...    
            'CoolSched',@(T) (cool*T),...
            'Generator',@(x,t) myRand(x,t),...
            'InitTemp',init,...
            'MaxConsRej',maxRejectAnneal,...
            'MaxSuccess',Inf,...
            'MaxTries',tries,...
            'StopTemp',stop,...
            'StopVal',0,...
            'Verbosity',1);

    extra.alpha = 1.9540; 
    % length of simulation
    extra.endSpan = 50;
    % mean for beta_1,gamma_1
    enhancer.mean = 50;
    % number of bins to simulate
    extra.numBinsRWdata = 26;
    extra.eggLength = linspace(lowerEgg, upperEgg, extra.numBinsRWdata);
    extra.ElapsedTime = 0:0.5:extra.endSpan;

    % defining r
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
    elseif distAtProx
        % distAtProx is the default unless model 'distal' and 'sifd' are given in the original
        % arguments.
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
       enhancer.rFull1 =  [0   0   67.0257   60.0482   57.4959   51.5140   53.9837   53.1280 ...
           60.4217   90.3269  104.0844  114.0616  119.2144  101.4468   87.6735   79.4335 ...
           68.7479   63.0808   86.9607   0   0   0   0   0 ...
           0   0];
       enhancer.rFull2 = [0    0    0    0    0     ...
           40.4270   47.7268   56.6218   57.3002  72.3034  95.2149  123.3258  134.0092  ... 
           132.5551  121.6160 102.2820   87.8931   71.4522   65.2000   58.7315 ...
           56.1441   50.0337   0   0   0  0];
    % for a duplicated distal case we want the r of the distal at prox followed by the 
    % r of the endogenous distal
    elseif duplicatedAnneal && ~prox
        enhancer.rFull1 = [0   0   0   0   0   0   0   0   58.8609 ...
           74.4993   86.2019  111.3550  141.3504  131.6084  127.4489  115.2295   98.5590 ...
           95.3566   81.3611   71.7599   67.8933   59.8787   60.7270   59.1140   62.3468 ... 
           64.2655 ];
        enhancer.rFull2 = [0    0    0    0    0     ...
           40.4270   47.7268   56.6218   57.3002  72.3034  95.2149  123.3258  134.0092  ... 
           132.5551  121.6160 102.2820   87.8931   71.4522   65.2000   58.7315 ...
           56.1441   50.0337   0   0   0  0];
    end

    yGlobal = [ySize/max(ySize) yIF/max(yIF) yFreq/max(yFreq) yDur/max(yDur)];
    yMaxGlobal =  [max(ySize) max(yIF) max(yFreq) max(yDur)];

    sizeIndex = 1;
    IFIndex = 2;
    freqIndex = 3;
    durIndex = 4;

    % we don't anneal to correlation in double, only 
    % to burst properties
    if duplicatedAnneal || shadow
        xIC = [k_onIC1, k_offIC1, k_onIC2, k_offIC2];
        kon1Index = 1;
        koff1Index = 2;
        kon2Index = 3;
        koff2Index = 4;
    elseif withoutBursts || fixedTF
        xIC = [scalingIC, sdIC, degIC, k_onIC, k_offIC]; 
    % fixed TF or just regular single enhancer cases
    else
        xIC = [nIC, scalingIC, sdIC, degIC, k_onIC, k_offIC];
    end 

    switch strat
        % fitting only to size frequency and duration IF
        case 'sifd'  
            withC = false;
        % fitting to size frequency,duration, IF
        % and correlation
        case 'sifdc'  
            yGlobal = [yGlobal yCorr];    
            withC = false;
        otherwise  
            disp('You gave a nonexistent strategy.')    
            return
    end

    f = @paramOptimizer;
    g = @(p)f(p);

    if bootstrapping
    
        % so we know which parameter had the most effect
        % 1 = kon1, 2 = koff1, 3 = kon2, 4 = koff2
        paramString = [kon1Index, koff1Index, kon2Index, koff2Index];

        % changing either kon1, kon2, koff1, or koff2, one at a time
        for bootIndex = 1:length(xIC) 

            [storageTemp, outTemp, energies, iters, times] =  anneal(g, xIC(bootIndex), def);
            storage(bootIndex) = storageTemp;
            out(bootIndex) = outTemp;

        end 

        % sorting to see which one was able to achieve the best fitting
        [~, outOrder] = sort(out);
        paramString = paramString(outOrder);
        % note that storage and out keep their original order

    else

        [storage, out, energies, iters, times] =  anneal(g, xIC, def);

        % making storage the right size for file writing
        if withoutBursts || fixedTF
            storage = [1, storage];
        end  

    end

    % writing to file
    format = '%f, ';

    if ~(duplicatedAnneal || shadow)

        % predefining names for the files based on parameter, model type, strategy, initial temperature
        % and random ICs
        scalingFile = strcat('scaling',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        degFile = strcat('deg',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        sdFile = strcat('sd',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        konFile = strcat('kon',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        koffFile = strcat('koff',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        nFile = strcat('n',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        cFile = strcat('c',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        energiesFile = strcat('energies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        timesFile = strcat('times',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        itersFile = strcat('iters',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
        finalEnergyFile = strcat('finalEnergies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');


        if ~withoutBursts || ~fixedTF
            % writing to file
            fid = fopen(nFile,'a');
            fprintf(fid,format,storage(nIndex));
            fclose all;
        end

        fid = fopen(scalingFile,'a');
        fprintf(fid,format,storage(scalingIndex));
        fclose all;

        % writing to file
        fid = fopen(sdFile,'a');
        fprintf(fid,format,storage(sdIndex));
        fclose all;

        % writing to file
        fid = fopen(degFile,'a');
        fprintf(fid,format,storage(degIndex));
        fclose all;

        % writing to file
        fid = fopen(konFile,'a');
        fprintf(fid,format,storage(konIndex));
        fclose all;

        % writing to file
        fid = fopen(koffFile,'a');
        fprintf(fid,format,storage(koffIndex));
        fclose all;

        if withC 
            % writing to file
            fid = fopen(cFile,'a');
            fprintf(fid,format,storage(cIndex));
            fclose all;
        end

        fid = fopen(energiesFile,'a');
        fprintf(fid,format,energies);
        fprintf(fid,'\n','');
        fclose all;

        fid = fopen(finalEnergyFile,'a');
        fprintf(fid,format,energies(end));
        fprintf(fid,'\n','');
        fclose all;

        fid = fopen(itersFile,'a');
        fprintf(fid,format,iters);
        fprintf(fid,'\n','');
        fclose all;


        fid = fopen(timesFile,'a');
        fprintf(fid,format,times);
        fprintf(fid,'\n','');
        fclose all;

    else

        if bootstrapping

            outFile = strcat('outBoot',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            storageFile = strcat('storageBoot',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            paramsFile = strcat('paramsBoot',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

            fid = fopen(outFile,'a');
            fprintf(fid,format,out);
            fclose all;

            fid = fopen(storageFile,'a');
            fprintf(fid,format,storage);
            fclose all;

            fid = fopen(paramsFile,'a');
            fprintf(fid,format,paramString);
            fclose all;

        else

            konFile1 = strcat('kon1',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            koffFile1 = strcat('koff1',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            konFile2 = strcat('kon2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            koffFile2 = strcat('koff2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            energiesFile = strcat('energies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            timesFile = strcat('times',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            itersFile = strcat('iters',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
            finalEnergyFile = strcat('finalEnergies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

            fid = fopen(konFile1,'a');
            fprintf(fid,format,storage(kon1Index));
            fclose all;

            fid = fopen(koffFile1,'a');
            fprintf(fid,format,storage(koff1Index));
            fclose all;

            fid = fopen(konFile2,'a');
            fprintf(fid,format,storage(kon2Index));
            fclose all;

            fid = fopen(koffFile2,'a');
            fprintf(fid,format,storage(koff2Index));
            fclose all;

            fid = fopen(energiesFile,'a');
            fprintf(fid,format,energies);
            fprintf(fid,'\n','');
            fclose all;

            fid = fopen(finalEnergyFile,'a');
            fprintf(fid,format,energies(end));
            fprintf(fid,'\n','');
            fclose all;

            fid = fopen(itersFile,'a');
            fprintf(fid,format,iters);
            fprintf(fid,'\n','');
            fclose all;


            fid = fopen(timesFile,'a');
            fprintf(fid,format,times);
            fprintf(fid,'\n','');
            fclose all;

        end
        
    end

    toc

    % Prepares the error function for anneal.m
    % --Input: 
    % x: the vector of parameters to be fitted
    % --Output:
    % error: the error of the data given by simulating
    % with x
    % Alvaro Fletcher
    function error = paramOptimizer(x)

        switch strat
            case 'sifd'

                if duplicatedAnneal
                    enhancer.tf = nIC;
                    enhancer.scaling = scalingIC; 
                    enhancer.sd = sdIC;
                    enhancer.degTF = degIC;

                    if bootstrapping 

                        % default for bootstrapping     
                        enhancer.k_on1 = xIC(kon1Index); 
                        enhancer.k_off1 =  xIC(koff1Index);   
                        enhancer.k_on2 = xIC(kon2Index); 
                        enhancer.k_off2 =  xIC(koff2Index); 

                        % note that when bootstrapping, x is a scalar with 
                        % the parameter to be fitted instead the usual 
                        % vector of parameters
                        switch bootIndex

                            case kon1Index
                                enhancer.k_on1 = x;
                            case koff1Index
                                enhancer.k_off1 =  x;   
                            case kon2Index
                                enhancer.k_on2 = x; 
                            case koff2Index
                                enhancer.k_off2 =  x; 
                            otherwise
                                disp('Something is not working with bootstrapping')

                        end

                    % if not bootstrapping then do regular case in double annealing
                    else
                        
                        enhancer.k_on1 = x(kon1Index); 
                        enhancer.k_off1 =  x(koff1Index);   
                        enhancer.k_on2 = x(kon2Index); 
                        enhancer.k_off2 =  x(koff2Index); 

                    end

                elseif shadow
                    enhancer.tf1 = nIC1; 
                    enhancer.scaling1 = scalingIC1; 
                    enhancer.sd1 = sdIC1;
                    enhancer.degTF1 = degIC1;    

                    enhancer.tf2 = nIC2;
                    enhancer.scaling2 = scalingIC2; 
                    enhancer.sd2 = sdIC2;
                    enhancer.degTF2 = degIC2; 

                    if bootstrapping

                        % default for bootstrapping     
                        enhancer.k_on1 = xIC(kon1Index); 
                        enhancer.k_off1 =  xIC(koff1Index);   
                        enhancer.k_on2 = xIC(kon2Index); 
                        enhancer.k_off2 =  xIC(koff2Index); 

                        switch bootIndex

                            case kon1Index
                                enhancer.k_on1 = x; 
                            case koff1Index
                                enhancer.k_off1 =  x;   
                            case kon2Index
                                enhancer.k_on2 = x; 
                            case koff2Index
                                enhancer.k_off2 =  x; 
                            otherwise
                                disp('Something is not working with bootstrapping')

                        end

                    else

                        enhancer.k_on1 = x(kon1Index); 
                        enhancer.k_off1 =  x(koff1Index);   
                        enhancer.k_on2 = x(kon2Index); 
                        enhancer.k_off2 =  x(koff2Index); 

                    end

                % single enhancer cases
                else
                    enhancer.tf = x(1);
                    enhancer.scaling = x(2); 
                    enhancer.sd = x(3);
                    enhancer.degTF = x(4);      
                    enhancer.k_on = x(5); 
                    enhancer.k_off =  x(6); 
                    
                end     

                simData = mainfunction(enhancer, extra);
                sizes =  simData.sizes/yMaxGlobal(sizeIndex);
                IF = simData.IF/yMaxGlobal(IFIndex);
                freqs = simData.freqs/yMaxGlobal(freqIndex);
                durs = simData.durs/yMaxGlobal(durIndex);
                simVec = [sizes, IF, freqs, durs];

            case 'sifdc' 

                % without is always assumed n = 1 but in the case of fixed the 
                % gaussian is used as the fixed number of TFs which means we don't 
                % need the value of n.
                if withoutBursts || fixedTF
                    enhancer.scaling = x(1); 
                    enhancer.sd = x(2);
                    enhancer.degTF = x(3);      
                    enhancer.k_on = x(4); 
                    enhancer.k_off =  x(5);                  
                else
                    % note that for fixed TF concentrations
                    % we have that enhancer.tf represents the
                    enhancer.tf = x(1);
                    enhancer.scaling = x(2); 
                    enhancer.sd = x(3);
                    enhancer.degTF = x(4);      
                    enhancer.k_on = x(5); 
                    enhancer.k_off =  x(6); 
                end

                simData = mainfunction(enhancer, extra);
                sizes =  simData.sizes/yMaxGlobal(sizeIndex);
                IF = simData.IF/yMaxGlobal(IFIndex);
                freqs = simData.freqs/yMaxGlobal(freqIndex);
                durs = simData.durs/yMaxGlobal(durIndex);
                corrs = simData.corr;   
                simVec = [sizes, IF, freqs, durs, corrs];
 
        end 

        error =  SSEcalculator(simVec, yGlobal);  

    end

    % ANNEAL  Minimizes a function with the method of simulated annealing
    % (Kirkpatrick et al., 1983)
    %   joachim.vandekerckhove@psy.kuleuven.be
    %   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $
    function [minAnneal,fval,energiesAnneal,total,timesAnneal] = anneal(loss, parent, def)

        options=def;

        % main settings
        newsol = options.Generator;      % neighborhood space function
        Tinit = options.InitTemp;        % initial temp
        minT = options.StopTemp;         % stopping temp
        coolAnneal = options.CoolSched;        % annealing schedule
        minF = options.StopVal;
        max_consec_rejections = options.MaxConsRej;
        max_try = options.MaxTries;
        max_success = options.MaxSuccess;
        boltzmannC = 1;                           % boltzmann constant
        
        % counters etc
        itry = 0;
        success = 0;
        finished = 0;
        consec = 0;
        Tanneal = Tinit;
        annealTic = tic;
        initenergy = loss(parent);
        oldenergy = initenergy;

        energiesAnneal = [];
        timesAnneal = [];
        
        total = 0;
        fprintf(1,'\n  Tanneal = %7.5f, loss = %10.5f\n',Tanneal,oldenergy); 

        while ~finished
            energiesAnneal = [energiesAnneal oldenergy];
            timesAnneal = [timesAnneal toc(annealTic)];
            annealTic = tic;
            itry = itry + 1; % just an iteration counter
            current = parent; 
            % stop the program if it takes more 
            % than 24 hours
            if toc > timeLim
                minAnneal = parent;
                fval = oldenergy;

                fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
                fprintf(1, '  Final temperature:       \t%g\n', Tanneal);
                fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
                fprintf(1, '  Number of function calls:\t%i\n', total);
                fprintf(1, '  Total final loss:        \t%g\n', fval);
               
                return
            end
            
            % % Stop / decrement Tanneal criteria
            if itry >= max_try || success >= max_success
                if Tanneal < minT || consec >= max_consec_rejections
                    finished = 1;
                    total = total + itry;
                    break;
                else
                    Tanneal = coolAnneal(Tanneal); % decrease Tanneal according to cooling schedule
                    fprintf(1,'  Tanneal = %7.5f, loss = %10.5f\n',Tanneal,oldenergy);
                    total = total + itry;
                    itry = 1;
                    success = 1;
                end
            end
            
            newparam = newsol(current, Tanneal);
            newenergy = loss(newparam);
            
            if (newenergy < minF)
                parent = newparam; 
                oldenergy = newenergy;
                break
            end
            
            if (oldenergy-newenergy > 0)
                parent = newparam;
                oldenergy = newenergy;
                success = success + 1;
                consec = 0;
            else
                if (rand < exp( (oldenergy-newenergy)/(boltzmannC*Tanneal) ))
                    parent = newparam;
                    oldenergy = newenergy; 
                    % it's counted as a success but also 
                    % a failure since consec is not updated
                    success = success + 1;
                else
                    consec = consec + 1;
                end
            end
        end

        minAnneal = parent;
        fval = oldenergy;

        fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
        fprintf(1, '  Final temperature:       \t%g\n', Tanneal);
        fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
        fprintf(1, '  Number of function calls:\t%i\n', total);
        fprintf(1, '  Total final loss:        \t%g\n', fval);
        

    end

    % Updates parameters for annealing routine
    % -- Input:
    % myParams: vector of parameters
    % currentTemp: the current temperature
    % -- Output: 
    % out: the updated parameter set
    % Alvaro Fletcher
    function out = myRand(myParams, currentTemp) 

        myBool = true;  
        % lower bound tolerance
        tol = 0;
        % lower bound for TFs 
        tolTf = 1;

        while myBool

            tempOut = myParams;
            % sample from normal dist centered at 1 with 
            % standard dev multiplied by current temperature and multiply to a randomly chosen parameter
            newMean = 1;
            newSigma = currentTemp;
            % choose a random parameter and add sampled value
            holderParams = (randperm(length(myParams)) == length(myParams)) * (newMean + randn * newSigma);
            % finding index of parameter to update and value by which it will be multiplied
            holderVal = nonzeros(holderParams);
            holderIndex = find(holderParams);
            tempOut(holderIndex) = myParams(holderIndex) * holderVal;

            % in these cases we have no n
            if duplicatedAnneal || shadow || withoutBursts || fixedTF

                if  all(tempOut > tol)

                    out = tempOut;
                    myBool = false;
                    
                end

            else

                % the 1st element is the number of TFs in a burst 
                % and we force this number to always be an integer 
                % and greater than 0.This is also the case when we are 
                % fitting to fixed concentrations of TF.
                tempOut(1) = floor(tempOut(1));

                if tempOut(1) > tolTf && all(tempOut(2:end) > tol)

                    out = tempOut;
                    myBool = false;
                    
                end

            end

        end
        
    end

    % Calculates the sum of squared errors between 
    % two vectors vec1 and vec2.
    % --Inputs:
    % vec1: the simulated data
    % vec2: the observations
    % --Outputs:
    % error: the sum of squared errors
    % Alvaro Fletcher
    function error = SSEcalculator(vec1, vec2)

        % stores total error
        error = 0;

        % checking to see vectors are the same length
        if length(vec1) ~= length(vec2)
            disp('Your simulated vector and your experimental vector are not the same length');
        end

        % the length of vec1 and vec2 should be the same
        for sseIndex = 1:length(vec2)
            
            if ~(isnan(vec2(sseIndex)) || isnan(vec1(sseIndex)))
                
                tempError = (vec2(sseIndex) - vec1(sseIndex))^2;
                error = error + tempError;
                
            else
                
                % for our purposes, nan means that there was no data recorded
                % therefore we just set the error to the square of whatever
                % quantity is not nan
                if isnan(vec2(sseIndex)) && ~isnan(vec1(sseIndex))
                    
                    tempError = vec1(sseIndex)^2;
                    error = error + tempError;

                elseif isnan(vec1(sseIndex)) && ~isnan(vec2(sseIndex))
                    
                    tempError = vec2(sseIndex)^2;
                    error = error + tempError;
        
                else
                    
                    tempError = 0;
                    error = error + tempError;
            
                end
                
            end
            
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
    % simData: struct containing burst properties and/or allele
    % correlation of simulated data
    % Alvaro Fletcher
    function simData = mainfunction(enhancer, extra)

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
        
        if shadow
            [prodTF1, prodTF2] = TFproductionShadow(enhancer, extra);
            prodTF = zeros(1,26);
        else
            prodTF1 = zeros(1,26);
            prodTF2 = zeros(1,26);
            prodTF = TFproduction(enhancer, extra);
        end

        % shadow case
        if shadow || (duplicatedAnneal && ~prox)
            rateRNA = zeros(1,26);
            rateRNA1 = enhancer.rFull1;
            rateRNA2 = enhancer.rFull2;
        else
            rateRNA1 = zeros(1,26);
            rateRNA2 = zeros(1,26);   
            rateRNA = enhancer.rFull;
        end

        numSpotsPar = extra.numSpots;

        % holds the mRNA values for SpotDiff
        tempSimR = zeros(holderEggLength, length(extra.ElapsedTime), holderNumSpots);
        tempCorr = zeros(1,holderEggLength);
        tempIF = zeros(1,holderEggLength);

        % going over different bins in the egg
        parfor i = 1:holderEggLength

            production = prodTF(i);
            myParR = rateRNA(i);
      
            % limits the number of iterations for the Gillespie algorithm.
            % However, we want all simulations to be limited by the time span of
            % nc14 and just choose a very large number for maxIterations
            maxIterations = 10000000;

            % this is only for single enhancer cases that fit to all burst properties
            if ~(duplicatedAnneal || shadow)

                % the locations of the mRNA and transcription factor for our network
                RNAloc = 3;

                parTempSimR = zeros(length(holderElapsedTime), holderNumSpots);
                parTempIF = zeros(1, holderNumSpots);

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

                end

                tempSimR(i,:,:) = parTempSimR; 
                tempIF(i) = nanmean(parTempIF);

            end


            % case for fitting to correlation or double enhancer systems since 
            % we are assuming that both enhancers can bind at once allowing us 
            % to get the correlation "for free"
            if corrBool || duplicatedAnneal || shadow

                % picking the production rate of TF for shadow; otherwise it just remains
                % the same as that given by the single case
                production1 = prodTF1(i);
                myParR1 = rateRNA1(i);
                production2 = prodTF2(i);
                myParR2 = rateRNA2(i);

                % note that specifying "prox" and "dist" does not necessarily mean 
                % the proximal enhancer or distal enhancer but rather whichever 
                % enhancers are found in the proximal and distal locations
                RNAlocProx = 5;
                RNAlocDist = 8;

                if duplicatedAnneal || shadow
                    parTempSimR = zeros(length(holderElapsedTime), holderNumSpots);
                    parTempIF = zeros(1, holderNumSpots);
                end

                parTempCorr = zeros(1, holderNumSpots);
                for spot = 1:numSpotsPar

                    if shadow
                        CRN = feval(fcnShadow, enhancer, extra, production1, ...
                            production2, myParR1, myParR2);  
                    % duplicated case or correlation
                    else
                        if fixedTF
                            CRN = feval(fcnFixedTFDoub, enhancer, extra, production, myParR);        
                        else
                            if prox
                                CRN = feval(fcnDupEnhancer, enhancer, extra, production, myParR, myParR);
                            % we have r values for distal and edistal
                            else
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

                    if duplicatedAnneal || shadow
                        proxR = expN(:,RNAlocProx);
                        distR = expN(:,RNAlocDist);
                        R = proxR + distR;

                        parTempSimR(:,spot) = R;
                        parTempIF(spot) = trapz(holderElapsedTime, R);

                    else 

                        % the RNA locations for the correlation models
                        % that include R1 and R2 (RNA produced by enhancer 1
                        % and enhancer 2 respectively) are selected per 
                        % choice of model at the beginning
                        proxR = expN(:,RNAlocProx);
                        distR = expN(:,RNAlocDist);

                        % storing the correlation between R1 and R2 
                        corrholder = corrcoef(proxR,distR);
                        parTempCorr(spot) = corrholder(1,2);

                    end
                    
                end

                % storing all mRNA data and IF averages into our main 
                % structure
                if duplicatedAnneal || shadow
                     tempSimR(i,:,:) = parTempSimR; 
                     tempIF(i) = nanmean(parTempIF);
                else
                     tempCorr(i) = nanmean(parTempCorr);
                end

            end    

        end

        simData.corr = tempCorr;
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

                simData.sizes(end + 1) = ... 
                      nanmean([BurstProperties(currentAPBinVals).BurstSize]); 
                simData.durs(end + 1) = ...
                      nanmean([BurstProperties(currentAPBinVals).Duration]);   
                simData.freqs(end + 1) = ...
                       nanmean([BurstProperties(currentAPBinVals).Frequency]); 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);

            else 

                simData.sizes(end + 1) = nan;       
                simData.durs(end + 1) = nan;
                simData.freqs(end + 1) = nan; 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);
           
            end 

        end

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

        % we need to convert the number from our TF 
        % production gaussian curve to an integer
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

        CRN.target=reaction(1:reactionsNum,speciesNum + 1:2 * speciesNum); % target complexes

        % default parameters - any vector with reactionsNum entries

        % the r's will be the same for correlation purposes and for 2x proximal.
        % however, they are different in 2x distal due to having the data to calculate
        % r of edistal vs. distal 
        if duplicatedAnneal
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

        % number of TFs per burst
        n1 = params.tf1;
        n2 = params.tf2;

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
        source=CRN.source ;  % source is for source
        target=CRN.target; % target is for target
        paramsGE=CRN.params ; %  all the parameters in the system
        initGE=CRN.init;
        tspan=CRN.tspan;

        %default maximum number of iterations
        if nargin==1    
            maxiterations=500*nGE*mGE; 
        end

        % make sure initGE is aGE column vector
        temp = size(initGE);

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
    % 
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