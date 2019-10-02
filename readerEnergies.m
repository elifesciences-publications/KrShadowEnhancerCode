% Sorts simulated annealing results and then
% plots information on the best parameter sets.
% -- Inputs: 
% myDir: the name of the directory to be read (already processed
% by cleaner.sh or cleanerS.sh)
% stratDesired: 1 for the best parameter set, 2 for the second
% best, and so on.
% og: the fitted kon,koff values for the single
% enhancer cases to be compared with double enhancer 
% values (optional)
% -- Outputs: 
% bestVals: the parameter set resulting from the choice
% of stratDesired
% bestParams: the names of parameters in bestVals.
% Last updated: 9/01/2019
% Alvaro Fletcher
function [bestVals,bestParams] = readerEnergies(myDir, stratDesired, og)

	% true for outputting plots 
    boolPlot = true;
    % true to see differences in quality among
    % many temperatures
    temperDist = true;

    % top n strategies to be kept for plots 
    topStrats = 10;

	% thresholds under which we will visualize the parameters sets or 
	% the top n strategies that will be used (uncomment below to use
	% this function)
	proxThreshold = 10; 
	distThreshold = 10;
	% for shadow
	shadowThreshold = 10;

	if nargin == 2
		ogKoff1 = 0;
		ogKoff2 = 0;
		ogKon1 = 0;
		ogKon2 = 0;
    else
        % original values for kon and koff of single enhancer enhancer systems
        % in this case 1 = proximal and 2 = distal. 
        ogKoff1 = og.ogKoff1;
        ogKoff2 = og.ogKoff2;
        ogKon1 = og.ogKon1;
        ogKon2 = og.ogKon2;
	end

	% to prevent the legend from updating everytime we 
	% add a new data point
	set(0,'DefaultLegendAutoUpdate','off')

	% true for annealing to double models
	doub = false;
	% true if you are analyzing bootstrapping results
	boot = false;
	% true if using random ICs
	random = false;
	% true if n=1 always (without bursts)
	without = false;
	% true if doing fixed TF 
	fixedTF = false;

	switch myDir
        case 'fixed'
            color1 = [238 123 23]./255; % orange prox
            color2 =  [1 64 172]./255; % blue distal 
            fixedTF = true;
        case 'without'
        	without = true;
        	color1 = [238 123 23]./255; % orange prox
            color2 =  [1 64 172]./255; % blue distal 
        case 'sifd'
        	color1 = [238 123 23]./255;
            color2 = [208 76 240] ./255; % violet
        case 'sifdc'
	        color1 = [238 123 23]./255;
	        color2 =  [1 64 172]./255; % blue distal
        case 'sifdr'
        	color1 = [238 123 23]./255;
            color2 = [208 76 240] ./255; % violet 
            random = true;
        case 'sifdcr'
        	color1 = [238 123 23]./255;
            color2 =  [1 64 172]./255; % blue distal 
            random = true;
        case 'boot'
        	boot = true;
        	doub = true;
        	color1 = [215, 183, 58]./255; % prox yellow
            color2 =  [73 184 253] ./255;  % light blue
            color3 = [52 119 71]./255;  % green
        case {'double','doublec','doubles'}
        	doub = true;
        	color1 = [215, 183, 58]./255; % prox yellow
            color2 =  [73 184 253] ./255;  % light blue
            color3 = [52 119 71]./255;  % green
         case {'hot','hotdiff'} 
            color1 =  [255 50 23]./255; % pink 
            color2 =   [116 64 172]./255; % purple
        case {'cold','colddiff'}
            color1 = [200 123 185]./255; % darker pink
            color2 = [0 0 255]./255; % blue
        % copy the proximal files from sifd to this dir
        % before running readerEnergies
        case 'sifds'
        	 color1 = [238 123 23]./255; % orange
        	 color2 =  [1 64 172]./255; % blue distal
        case 'sifdrs'
        	 color1 = [238 123 23]./255; % orange
        	 color2 =  [1 64 172]./255; % blue distal
        	 random = true;
        otherwise
        	color1 = [238 123 23]./255; % orange prox
            color2 =  [1 64 172]./255; % blue distal
    end

    % changing to directory with two subdirs: finalenergies and 
    % param values
	cd(myDir)

	if doub
		numParams = 4;
		energyThreshold = [proxThreshold, distThreshold, shadowThreshold];
	    models = {'proximal','distal','shadow'};
		energyFolders = {'proxEnergies','distEnergies', 'shadEnergies'};
	elseif without || fixedTF
		numParams = 5;
		energyThreshold = [proxThreshold, distThreshold];
		models = {'proximal','distal'};
		energyFolders = {'proxEnergies','distEnergies'}; 
	else
		numParams = 6;
		energyThreshold = [proxThreshold, distThreshold];
		models = {'proximal','distal'};
		energyFolders = {'proxEnergies','distEnergies'}; 
	end

	% keeps track of what element in energyFolders we are in 
	energyCounter = 1;
	% store all the parameters values for the proximal and distal enhancer
	% in a cell array
	data = cell(length(models),numParams);

    if ~boot

		% stores information of the lowest energy values
		bestParams = {};
		bestVals = [];
		bestCounter = 1;

		for enhancer = models

			enhancerMod = char(enhancer);

			% changing into directory with the
			% energies
			cd(energyFolders{energyCounter})
            
			% choosing the adequate color
			if energyCounter == 1
				tempColor = color1;
			elseif energyCounter == 2
				tempColor = color2;
			else 
				tempColor = color3;
			end

			% getting data for the current directory with energies
			currentDir = dir;

			% initializing
			stratHolder = [];
			energyHolder = [];
			temperHolder = [];

			% doing the energies
			% first elements are '.' and '..' used for navigation. 
			for currentFile = 3:length(currentDir)
				% getting the current file name
				fileName = currentDir(currentFile).name; 
				% getting the strategy name and the temperature along 
				% with whether the annealing was random or not
				holderString = extractBetween(fileName,enhancerMod,'.txt');
		 		holderStringChar = char(holderString);
				% reads all the parameters for this strategy
				tempHolder = readmatrix(fileName);
				% removes empty column at the end due to formatting
				tempHolder(:,end) = [];

				% going over all the energies in this file
			    for j = 1:length(tempHolder) 
					 tempHolder = tempHolder';
					 % storing the energies in one array
					 energyHolder = [energyHolder, tempHolder(j)];
					 % storing the corresponding strategy of this energy
					 % along with its position within the file so that 
					 % we can retrieve it later if necessary
					 stratHolder = [stratHolder,strcat(holderString,'ZZZ',num2str(j))]; 
					 % stores the corresponding temperature so that we can sort later
					 % (the last character is just whether the run was random or not)
					 temperHolder = [temperHolder, str2double(holderStringChar(end - 1))];
		        end   
		      
		    end



			% sorting in ascending order
			[energyHolder, energyOrder] = sort(energyHolder);
			stratHolder = stratHolder(energyOrder);
			temperHolder = temperHolder(energyOrder);


			if boolPlot && temperDist
				figure(energyCounter + 500)
				% plotting energy-temperature distributions
				scatter(1:length(temperHolder),temperHolder,[],tempColor)
	            ylabel('Temperature','fontweight','bold','fontsize',10, 'FontName', 'Arial')
	            xlabel('Iteration Number','fontweight','bold','fontsize',10, 'FontName', 'Arial')
	            title('Temperatures sorted by SSE')
	            ylim([0 6])
	            hold on 
	        end
		    
		    myDispStr = ['Strat desired in' ' ' enhancerMod  ' ' 'is' ' ' ...
		        num2str(energyHolder(stratDesired)) ' ' 'from' ' ' ...
		        char(stratHolder(stratDesired))];
		    disp(myDispStr);

            if boolPlot 
                % plotting all energies in a scatter plot
                figure(1)
                % random ICs will be printed in diamonds
                if random 
                     scatter(1:length(energyHolder),energyHolder,[],tempColor,'d')
                else
                    scatter(1:length(energyHolder),energyHolder,[],tempColor)
                end
                xlabel('Iteration number','fontweight','bold','fontsize',10, 'FontName', 'Arial')
                ylabel('Energy','fontweight','bold','fontsize',10, 'FontName', 'Arial')
                title(['Sorted SSE for' ' ' myDir])
                hold on 
            end

		    % plotting all energies in a histogram
            if boolPlot
            
                figure(energyCounter + 200)
                if random
                    myHist = histogram(energyHolder);
                    % random histograms will have different color
                    myHist.FaceColor = [128 128 128]/255; % grey
                else
                    myHist = histogram(energyHolder,10);
                    myHist.FaceColor = tempColor;
                end

                xlabel('Energies (SSE)','fontweight','bold','fontsize',10, 'FontName', 'Arial')
                ylabel('Counts','fontweight','bold','fontsize',10, 'FontName', 'Arial')
                % normalizing histogram and making bins same 
                % size for comparing random vs nrandom ICs
                % myHist.Normalization = 'probability';
                myHist.BinWidth = 0.25;
                hold on

				%title('Histogram of Energies')
				xlim([0 max(energyThreshold) + 2])

           
            end 


            % keeping only top ten parameter sets
			energyHolder = energyHolder(1:topStrats);
			stratHolder = stratHolder(1:topStrats);

		    % getting rid of energies that are below the threshold 

			% badIndices = find(energyHolder > energyThreshold(energyCounter));
			% energyHolder(badIndices) = []; 
			% stratHolder(badIndices) = [];

			% going back to dir with final energies and param subdir
			cd('..')

			% going to subdir with parameters
			cd(enhancerMod)

			% getting data for the current directory
			currentDir = dir; 

			% hold all the parameter names in the same index as their value
			% in allVals
			allParams = {};
			allVals = [];

			% counter for allParams so that we can add values
			% simultaneously with allVals
			tempCounter = 1;

		 	% iterating over all param files
			for currentFile = 3:length(currentDir)

			  % getting the current file name
			  fileName = currentDir(currentFile).name; 
			  % hold the strategy used to obtain these parameters
			  stratString = extractBetween(fileName, enhancerMod,'.txt');
			  % holds the name of the parameters in the current file
			  paramString = extractBefore(fileName, enhancerMod);

			  % this part is very inefficient for large data sets, but
			  % in our case we have only a few kb of text

			  % trying to find a match in the winner strats
			    for j = 1:length(stratHolder)
			  		% comparing strats to find a match
			  		if startsWith(stratHolder{j},char(stratString)) 

			  			% note that the information which was the parameter
			  			% corresponding the to best run is already included in 
			  			% stratHolder{1} after the 'ZZZ'
			  			temp = extractAfter(stratHolder{j},'ZZZ');
			  			% picking which the number of the run 
			  			% to extract from this file
			  			numRun = str2double(temp);

			  			% reading matrix in file and picking 
			  			% the right value
						valHolder = readmatrix(fileName);
						valHolder(:,end) = [];

			  			allVals = [allVals, valHolder(numRun)];
			  			allParams{tempCounter} = char(paramString);
			  			tempCounter = tempCounter + 1;	

			  		    % if we found the best strategy. note that 
			  		    % stratholder is sorted already.
			  		    if j == stratDesired 

				  			bestVals = [bestVals, valHolder(numRun)];
				  			bestParams{bestCounter} = char(paramString);
				  			bestCounter = bestCounter + 1;

			  		    end

			  		end

			    end

			end

			% getting indices. Consecutive values in allParams are the same since the files
			% are read alphabetically. ia is an array containing 
			% where the consecutive repeats of one values switch to 
			% those of another
			[uniqueParams, ia, ~] = unique(allParams);

			% storing all the values for each parameter 
			% in different columns of the cell and in the 
			% the current (each row is either proximal,
			% distal or shadow
			for i = 1:length(ia) - 1
				data{energyCounter, i} = allVals(ia(i):ia(i+1) - 1);
			end

			% doing the last value
			data{energyCounter, length(ia)} = allVals(ia(end):end);

			% going back to directory with final energies and parameters
			cd('..') 

			energyCounter = energyCounter + 1;

		end

		% going back to the top level directory where this text
		% file is located

		cd('..')

		% renaming some parameters to their real names
		uniqueParams = strrep(uniqueParams,'deg', [char(946),'_{-1}']);
        uniqueParams = strrep(uniqueParams,'sd', char(963));
        uniqueParams = strrep(uniqueParams,'scaling', 'c');
		
		if doub
			% the parameter names for the x axis of the boxplot
			xlab={char(uniqueParams(1)),char(uniqueParams(2)), ...
			char(uniqueParams(3)),char(uniqueParams(4))};
			col=[color3, 200/255;
			color2, 200/255;
			color1, 200/255];
		elseif without || fixedTF
			xlab={char(uniqueParams(1)),char(uniqueParams(2)),char(uniqueParams(3)),char(uniqueParams(4)), ...
			char(uniqueParams(5))};  
			col=[color2, 200/255;
			color1, 200/255];	
		else
			xlab={char(uniqueParams(1)),char(uniqueParams(2)),char(uniqueParams(3)),char(uniqueParams(4)), ...
			char(uniqueParams(5)),char(uniqueParams(6))};  
			col=[color2, 200/255;
			color1, 200/255];
		end

		% building matrices for correlation heatmap 
		% different matrices are necessary since our thresholds
		% are distinct for different systems
		for k = 1:numParams
			myPairMat1(:,k) = data{1,k};
		end

		for k = 1:numParams
			myPairMat2(:,k) = data{2,k};
		end  

		if doub
			for k = 1:numParams
				myPairMat3(:,k) = data{3,k};
			end  
		end

		% note that each column corresponds to a parameter while
		% each row corresponds to a particular parameter set 
		% resulting from annealing. The pairwise correlation
		% between all parameters is then calculated.

		allCorrs1 = corr(myPairMat1);
		allCorrs2 = corr(myPairMat2);
		if doub
			allCorrs3 = corr(myPairMat3);
		end

		% standardizing the parameter matrix
		zVals1 = zscore(myPairMat1);
		zVals2 = zscore(myPairMat2);
		if doub
			zVals3 = zscore(myPairMat3);
		end

		% doing pca on all results
		[coeff1,score1,~] = pca(zVals1);
		[coeff2,score2,~] = pca(zVals2);
		if doub
			[coeff3,score3,~] = pca(zVals3);
        end
        
        if boolPlot 
            % plotting correlation
            cpCounter = 400; 
            figure(cpCounter)
            imagesc(allCorrs1)
            colorbar 
            title([myDir ' ' 'pairwise correlation for proximal'])
            if doub		
                set(gca,'YTick',[1:4])
                set(gca,'XTick',[1:4])
            elseif without || fixedTF
            	set(gca,'YTick',[1:5])
                set(gca,'XTick',[1:5])
            end
            set(gca,'YTickLabel',xlab)
            set(gca,'XTickLabel',xlab)


            cpCounter = cpCounter + 1;
            figure(cpCounter)
            imagesc(allCorrs2)
            colorbar
            title([myDir ' ' 'pairwise correlation for distal'])
            if doub		
                set(gca,'YTick',[1:4])
                set(gca,'XTick',[1:4])
            elseif without || fixedTF
            	set(gca,'YTick',[1:5])
                set(gca,'XTick',[1:5])
            end
            set(gca,'YTickLabel',xlab)
            set(gca,'XTickLabel',xlab)

            if doub
                cpCounter = cpCounter + 1;
                figure(cpCounter)
                imagesc(allCorrs3)
                colorbar
                title([myDir ' ' 'pairwise correlation for shadow'])
                set(gca,'YTickLabel',xlab)
                set(gca,'YTick',[1:4])
                set(gca,'XTickLabel',xlab)
                set(gca,'XTick',[1:4])
            end

            % plotting pca
            cpCounter = cpCounter + 1;
            figure(cpCounter)
            biplot(coeff1(:,1:2),'scores',score1(:,1:2), ...
                'Color',color1,'Marker','o','varlabels',xlab);
            title([myDir ' ' 'PCA for proximal'])

            cpCounter = cpCounter + 1;
            figure(cpCounter)
            biplot(coeff2(:,1:2),'scores',score2(:,1:2), ...
                'Color',color2,'Marker','o','varlabels',xlab);
            title([myDir ' ' 'PCA for distal'])

            if doub
                cpCounter = cpCounter + 1;
                figure(cpCounter)
                biplot(coeff3(:,1:2),'scores',score3(:,1:2), ...
                    'Color',color3,'Marker','o','varlabels',xlab);
                title([myDir ' ' 'PCA for shadow'])
            end

            % plotting table with mean and cv's of each parameter
            cpCounter = cpCounter + 1;
            figure(cpCounter)
            tableMat1 = [median(myPairMat1); std(myPairMat1)./mean(myPairMat1)]
            uit1 = uitable('Data',tableMat1');
            uit1.RowName = xlab;
            uit1.ColumnName = {'Median','CV'};

            cpCounter = cpCounter + 1;
            figure(cpCounter)
            tableMat2 = [median(myPairMat2); std(myPairMat2)./mean(myPairMat2)]
            uit2 = uitable('Data',tableMat2');
            uit2.RowName = xlab;
            uit2.ColumnName = {'Median','CV'};

            if doub
                cpCounter = cpCounter + 1;
                figure(cpCounter)
                tableMat3 = [median(myPairMat3); std(myPairMat3)./mean(myPairMat3)]
                uit3 = uitable('Data',tableMat3');
                uit3.RowName = xlab;
                uit3.ColumnName = {'Median','CV'};
            end

            % adding legends to the scatter plot of energies and to histograms
            figure(1)

            title('Sorted Energies')
            ylim([0 max(energyThreshold) + 2])

            figure(3)

            if doub	
                boxplot2(data',xlab,{'proximal','distal','shadow'},col');
            else 		
                boxplot2(data',xlab,{'proximal','distal'},col') 
            end
            hold on

            if doub && length(models) == 3

                 plot(1.25,ogKoff1,'g*')
                 plot(1.5,ogKoff2,'g*')
                 plot(1.75,ogKoff1,'g*')

                 plot(2.25,ogKoff1,'g*')
                 plot(2.5,ogKoff2,'g*')
                 plot(2.75,ogKoff2,'g*')

                 plot(3.25,ogKon1,'g*')
                 plot(3.5,ogKon2,'g*')
                 plot(3.75,ogKon1,'g*')

                 plot(4.25,ogKon1,'g*')
                 plot(4.5,ogKon2,'g*')
                 plot(4.75,ogKon2,'g*')

            elseif doub && length(models) == 2

                plot(1.25,ogKoff1,'g*') 
                plot(1.5,ogKoff2,'g*')

                plot(2,ogKoff1,'g*')
                plot(2.25,ogKoff2,'g*')

                plot(2.75,ogKon1,'g*')
                plot(3,ogKon2,'g*')

                plot(3.5,ogKon1,'g*')
                plot(3.75,ogKon2,'g*')

            end
            
        end
        
	% doing the bootstrapping case
	else 

		modelCounter = 1;
		plotTracker = 300;
		folderPrefixes = {'prox','dist','shad'};
		paramNames = {'prox','dist','shad'};

		% doing proximal, distal, or shadow (one at a time)
		for enhancer = models

			enhancerMod = char(enhancer);
			cd(enhancerMod)


		    endingFiles = {'Params','Storage','Out'};
			endingCounter = 1;


			for endingIndex = endingFiles

				clear dataOther

			    endingIndex = char(endingIndex);

				% moving into the folder with param nums
				folderString = [folderPrefixes{modelCounter},endingIndex];
				cd(folderString)
				% getting data for the current directory
				currentDir = dir; 

				% keeps track of the current column
				colCounter = 1;

				% doing either the param files, out files, or storage files
				% first elements are '.' and '..' used for navigation. 
				for currentFile = 3:length(currentDir)
					% getting the current file name
					fileName = currentDir(currentFile).name; 
					% reads all the parameters for this strategy
					tempHolder = readmatrix(fileName);
					% removes empty column at the end due to formatting
					tempHolder(:,end) = [];


					% resets every time it reaches 4. each annealing run 
					% of this method corresponds to 4 consecutive values
					% in each folder
					quadCounter = 1;
					% assigning the strategy for the current energies
					for j = 1:length(tempHolder) 
						 dataOther(quadCounter,colCounter) = tempHolder(j);
						 quadCounter = quadCounter + 1;
						 if quadCounter == 5
						 	quadCounter = 1;
						 	colCounter = colCounter + 1;
						 end

			        end   

			    end

			    % choosing color and titles for plots
			    switch enhancerMod
			    	case 'proximal'
			    		tempColor = color1;
			    		tempString = '2x Proximal';
			    	case 'distal'
			    		tempColor = color2;
			    		tempString = '2x Distal';
			    	case 'shadow'
			    		tempColor = color3;
			    		tempString = 'Shadow';
			    end

			    % 1 = kon1, 2 = koff1, 3 = kon2, 4 = koff2
			    switch endingIndex
			    	% generating 1st,2nd...place histograms
			    	case 'Params'
			    		for k = 1:numParams
				    		figure(plotTracker)
				    		plotTracker = plotTracker + 1;


				    		dataC = categorical(dataOther(k,:),[1 2 3 4], ...
				    			{'kon1','koff1','kon2','koff2'});
					    	myHist = histogram(dataC);
					    	myHist.FaceColor = tempColor;

					    	title([tempString ' ' num2str(k) ' ' 'st/nd/rd/th place distribution'])
					    end
					% storing data for box plot
			    	case 'Storage'
						for k = 1:numParams
						    % filling up rows of data
						    data{modelCounter, k} = dataOther(k,:);
						end
					% generating histograms and scatter of energies
			    	case 'Out'  	
			    		for k = 1:numParams
			    		 	figure(plotTracker)
					    	scatter(1:length(dataOther),sort(dataOther(k,:)))
					    	hold on
					    end
					    legend('kon1','koff1','kon2','koff2')
					    title([tempString ' ' 'sorted SSE for each parameter'])
					    plotTracker = plotTracker + 1;

				end

				endingCounter = endingCounter + 1;

				% going back to directory with out,params and storage
		   		cd('..')

			 end
	
			% going back to directory with enhancer systems
		    cd('..')

		    % keeps track of which enhancer we are on
		    modelCounter = modelCounter + 1;

		end

		% going back to top level dir where boxplot function 
		% is located
		cd('..')

		 % generating box plot
		 % 1 = kon1, 2 = koff1, 3 = kon2, 4 = koff2
		figure(plotTracker)
		% the parameter names for the x axis of the boxplot
		xlab={'kon1','koff1','kon2','koff2'}; 	
		col=[color3, 200/255;
		color2, 200/255;
		color1, 200/255];

		boxplot2(data',xlab,{'proximal','distal','shadow'},col');

		hold on

		plot(1.25,ogKon1,'g*')
		plot(1.5,ogKon2,'g*')
		plot(1.75,ogKon1,'g*')

		plot(2.25,ogKoff1,'g*')
		plot(2.5,ogKoff2,'g*')
		plot(2.75,ogKoff1,'g*')

		plot(3.25,ogKon1,'g*')
		plot(3.5,ogKon2,'g*')
		plot(3.75,ogKon2,'g*')

		plot(4.25,ogKoff1,'g*')
		plot(4.5,ogKoff2,'g*')
		plot(4.75,ogKoff2,'g*')

end


