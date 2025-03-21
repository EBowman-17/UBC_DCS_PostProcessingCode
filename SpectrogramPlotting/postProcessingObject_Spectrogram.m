classdef postProcessingObject_Spectrogram < handle & postProcessingObject
    %% This class is a sub class of the post processing object, 
    %  its purpose is for creating spectrogram plot data without 
    %  modifying its super class.
    
    properties

        % Number of data sets included in avgData
        numAvg;
        % Time axis lab time, on the scale of 10's of seconds
        labTime;
        % Ablation event repetition rate
        ablaFreq = 20;
        % Repetition rate data from both combs, and the third column is a
        % time vector in lab time
        repRates;
        
    end
    
    methods
    
        %% Constructor
        % The constructor is taken from postProcessingObject.m
        function obj = postProcessingObject_Spectrogram(sigLoc,refLoc,metaLoc)
            % Call the same constructor as postProcessingObject
            obj@postProcessingObject(sigLoc,refLoc,metaLoc);
        end
        
        
        %% Loading rep rate data
        % Rep rate data should be saved in the format: [comb2 data,comb1
        % data, time data]
        function obj = setRepRateData(obj,path)
           obj.repRates = importdata(path);
        end
        
        %% Averaging and calculating absorbance
        
        % Average the data in bins the size of numAvg, then consolidate 
        % then calculate absorbance. Also create a lab time vector with the
        % same amount of data points as the absorbance vector based off the
        % ablation laser frequency. If numAvg does not fit evenly into the
        % number of data traces, the last traces that don't form a complete
        % set will be discarded.
        function obj = freqAvg(obj,numAvg)
            obj.numAvg = numAvg;
            numSets = idivide(obj.num,int16(numAvg));
            tempSig = reshape(obj.sigSpec(:,1:numSets*numAvg),length(obj.sigSpec(:,1)),numAvg,numSets);
            tempRef = reshape(obj.refSpec(:,1:numSets*numAvg),length(obj.refSpec(:,1)),numAvg,numSets);
            tempSig = mean(tempSig,2);
            tempRef = mean(tempRef,2);
            obj.sigSpec = reshape(tempSig,length(obj.sigSpec(:,1)),numSets);
            obj.refSpec = reshape(tempRef,length(obj.refSpec(:,1)),numSets);
            obj.absorbance = log10(obj.refSpec./obj.sigSpec);
            obj.labTime = (0:numSets-1)*numAvg/obj.ablaFreq;
            
            % Update numSpec
            obj.numSpec = numSets;
        end
        
        %% Absorbance plotting
        
        % Show all the absorbance plots
        function plotAbsorbance(obj)
            % Create figure
            figure;
            index = 1; % Start with the first data set in the analysis range    
            % Micro second time scale
            fMHz = obj.freq;
            
            % Update figure y data
            plot(fMHz,obj.absorbance(:,index));
            title(strcat('Absorbance data from ',string(obj.numAvg),' averages. (',string(index),')'));
            xlabel('Frequency (THz)');
            ylabel('Absorbance');
            
            plotting = true;
            while plotting
                
                % Wait for user keypress
                if waitforbuttonpress
                    key = get(gcf, 'CurrentCharacter');

                    % Handle key presses
                    if key == 28 % Left arrow
                        index = max(index - 1, 1);
                        clf;
                        % Update figure y data
                        plot(fMHz,obj.absorbance(:,index));
                        title(strcat('Absorbance data from ',string(obj.numAvg),' averages. (',string(index),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Absorbance');
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.numSpec);
                        clf;
                        % Update figure y data
                        plot(fMHz,obj.absorbance(:,index));
                        title(strcat('Absorbance data from ',string(obj.numAvg),' averages. (',string(index),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Absorbance');
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end
        
        % Create a waterfall plot from the avg'd absorption data, range is
        % the frequency range to plot the absorption data.
        function plotWaterfall(obj,range)
            indices = find(obj.freq>min(range) & obj.freq<max(range));
  
            figure()
            waterfall(obj.labTime,obj.freq(indices(1):indices(end)),obj.absorbance(indices(1):indices(end),:))
            xlabel('Lab Time (s)')
            ylabel('Frequency (Hz)')
            zlabel('Absorbance')
            title('Absorbance Peak Vs. Lab Time')
        end
        
        function plotHeatmap(obj,range)
            indices = find(obj.freq>min(range) & obj.freq<max(range));
            centerFreq = obj.freq(indices(1)+ round((indices(end)-indices(1))/2))/1e9;
  
            figure()
            set(gcf,'Color','white')
            hold on
            title('Absorbance Peak Vs. Lab Time')
            xlabel('Lab Time (s)')
            
            yyaxis left
            imagesc(obj.labTime,obj.freq(indices(1):indices(end))/1e9-centerFreq,obj.absorbance(indices(1):indices(end),:))
            % Plot a line through the maximum of each peak
            [~,I] = max(obj.absorbance(indices(1):indices(end),:));
            I = I + indices(1) - 1;
            plot(obj.labTime,obj.freq(I)/1e9-centerFreq,'color','k')
            ylabel('Absorption Peak Deviation (GHz)')
            zlabel('Absorbance')
            ylim([obj.freq(indices(1)) obj.freq(indices(end))]/1e9-centerFreq);
            
            yyaxis right
            deviation = abs(obj.repRates(:,1)-obj.repRates(:,2))-mean(abs(obj.repRates(:,1)-obj.repRates(:,2)));
            plot(obj.repRates(:,3),deviation,'color','r')
            ylabel('Delta Frep Deviation (Hz)','color','r')
            ax = gca;
            ax.YColor = 'r';
            ylim([-0.25 0.2]);
            
            hold off
            
            %obj.freq(I(1:end-1))-obj.freq(I(2:end))
        end
        
        % Plot repetition rate data along with delta frep data
        function plotRepRates(obj)
           figure()
           hold on
           title("Comb Repetition Rates")
           xlabel("Time (s)")
           yyaxis left
           plot(obj.repRates(:,3),obj.repRates(:,2)-mean(obj.repRates(:,2)))
           plot(obj.repRates(:,3),obj.repRates(:,1)-mean(obj.repRates(:,1)))
           
           yyaxis right
           plot(obj.repRates(:,3),abs(obj.repRates(:,1)-obj.repRates(:,2))-mean(abs(obj.repRates(:,1)-obj.repRates(:,2))))
           
           legend("Frep Deviation C1","Frep Deviation C2","Delta Frep Deviation")
        end
    
        % Plot predicted peak displacement vs. actual peak displacement 
        function plotDisplacement(obj, range)
            
            indices = find(obj.freq>min(range) & obj.freq<max(range));
  
            figure()
            set(gcf,'Color','white')
            hold on
            title('Absorbance Peak Displacement Vs. Predicted Displacement')
            xlabel('Lab Time (s)')
            ylabel('Absorption Peak Deviation')
            % Plot a line through the maximum of each peak
            [~,I] = max(obj.absorbance(indices(1):indices(end),:));
            I = I + indices(1) - 1;
            plot(obj.labTime,obj.freq(I)/1e9-obj.freq(I(1))/1e9)
            
            % Generate the predicted displacement
            k = 5710001;
            n = 2910606-8;
            
            deltaFrep = abs(obj.repRates(:,1)-obj.repRates(:,2));
            RFDrift = (deltaFrep-deltaFrep(1))*(k-n*2)*obj.frep/obj.dFrep;
            plot(obj.repRates(:,3),RFDrift/1e9)
            
            legend("Data","Predicted")
            hold off
        end
    end
end