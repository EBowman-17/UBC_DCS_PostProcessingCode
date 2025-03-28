classdef postProcessingObject_SNR < handle & postProcessingObject
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
        
        % SNR data where each consecutive point contains numAvg more data
        % points than the next
        SNRDat;
        
        % Width of absorption line, this data goes with SNRDat
        peakWidth;
        
        % Chi^2 of peak fitting
        chi2;
        
        % Vector containing the number of averages per data point in SNRDat
        SNRDatNum;
        
    end
    
    methods
    
        %% Constructor
        % The constructor is taken from postProcessingObject.m
        function obj = postProcessingObject_SNR(sigLoc,refLoc,metaLoc)
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
        
        %% SNR finding
        
        % This function takes the peak location and floor location given in
        % THz and fits a voigt function to the peak to get its amplitude
        % then calculates the noise level for SNR calculations from the STD
        % of the floor range. If a 4th parameter is passed, it saves the
        % chi^2 of the fit, the peak width, and the SNR at all different
        % numbers of averages.
        function obj = SNRFinding(obj,peakRange,floorRange,peakProm,saveName)
            % Find the index corresponding to the peak in peakRange
            tempMeanSig = mean(obj.sigSpec,2);
            tempMeanRef = mean(obj.refSpec,2);
            tempAbsorb = log10(tempMeanRef./tempMeanSig);
            
            % Peak range
            peakStart = find(min(peakRange)<obj.freq & max(peakRange)>obj.freq,1);
            peakEnd = find(min(peakRange)<obj.freq & max(peakRange)>obj.freq,1,'last');
            
            % Find the psd of the looking range
            psd = 10*log10(mean(tempMeanRef(peakStart:peakEnd))*1000);
            
            % Find maximum peak height, also make sure a peak can be fitted
            % for
            [maxPeak,~,chiMaster] = voigtFit(tempAbsorb(peakStart:peakEnd)',obj.freq(peakStart:peakEnd));
            
            if chiMaster < 0.1
                obj.SNRDat = zeros(obj.numSpec,1);
                obj.SNRDatNum = zeros(obj.numSpec,1);
                obj.peakWidth = zeros(obj.numSpec,1);
                obj.chi2 = zeros(obj.numSpec,1);
                
                % Floor range
                rangeStart = find(min(floorRange)<obj.freq & max(floorRange)>obj.freq,1);
                rangeEnd = find(min(floorRange)<obj.freq & max(floorRange)>obj.freq,1,'last');

                for i = 1:obj.numSpec

                   j = double(i);
                   tempMeanSig = (j-1)*obj.numAvg*tempMeanSig/(j*obj.numAvg)+obj.numAvg*obj.sigSpec(:,j)/(j*obj.numAvg);
                   tempMeanRef = (j-1)*obj.numAvg*tempMeanRef/(j*obj.numAvg)+obj.numAvg*obj.refSpec(:,j)/(j*obj.numAvg);
                   tempAbsorb = log10(tempMeanRef./tempMeanSig);
                   
                   [maxi,fwhm,chi] = voigtFit(tempAbsorb(peakStart:peakEnd)',obj.freq(peakStart:peakEnd));
                   
                   obj.SNRDat(i) = maxi/maxPeak/std(tempAbsorb(rangeStart:rangeEnd));
                   obj.peakWidth(i) = fwhm;
                   obj.chi2(i) = chi;
                   obj.SNRDatNum(i) = i*obj.numAvg;
                end
            else
                disp("Peak Fitting Routine did not find a good fit, try narrowing your peak window")
            end
            
            % If a save location is given, save the data
            if nargin == 5
                SNRDat = obj.SNRDat;
                SNRDatNum = obj.SNRDatNum;
                peakWidth = obj.peakWidth;
                chi2 = obj.chi2;
                save(saveName,"SNRDat","SNRDatNum","peakRange","floorRange","psd","peakWidth","chi2");
            end
        end
        
        %% Absorbance plotting
        
        % Show all the absorbance plots
        function plotAbsorbance(obj)
            % Create figure
            figure;
            index = 1; % Start with the first data set in the analysis range    
            % Micro second time scale
            fMHz = obj.freq/1e12;
            
            % Update figure y data
            plot(fMHz,obj.absorbance(:,index));
            title(strcat("Absorbance data from ",string(obj.numAvg)," averages. (",string(index),")"));
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
                        title(strcat("Absorbance data from ",string(obj.numAvg)," averages. (",string(index),")"));
                        xlabel('Frequency (THz)');
                        ylabel('Absorbance');
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.numSpec);
                        clf;
                        % Update figure y data
                        plot(fMHz,obj.absorbance(:,index));
                        title(strcat("Absorbance data from ",string(obj.numAvg)," averages. (",string(index),")"));
                        xlabel('Frequency (THz)');
                        ylabel('Absorbance');
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end
        
        
    end
end