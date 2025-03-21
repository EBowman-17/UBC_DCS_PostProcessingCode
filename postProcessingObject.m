classdef postProcessingObject < handle
    % Post Processing object
    
    properties
        % Reference, Signal, and Metadata file locations
        sigLoc = ""
        refLoc = ""
        metaLoc = ""
        
        % Range of data points to analyze, default is the full range, this
        % can be changed with the mutator: setAnaRange
        anaRange;
        
        % Number of data traces in the current analysis window, for time
        % domain data
        num;
        
        % Number of data traces in the current analysis window, for
        % frequency domain data
        numSpec;
        
        % Data acquisition clock frequency, used for generating a frequency
        % axis default is 200e6 samples/s, modify it with the daqFreq
        % Mutator: setDaqFreq
        daqFreq = 200e6;

        % For peak finding in the time domain, this is the minimum value
        % that constitutes a peak, default is 5000, can be set with the
        % mutator: setPeakProm
        peakProm = 5000;
        
        % The signal, reference data, and the frequency and time axis
        sig;
        ref;
        % The signal amplitude spectrum and reference amplitude spectrum
        sigSpec;
        refSpec;
        
        % Absorbance data, should only be filled if averaging is performed
        % and absorbance is calculated
        absorbance;
        
        % Frep and Delta Frep
        frep;
        dFrep;

        % The frequency and time axis in us and Hz
        t;
        freq;
        
        % Rectangual Apodization window length us 
        % Mutator is : setApoLen
        apoLen;
        
        % Zero padding, this adds the specified amount of zeros to the data
        % in us to both sides of the interferogram.
        padNum;
        
        % Appodization window, this will be blackmanHarris nearly always
        window;
        
    end
    
    methods
%% Object constructor
        %Instantiate an object, load variables
        function obj = postProcessingObject(sigPath,refPath,metaPath)
            
            if nargin == 3
                % Signal, Reference, and MetaData paths and loading data
                obj.sigLoc = sigPath;
                obj.refLoc = refPath;
                obj.metaLoc = metaPath;

                load(obj.sigLoc);
                obj.sig = sig;
                load(obj.refLoc);
                obj.ref = ref;
                load(obj.metaLoc,'delta_frep');
                obj.dFrep = delta_frep;
                load(obj.metaLoc,'comb1_frep');
                obj.frep = comb1_frep;

                % Analysis range of data, default is all the data
                obj.anaRange = 1:length(obj.sig(1,:));

                % Number of data traces currently in the analysis object
                obj.num = length(obj.sig(1,:));

                % Create a time axis for reference in the data in us
                obj.t = 1e6*(0:length(obj.sig(:,1))-1)./obj.daqFreq;
            end

        end

%% Time Domain Plotting
        % Plot the time domain signal data in its current state default is it
        % plots every time domain trace.
        % Also labels peak positions, and ablation delays.
        function plotTimeDomSig(obj)

            % Create figure
            figure();
            index = 1; % Start with the first data set in the analysis range  
            % Micro second time scale
            tus = obj.t;
            
            % Plot current data
            plot(tus,obj.sig(:,index));
            % Find peaks and label them
            [pks,locs] = findpeaks(obj.sig(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
            text(locs+1,pks,string(locs))
            % If two peaks are detected, display the delta
            if length(locs)>=2
                text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
            end
            title(strcat('Time Domain, current signal data (',string(obj.anaRange(index)),')'));
            xlabel('Time (us)');
            ylabel('Amplitude (A.U.)');

            plotting = true;
            while plotting

                % Wait for user keypress
                if waitforbuttonpress
                    key = get(gcf, 'CurrentCharacter');

                    % Handle key presses
                    if key == 28 % Left arrow
                        index = max(index - 1, 1);
                        clf;
                        plot(tus,obj.sig(:,index));
                        title(strcat('Time Domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        % Find peaks and label them
                        [pks,locs] = findpeaks(obj.sig(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
                        text(locs+1,pks,string(locs))
                        % If two peaks are detected, display the delta
                        if length(locs)>=2
                            text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
                        end
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.num);
                        clf;
                        plot(tus,obj.sig(:,index));
                        title(strcat('Time Domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        % Find peaks and label them
                        [pks,locs] = findpeaks(obj.sig(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
                        text(locs+1,pks,string(locs))
                        % If two peaks are detected, display the delta
                        if length(locs)>=2
                            text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
                        end
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end

        % Plot the time domain reference data in its current state default is it
        % plots every time domain trace.
        % Also labels peak positions
        function plotTimeDomRef(obj)

                       % Create figure
            figure();
            index = 1; % Start with the first data set in the analysis range  
            % Micro second time scale
            tus = obj.t;
            
            % Plot current data
            plot(tus,obj.ref(:,index));
            % Find peaks and label them
            [pks,locs] = findpeaks(obj.ref(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
            text(locs+1,pks,string(locs))
            % If two peaks are detected, display the delta
            if length(locs)>=2
                text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
            end
            title(strcat('Time Domain, current reference data (',string(obj.anaRange(index)),')'));
            xlabel('Time (us)');
            ylabel('Amplitude (A.U.)');

            plotting = true;
            while plotting

                % Wait for user keypress
                if waitforbuttonpress
                    key = get(gcf, 'CurrentCharacter');

                    % Handle key presses
                    if key == 28 % Left arrow
                        index = max(index - 1, 1);
                        clf;
                        plot(tus,obj.ref(:,index));
                        title(strcat('Time Domain, current reference data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        % Find peaks and label them
                        [pks,locs] = findpeaks(obj.ref(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
                        text(locs+1,pks,string(locs))
                        % If two peaks are detected, display the delta
                        if length(locs)>=2
                            text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
                        end
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.num);
                        clf;
                        plot(tus,obj.ref(:,index));
                        title(strcat('Time Domain, current reference data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        % Find peaks and label them
                        [pks,locs] = findpeaks(obj.ref(:,index),tus,'MinPeakProminence',obj.peakProm,'MinPeakDistance',1);
                        text(locs+1,pks,string(locs))
                        % If two peaks are detected, display the delta
                        if length(locs)>=2
                            text(0.72,0.94,append("Delta = ",string(locs(2)-locs(1))),'units','normalized')
                        end
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end

        %% Frequency Domain Plotting
        function plotFreqDomSig(obj)
            % Create figure
            figure;
            index = 1; % Start with the first data set in the analysis range    
            % Micro second time scale
            fTHz = obj.freq*obj.dFrep/obj.frep/1e6;
            
            % Update figure y data
            plot(fTHz,obj.sigSpec(:,index));
            title(strcat('Frequency domain, current signal data (',string(obj.anaRange(index)),')'));
            xlabel('Frequency (MHz)');
            ylabel('Power Spectral Density (W/Hz)');
            
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
                        plot(fTHz,obj.sigSpec(:,index));
                        title(strcat('Frequency domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Power Spectral Density (W/Hz)');
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.numSpec);
                        clf;
                        % Update figure y data
                        plot(fTHz,obj.sigSpec(:,index));
                        title(strcat('Frequency domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Power Spectral Density (W/Hz)');
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end
        
        function plotFreqDomRef(obj)
            % Create figure
            figure;
            index = 1; % Start with the first data set in the analysis range    
            % Micro second time scale
            fTHz = obj.freq*obj.dFrep/obj.frep/1e6;
            
            % Update figure y data
            plot(fTHz,obj.sigSpec(:,index));
            title(strcat('Frequency domain, current reference data (',string(obj.anaRange(index)),')'));
            xlabel('Frequency (MHz)');
            ylabel('Power Spectral Density (W/Hz)');
            
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
                        plot(fTHz,obj.sigSpec(:,index));
                        title(strcat('Frequency domain, current reference data (',string(obj.anaRange(index)),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Power Spectral Density (W/Hz)');
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.numSpec);
                        clf;
                        % Update figure y data
                        plot(fTHz,obj.sigSpec(:,index));
                        title(strcat('Frequency domain, current reference data (',string(obj.anaRange(index)),')'));
                        xlabel('Frequency (MHz)');
                        ylabel('Power Spectral Density (W/Hz)');
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end
                end
            end
        end
        
        % Plot absorbance data
        function plotAbsorbance(obj)
            figure;
            plot(obj.freq/1e12,obj.absorbance)
            title('Absorbance of averaged data')
            xlabel('Frequency (THz)')
            ylabel('Absorbance')
        end
        %% Apodization windowing 

        % Displays the locations where the data would be truncated when the
        % ablation spike removal is applied. "len" is the rectangular
        % apodization window length, this function can be called without
        % specifying len if apoLen was specified with its mutator. Calling
        % this function sets apoLen to the value of len in us. Note that
        % this will not work if the apodization window is large enough to
        % contain the ablation spike as well.
        function obj = showApoWindow(obj,len)
            
            if nargin() == 2
                obj.apoLen = len; %In microseconds
            end

            numInd = round(0.5e-6*obj.apoLen*obj.daqFreq);
            centInd = round(length(obj.sig(:,1))/2);
            
            % Create figure
            figure;
            
            index = 1; % Start with the first data set in the analysis range    
            % Micro second time scale
            tus = obj.t;
            
            % Plot current data
            plot(tus,obj.sig(:,index));
            title(strcat('Time domain, current signal data (',string(obj.anaRange(index)),')'));
            xlabel('Time (us)');
            ylabel('Amplitude (A.U.)');
            hold on
            % Find peaks and label them
            [pks,locs] = findpeaks(obj.sig(centInd-numInd:centInd+numInd,index),'MinPeakProminence',obj.peakProm,'MinPeakDistance',200);
            locs = locs+centInd-numInd-1;
            
            % Show the apodization window
            % Find the window positions
            xPos1 = [tus(locs(1))-obj.apoLen;tus(locs(1))-obj.apoLen];
            xPos2 = [tus(locs(1))+obj.apoLen;tus(locs(1))+obj.apoLen];
            yPos1 = [-pks(1),pks(1)];
            yPos2 = [-pks(1),pks(1)];
            plot(xPos1,yPos1,'color','r')
            plot(xPos2,yPos2,'color','r')  
            hold off
            
            plotting = true;
            while plotting
                
                % Wait for user keypress
                if waitforbuttonpress
                    key = get(gcf, 'CurrentCharacter');

                    [pks,locs] = findpeaks(obj.sig(centInd-numInd:centInd+numInd,index),'MinPeakProminence',obj.peakProm,'MinPeakDistance',200);
                    locs = locs+centInd-numInd-1;
                    % Find the window positions
                    xPos1 = [tus(locs(1))-obj.apoLen;tus(locs(1))-obj.apoLen];
                    xPos2 = [tus(locs(1))+obj.apoLen;tus(locs(1))+obj.apoLen];
                    yPos1 = [-pks(1),pks(1)];
                    yPos2 = [-pks(1),pks(1)];
                    
                    % Handle key presses
                    if key == 28 % Left arrow
                        index = max(index - 1, 1);
                        clf;
                        % Plot current data
                        plot(tus,obj.sig(:,index));
                        title(strcat('Time domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        hold on
                        % Show the apodization window
                        plot(xPos1,yPos1,'color','r')
                        plot(xPos2,yPos2,'color','r')  
                        hold off
                    elseif key == 29 % Right arrow
                        index = min(index + 1, obj.num);
                        clf;
                        % Plot current data
                        plot(tus,obj.sig(:,index));
                        title(strcat('Time domain, current signal data (',string(obj.anaRange(index)),')'));
                        xlabel('Time (us)');
                        ylabel('Amplitude (A.U.)');
                        hold on
                        % Show the apodization window
                        plot(xPos1,yPos1,'color','r')
                        plot(xPos2,yPos2,'color','r')  
                        hold off
                    elseif key == 27 % Escape key to exit
                        close;
                        plotting = false;
                    end  
                end
            end
        end

% Apply a rectangular apodization windo to the data of length apoLen
% centered around the interferograms
        function obj = applyApoWindow(obj)
            % Find the interferogram peak, assuming there is two
            % peaks in the signal data, an ablation spike and an
            % interferogram spike
            
            % Index is half the index length of the apodization window
            index = round(1e-6*obj.apoLen*obj.daqFreq);
            tempSig = zeros(2*index+1,obj.num);
            tempRef = zeros(2*index+1,obj.num);
            
            % Only search the center apodization window of the array, the user has to
            % make sure the interferograms are centered
            numInd = round(0.5e-6*obj.apoLen*obj.daqFreq);
            centInd = round(length(obj.sig(:,1))/2);
            
            for i=1:obj.num
                [~,locs] = findpeaks(obj.sig(centInd-numInd:centInd+numInd,i),'MinPeakProminence',obj.peakProm,'MinPeakDistance',200);
                locs = locs+centInd-numInd-1;
                tempSig(:,i) = obj.sig(locs(1)-index:locs(1)+index,i);

                [~,locs] = findpeaks(obj.ref(centInd-numInd:centInd+numInd,i),'MinPeakProminence',obj.peakProm,'MinPeakDistance',200);
                locs = locs+centInd-numInd-1;
                tempRef(:,i) = obj.ref(locs(1)-index:locs(1)+index,i);
            end

            obj.sig = tempSig;
            obj.ref = tempRef;

            % Update time vector
            obj.t = obj.t(1:2*index+1);
        end

        % Apply a Blackman Harris apodization window to the data,
        % applyApoWindow needs to be called first
        function obj = blackmanHarris(obj)
            % Generate the window
            obj.window = blackmanharris(2*round(1e-6*obj.apoLen*obj.daqFreq)+1);
            
            % Apply the window
            obj.sig = obj.sig.*obj.window;
            obj.ref = obj.ref.*obj.window;
        end
        %% Remove cubic spline background from the data
        
        % Signal and reference data get subtracted by the data after a cubic
        % smoothing spline is applied. Cubic Smoothing Spline, is used for 
        % smoothing data while preserving trends.
        function obj = removeSplineBackground(obj)
            for i=1:obj.num
                spline = csaps(obj.t,obj.sig(:,i),0.8);
                obj.sig(:,i) = obj.sig(:,i) - fnval(spline,obj.t)';
            end
        end
        
        %% Zero Padding
        % Pad both sides of the time domain signal with the number of zeros
        % specified in us, makes sure the number of points is a power of 2
        % for fft speed
        function obj = zeroPad(obj,numZero)
            
            nZ = 1e-6*numZero*obj.daqFreq;
            n = length(obj.sig(:,1));
            
            p = nextpow2(n+nZ);
            nZ = round((2^p-n)/2);
            obj.sig = padarray(obj.sig,[nZ 0]);
            obj.ref = padarray(obj.ref,[nZ 0]);
            % Make a new time array for plotting reasons
            obj.t = 1e6*(0:length(obj.sig(:,1))-1)./obj.daqFreq;
            % Save the amount of zeros added per side in us of data
            obj.padNum = nZ/obj.daqFreq*1e6;
        end
        %% FFT and MERTZ phase
        
        % This function performs an fft and applies Mertz phase correction
        % on the signal data and reference data. It also generates a
        % frequency vector to go with the data. Note the Mertz phase data
        % is taken over the whole data range as recommended by: "Fourier
        % Transform Infrared Spectrometry" By: Peter R. Griffiths, James A.
        % de Haseth by default. If a second argument is passed, it should
        % be the Mertz phase window size in us.
        function obj = fftAndMertz(obj,mertzWind)
            sigLen = ceil(length(obj.sig(:,1))/4);
            obj.sigSpec = fft(obj.sig,[],1);
            obj.sigSpec = obj.sigSpec(1:sigLen,:);
            obj.refSpec = fft(obj.ref,[],1);
            obj.refSpec = obj.refSpec(1:sigLen,:);
            
            if nargin == 1
                % Apply Mertz phase correction
                obj.sigSpec = real(obj.sigSpec).^2./abs(obj.sigSpec) + imag(obj.sigSpec).^2./abs(obj.sigSpec);
                obj.refSpec = real(obj.refSpec).^2./abs(obj.refSpec) + imag(obj.refSpec).^2./abs(obj.refSpec); 
            
            elseif nargin == 2
                mertzSig = zeros(size(obj.sig));
                mertzRef = zeros(size(obj.sig));
                
                HalfWindSize = round(0.5e-6*mertzWind*obj.daqFreq);
                centInd = round(length(obj.sig(:,1))/2);
                mertzSig(centInd-HalfWindSize:centInd+HalfWindSize,:) = obj.sig(centInd-HalfWindSize:centInd+HalfWindSize,:);
                mertzRef(centInd-HalfWindSize:centInd+HalfWindSize,:) = obj.ref(centInd-HalfWindSize:centInd+HalfWindSize,:);
                
                mertzSig = fft(mertzSig,[],1);
                mertzSig = mertzSig(1:sigLen,:);
                mertzRef = fft(mertzRef,[],1);
                mertzRef = mertzRef(1:sigLen,:);
                
                obj.sigSpec = real(obj.sigSpec).*real(mertzSig)./abs(mertzSig) + imag(obj.sigSpec).*imag(mertzSig)./abs(mertzSig);
                obj.refSpec = real(obj.refSpec).*real(mertzRef)./abs(mertzRef) + imag(obj.refSpec).*imag(mertzRef)./abs(mertzRef);
            end
            
            % Scale appropriately to get Watts with a 50 Ohm load
            % assumption. This will give power spectral density, which you 
            % get by squaring the amplitude of the signal, then dividing by
            % ENBW of the apodization window. 
            
            % *******************************************
            % Need to convert GaGe Card
            % amplitude to volts for this to make sense.
            % *******************************************
            
%             scaleFact = sum(obj.window);
%             ENBW = obj.daqFreq*sum(abs(obj.window.^2))./abs(sum(obj.window)).^2;
%             
%             obj.sigSpec = 25*(obj.sigSpec/scaleFact).^2./ENBW;
%             obj.refSpec = 25*(obj.refSpec/scaleFact).^2./ENBW;
            
            len = length(obj.sig(:,1));
            freqBin = obj.daqFreq/len;
            % Create frequency axis in Hz
            obj.freq = (0:length(obj.sigSpec)-1)*freqBin*obj.frep/obj.dFrep;
            obj.numSpec = obj.num;
        end
        
        %% Frequency Domain Averaging and Absorbance Calculation
        % Average frequency domain data, this reduces sigSpec and refSpec
        % to have 1 column and changes num to have a value of 
        function obj = freqAvg(obj)
            obj.sigSpec = mean(obj.sigSpec,2);
            obj.refSpec = mean(obj.refSpec,2);
            obj.numSpec = 1;
            obj.absorbance = log10(obj.refSpec./obj.sigSpec);
        end
        
        
        %% Accessors and Mutators
        % Rectangular apodization window size mutator
        function obj = setApoLen(obj,apoLen)
            obj.apoLen = apoLen;
        end

        % Set analysis range, this throws away the data outside of the
        % analysis range, data range should take the form [start;end]
        function obj = setAnaRange(obj,anaRange)
            obj.anaRange = anaRange(1):anaRange(2);
            obj.sig = obj.sig(:,anaRange(1):anaRange(2));
            obj.ref = obj.ref(:,anaRange(1):anaRange(2));
            obj.num = anaRange(2)-anaRange(1)+1;
        end
        
        % Set the daq frequency in samples/s
        function obj = setDaqFreq(obj,daqFreq)
            obj.daqFreq = daqFreq;
        end
        
        % Set the peak prominence for peak finding
        function obj = setPeakProm(obj,peakProm)
            obj.peakProm = peakProm;
        end
        
        %% Custom saving, saving the whole object does not work, it's too much data
        function saveData(obj,folderName)
            
            analysisRange = obj.anaRange;

            daqFreq = obj.daqFreq;

            sigSpectrum = obj.sigSpec;
            refSpectrum = obj.refSpec;

            absorbance = obj.absorbance;
            frep = obj.frep;
            dFrep = obj.dFrep;
            freq = obj.freq;
            apoLen = obj.apoLen;
            zeroPadNum = obj.padNum;
            apoWindow = obj.window;
            
            % Make directory for the analyzed data
            mkdir(fileparts(obj.sigLoc)+"\"+folderName);
            
            saveLoc = fileparts(obj.sigLoc)+"\"+folderName+"\"+"analyzed_Data"; 
            
            save(saveLoc,"analysisRange","daqFreq","sigSpectrum","refSpectrum"...
                ,"absorbance","frep","dFrep","freq","apoLen","zeroPadNum","apoWindow");
        end
    end
end

