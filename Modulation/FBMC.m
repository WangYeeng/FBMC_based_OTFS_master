classdef FBMC < handle
    % =====================================================================        
    % This MATLAB class represents an implementation of FBMC. The 
    % modulation parameters are initialized by the class contructor. 
    % The modulation of data symbols x and the demodulation of the
    % received samples r is then performed by the methods ".Modulation(x)"
    % and ".Demodulation(r)".
    % =====================================================================    
    % WangYing, wangyingstu@163.com//
    % https://github.com/WangYeeng/FBMC_based_OTFS_master
    % =====================================================================   
    
    properties (SetAccess = private)
        Method              % defines the modulation method (prototype filter)
        Nr                  % for dimensionless parameters
        PHY                 % for parameters with physical interpretation
        PrototypeFilter     % for prototype filter parameters
        Implementation      % implmentation relevent parameters
    end
    
    
    methods
        % Class constructor, define default values.
        function obj = FBMC(varargin)
            % Initialize parameters, set default values
            if numel(varargin) == 10
                obj.Nr.Subcarriers                      = varargin{1};      % Number of subcarriers
                obj.Nr.MCSymbols                        = varargin{2};      % Number FBMC symbols in time
                obj.PHY.SubcarrierSpacing               = varargin{3};      % Subcarrier spacing (Hz)
                obj.PHY.SamplingRate                    = varargin{4};      % Sampling rate (Samples/s)
                obj.PHY.IntermediateFrequency           = varargin{5};      % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
                obj.PHY.TransmitRealSignal              = varargin{6};      % Transmit real valued signal (sampling theorem must be fulfilled!)
                obj.Method                              = varargin{7};      % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
                obj.PrototypeFilter.OverlappingFactor   = varargin{8};      % Overlapping factor (also determines oversampling in the frequency domain)
                obj.Implementation.InitialPhaseShift    = varargin{9};      % Initial phase shift
                obj.Implementation.UsePolyphase         = varargin{10};     % Efficient IFFT implementation, true or false
            elseif numel(varargin) == 0
                % Default Values (can later be changed using FBMC.Set...)
                obj.Nr.Subcarriers                      = 12;
                obj.Nr.MCSymbols                        = 30;
                obj.PHY.SubcarrierSpacing               = 15e3;
                obj.PHY.SamplingRate                    = obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing;
                obj.PHY.IntermediateFrequency           = 0;
                obj.PHY.TransmitRealSignal              = false;
                obj.Method                              = 'Hermite-OQAM';
                obj.PrototypeFilter.OverlappingFactor   = 8;
                obj.Implementation.InitialPhaseShift    = 0;
                obj.Implementation.UsePolyphase         = true;
            else
                error('Number of input variables must be either 0 (default values) or 10');
            end
                        
            % calculate and set all dependent parameters
            obj.SetDependentParameters();       
        end
        
        function SetDependentParameters(obj)
            % method that sets all parameters which are dependent on other parameters
            
            % Check Parameters
            if mod(obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing),1)~=0
                obj.PHY.SubcarrierSpacing=obj.PHY.SamplingRate/(2*round(obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing)));
                disp('Sampling Rate divided by (Subcarrier spacing times 2) must be must be an integer!');
                disp(['Therefore, the subcarrier spacing is set to: ' int2str(obj.PHY.SubcarrierSpacing) 'Hz']);
            end
            
            if mod(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing,1)~=0
                obj.PHY.IntermediateFrequency = round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing)*obj.PHY.SubcarrierSpacing;
                disp('The intermediate frequency must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the intermediate frequency is set to ' int2str(obj.PHY.IntermediateFrequency) 'Hz']);
            end
            
            if (obj.PHY.SamplingRate<obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)
                error('Sampling Rate must be higher: at least Number of Subcarriers times Subcarrier Spacing');
            end
            
            % dependent parameters
            obj.PHY.dt = 1/obj.PHY.SamplingRate;
            
            % Different Prototype Filters and OQAM or QAM
            switch obj.Method
                case 'Hermite-OQAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_Hermite(obj.PHY.TimeSpacing*2,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor/2);
                case 'Hermite-QAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)*2;
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor*4;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_Hermite(obj.PHY.TimeSpacing,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor);
                case 'Hermite&Nuttall-OQAM'
                     obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor; 
                    ProFilter = PrototypeFilter_Hermite(obj.PHY.TimeSpacing*2,obj.PHY.dt,2/2);
                    FFTSize = round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing));    
                    PrototypeFilter_New = zeros(size(ProFilter));   
                    O = 1.5;
                    PrototypeFilter_Temp = ProFilter(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1));
                    PrototypeFilter_Temp = PrototypeFilter_Temp.*nuttallwin(length(PrototypeFilter_Temp)); 
                    PrototypeFilter_New(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1)) = PrototypeFilter_Temp;                 
                    obj.PrototypeFilter.TimeDomain = zeros(obj.PrototypeFilter.OverlappingFactor*2*obj.Implementation.TimeSpacing,1);
                    obj.PrototypeFilter.TimeDomain(end/2+((-obj.Implementation.TimeSpacing*2+1):(obj.Implementation.TimeSpacing*2))) = PrototypeFilter_New;
                    obj.PrototypeFilter.TimeDomain = obj.PrototypeFilter.TimeDomain/sqrt(sum(abs(obj.PrototypeFilter.TimeDomain).^2)/obj.PHY.SamplingRate);
                case 'RRC-OQAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_RootRaisedCosine(obj.PHY.TimeSpacing*2,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor/2);
                case 'RRC-QAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)*2;
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor*4;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_RootRaisedCosine(obj.PHY.TimeSpacing,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor);
                case 'PHYDYAS-OQAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_PHYDYAS(obj.PHY.TimeSpacing*2,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor/2);
                case 'PHYDYAS-QAM'
                    obj.Implementation.TimeSpacing      = obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)*2;
                    obj.PHY.TimeSpacing                 = obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor*4;
                    obj.PrototypeFilter.TimeDomain      = PrototypeFilter_PHYDYAS(obj.PHY.TimeSpacing,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor);
                case 'PSRO-QAM'
                    obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;
                    obj.PrototypeFilter.TimeDomain = PrototypeFilter_PositiveSymmetricalRollOff(obj.PHY.TimeSpacing,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor); 
                case 'PSRO-OQAM'
                    obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;
                    obj.PrototypeFilter.TimeDomain = PrototypeFilter_PositiveSymmetricalRollOff(obj.PHY.TimeSpacing*2,obj.PHY.dt,obj.PrototypeFilter.OverlappingFactor/2); 
                case 'tRRC-OQAM'
                    obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor;                     
                    FFTSize = round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing));    
                    obj.PrototypeFilter.TimeDomain = zeros(obj.PrototypeFilter.OverlappingFactor*FFTSize,1);  
                    obj.PrototypeFilter.TimeDomain(end/2+(-FFTSize/2:FFTSize/2-1)+1)=sqrt(cos(2*pi*(-FFTSize/2:FFTSize/2-1)/FFTSize)+1);   
                    obj.PrototypeFilter.TimeDomain = obj.PrototypeFilter.TimeDomain/sqrt(sum(abs(obj.PrototypeFilter.TimeDomain).^2)/obj.PHY.SamplingRate);
                case 'Hermite&Tukey-OQAM'
                     obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor; 
                    ProFilter = PrototypeFilter_Hermite(obj.PHY.TimeSpacing*2,obj.PHY.dt,2/2);
                    FFTSize = round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing));    
                    PrototypeFilter_New = zeros(size(ProFilter));   
                    O = 1.5;
                    WindowFactor = 1;
                    PrototypeFilter_Temp = ProFilter(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1));
                    PrototypeFilter_Temp = PrototypeFilter_Temp.*tukeywin(length(PrototypeFilter_Temp),WindowFactor); 
                    PrototypeFilter_New(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1)) = PrototypeFilter_Temp;                 
                    obj.PrototypeFilter.TimeDomain = zeros(obj.PrototypeFilter.OverlappingFactor*2*obj.Implementation.TimeSpacing,1);
                    obj.PrototypeFilter.TimeDomain(end/2+((-obj.Implementation.TimeSpacing*2+1):(obj.Implementation.TimeSpacing*2))) = PrototypeFilter_New;
                    obj.PrototypeFilter.TimeDomain = obj.PrototypeFilter.TimeDomain/sqrt(sum(abs(obj.PrototypeFilter.TimeDomain).^2)/obj.PHY.SamplingRate);
                case 'Hermite&Kaiser-OQAM'
                    obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(2*obj.PHY.SubcarrierSpacing);
                    obj.PHY.TimeSpacing =obj.Implementation.TimeSpacing*obj.PHY.dt;
                    obj.Implementation.FrequencySpacing = obj.PrototypeFilter.OverlappingFactor; 
                    ProFilter = PrototypeFilter_Hermite(obj.PHY.TimeSpacing*2,obj.PHY.dt,2/2);
                    FFTSize = round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing));    
                    PrototypeFilter_New = zeros(size(ProFilter));   
                    O = 1.5;
                    WindowFactor = 9;
                    PrototypeFilter_Temp = ProFilter(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1));
                    PrototypeFilter_Temp = PrototypeFilter_Temp.*kaiser(length(PrototypeFilter_Temp),WindowFactor);
                    PrototypeFilter_New(end/2+1+(-round(FFTSize/2*O):round(O*FFTSize/2)-1)) = PrototypeFilter_Temp;                 
                    obj.PrototypeFilter.TimeDomain = zeros(obj.PrototypeFilter.OverlappingFactor*2*obj.Implementation.TimeSpacing,1);
                    obj.PrototypeFilter.TimeDomain(end/2+((-obj.Implementation.TimeSpacing*2+1):(obj.Implementation.TimeSpacing*2))) = PrototypeFilter_New;
                    obj.PrototypeFilter.TimeDomain = obj.PrototypeFilter.TimeDomain/sqrt(sum(abs(obj.PrototypeFilter.TimeDomain).^2)/obj.PHY.SamplingRate);
                otherwise
                    error(['Method (prototype filter) "' obj.Method '" is not supported']);
            end
            obj.Nr.SamplesPrototypeFilter   = length(obj.PrototypeFilter.TimeDomain);
            obj.Nr.SamplesTotal             = obj.Nr.SamplesPrototypeFilter+(obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing;
            
            % We assume a symmetric filter (=> real) and set small values to zero (=> lower computational complexity)
            PrototypefilterFFT = obj.PHY.dt*(fft(circshift(obj.PrototypeFilter.TimeDomain,obj.Nr.SamplesPrototypeFilter/2)));
            PrototypefilterFFT(abs(PrototypefilterFFT)./PrototypefilterFFT(1)<10^-14)=0;
            obj.PrototypeFilter.FrequencyDomain = PrototypefilterFFT;
            
            % Prepare paramters for an efficient implementation
            % Phase shift to guarentee that the interference is purely imaginary
            [k,l] = meshgrid(0:obj.Nr.MCSymbols-1,0:obj.Nr.Subcarriers-1);
            obj.Implementation.PhaseShift = exp(1j*pi/2*(l+k))*exp(1j*obj.Implementation.InitialPhaseShift);
            
            % index for the time shift of different FBMC symbols
            IndexAfterIFFT =[ones(obj.Nr.SamplesPrototypeFilter,obj.Nr.MCSymbols);zeros(obj.Implementation.TimeSpacing*(obj.Nr.MCSymbols-1),obj.Nr.MCSymbols)];
            IndexNumberAfterIFFT = zeros(obj.Nr.SamplesPrototypeFilter,obj.Nr.MCSymbols);
            for i_k = 1:obj.Nr.MCSymbols
                IndexAfterIFFT(:,i_k) = circshift(IndexAfterIFFT(:,i_k),(i_k-1)*obj.Implementation.TimeSpacing);
                IndexNumberAfterIFFT(:,i_k) = find(IndexAfterIFFT(:,i_k));
            end
            obj.Implementation.IndexAfterIFFT = logical(IndexAfterIFFT);
            obj.Implementation.IndexNumberAfterIFFT = IndexNumberAfterIFFT; 
            
            % index for the polyphase implementation
            obj.Implementation.FFTSize               = round(obj.Nr.SamplesPrototypeFilter./obj.Implementation.FrequencySpacing);
            obj.Implementation.IntermediateFrequency = round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing);
            obj.Implementation.IndexPolyphaseMap     = logical(circshift([ones(obj.Nr.Subcarriers,obj.Nr.MCSymbols);...
                zeros(obj.Implementation.FFTSize-obj.Nr.Subcarriers,obj.Nr.MCSymbols)],...
                [obj.Implementation.IntermediateFrequency 1]));
            
            % Normalization factor so that the default power = 1
            obj.Implementation.NormalizationFactor = sqrt(obj.PHY.SamplingRate^2/obj.PHY.SubcarrierSpacing^2*obj.PHY.TimeSpacing/obj.Nr.Subcarriers);
        end
        
        % Set Functions
        function SetNrSubcarriers(varargin)
            % set the number of subcarriers
            
            obj = varargin{1};
            % set specific property
            obj.Nr.Subcarriers = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetNrMCSymbols(varargin)
            % set the number of symbols
            
            obj = varargin{1};
            % set specific property
            obj.Nr.MCSymbols = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetSubcarrierSpacing(varargin)
            % set the subcarrier spacing
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SubcarrierSpacing = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetSamplingRate(varargin)
            % set the sampling rate
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SamplingRate = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetIntermediateFrequency(varargin)
            % set intermediate frequency
            
            obj = varargin{1};
            % set specific property
            obj.PHY.IntermediateFrequency = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetTransmitRealSignal(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            obj.PHY.TransmitRealSignal = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetMethod(varargin)
            % set method (prototype filter)
            
            obj = varargin{1};
            % set specific property
            obj.Method = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetOverlappingFactor(varargin)
            % set overlapping factor
            
            obj = varargin{1};
            % set specific property
            obj.PrototypeFilter.OverlappingFactor = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetInitialPhaseShift(varargin)
            % set initial phase shift
            
            obj = varargin{1};
            % set specific property
            obj.Implementation.InitialPhaseShift = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        function SetUsePolyphase(varargin)
            % set polyphase filter indicator
            
            obj = varargin{1};
            % set specific property
            obj.Implementation.UsePolyphase = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
     
        
        % Modulation and Demodulation
        function TransmitSignal = Modulation(obj, DataSymbols)
            % modulates the data symbols. The input argument is a matrix
            % of size "Number of subcarriers" \times "Number of FBMC symbols"
            % which represents the transmit data symbols
            
            TransmitSignal = zeros(size(obj.Implementation.IndexAfterIFFT));
            if obj.Implementation.UsePolyphase
                DataSymbolsTemp = zeros(size(obj.Implementation.IndexPolyphaseMap));
                DataSymbolsTemp(obj.Implementation.IndexPolyphaseMap) = DataSymbols.*obj.Implementation.PhaseShift.*obj.Implementation.NormalizationFactor;
                if obj.PHY.TransmitRealSignal
                    DataSymbolsTemp = (DataSymbolsTemp+conj(DataSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
                end
                TransmitSignal(obj.Implementation.IndexAfterIFFT) = repmat(ifft(DataSymbolsTemp),[obj.Implementation.FrequencySpacing 1]).*repmat(obj.PrototypeFilter.TimeDomain,[1 obj.Nr.MCSymbols]);
                TransmitSignal = sum(TransmitSignal,2);
            else
                % Include also the design in the frequency domain because it provides an alternative understanding of FBMC, but is less efficient!
                % Note that this matrix could be precomputed!
                FrequencyGeneration = zeros(size(obj.PrototypeFilter.FrequencyDomain,1),obj.Nr.Subcarriers);
                for i_l=1:obj.Nr.Subcarriers
                    FrequencyGeneration(:,i_l) = circshift(obj.PrototypeFilter.FrequencyDomain, obj.Implementation.FrequencySpacing*(i_l-1+obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing));
                end
                % normalized to have power 1
                FrequencyGeneration = FrequencyGeneration*sqrt(1/obj.Nr.Subcarriers*obj.PHY.TimeSpacing)*obj.PHY.SamplingRate;
                FrequencyDomain = FrequencyGeneration*(DataSymbols.*obj.Implementation.PhaseShift);
                if obj.PHY.TransmitRealSignal
                    FrequencyDomain = (FrequencyDomain+conj(FrequencyDomain([1 end:-1:2],:)))/sqrt(2);
                end
                TransmitSignal(obj.Implementation.IndexAfterIFFT) = circshift(ifft(FrequencyDomain),[-obj.Nr.SamplesPrototypeFilter/2 0]);
                TransmitSignal = sum(TransmitSignal,2);
            end
        end
        
        function ReceivedSymbols = Demodulation(obj, ReceivedSignal)
            % demodulates the received time signal and returns a matrix of 
            % size "Number of subcarriers" \times "Number of FBMC symbols"
            % which represents the received symbols after demodulation but
            % before channel equalization
            ReceivedSignal_CorrespondingSamplesToFBMCSymbol = ReceivedSignal(obj.Implementation.IndexNumberAfterIFFT);
            if obj.Implementation.UsePolyphase
                % similar to the transmitter, just in reversed order
                FilteredReceivedSignal = ReceivedSignal_CorrespondingSamplesToFBMCSymbol.*repmat(obj.PrototypeFilter.TimeDomain,[1 obj.Nr.MCSymbols]);
                ReceivedSymbolsTemp = fft(sum(reshape(FilteredReceivedSignal,[size(obj.Implementation.IndexPolyphaseMap,1) obj.Implementation.FrequencySpacing obj.Nr.MCSymbols]),2));
                if obj.PHY.TransmitRealSignal
                    ReceivedSymbolsTemp = ReceivedSymbolsTemp *sqrt(2);
                end
                ReceivedSymbols = reshape(ReceivedSymbolsTemp(obj.Implementation.IndexPolyphaseMap),size(obj.Implementation.PhaseShift)).*conj(obj.Implementation.PhaseShift)/(obj.Implementation.NormalizationFactor*obj.PHY.SubcarrierSpacing);
            else
                % same as before
                FrequencyGeneration = zeros(size(obj.PrototypeFilter.FrequencyDomain,1),obj.Nr.Subcarriers);
                for i_l=1:obj.Nr.Subcarriers
                    FrequencyGeneration(:,i_l) = circshift(obj.PrototypeFilter.FrequencyDomain, obj.Implementation.FrequencySpacing*(i_l-1+obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing));
                end
                % normalized to have power 1
                FrequencyGeneration = FrequencyGeneration*sqrt(1/obj.Nr.Subcarriers*obj.PHY.TimeSpacing)*obj.PHY.SamplingRate;
                
                FrequencyDomain = fft(circshift(ReceivedSignal_CorrespondingSamplesToFBMCSymbol,[obj.Nr.SamplesPrototypeFilter/2 0]));
                ReceivedSymbols = (FrequencyGeneration'*FrequencyDomain).*conj(obj.Implementation.PhaseShift)*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)/(obj.PHY.SamplingRate^2*obj.PHY.TimeSpacing*obj.Implementation.FrequencySpacing);
            end
        end
        
        % Matrix Description
        function TXMatrix = GetTXMatrix(obj)
            % returns a matrix G so that s=G*x(:) is the same as 
            % s=obj.Modulation(x)
            
            TransmitRealSignal=false;
            if obj.PHY.TransmitRealSignal
                obj.PHY.TransmitRealSignal = false;
                TransmitRealSignal=true;
            end
            TXMatrix=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers*obj.Nr.MCSymbols);
            TXMatrixTemp=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers);
            x = zeros(obj.Nr.Subcarriers, obj.Nr.MCSymbols);
            for i_l= 1:obj.Nr.Subcarriers
                x(i_l)=1;
                TXMatrixTemp(:,i_l) = obj.Modulation(x);
                x(i_l)=0;
            end
            for i_k=1:obj.Nr.MCSymbols
                TXMatrix(:,(1:obj.Nr.Subcarriers)+(i_k-1)*obj.Nr.Subcarriers)=circshift(TXMatrixTemp,[(i_k-1)*obj.Implementation.TimeSpacing,0])*1j^(i_k-1);
            end
            if TransmitRealSignal
                obj.PHY.TransmitRealSignal = true;
                TXMatrix = real(TXMatrix)*sqrt(2);
            end
        end
        function RXMatrix = GetRXMatrix(obj)
            % returns a matrix Q so that y=Q*r is the same as 
            % y=reshape(obj.Demodulation(r),[],1)            
            
            if obj.PHY.TransmitRealSignal
                obj.PHY.TransmitRealSignal = false;
                RXMatrix = sqrt(2)*obj.GetTXMatrix'*(obj.Nr.Subcarriers/(obj.PHY.SamplingRate*obj.PHY.TimeSpacing));
                obj.PHY.TransmitRealSignal = true;
            else
                RXMatrix = obj.GetTXMatrix'*(obj.Nr.Subcarriers/(obj.PHY.SamplingRate*obj.PHY.TimeSpacing));
            end
        end
        function FBMCMatrix = GetFBMCMatrix(obj,FastCalculation)
            % returns a matrix D, so that y=D*x is the same as 
            % y=reshape(obj.Demodulation(obj.Modulation(x)),[],1)  
            
            if not(exist('FastCalculation','var'))
                FastCalculation = true;
            end
            
            if FastCalculation
                InterferenceMatrix = obj.GetInterferenceMatrix;
                [Symbol_i,Subcarrier_i] = meshgrid(1:obj.Nr.MCSymbols,1:obj.Nr.Subcarriers);
                LK = obj.Nr.Subcarriers*obj.Nr.MCSymbols;
                Index_DeltaSubcarrier = repmat(Subcarrier_i(:),1,LK)-repmat(Subcarrier_i(:)',LK,1);
                Index_DeltaSymbols = repmat(Symbol_i(:),1,LK)-repmat(Symbol_i(:)',LK,1);
                Index_Subcarrier = repmat(Subcarrier_i(:),1,LK)-1;
                FBMCMatrix=reshape(InterferenceMatrix(Index_DeltaSubcarrier(:)+obj.Nr.Subcarriers+(Index_DeltaSymbols(:)+obj.Nr.MCSymbols-1)*size(InterferenceMatrix,1)),LK,LK);
                if obj.Method(end-3)=='-'
                    % QAM, works only for an odd number of subcarriers, needs to be fixed!  Hower it does not really matter because values are too small to have any effect!
                    FBMCMatrix=FBMCMatrix.*exp(-1j*pi/2*(Index_DeltaSubcarrier+Index_DeltaSymbols)).*exp(-1j*pi*3/2*Index_DeltaSymbols.*Index_DeltaSubcarrier);
                else
                    % OQAM, this works
                    FBMCMatrix=FBMCMatrix.*exp(-1j*pi/2*(Index_DeltaSubcarrier+Index_DeltaSymbols)).*exp(-1j*2*pi*(obj.PHY.TimeSpacing*obj.PHY.SubcarrierSpacing)*Index_DeltaSymbols.*(Index_Subcarrier+Index_DeltaSubcarrier/2));
                end
            else
                % straightforward slow implementation
                FBMCMatrix = zeros(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                DataImpulse=zeros(obj.Nr.Subcarriers,obj.Nr.MCSymbols);
                for i_lk=1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                    DataImpulse(i_lk)=1;
                    FBMCMatrix(:,i_lk) = reshape(obj.Demodulation(obj.Modulation(DataImpulse)),[],1);
                    DataImpulse(i_lk)=0;
                end
            end
        end

        function InterferenceMatrix = GetInterferenceMatrix(obj)
            % returns a matrix which represents the imaginary interference
            % weights in FBMC   

            DataSymbols=zeros(obj.Nr.Subcarriers,obj.Nr.MCSymbols);
            DataSymbols(1,1)=1;
            Y11=reshape(obj.Demodulation(obj.Modulation(DataSymbols)),obj.Nr.Subcarriers,obj.Nr.MCSymbols);
            [k_all,l_all] = meshgrid(0:obj.Nr.MCSymbols-1,0:obj.Nr.Subcarriers-1);
            Y11=Y11.*(exp(1j*pi/2*(l_all+k_all)).*exp(-1j*pi*k_all.*(l_all/2)));
            InterferenceMatrix = [Y11(end:-1:2,end:-1:2),Y11(end:-1:2,1:end);Y11(:,end:-1:2),Y11(1:end,1:end)];
        end

        function Pn = GetSymbolNoisePower(obj,Pn_time)
            % returns the symbol noise power, that is, the noise power
            % after demodulation. The input argument is the noise power 
            % in the time domain. 
          
            Pn = (Pn_time*(obj.Nr.Subcarriers/(obj.PHY.SamplingRate*obj.PHY.TimeSpacing)));
        end
    end
end

function PrototypeFilter = PrototypeFilter_Hermite(T0,dt,OF)
% The pulse is orthogonal for a time-spacing of T=T_0 and a frequency-spacing of F=2/T_0
% Values taken from: R.Nissel et al. "ON PILOT-SYMBOL AIDED CHANNEL ESTIMATION IN FBMC-OQAM"
    t_filter=-(OF*T0):dt:(OF*T0-dt);
    D0=1/sqrt(T0)*HermiteH(0,sqrt(2*pi)*(t_filter./(T0/sqrt(2))))  .*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    D4=1/sqrt(T0)*HermiteH(4,sqrt(2*pi)*(t_filter./(T0/sqrt(2))))  .*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    D8=1/sqrt(T0)*HermiteH(8,sqrt(2*pi)*(t_filter./(T0/sqrt(2))))  .*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    D12=1/sqrt(T0)*HermiteH(12,sqrt(2*pi)*(t_filter./(T0/sqrt(2)))).*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    D16=1/sqrt(T0)*HermiteH(16,sqrt(2*pi)*(t_filter./(T0/sqrt(2)))).*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    D20=1/sqrt(T0)*HermiteH(20,sqrt(2*pi)*(t_filter./(T0/sqrt(2)))).*exp(-pi*(t_filter./(T0/sqrt(2))).^2);
    H0= 1.412692577;
    H4= -3.0145e-3;
    H8=-8.8041e-6;
    H12=-2.2611e-9;
    H16=-4.4570e-15;
    H20 = 1.8633e-16;
    PrototypeFilter=(D0.*H0+D4.*H4+D8.*H8+D12.*H12+D16.*H16+D20.*H20).';
    PrototypeFilter = PrototypeFilter/sqrt(sum(abs(PrototypeFilter).^2)*dt);
end

function PrototypeFilter = PrototypeFilter_RootRaisedCosine(T0,dt,OF)
    % The pulse is orthogonal for a time-spacing of T=T0 and a frequency-spacing of F=2/T0
    t_filter=-(OF*T0):dt:(OF*T0-dt);
    PrototypeFilter=(1/sqrt(T0)*(4*t_filter/T0.*cos(2*pi*t_filter/T0))./(pi*t_filter/T0.*(1-(4*t_filter/T0).^2)));
    PrototypeFilter(abs(t_filter)<10^-14) =1/sqrt(T0)*(4/pi);
    PrototypeFilter(abs(abs(t_filter)-T0/4)<10^-14)=1/sqrt(2*T0)*((1+2/pi)*sin(pi/4)+(1-2/pi)*cos(pi/4));
    PrototypeFilter=PrototypeFilter.';
    PrototypeFilter = PrototypeFilter/sqrt(sum(abs(PrototypeFilter).^2)*dt);
end

function PrototypeFilter = PrototypeFilter_PHYDYAS(T0,dt,OF)
    % The pulse is orthogonal for a time-spacing of T=T0 and a frequency-spacing of F=2/T0
    t_filter=-(OF*T0):dt:(OF*T0-dt);
    switch OF*2
        case 2
            H = sqrt(2)/2;
        case 3
            H = [0.91143783 0.41143783];
        case 4
            H = [0.97195983 sqrt(2)/2 0.23514695];
        case 5
            H = [0.99184131 0.86541624 0.50105361 0.12747868];
        case 6
            H = [0.99818572 0.94838678 sqrt(2)/2 0.31711593 0.06021021];
        case 7
            H = [0.99938080 0.97838560 0.84390076 0.53649931 0.20678881 0.03518546];
        case 8
            H = [0.99932588 0.98203168 0.89425129 sqrt(2)/2 0.44756522 0.18871614 0.03671221];
        otherwise
            error('Oversampling factor must be an integer between 1 and 8 for OQAM or betwen 1 and 4 for QAM');
    end
    PrototypeFilter = 1+2*sum(repmat(H,length(t_filter),1).*cos(2*pi*repmat(t_filter',1,length(H)).*repmat(1:length(H),length(t_filter),1)/((length(H)+1)*T0)),2);
    PrototypeFilter = PrototypeFilter/sqrt(sum(abs(PrototypeFilter).^2)*dt);
end


function Hermite = HermiteH(n,x)
    % Hermite polynomials (obtained by Mathematica, "ToMatlab[HermiteH[n, x]]")
    if n==0
        Hermite=ones(size(x));
    elseif n==4
        Hermite=12+(-48).*x.^2+16.*x.^4;
    elseif n==8
        Hermite = 1680+(-13440).*x.^2+13440.*x.^4+(-3584).*x.^6+256.*x.^8;
    elseif n==12
        Hermite = 665280+(-7983360).*x.^2+13305600.*x.^4+(-7096320).*x.^6+1520640.* ...
            x.^8+(-135168).*x.^10+4096.*x.^12;
    elseif n==16
        Hermite = 518918400+(-8302694400).*x.^2+19372953600.*x.^4+(-15498362880).* ...
            x.^6+5535129600.*x.^8+(-984023040).*x.^10+89456640.*x.^12+( ...
            -3932160).*x.^14+65536.*x.^16;
    elseif n==20
        Hermite = 670442572800+(-13408851456000).*x.^2+40226554368000.*x.^4+( ...
            -42908324659200).*x.^6+21454162329600.*x.^8+(-5721109954560).* ...
            x.^10+866834841600.*x.^12+(-76205260800).*x.^14+3810263040.*x.^16+ ...
            (-99614720).*x.^18+1048576.*x.^20;
    end
end

function [PrototypeFilter] = PrototypeFilter_PositiveSymmetricalRollOff(T0,dt,OF)
% The pulse is orthogonal for a time-spacing of T=T0 and a frequency-spacing of F=2/T0
t_filter=-(OF*T0):dt:(OF*T0-dt);
switch 2*OF
    case 2
        p_PositiveSymmetricalRollOff_O2 = @(t,T0) (...
            sqrt(2)/2 * sinc((1+2*t/T0*2)) +...
            1 *         sinc((0+2*t/T0*2)) +...
            sqrt(2)/2 * sinc((1-2*t/T0*2)) ...
            )*sqrt(2/T0);            
        PrototypeFilter = p_PositiveSymmetricalRollOff_O2(t_filter.',T0);
    case 3
        p_PositiveSymmetricalRollOff_O3 = @(t,T0) (...
           0.41143783 * sinc((2+3*t/T0*2)) + ...
           0.91143783 * sinc((1+3*t/T0*2)) + ...
           1 *          sinc((0+3*t/T0*2)) + ...
           0.91143783 * sinc((1-3*t/T0*2)) + ...
           0.41143783 * sinc((2-3*t/T0*2)) ...    
           )*sqrt(2/T0);
        PrototypeFilter = p_PositiveSymmetricalRollOff_O3(t_filter.',T0);
    case 4
        p_PositiveSymmetricalRollOff_O4 = @(t,T0) (...
           0.23514695 * sinc((3+4*2*t/T0)) + ...
           sqrt(2)/2  * sinc((2+4*2*t/T0)) + ...
           0.97195983 * sinc((1+4*2*t/T0)) + ...
           1 *          sinc((0+4*2*t/T0)) + ...
           0.97195983 * sinc((1-4*2*t/T0)) + ...
           sqrt(2)/2  * sinc((2-4*2*t/T0)) + ...
           0.23514695 * sinc((3-4*2*t/T0))   ...    
           )*sqrt(2/T0);
        PrototypeFilter = p_PositiveSymmetricalRollOff_O4(t_filter.',T0);
    case 5
         p_PositiveSymmetricalRollOff_O5 = @(t,T0) (...
           0.12747868 * sinc((4+5*t/T0*2)) +...
           0.50105361 * sinc((3+5*t/T0*2)) + ...
           0.86541624 * sinc((2+5*t/T0*2)) + ...
           0.99184131 * sinc((1+5*t/T0*2)) + ...
           1 *          sinc((0+5*t/T0*2)) + ...
           0.99184131 * sinc((1-5*t/T0*2)) + ...
           0.86541624 * sinc((2-5*t/T0*2)) + ...
           0.50105361 * sinc((3-5*t/T0*2)) + ...
           0.12747868 * sinc((4-5*t/T0*2))   ... 
           )*sqrt(2/T0);
        PrototypeFilter = p_PositiveSymmetricalRollOff_O5(t_filter.',T0);
    case 6
        p_PositiveSymmetricalRollOff_O6 = @(t,T0) (...
           0.06021021 * sinc((5+6*t/T0*2)) + ... 
           0.31711593 * sinc((4+6*t/T0*2)) + ...
           sqrt(2)/2  * sinc((3+6*t/T0*2)) + ...
           0.94838678 * sinc((2+6*t/T0*2)) + ...
           0.99818572 * sinc((1+6*t/T0*2)) + ...
           1 *          sinc((0+6*t/T0*2)) + ...
           0.99818572 * sinc((1-6*t/T0*2)) + ...
           0.94838678 * sinc((2-6*t/T0*2)) + ...
           sqrt(2)/2  * sinc((3-6*t/T0*2)) + ...
           0.31711593 * sinc((4-6*t/T0*2)) + ...
           0.06021021 * sinc((5-6*t/T0*2))   ... 
           )*sqrt(2/T0);
        PrototypeFilter = p_PositiveSymmetricalRollOff_O6(t_filter.',T0);
    case 7
        p_PositiveSymmetricalRollOff_O7 = @(t,T0) (...
           0.03518546 * sinc((6+7*t/T0*2)) + ...
           0.20678881 * sinc((5+7*t/T0*2)) + ...
           0.53649931 * sinc((4+7*t/T0*2)) + ...
           0.84390076 * sinc((3+7*t/T0*2)) + ...
           0.97838560 * sinc((2+7*t/T0*2)) + ...
           0.99938080 * sinc((1+7*t/T0*2)) + ...
           1 *          sinc((0+7*t/T0*2)) + ...
           0.99938080 * sinc((1-7*t/T0*2)) + ...
           0.97838560 * sinc((2-7*t/T0*2)) + ...
           0.84390076 * sinc((3-7*t/T0*2)) + ...
           0.53649931 * sinc((4-7*t/T0*2)) + ...
           0.20678881 * sinc((5-7*t/T0*2)) + ...
           0.03518546 * sinc((6-7*t/T0*2))   ... 
           )*sqrt(2/T0);
        PrototypeFilter = p_PositiveSymmetricalRollOff_O7(t_filter.',T0);
    case 8
       p_PositiveSymmetricalRollOff_O8 = @(t,T0) (...
          0.03671221 * sinc((7+8*t/T0*2)) + ...    
          0.18871614 * sinc((6+8*t/T0*2)) + ...    
          0.44756522 * sinc((5+8*t/T0*2)) + ...    
          sqrt(2)/2  * sinc((4+8*t/T0*2)) + ...    
          0.89425129 * sinc((3+8*t/T0*2)) + ...
          0.98203168 * sinc((2+8*t/T0*2)) + ...
          0.99932588 * sinc((1+8*t/T0*2)) + ...
          1 *          sinc((0+8*t/T0*2)) + ...
          0.99932588 * sinc((1-8*t/T0*2)) + ...
          0.98203168 * sinc((2-8*t/T0*2)) + ...
          0.89425129 * sinc((3-8*t/T0*2)) + ...
          sqrt(2)/2  * sinc((4-8*t/T0*2)) + ...
          0.44756522 * sinc((5-8*t/T0*2)) + ...
          0.18871614 * sinc((6-8*t/T0*2)) + ...
          0.03671221 * sinc((7-8*t/T0*2))   ...    
          )*sqrt(2/T0);
     PrototypeFilter = p_PositiveSymmetricalRollOff_O8(t_filter.',T0);
    otherwise
     error('Oversampling factor must be an integer between 1 and 8 for OQAM or betwen 1 and 4 for QAM');
end
end

