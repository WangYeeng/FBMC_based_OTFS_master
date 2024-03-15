classdef SignalConstellation < handle 
    % =====================================================================        
    % This MATLAB class represents a QAM or PAM signal constellation. The 
    % parameters are initialized by the class contructor. 
    % Afterwards we can transform a bit stream into symbols by 
    % ".Bit2Symbol()" and back from symbol to bit by ".Symbol2Bit".
    % =====================================================================    
    % WangYing, wangyingstu@163.com//
    % https://github.com/WangYeeng/FBMC_based_OTFS_master
    % =====================================================================     

  properties (SetAccess = private)
       Method
       ModulationOrder
       BitMapping
       SymbolMapping
       Implementation
  end
  
  methods
  	function obj = SignalConstellation(ModulationOrder,Method)
        % Generates a QAM or PAM signal constellation object with 
        % corresponding Gray-coded bit mapping. The first argument 
        % represents the modulation order and the second argument the
        % method, either 'QAM' or 'PAM'. For example (4,'QAM'), (256,'QAM')
        % or (2,'PAM'), (16,'PAM')
        
        obj.ModulationOrder = ModulationOrder;
        obj.Method          = Method;
        
        % Gray coded bitmapping
        if strcmp( obj.Method,'QAM' )
            BitMappingAtom = [ones(sqrt(obj.ModulationOrder)/2,1);zeros(sqrt(obj.ModulationOrder)/2,1)];
            for i_temp = 2:log2(sqrt(obj.ModulationOrder))
                BinaryTemp = BitMappingAtom(1:2:end,i_temp-1);
                BitMappingAtom(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            IQ = 2*(1:sqrt(obj.ModulationOrder))-sqrt(obj.ModulationOrder)-1;
            [I_rep,Q_rep]=meshgrid(IQ,IQ);
            obj.SymbolMapping = I_rep(:)+1i*Q_rep(:);
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
            obj.BitMapping = false(obj.ModulationOrder,log2(obj.ModulationOrder));
            for x_IQ = IQ
                obj.BitMapping(I_rep(:)==x_IQ,2:2:end) = BitMappingAtom;
                obj.BitMapping(Q_rep(:)==x_IQ,1:2:end) = BitMappingAtom;
            end
        elseif strcmp(obj.Method,'PAM')
            obj.BitMapping = [ones(obj.ModulationOrder/2,1);zeros(obj.ModulationOrder/2,1)];
            for i_temp = 2:log2(obj.ModulationOrder)
                BinaryTemp = obj.BitMapping(1:2:end,i_temp-1);
                obj.BitMapping(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            obj.SymbolMapping = (2*(1:obj.ModulationOrder)-obj.ModulationOrder-1).';
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
        else
           error('Signal constellation method must be QAM or PAM!');
        end
        
        % Determine the underlying symbol alphabet and the corresponding
        % bit mapping
        [~,SortOrder] = sort(bi2de(obj.BitMapping),'ascend');
        obj.SymbolMapping = obj.SymbolMapping(SortOrder);
        obj.BitMapping = obj.BitMapping(SortOrder,:);
        
        % For the LLR detection, we determine all data symbols which have a
        % bit value of one (zero) at a certain position
        obj.Implementation.DataSymbolsBitvalueOne   = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueOne   = reshape(obj.Implementation.DataSymbolsBitvalueOne(logical(obj.BitMapping)),obj.ModulationOrder/2,[]);
        obj.Implementation.DataSymbolsBitvalueZero  = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueZero  = reshape(obj.Implementation.DataSymbolsBitvalueZero(not(logical(obj.BitMapping))),obj.ModulationOrder/2,[]);       
    end   
    
    function DataSymbols = Bit2Symbol(obj,BinaryStream)
        % Maps a bit stream to the correpsonding symbol alphabet
        tmpSize = size(BinaryStream);
        DataSymbols = obj.SymbolMapping( bi2de(reshape(BinaryStream(:),log2(obj.ModulationOrder),[])')+1 );
        DataSymbols = reshape( DataSymbols, tmpSize(1)/log2(obj.ModulationOrder), tmpSize(2) );
    end

    function EstimatedBitStream = Symbol2Bit(obj,EstimatedDataSymbols)
        % Maps symbols (nearest neighbor detection) to the corresponding
        % bit stream
        EstimatedDataSymbols = EstimatedDataSymbols(:);
        
        [~,b] = min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        EstimatedBitStream = obj.BitMapping(b(:),:).';
        EstimatedBitStream = EstimatedBitStream(:);         
    end

    function LLR = LLR_AWGN(obj,y,Pn)
        % Calculates the LLR for an AWGN channel, that is, y=x+n with y
        % denoting the received data symbol, x the transmitted data symbol
        % and n the Gaussian distributed noise with power Pn
    
        if numel(Pn)>1
            PnRepeated = reshape(repmat(Pn.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        else 
            PnRepeated = Pn;
        end
  
        ReceivedDataSymbolsRepeated     = reshape(repmat(y.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        DataSymbolsBitvalueOneRepeated  = repmat(obj.Implementation.DataSymbolsBitvalueOne,1,length(y));
        DataSymbolsBitvalueZeroRepeated = repmat(obj.Implementation.DataSymbolsBitvalueZero,1,length(y));
        
        LLR =   log(sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueOneRepeated).^2./PnRepeated),1)./...
                    sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueZeroRepeated).^2./PnRepeated),1)).';
        LLR(LLR==Inf)=10^10;
        LLR(LLR==-Inf)=-10^10;     
    end
    
    function QuantizedDataSymbols = SymbolQuantization(obj,EstimatedDataSymbols)
        % Performs quantization of the received symbols, that is nearest
        % neighbor detection
        EstimatedDataSymbols = EstimatedDataSymbols(:);

        [~,b] =min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        QuantizedDataSymbols = obj.SymbolMapping(b(:),:).';
        QuantizedDataSymbols = QuantizedDataSymbols(:);
    end      
  end
end
      