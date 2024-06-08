classdef MemoParityList < handle % !!
	
	properties (Access = public)
		entries;
	end
	
	methods (Access = public)
		function this = MemoParityList(quant)
			this.entries = repmat( htcurve.MemoParity(), quant, 1 );
		end
	end
	
end