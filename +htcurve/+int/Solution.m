classdef Solution % not handle
	
	properties (Access = public)
		% This indicates what problem it's a solution to
		memo (1,1) htcurve.int.SolutionMemo = htcurve.int.SolutionMemo();
		
		% This describes the solution
		solnIsDeferred (1,1) logical = false; % This controls how to interpret the following properties
		% If not deferred:
		order uint32 = uint32.empty(0,0); % starting from 0, incrementing up to max number of pixels less one
		% If deferred:
		numSubSolns (1,1) uint8 = 0;
		subSolutions (:,1) htcurve.int.SolutionPiece = htcurve.int.SolutionPiece.empty(0,1);
		reconstruction uint8 = uint16.empty(0,0); % will show how to reassemble the pieces into the whole region (1-based indices)
	end
	
end