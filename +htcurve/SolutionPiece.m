classdef SolutionPiece % not handle
	
	% This class instance represents one piece of a Solution, and indicates
	% how to interpret in that grander context
	
	properties (Access = public)
		orderOffset           (1,1) uint32 = 0;
		subRegionWindow       (1,4) uint32 = zeros(1,4,'uint32'); % [first row, last row, first col, last col] (inclusive)
		transformToSubProblem (1,1) htcurve.Transform    = htcurve.Transform();
		subProblemMemo        (1,1) htcurve.SolutionMemo = htcurve.SolutionMemo();
	end
	% NOTE: the transformToSubProblem is intended for use *after*
	% truncating things to the subRegionWindow (and becoming relative to
	% it)
	
end