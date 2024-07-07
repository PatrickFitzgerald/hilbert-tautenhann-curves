classdef MemoParity < htcurve.MemoBase % not handle
	
	properties (Access = public)
		isOddH (1,1) logical = false;
		isOddW (1,1) logical = false;
	end
	
	% Satisfying MemoBase inheritance
	methods (Access = public, Hidden) % not worth the hassle of giving the CurveGenerator access too...
		function parityH = getParityH(this)
			parityH = this.isOddH; % isOddH is true=1 when mod(H,2)=1
		end
		function parityW = getParityW(this)
			parityW = this.isOddW; % isOddW is true=1 when mod(W,2)=1
		end
		function effSizeH = getSizeH(this)
			effSizeH = this.isOddH + 2; % if false, show size 2. if true, show size 3
		end
		function effSizeW = getSizeW(this)
			effSizeW = this.isOddW + 2; % if false, show size 2. if true, show size 3
		end
	end
	methods (Access = public)
		% This assesses the validity of the memo, i.e. determines if it
		% could have a solution. Throws an error if no outputs are
		% extracted. strictness (optional) can be either 'any' (default) to
		% return whether the memo *may* have a solution, or 'guaranteed' to
		% return whether the memo definitely has a solution.
		function [success,errIDs,errMessages] = hasSolution(this,strictness)
			
			if ~exist('strictness','var')
				strictness = 'any';
			end
			
			% Defer to the superclass's definition of parity validity
			[success,errIDs,errMessages,isSpecial,~] = this.hasSolution_parity();
			% We don't know the actual dimensions, so we don't care about
			% the special cases.
			switch lower(strictness)
				case 'any'
					% leave success as is
				case {'guarantee','guaranteed'}
					if isSpecial
						success = false;
						errIDs{end+1,1} = 'MemoParity:solnNotCertain';
						errMessages{end+1,1} = 'This memo falls into a special case which (without more information) does have a guaranteed solution.';
					end
				otherwise
					error('Unexpected value for "strictness" argument. Use either ''any'' or ''guaranteed''.')
			end
			
			% If no outputs are used and there's an issue, throw a formal
			% error.
			if nargout == 0 && ~success
				% We'll just throw the first error
				error( errIDs{1}, errMessages{1} );
			end
			
		end
		
		% This checks the symmetry of the proposed splitting solution. This
		% determines whether the path is symmetric and whether the parity
		% of the splitting is symmetric. Overall diagonal transits are
		% symmetric if they are invariant to 180 deg rotation and path
		% reversal. Overall left or right transits are symmetric if they
		% are invariant to an appropriate up-down/left-right flip and path
		% reversal.
		function [isSymPath,isSymParity] = checkSymmetry(thisVec,orderArray)
						
			% Validate the inputs and consistency of concatenating based
			% according to orderArray. Extract the concatenated sizes.
			[subParitiesRow,subParitiesCol] = checkAndExtractSizes(thisVec,orderArray,'parity');
			
			% The start and stop corners are specified by the first and
			% last memos present.
			startCorner = thisVec( 1 ).start;
			stopCorner  = thisVec(end).stop;
			% checkAndExtractSizes() called above confirms that these memos
			% are in their corresponding corners of orderArray.
			
			% Determine what sort of symmetry we should hope for. It breaks
			% down into left/right and diagonal cases.
			switch mod( stopCorner - startCorner, 4)
				case 0
					error('The start and stop positions cannot be the same');
				case 2 % diagonal
					% We'll consider rotating the solution by 180 degrees
					% and reversing the resulting path.
					
					% For the orderArray:
					orderArrayRot = rot90(orderArray,2);
					orderArrayConj = numel(orderArray) - orderArrayRot + 1;
					% For the parities:
					subParitiesRowConj = flip(subParitiesRow);
					subParitiesColConj = flip(subParitiesCol);
					% And for the sub-solutions (only accounts for
					% rotation):
					cornerMapping = [2,3,0,1]; % corresponds to initial corner [0,1,2,3]
				otherwise % left/right
					% We'll consider flipping the solution horizontally or
					% vertically so the start and stop positions just swap
					% places, and then we'll reverse the resulting path.
					
					% If the start and stop positions are both in the top,
					% or both in the bottom, then we will flip left-right.
					% If they transition from top to bottom (or vis versa),
					% then we'll flip up-down.
					topBottomByCorner = [0,0,1,1]; % 0=top,1=bottom. corresponds to corners [0,1,2,3]
					if topBottomByCorner(startCorner+1)==topBottomByCorner(stopCorner+1)
						% Flip left-right
						
						% For the orderArray:
						orderArrayFlip = fliplr(orderArray);
						% For the parities:
						subParitiesRowConj =      subParitiesRow;
						subParitiesColConj = flip(subParitiesCol);
						% And for the sub-solutions (only accounts for
						% rotation):
						cornerMapping = [1,0,3,2]; % corresponds to initial corner [0,1,2,3]
					else
						% Flip up-down
						
						% For the orderArray:
						orderArrayFlip = flipud(orderArray);
						% For the parities:
						subParitiesRowConj = flip(subParitiesRow);
						subParitiesColConj =      subParitiesCol;
						% And for the sub-solutions (only accounts for
						% rotation):
						cornerMapping = [3,2,1,0]; % corresponds to initial corner [0,1,2,3]
					end
					orderArrayConj = numel(orderArray) - orderArrayFlip + 1;
			end
			% Both of these states keep rows as rows and columns as
			% columns, so the overall shape is maintained.
			
			% We can now test the equivalence between the flipped/rotated +
			% reversed solution and the original solution.
			% For the path, this requires that the overall order matches
			% its conjugate, and also that the flipped sub-solutions match
			% the reverse order of the original sub-solutions.
			isSymPath = isequal(orderArray,orderArrayConj);
			if isSymPath
				starts = [thisVec.start];
				stops  = [thisVec.stop];
				mappedStarts = cornerMapping(starts+1); % "mapped" = "rotated" or "flipped"
				mappedStops  = cornerMapping(stops +1);
				% Reversing the path involves swapping the start and stop
				% corners, and flipping the order of sub-solutions.
				isSymPath = ...
					all(starts == flip(mappedStops )) && ...
					all(stops  == flip(mappedStarts));
			end
			% For the parity symmetry, we only need to compare the
			% sub-parities before and after the transformation.
			isSymParity = ...
				isequal(subParitiesRow,subParitiesRowConj) && ...
				isequal(subParitiesCol,subParitiesColConj);
			
		end
	end
	
end