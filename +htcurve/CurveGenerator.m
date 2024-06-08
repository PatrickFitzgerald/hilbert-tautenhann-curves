classdef CurveGenerator < handle
	
	% This is a list of the gathered solutions, which we can build off as
	% we do more processing.
	properties (GetAccess = public, SetAccess = private)
		storedSolutions (:,1) htcurve.Solution = htcurve.Solution.empty(0,1);
	end
	
	methods (Access = public)
		% Run this method to ensure the requested size has been solved.
		% This is valuable to pre-process the curve generation without
		% specifically trying to *use* it for anything.
		function forceSolve(this,H,W)
			
			% We're not going to do anything beyond running the internal
			% force-solve function.
			this.internalForceSolve(H,W);
			
		end
		% This looks up the order of the requested coordinate without
		% needing to evaluate a full image of orders. This should be quite
		% efficient
		function order = lookupPoint(this,H,W,r,c)
			
			% So the request can be delivered, make sure we've solved this
			% case.
			[memo,transform] = this.internalForceSolve(H,W);
			% This will implicitly verify the values of H,W too
			
			% Ensure that the row/column request is reasonable
			assert(isscalar(r) && isnumeric(r) && isfinite(r) && imag(r)==0 && 0<r && r<=H && mod(r,1)==0,...
				'The row value needs to be an integer inside [1,H] (inclusive)');
			assert(isscalar(c) && isnumeric(c) && isfinite(c) && imag(c)==0 && 0<c && c<=W && mod(c,1)==0,...
				'The column value needs to be an integer inside [1,W] (inclusive)');
			
			% We can now use this info to lookup the point. We can be
			% confident everything we need exists.
			order = this.internalLookupPoint(memo,transform,uint32(r),uint32(c));
			
		end
		% This method extracts an image of the pixel orders over the
		% requested region.
% NOT DONE:
		function orderImage = lookupAll(this,H,W)
			
			% So the request can be delivered, make sure we've solved this
			% case.
			[memo,transform] = this.internalForceSolve(H,W);
			
% TODO:
% NOT DONE: propagate request to populate `order`. we may want to offer
% deleting these (if ~solnIsDeferred) in case it needs to much data. If we
% serialize solutions, these should also be discarded, since they can be
% reconstructed.
			
		end
	end
	methods (Access = private)
		% This is used to convert a user's request to our baseline
		% assumption on start and stop conditions
		function [memoStd,transform] = convertUserInput(~,height,width)
			
			% Make a default memo of the desired size (not standardized)
			memo = htcurve.SolutionMemo();
			memo.height = height;
			memo.width  = width;
			memo.start  = 0; % default assumption
			% We'll choose the ending location to be along the longer
			% dimension, if possible
			if height > width
				memo.stop = 3;
			elseif width < height
				memo.stop = 1;
			else % height == width
				% Tie-breaker. Just choose left
				memo.stop = 1;
			end
			% Check if this definition is valid. If not, we'll defer to the
			% diagonal solution, which is valid for all cases where these
			% selections are not.
			[~,errIDs,~] = memo.hasSolution();
			expectedIDs = {'MemoBase:invalidParity','SolutionMemo:specialNoSolution'};
			if any(ismember(expectedIDs,errIDs))
				memo.stop = 2; % diagonal
			end
			% Make sure nothing else went wrong
			if any(~ismember(errIDs,expectedIDs))
				list = sprintf('\t%s\n',errIDs{:});
				error('Something went wrong: caught these errors\n%s\n(can only tolerate %s or %s)',list,expectedIDs{1},expectedIDs{2})
			end
			
			% We should be good now
			[memoStd,transform] = memo.standardize(); % This will throw an error if something is still off
			
		end
		% This overwrites the storedSolutions with the handful of
		% hard-coded base cases.
		function defineBaseSolutions(this)
			
			% Make sure we're starting from a consistent definition
			this.storedSolutions = htcurve.Solution.empty(0,1);
			% Add each hardcoded case in turn
			
			% We'll define everything in its standardized representation.
			% We'll build from a prototype object
			proto = htcurve.Solution();
			proto.memo.height = 0; % to be tailored below
			proto.memo.width  = 0; % to be tailored below
			proto.memo.start  = 0; % standard start
			proto.memo.stop   = 1; % left - to be tailored below
			proto.solnIsDeferred = false; % none of these base solutions are deferred
			proto.order = nan; % not yet populated
			proto.subSolutions = htcurve.SolutionPiece.empty(0,1); % not used for non-deferred solutions
			proto.reconstruction = double.empty(0,0); % not used for non-deferred solutions
			
			% 1x1
			k = 1; % 1x1 - ending left
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 1;
			this.storedSolutions(k,1).memo.width  = 1;
			this.storedSolutions(k,1).memo.stop   = 1; % left
			this.storedSolutions(k,1).order = 0; % 1x1 image of orders
			% Include the other ending cases for this size
			k = k + 1; % 1x1 - ending diagonal
			this.storedSolutions(k,1) = this.storedSolutions(k-1,1); % copy
			this.storedSolutions(k,1).memo.stop = 2; % diagonal
			k = k + 1; % 1x1 - ending right
			this.storedSolutions(k,1) = this.storedSolutions(k-1,1); % copy
			this.storedSolutions(k,1).memo.stop = 3; % right
			
			% 2x1
			k = k + 1; % 2x1 - ending diagonal
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 2;
			this.storedSolutions(k,1).memo.width  = 1;
			this.storedSolutions(k,1).memo.stop   = 2; % diagonal
			this.storedSolutions(k,1).order = [0;1]; % 2x1 image of orders
			k = k + 1; % 2x1 - ending right
			this.storedSolutions(k,1) = this.storedSolutions(k-1,1); % copy
			this.storedSolutions(k,1).memo.stop = 3; % right
			
			% 2x2
			k = k + 1; % 2x2 - ending left
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 2;
			this.storedSolutions(k,1).memo.width  = 2;
			this.storedSolutions(k,1).memo.stop   = 1; % left
			this.storedSolutions(k,1).order = [0,3;1,2]; % 2x2 image of orders
			
			% 3x1
			k = k + 1; % 3x1 - ending diagonal
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 3;
			this.storedSolutions(k,1).memo.width  = 1;
			this.storedSolutions(k,1).memo.stop   = 2; % diagonal
			this.storedSolutions(k,1).order = [0;1;2]; % 3x1 image of orders
			k = k + 1; % 3x1 - ending right
			this.storedSolutions(k,1) = this.storedSolutions(k-1,1); % copy
			this.storedSolutions(k,1).memo.stop = 3; % right
			
			% 3x2
			k = k + 1; % 3x2 - ending left
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 3;
			this.storedSolutions(k,1).memo.width  = 2;
			this.storedSolutions(k,1).memo.stop   = 1; % left
			this.storedSolutions(k,1).order = [0,5;1,4;2,3]; % 3x2 image of orders
			k = k + 1; % 3x2 - ending diagonal
			this.storedSolutions(k,1) = proto; % copy
			this.storedSolutions(k,1).memo.height = 3;
			this.storedSolutions(k,1).memo.width  = 2;
			this.storedSolutions(k,1).memo.stop   = 2; % diagional
			this.storedSolutions(k,1).order = [0,1;3,2;4,5]; % 3x2 image of orders
			
		end
		% This function converts the user specification to our standard
		% internal representation, and then kicks off the recursive solver
		% to ensure the full solution exists (or make that so).
		function [memo,transform] = internalForceSolve(this,height,width)
			
			% Convert the user's inputs to something more compatible with
			% our framework. We'll return them too, since they're useful
			% elsewhere.
			[memo,transform] = this.convertUserInput(height,width);
			
			% If we haven't started anything, make sure the baseline
			% solutions are present, so the algorithm can actually
			% terminate.
			if isempty(this.storedSolutions)
				this.defineBaseSolutions();
			end
			% Now kick off the recursive solver, which will terminate
			% whenever everything has been solved. If everything was
			% already solved, then it will terminate immediately.
			this.recursiveSolver(memo);
			
		end
		% This looks up the requested coordinate based on our internal
		% format. This function assumes that the curve has already been
		% solved. The r,c is interpreted in the parent units, and the
		% transform to the standardized memo's space has not yet been
		% applied.
		function order = internalLookupPoint(this,memo,transform,r,c)
			
			% Transform the query point to the standardized memo's
			% coordinates.
			[rs,cs] = transform.transformPoint(r,c); % s = standardized
			
			% Find the relevant item for this memo
			[~,solnIndex] = this.lookupMemo(memo); % assume we're good.
			% Get a copy of this object
			relevantSoln = this.storedSolutions(solnIndex);
			
			% Determine how to proceed
			if relevantSoln.solnIsDeferred
				% There's at least 2 sub-solutions. We need to determine
				% which sub-solution is relevant.
				numSub = double(relevantSoln.numSubSolns);
				for k = 1:numSub
					window = relevantSoln.subSolutions(k).subRegionWindow; % [first row, last row, first col, last col] (inclusive)
					if window(1) <= rs && rs <= window(2) && window(3) <= cs && cs <= window(4)
						% We fall inside this window. Convert the
						% row/column to be relative to this smaller window.
						rst = (rs + 1) - window(1); % t = truncated
						cst = (cs + 1) - window(3);
						% Recur!
						subOrder = this.internalLookupPoint(...
							relevantSoln.subSolutions(k).subProblemMemo,... % use sub-problem's info on linking
							relevantSoln.subSolutions(k).transformToSubProblem,... % to next lower level
							rst,cst);
						% The order returned therein will be local to that
						% smaller scope, so we'll re-contextualize it by
						% adding the appropriate offset.
						order = subOrder + relevantSoln.subSolutions(k).orderOffset;
						return
					end
				end
				% If we're still running, something went wrong
				error('Point lookup process entered an invalid state...');
			else % not deferred, actual solution right here.
				% Just index into the matrix of orders
				order = relevantSoln.order(rs,cs);
			end
			
		end
		% This finds a match of the memo (assumed already standardized) to
		% the list of available solutions already generated. wasFound
		% (boolean) indicates if the solution was found, and if it was
		% found, then solnIndex will be meaningfully populated with that
		% solution's index in the list. Note that other operations may
		% reorder the solution list, so this index may become invalid if
		% not used immediately.
		function [wasFound,solnIndex] = lookupMemo(this,memo)
			
			% Default results
			wasFound = false;
			solnIndex = nan;
			% Search for a matching memo
			for k = 1:numel(this.storedSolutions)
				if isequal( this.storedSolutions(k).memo, memo )
					wasFound = true;
					solnIndex = k;
					return
				end
			end
% binary search? periodically reorder data for efficiency?
% TODO: sort storedSolutions so it can be searched more efficiently?
			
		end
		% This chooses how to break the problem into smaller pieces, and
		% recursively solves those subproblems. memo should be a
		% standardized memo.
		function recursiveSolver(this,memo)
			
			% Determine if this has already been solved
			if this.lookupMemo(memo) % wasFound is condition
				% Terminate if solved
				return
			end
			% If we're still running, then we need to solve it.
			
			
			% For the time being, I'll implement something quite close to
			% Tautenhahn's approach
			
			
			
			
			
			
		end
	end
	
end