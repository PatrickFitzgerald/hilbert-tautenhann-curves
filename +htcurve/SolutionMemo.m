classdef SolutionMemo < htcurve.MemoBase % not handle
	
	properties (Access = public)
		height (1,1) uint32 = 0; % > 0 is valid
		width  (1,1) uint32 = 0; % > 0 is valid
	end
	
	methods (Access = public)
		% This assesses the validity of the memo, i.e. determines if it has
		% a solution. Throws an error if no outputs are extracted.
		function [success,errIDs,errMessages] = hasSolution(this)
			
			% Defer to the superclass's definition of parity validity
			[success,errIDs,errMessages,isSpecial,specialCase] = this.hasSolution_parity();
			% If we're not in a special edge case, then parity arguments
			% are sufficient.
			
			if this.height == 0 || this.width == 0
				success = false;
				errIDs{end+1,1} = 'SolutionMemo:zeroSize';
				errMessages{end+1,1} = 'Invalid memo: the height and/or width is zero';
			end
			
			% If we *are* in a special case, we need to check a bit more
			if isSpecial
				switch specialCase % this label is what can't be 1 (unless the other dimension is also 1)
					case 'H'
						cantBe1 = this.height;
						unlessAlso1 = this.width;
					case 'W'
						cantBe1 = this.width;
						unlessAlso1 = this.height;
					otherwise % including 'unused'
						error('Unexpected specialCase');
				end
				
				if unlessAlso1==1 % see if the "unless" case saves the day
					% No problem
				else
					if cantBe1 == 1
						success = false;
						errIDs{end+1,1} = 'SolutionMemo:specialNoSolution';
						errMessages{end+1,1} = 'This memo falls into a special case which disallows some singleton dimension regions';
					end
				end
			end
			
			% If no outputs are used and there's an issue, throw a formal
			% error.
			if nargout == 0 && ~success
				% We'll just throw the first error
				error( errIDs{1}, errMessages{1} );
			end
			
		end
		% This removes non-substantive degrees of freedom (rotation,
		% flipping) from the memo so that memo comparison is maximally
		% likely to find a match. The `transform` indicates what
		% rotation/flipping was applied to convert to this standardized
		% form.
		function [this_standard,transform] = standardize(this)
			
			% Do some error checking to make sure this is a reasonable
			% object to operate on.
			this.hasSolution(); % no outputs -> throw error on invalid
			
			% STANDARDIZATION PROCEDURE
			%  1. Rotate domain so starting position (start) is top-left (0)
			%  2. Transpose domain so width <= height
			%  3. Tie-breaker: transpose domain so ending position is left or diagonal
			
			% Start from an identity transform
			transform = htcurve.Transform();
			% Copy memo's size
			transform.height = this.height;
			transform.width  = this.width;
			
			% Apply a rotation aimed solely at placing the start at 0.
			transform.rotation = uint8( mod(-this.start,4) ); % negative, mod 4
			
			% Test it out
			this_trial = transform.transformMemo(this);
			% Ensure that worked as desired
			assert( this_trial.start == 0, 'Rotation standardization failed.');
			
			% Determine what flipping (transpose) is necessary
			if this_trial.width < this_trial.height
				% Meets requirement, no tie-breaker necessary
				% leave as is
			elseif this_trial.height < this_trial.width
				% Strictly does not meet requirement, flip; no tie-breaker
				% necessary
				transform.flip = true;
			else % this_trial.width == this_trial.height
				% Tie-breaker needed: if the stop position is in position
				% 3, we'll need to flip to get it into position 1.
				% Otherwise no change needed
				if this_trial.stop == 3
					transform.flip = true;
				end
			end
			
			% Generate the final version, which is now standardized
			this_standard = transform.transformMemo(this);
			% Verify this worked as intended
			assert( this_standard.start == 0, 'Rotation standardization failed.');
			assert( this_standard.width < this_standard.height || any(this_standard.stop==[1,2]),...
				'Flipping standardization failed.')
			
		end
	end
	
	methods (Access = public, Hidden) % not worth the hassle of giving the CurveGenerator access too...
		function parityH = getParityH(this)
			parityH = mod(this.height,2);
		end
		function parityW = getParityW(this)
			parityW = mod(this.width,2);
		end
		function effSizeH = getSizeH(this)
			effSizeH = this.height;
		end
		function effSizeW = getSizeW(this)
			effSizeW = this.width;
		end
	end
	
end