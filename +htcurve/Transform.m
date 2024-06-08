classdef Transform % not handle
	
	properties (Access = public)
		% Input coordinate system
		height   (1,1) uint32 = 0; % > 0
		width    (1,1) uint32 = 0; % > 0
		% How to transform to the output coordinate system
		rotation (1,1) uint8   = 0; % 0,1,2,3 - see SolutionMemo's Direction enums. 1 = 90deg CW rotation
		flip     (1,1) logical = false;
	end
	
	methods (Access = public)
		% This transforms a solution memo instance by the current transform
		% instance
		function memo_rot = transformMemo(this,memo)
			
			% Ensure everything is compatible
			assert( this.height==memo.height && this.width==memo.width,...
				'Solution memo and transform have inconsistent sizes.');
			
			% Prepare the output (representing an identity transform)
			memo_rot = memo;
			
			% Apply each 90deg CW rotation separately
			for rot = 1:this.rotation
				% Swap the width and height
				[memo_rot.width,memo_rot.height] = deal(memo_rot.height,memo_rot.width);
				
				% Rotate the start/stop coordinates
				memo_rot.start = mod( memo_rot.start + 1, 4 );
				memo_rot.stop  = mod( memo_rot.stop  + 1, 4 );
			end
			
			% Apply the flip, if present
			if this.flip
				% Again, flip the sizes
				[memo_rot.width,memo_rot.height] = deal(memo_rot.height,memo_rot.width);
				% Flip the start,stop coordinates over the diagonal
				% 0 and 2 are unchanged. 1,3 swap. This is equivalent to
				% negation mod 4
				memo_rot.start = mod( -memo_rot.start, 4 );
				memo_rot.stop  = mod( -memo_rot.stop , 4 );
			end
			
		end
		% This applies the transform instance to the row/column coordinate,
		% following the transform that would happen to the whole space.
		function [rt,ct] = transformPoint(this,r,c)
			
			% Create temporary height and width that we can modify as we
			% rotate
			height_ = this.height;
			width_  = this.width;
			
			% Accommodate a zero rotation case too
			rt = r;
			ct = c; % will be overwritten if we have any rotation
			
			% Apply each 90deg CW rotation separately
			for rot = 1:this.rotation
				% Swap the live width and height
				[width_,height_] = deal(height_,width_);
				
				% Rotate the coordinates
				rt = c;
				ct = width_ - r + 1;
				
				% Enable compounding rotation
				r = rt;
				c = ct;
			end
			
			% Apply the flip, if present
			if this.flip
				% [width_,height_] = deal(height_,width_);
				[ct,rt] = deal(rt,ct);
			end
			
		end
		% This applies the *inverse* transform on an image
		function imit = invTransformImage(this,im)
			
			% We'll perform the inverse of the flip first (also a flip),
			% and then the inverse of the rotation (also a rotation).
			
			% rot90 already does a CCW rotation, so no negation needed
			if this.flip
				imit = rot90(transpose(im), double(this.rotation));
			else
				imit = rot90(          im , double(this.rotation));
			end
			
		end
	end
	
end