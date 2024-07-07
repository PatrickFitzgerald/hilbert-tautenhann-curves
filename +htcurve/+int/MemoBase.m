classdef (Abstract) MemoBase % not handle
	
	properties (Access = public)
		start (1,1) int8 = 0; % See direction enums
		stop  (1,1) int8 = 0; % ^^^^
	end
	% Direction enums:
	%    0 - (1,1) - top left
	%    1 - (1,W) - top right
	%    2 - (H,W) - bottom right
	%    3 - (H,1) - bottom left
	% Arithemtic mod 4 is isometric to rotations (+1 = 90deg CW)
	
	
	methods (Access = public)
		function plot(thisVec,orderArray)
			
			% Validate inputs and extract sizes along rows/columns
			[sizesRow,sizesCol] = checkAndExtractSizes(thisVec,orderArray,'effSize');
			
			% If these are empty, we have nothing to plot.
			if isempty(sizesRow) || isempty(sizesCol)
				return
			end
			quant = numel(orderArray);
			
			% Repmat it to directly correspond to orderArray.
			sizesHmat = repmat(sizesRow,1,numel(sizesCol));
			sizesWmat = repmat(sizesCol,numel(sizesRow),1);
			
			% Extract the start/stop corner
			startCorner = arrayfun( @(t) t.start, thisVec(:));
			stopCorner  = arrayfun( @(t) t.stop,  thisVec(:));
			
			% Plan what to plot
			offsetH = zeros(size(sizesHmat));
			offsetW = zeros(size(sizesWmat));
			offsetH(2:end,:) = cumsum(sizesHmat(1:end-1,:),1);
			offsetW(:,2:end) = cumsum(sizesWmat(:,1:end-1),2);
			% Plan how we'll enter and exist each region
			isRight = [0,1,1,0]; % corner 0,1,2,3
			isDown  = [0,0,1,1]; % corner 0,1,2,3
			R = @(x) reshape(x,size(offsetH));
			startCoordR = offsetH + R(isDown (startCorner(orderArray)+1)).*(sizesHmat-1) + 1;
			stopCoordR  = offsetH + R(isDown ( stopCorner(orderArray)+1)).*(sizesHmat-1) + 1;
			startCoordC = offsetW + R(isRight(startCorner(orderArray)+1)).*(sizesWmat-1) + 1;
			stopCoordC  = offsetW + R(isRight( stopCorner(orderArray)+1)).*(sizesWmat-1) + 1;
			% Sort these based on order. orderArray is the most dominant
			% effect, but start comes before stop
			coordList = [...
				orderArray(:),ones(quant,1)*1,startCoordR(:),startCoordC(:);...
				orderArray(:),ones(quant,1)*2, stopCoordR(:), stopCoordC(:);...
			];
			coordList = sortrows(coordList,[1,2]);
			coordListRC = coordList(:,3:4);
			
			% Plot
			hold(gca,'on');
			eps_ = 0.125;
			for k = 1:quant
				rectangle(...
					'Position',[ offsetW(k)+0.5+eps_, offsetH(k)+0.5+eps_, sizesWmat(k)-2*eps_, sizesHmat(k)-2*eps_ ],...
					'FaceColor',[1,1,1] * 0.8)
			end
			plot(coordListRC(:,2),coordListRC(:,1),'k'); % main path
			plot(coordListRC(1:2:end,2),coordListRC(1:2:end,1),'og'); % green empty circles for starts
			plot(coordListRC(2:2:end,2),coordListRC(2:2:end,1),'sr'); % red empty squares for stops
			plot(coordListRC(1,2),coordListRC(1,1),'og','MarkerFaceColor','g'); % green filled circle for first start
			plot(coordListRC(end,2),coordListRC(end,1),'sr','MarkerFaceColor','r'); % red filled square for last stop
			set(gca,'YDir','reverse');
			daspect(gca,[1,1,1]);
			box(gca,'on')
			s1 = sum(sizesHmat,1); s1 = s1(1);
			s2 = sum(sizesWmat,2); s2 = s2(1);
			set(gca,'XTick',1:s2,'YTick',1:s1,'XTickLabels','','YTickLabels','');
			grid(gca,'on')
			axis([0.5-eps_,s2+0.5+eps_,0.5-eps_,s1+0.5+eps_]);
			
		end
	end
	
	methods (Access = public, Hidden) % not worth the hassle of giving the CurveGenerator access too...
		% These functions extract the parity (0 for even, 1 for odd) of the
		% corresponding height (H) and width (W). Values returned should be
		% logical scalars.
		parityH = getParityH(this);
		parityW = getParityW(this);
		% These functions extract the effective size of the region they
		% represent. These are only used for plotting and the
		% checkAndExtractSizes method.
		effSizeH = getSizeH(this);
		effSizeW = getSizeW(this);
	end
	methods (Access = public, Abstract)
		% This assesses the validity of the memo, i.e. determines if it
		% could have a solution. Throws an error if no outputs are
		% extracted.
		[success,errIDs,errMessages] = hasSolution(this);
	end
	
	methods (Access = protected)
		% This is a helper function for subclasses to determine whether
		% this is precluded from having solutions on grounds of parity.
		% Will not intentionall throw error messages. isSpecial is true
		% when both dimensions have odd parity; if true, specialCase
		% indicates which dimension can't be length 1 (unless the other
		% dimension is also 1).
		function [success,errIDs,errMessages,isSpecial,specialCase] = hasSolution_parity(this)
			
			% Default to success
			success = true;
			errIDs = cell(0,1);
			errMessages = cell(0,1);
			
			% Do some basic error checking
			if mod(this.start - this.stop,4) == 0 % start==stop
				success = false;
				errIDs{end+1,1} = 'MemoBase:startEqStop';
				errMessages{end+1,1} = 'Invalid memo: start/stop positions cannot be the same';
			end
			if this.start~=mod(this.start,4) || this.stop~=mod(this.stop,4) % invalid start/stop values
				success = false;
				errIDs{end+1,1} = 'MemoBase:startStopOutOfBounds';
				errMessages{end+1,1} = 'Invalid memo: start/stop positions should be one of 0,1,2,3';
			end
			
			% Extract the parity of this memo
			parityH = this.getParityH();
			parityW = this.getParityW();
			
			
			% Default to not special
			isSpecial = false;
			specialCase = 'unused';
			
			
			% Check parity
			parityViolated = false;
			% We'll do some accounting so the start is effectively 0
			paritiesMN = [parityH,parityW]; % start from M=H, N=W
			dimLabelMN = 'HW'; % start from M=H, N=W
			if mod(this.start,2) == 1 % each 90deg rotation needed requires flipping the M,N associations
				paritiesMN = flip(paritiesMN);
				dimLabelMN = flip(dimLabelMN);
			end
			direction = mod( this.stop - this.start, 4);
			% Determine whether the memo's parity violates our constraints
			parityM = paritiesMN(1); % for readability
			parityN = paritiesMN(2); % ^^^^
			switch direction
				case 0 % equal, error will be thrown for previous reason
					% do nothing
				case 1 % left turn
					parityViolated = parityM==0 && parityN==1;
					if parityM==1 && parityN==1 % "*L" condition
						isSpecial = true;
						specialCase = dimLabelMN(2); % N can't be 1 (unless M is too)
					end
				case 2 % diagonal
					parityViolated = parityM==0 && parityN==0;
				case 3 % right turn
					parityViolated = parityM==1 && parityN==0;
					if parityM==1 && parityN==1 % "*R" condition
						isSpecial = true;
						specialCase = dimLabelMN(1); % M can't be 1 (unless N is too)
					end
			end
			if parityViolated
				success = false;
				errIDs{end+1,1} = 'MemoBase:invalidParity';
				errMessages{end+1,1} = 'Invalid memo: its parity indicates it is not solvable';
			end
			
		end
		
		% This function validates the array of memos is consistent with the
		% orderArray specification. This checks quantities are
		% interconsistent, and that orderArray is valid. Also confirms that
		% the size of each memo is consistent with how it is effectively
		% getting concatenated by orderArray. This function returns the
		% sizes along the rows/columns of that concatenated representation.
		%    sizeFlag = 'parity' or 'effSize' (default)
		function [sizesRow,sizesCol] = checkAndExtractSizes(thisVec,orderArray,sizeFlag)
			
			% Check the number of elements of thisVec
			quant = numel(orderArray);
			assert( numel(thisVec)==quant, 'Number of elements in the Memo array and in the orderArray must match');
			if quant == 0
				sizesRow = double.empty(0,1);
				sizesCol = double.empty(1,0);
				return
			end
			assert( min(orderArray(:))==1 && max(orderArray(:))==quant && all(diff(sort(orderArray(:)))==1),...
				'Invalid orderArray matrix.');
			
			
			% Check the orderArray specifies a contiguous order
			orderListRC = htcurve.int.MemoBase.orderArrayToList(orderArray);
			absDeltas = abs(diff(orderListRC,[],1));
			assert( all(all(sort(absDeltas,2) == [0,1])), ... all deltas should change one dimension by +/-1, and the other dimension should be unchanged
				'The orderArray specified does not specify a contiguous path.');
			
			
			% We can also check that the orderArray starts and stops in
			% reasonable corner. Make sure those respective memos are
			% consistent (i.e. if the big picture says it starts in corner
			% 3, then the first memo better start in its local corner 3
			% too)
			cornerOrders = orderArray([1,end],[1,end]);
			cornerOrders = cornerOrders([1,3,4,2]); % reorder to correspond to "Direction enums" above
			assert(...
				cornerOrders(thisVec(  1  ).start+1) == 1     && ...
				cornerOrders(thisVec(quant).stop +1) == quant,...
				'The start/stop corner of the first/last memo is not consistent with the start/stop corner of the ocncatenated specification.' );
			
			
			% Handle the optional behavior of sizeFlag
			if ~exist('sizeFlag','var') || isequal(sizeFlag,[])
				sizeFlag = 'effSize';
			end
			assert( ischar(sizeFlag)||isstring(sizeFlag), 'The sizeFlag needs to be a char/string' );
			sizeFlag = lower( char(sizeFlag) );
			switch sizeFlag
				case 'parity'
					funcH = @getParityH;
					funcW = @getParityW;
				case 'effsize'
					funcH = @getSizeH;
					funcW = @getSizeW;
				otherwise
					error('The sizeFlag must either be ''parity'' or ''effSize''.');
			end
			
			
			% Extract the effective sizes of each item
			sizesH = arrayfun( funcH, thisVec(:) ); % column vector
			sizesW = arrayfun( funcW, thisVec(:) ); % ^^^
			
			% Arrange these according to orderArray
			sizesHmat = reshape( sizesH(orderArray), size(orderArray)); % matrix indexing can behave unreliably if one dimension is singleton
			sizesWmat = reshape( sizesW(orderArray), size(orderArray)); % ^^^
			% Extract a representative set of the sizes along rows/columns
			sizesRow = sizesHmat(:,1);
			sizesCol = sizesWmat(1,:);
			% Confirm they are consistent with rectangularized plotting
			assert( all(all(sizesRow == sizesHmat,2)) && ... 
					all(all(sizesCol == sizesWmat,1)), 'Sizes incompatible for naive rectangular concatenation');
		end
	end
	methods (Access = public, Static)
		% These functions help convert between different ways of defining
		% an order through some rectangular region
		% 
		% e.g.
		%     orderArray:                 1,1
		%         1 4                     2,1
		%         2 3       orderListRC:  2,2
		%                                 1,2
		% 
		function orderListRC = orderArrayToList(orderArray)
			sz = size(orderArray);
			[R,C] = ndgrid(1:sz(1),1:sz(2));
			linearlizedOrder = [ orderArray(:), R(:), C(:) ]; % all linearize in the same way, so association maintained
			linearlizedOrder = sortrows(linearlizedOrder,1); % no tie breakers needed
			orderListRC = linearlizedOrder(:,2:3); % only need the row,column ordering. order is now implied
		end
		function orderArray = orderListToArray(orderListRC)
			orderArray = accumarray(...
				orderListRC,...
				(1:size(orderListRC,1)).' ...
			);
		end
	end
	
end