
% These splittingCases assume we're starting at start=0
sol_12 = [1,2];
sol_13 = [1,2,3];
sol_22 = [1,2;4,3];
sol_23a = [1,4,5;2,3,6];
sol_23b = [1,2,3;6,5,4];
sol_33a = [1,2,9;4,3,8;5,6,7];
sol_33b = [1,2,3;6,5,4;7,8,9];
splittingCases = {
%   #H,#W, ordering (matrix), ordering list [POPULATED AFTER] 
	1, 2, {sol_12}; % half
	2, 1, {sol_12.'}; % half
	1, 3, {sol_13}; % thirds
	3, 1, {sol_13.'}; % thirds
	2, 2, {sol_22, sol_22.'}; % half twice
	2, 3, {sol_23a, sol_23b}; % half-thirds
	3, 2, {sol_23a.',sol_23b.'}; % half-thirds
	3, 3, {sol_33a,sol_33a.',sol_33b,sol_33b.'};% thirds twice
};
numSplitCases = size(splittingCases,1);

%% Populate the ordering list (column 4)
for splitInd = 1:numSplitCases
	splittingCases{splitInd,4} = {}; % prepare where we're storing these items to also be a cell array
	for variantInd = 1:numel(splittingCases{splitInd,3})
		% Convert each order array to an order list [r1,c1; r2,c2; ...]
		splittingCases{splitInd,4}{variantInd} = ...
			htcurve.int.MemoBase.orderArrayToList(...
				splittingCases{splitInd,3}{variantInd}...
			);
	end
end

%%
% Ensure these were defined correctly (error catching and visual
% inspection)
close all
for splitInd = 1:numSplitCases
	generateAllSplits( splittingCases(splitInd,:) );
end

%%
close all

%% Group the splitting cases into different classes of splits
% Classes of splitting
%  1. Splitting along just the vertical dimension
%  2. Splitting along just the horizontal dimension
%  3. Splitting along both dimensions simultaneously
splitsAlongVert = cellfun( @(s)s~=1, splittingCases(:,1));
splitsAlongHorz = cellfun( @(s)s~=1, splittingCases(:,2));
% Gather indices for these different splitting cases, so we may loop over
% them together later
splittingClasses = cell(3,1);
splittingClasses{1} = find( splitsAlongVert & ~splitsAlongHorz);
splittingClasses{2} = find(~splitsAlongVert &  splitsAlongHorz);
splittingClasses{3} = find( splitsAlongVert &  splitsAlongHorz);

%% Downselect solutions

% We'll prepare a place to store all our solutions
%   dim 1: isOddH:     {false,true}
%   dim 2: isOddW:     {false,true}
%   dim 3: stop:       {1,2,3} conditional on starting on 0
%   dim 4: splitClass: {1,2,3}
splittings = cell(2,2,3,3);
% This will have empty contents wherever that case is impossible.
uiwait(msgbox('When each figure pops up, uncheck any splittings which are redundant to other checked splittings. Close the figure when you''re happy with your selection.','Instructions'));

original = htcurve.int.MemoParity();
original.start = 0; % fixed
for isOddH = [false,true]
	for isOddW = [false,true]
		for stop = 1:3
			original.isOddH = isOddH;
			original.isOddW = isOddW;
			original.stop = stop;
			
			% If the original is impossible, skip trying to split it
			if ~original.hasSolution()
				continue
			end
			
			for splitClass = 1:3
				
				subOrderings = cell(0,2);
				for splitInd = splittingClasses{splitClass}.'
					tempOrderings = generateAllSplits( splittingCases(splitInd,:), original );
					
					% Accumulate these solutions
					subOrderings = [ subOrderings; tempOrderings ]; %#ok<AGROW>
				end
				
				% If this split class is empty (e.g. trying to split in an
				% invalid direction), no need to prune.
				if isempty(subOrderings)
					continue
				end
				
				% Simplify these solutions down as much as possible. We
				% don't need many representatives from each size/splitting
				% class.
				[subOrderings,goodnessMetrics] = pruneSolutions(subOrderings,original);
				
				numSol = size(subOrderings,1);
				switch numSol
					case 0 % no solutions found, even though we solutions provided...
						error('Solution pruning failed');
					case 1 % exactly one. nothing fancy to do
						% no op
					otherwise % 2 or more solutions
						% Let the user downselect to remove symmetry, if
						% necessary
						
						% Plan how to format the figure
						numRows = floor(sqrt(numSol));
						numCols = ceil(numSol ./ numRows);
						f = figure('MenuBar','none','ToolBar','none',...
							'CloseRequestFcn',@recordAndClose);
						checkBoxes = gobjects(numSol,1);
						for k = 1:numSol
							ax = subplot(numRows,numCols,k);
							plot(subOrderings{k,2},subOrderings{k,1});
							checkBoxes(k) = uicontrol('Parent',f,'Style','checkbox',...
								'Value',true,'String','Keep?',...
								'Units','normalized'...
							);
							checkBoxes(k).Position(1:2) = ax.Position(1:2);
						end
						f.UserData = checkBoxes; % unambiguous, self-managed copy
						
						% Wait for the figure to close
						uiwait(f);
						% Retrieve the last copy of the checkbox settings
						selections = recordAndClose('retrieve');
						
						% Downselect and proceed
						subOrderings    = subOrderings(    logical(selections), :);
						goodnessMetrics = goodnessMetrics( logical(selections), :);
				end
				numSol = size(subOrderings,1);
				
				% Store answers and move to the next case
				splittings{isOddH+1,isOddW+1,stop,splitClass} = [ subOrderings, mat2cell(goodnessMetrics,ones(numSol,1),size(goodnessMetrics,2)) ];
				
			end
		end
	end
end

%% Package and save splittings

% We'll have the main set of solutions contain less dependence on cells,
% just to remove the overhead there. That necessitates some mapping from
% problem scenario to a set of solutions in a global list.
e = cell(1,0);
linearized = struct(...
	'numH',e,...          1 x 1
	'numW',e,...          1 x 1
	'isOddH',e,...        numH x 1
	'isOddW',e,...        1 x numW
	'subOrderings',e,...  numH x numW
	'subStarts',e,...     numH x numW
	'subStops',e ...      numH x numW
);
lookupOffset = nan(2,2,3); % Collapse on splitting class dimension
lookupCount  = nan(2,2,3);

for isOddH = [false,true]
	for isOddW = [false,true]
		for stop = 1:3
			% We cared about split class when we were pruning solutions,
			% but when we're using these splittings to prepare
			% sub-solutions, they are equivalent. We'll bundle them
			% together.
			splitClass = 1:3;
			
			origProbInds = arrayfun(@(sc) sub2ind([2,2,3,3],isOddH+1,isOddW+1,stop,sc), splitClass);
			testProbInds = sub2ind([2,2,3],  isOddH+1,isOddW+1,stop);
			
			% Group the solutions from separate splitting classes
			splittings_ = cat(1,splittings{origProbInds});
			
			% If this scenario has no solutions, then it was an invalid
			% problem, so we will skip formatting it. This will leave
			% the lookup terms nan, which will make downstream errors
			% more obvious.
			numSol = size(splittings_,1);
			if numSol == 0
				continue
			end
			
			% If the problem was valid, record where to find the
			% solutions.
			o = numel(linearized); % before we've appended. o = offset
			lookupOffset(testProbInds) = o;
			lookupCount (testProbInds) = numSol;
			
			% Add those solutions, where the solutions with the greatest
			% goodness appear in the list first
			[~,order] = sortrows( cell2mat(splittings_(:,3)), 'descend' );
			for k = 1:numSol % output index
				
				temp1 = splittings_{order(k),1};
				temp2 = splittings_{order(k),2};
				
				% Store size
				linearized(o+k).numH = size(temp1,1);
				linearized(o+k).numW = size(temp1,2);
				
				
				% Extract a dimensionally projected version of the
				% parity
				
				reshapedParityMemos = temp2(temp1);
				% Matlab has some discrepancies for the orientation of
				% arrays indexed with N (N>1) dimensional arrays vs 1D
				% arrays...
				if ~isequal(size(reshapedParityMemos),size(temp1))
					reshapedParityMemos = reshape(reshapedParityMemos,size(temp1));
				end
				
				% I will defer error catching on inconsistency upon
				% projecting onto each dimension... This should have
				% been elsewhere
				linearized(o+k).isOddH = arrayfun(@(pm)pm.isOddH, reshapedParityMemos(:,1));
				linearized(o+k).isOddW = arrayfun(@(pm)pm.isOddW, reshapedParityMemos(1,:));
				
				
				% Store the sub-solution ordering
				linearized(o+k).subOrderings = temp1;
				
				
				% At this point, the only uncaptured info in the solved
				% MemoParity objects is the start/stop position. Gather
				% those.
				linearized(o+k).subStarts = arrayfun(@(pm)pm.start, reshapedParityMemos);
				linearized(o+k).subStops  = arrayfun(@(pm)pm.stop,  reshapedParityMemos);
				
			end
			
		end
	end
end

% Save this to file, to be accessed by the solver
save( fullfile( htcurve.getPackagePath(),'+int','splittings.mat'),...
	'linearized','lookupOffset','lookupCount');

%% Support functions

% No output arguments will only plot the splittingCases. The subOrderings
% will be a cell array. The original will be a MemoParity object,
% indicating what we're trying to split
function subOrderings = generateAllSplits(splittingCase,original)
	
	s1 = splittingCase{1};
	s2 = splittingCase{2};
	
	% We'll construct this list of objects once to reduce overhead. It is a
	% handle object so, so only one copy of it will exist. Given that as we
	% recur later on, prior selections will remain intact until we're done
	% testing, we don't need to worry about reusing the memory.
	list = htcurve.int.MemoParityList(s1*s2);
	
	% We'll gather valid suborderings as we find them
	subOrderings = cell(0,2);
	
	
	subSolns = splittingCase{4};
	for ssInd = 1:numel(subSolns)
		orderList = subSolns{ssInd};
		
		% Error checking and visualization
		orderArray = htcurve.int.MemoBase.orderListToArray( orderList ); % reconstruct, for error checking and visualization
		assert( isequal(size(orderArray),[s1,s2]), 'Size mismatch' );
		if nargout == 0
			figure;
			imagesc(orderArray);
			daspect([1,1,1]);
			colormap gray
			title(sprintf('%u,%u - solution %u',s1,s2,ssInd));
			set(gca,'XTick',[],'YTick',[])
			continue
		end
		
		% Determine whether this sub-solution ends at the correct place.
		assert( original.start==0, 'Start position not 0...') % assumed format. needlessly complicated to handle general case
		switch original.stop
			case 1 % left (top right corner)
				lastExpected = [1,s2];
			case 2 % diagonal (bottom right corner)
				lastExpected = [s1,s2];
			case 3 % right (bottom left corner)
				lastExpected = [s1,1];
		end
		% If the last coordinate in the order list is not compatible with
		% where we need to end, skip it.
		if ~isequal(lastExpected,orderList(end,:))
			continue
		end
		
		% Generate a parity breakdown for the rows/cols, and determine if
		% it's compatible with the original. Discard if not
		paritiesRows = getListOfParities(s1);
		paritiesCols = getListOfParities(s2);
		paritiesRows = paritiesRows( mod(sum(paritiesRows,2),2) == original.isOddH, : );
		paritiesCols = paritiesCols( mod(sum(paritiesCols,2),2) == original.isOddW, : );
		for pIndRows = 1:size(paritiesRows,1)
			for pIndCols = 1:size(paritiesCols,1)
				% All combinations here are valid
				
				% Assign these parities
				for k = 1:s1*s2
					r = orderList(k,1);
					c = orderList(k,2);
					list.entries(k).isOddH = paritiesRows(pIndRows,r);
					list.entries(k).isOddW = paritiesCols(pIndCols,c);
				end
				
				% Apply the boundary conditions
				list.entries( 1 ).start = original.start; % should be 0, but meh
				list.entries(end).stop  = original.stop;
				
				% Generate all solutions 
				[foundSolutions,solutionList] = depthFirstMakeSolutions(2,list,orderList); % start processing at depth=2
				if foundSolutions
					for k = 1:numel(solutionList)
						subOrderings(end+1,:) = {orderArray,solutionList{k}}; %#ok<AGROW>
					end
				end
				
			end
		end
	end
end
function listOfParities = getListOfParities(numVars)
	if numVars == 1
		listOfParities = [false;true];
	else
		subList = getListOfParities(numVars-1);
		sz = size(subList,1);
		listOfParities = [ false(sz,1), subList; true(sz,1), subList ];
	end
end
% In each scope of this function, we guarantee
%  * the first "depth" "start" values are defined and will not be modified
%    *internally*,
%  * the first "depth"-1 "stop" values are defined and will not be modified
%    *internally*, and
%  * the last "stop" value is defined and will not be modified
% These *internal* remarks paired with serial operation guarantee we aren't
% messing things up when we don't intend to
function [foundSolutions,solutionList] = depthFirstMakeSolutions(depth,list,orderList)
	
	foundSolutions = false;
	solutionList = cell(0,1);
	
	isTerminal = depth == numel(list.entries);
	
	% Determine the direction we were moving to reach our current
	% position. For each case, determine which corners are valid to
	% perform that action. For example, going to the right, we could
	% pass through corner 1 and pop out at corner 0, or we could pass
	% through corner 2 and pop out at corner 3.
	delta = diff(orderList(depth+[-1,0],:),[],1); % [row delta, col delta]
	if isequal(delta,[0,1]) % right
		allowed = int8([ 1,0; 2,3 ]); % each row is a separate option
	elseif isequal(delta,[1,0]) % down
		allowed = int8([ 3,0; 2,1 ]);
	elseif isequal(delta,[0,-1]) % left
		allowed = int8([ 0,1; 3,2 ]);
	elseif isequal(delta,[-1,0]) % up
		allowed = int8([ 0,3; 1,2 ]);
	else % should not happen
		error('Something went wrong: a path was not unit-distance contiguous');
	end
	
	% Try out each of these cases
	for choice = 1:2
		% Apply the choice
		list.entries(depth-1).stop  = allowed(choice,1);
		list.entries(depth  ).start = allowed(choice,2);
		
		% Determine whether this choice is already invalid.
		if ~list.entries(depth-1).hasSolution()
			continue
		end
		% If we're still running, things appear valid
		
		% If we're at the end of the path, we can also check the last item
		% and confirm that one's good too
		if isTerminal
			if ~list.entries(end).hasSolution()
				continue
			end
			% If it's still good, then we can gather the answer from the
			% current state of the list
			foundSolutions = true;
			solutionList{end+1,1} = list.entries; %#ok<AGROW> % makes a copy since the "entries" is an array of value-type variables (not handle type)
		else % not terminal
			% We'll recur and see what remaining choices lead to valid
			% solutions
			[fS,sL] = depthFirstMakeSolutions(depth+1,list,orderList); % almost all of the incremental knowledge is stored in "list"
			if fS % only bother appending solutions if they exist
				foundSolutions = true;
				solutionList(end+(1:numel(sL)),1) = sL;
			end
		end
	end
end

function [subOrderingsBest,goodnessBest] = pruneSolutions(subOrderings, originalMemo)
	
	% Determine which of these sub-orderings have guaranteed solutions
	% (i.e. avoid column/row sub-sizes)
	isGuaranteed = cellfun( @(mpVec) all(arrayfun( @(mp) mp.hasSolution('guaranteed'), mpVec )), subOrderings(:,2) );
	
	% Check the symmetry of these solutions
	[isSymPath,isSymParity] = cellfun( @checkSymmetry, subOrderings(:,2), subOrderings(:,1) );
	
	% Check how uniformly different directions are sub-transit directions
	% are used in the full split path. If only one type is used, this
	% "variance" will return a 0. If there's a mix, it will be nonzero.
	% Also measure consistency with the parent memo's direction (diagonal
	% or not diagonal).
	[directionVariances,directionConsistencies] = cellfun(@(subMemos) measureDirectionMetrics(subMemos,originalMemo), subOrderings(:,2) );
	
	% Measure how uniformly the sub-solution parities are. Also measure the
	% portion of consistent parities with the parent parity.
	[parityVariances,parityConsistencies,paritySplittingQuality] = cellfun(@(subMemos,orderArray) measureParityMetrics(subMemos,originalMemo,orderArray), subOrderings(:,2), subOrderings(:,1) );
	
	% Summarize these qualities (weighting them differently). Here's the
	% priority: (#1 is most important)
	%   1. symmetric path                          MOST OBVIOUS "SYMMETRY"
	%   2. low direction variance                  INTER-CONSISTENCY
	%   3. low parity variance                     INTER-CONSISTENCY
	%   4. high direction consistency with parent  RECURSION CONSISTENCY
	%   5. high parity consistency with parent     RECURSION CONSISTENCY
	%   6. split size and parity consistency       RECURSION CONSISTENCY
	%   7. symmetric splitting parity              NICE TO HAVE
	%   8. solution guaranteed                     NICE TO HAVE (until we require it...)
	goodness = [
		isSymPath,...
		-directionVariances,...
		-parityVariances,...
		directionConsistencies,...
		parityConsistencies,...
		paritySplittingQuality,...
		isSymParity,...
		isGuaranteed...
	];
	meaningfulMetrics = 1:7;
	% Define what the best value could be for each column, which will help
	% with debugging
	maxGoodness = [
		1.0,... % isSymPath = true
		0.0,... % -directionVariances = 0.0
		0.0,... % -parityVariances = 0.0
		1.0,... % directionConsistencies = 1.0
		1.0,... % parityConsistencies = 1.0
		1.0,... % paritySplittingQuality = 1.0
		1.0,... % isSymParity = true
		1.0 ... % isGuaranteed = true
	];
	
	[~,order] = sortrows(goodness,'descend'); % put the best ones first
	
	% Gather a list of all solutions which are at the best available
	% goodness.
	selection = find(all(goodness == goodness(order(1),:),2));
	
	% This prior goodness was measured without *requiring* the solution is
	% guaranteed. If none of the selected values were guaranteed, then we
	% need to add some extra solutions to our selection so that
	% collectively there is always a guaranteed solution available.
	if ~any(isGuaranteed(selection))
		% Verify that a solution does exist (even if it's bad)
		assert(any(isGuaranteed), 'There were no solutions which were guaranteed...');
		% Prevent non-guaranteed solutions from being selected.
		tempGoodness = goodness;
		tempGoodness(~isGuaranteed,:) = -inf;
		% Solve in the same way as before
		[~,order] = sortrows(tempGoodness,'descend');
		guaranteedSelection = find(all(tempGoodness == tempGoodness(order(1),:),2));
		% Concatenate the selections (the non-guaranteed solutions go
		% first, since they are preferred.
		selection = [selection; guaranteedSelection ];
		% Since the old selection and the new selection are perfecly
		% partitioned by their isGuaranteed status, there's no possibility
		% of selections being redundant between the two parts.
	end
% 	close all
% 	for k = selection.'
% 		figure;
% 		plot(subOrderings{k,2},subOrderings{k,1});
% 		title(sprintf('is guaranteed: %u',isGuaranteed(k)));
% 	end
	
	% This is more for debugging
	if ~any(all(goodness==maxGoodness,2),1)
		disp('nothing perfect')
	end
	
	subOrderingsBest = subOrderings(selection,:);
	goodnessBest = goodness(selection,meaningfulMetrics);
	
end
% This measures the variance of diagonal/non-diagonal sub-crossings (zero
% if they're all the same, biggest if they're evenly mixed) and the
% consistency of diagonal/non-diagonal sub-crossings (relative to parent
% crossing, 1.0 for perfectly consistent with parent, lower if not).
function [directionVariance,directionConsistency] = measureDirectionMetrics(subMemos,parentMemo)
	
	% Rotating a diagonal still leaves a diagonal, which feels effectively
	% unchanged, and therefore worth grouping. Similarly, left and right
	% directions are equivalent up to rotations, so we'll group those too.
	% This leaves just distinguishing diagonal vs non-diagonal
	
	starts = [subMemos.start];
	stops  = [subMemos.stop];
	isDiag = mod(stops-starts,4) == 2;
	isDiagOrig = mod(parentMemo.stop-parentMemo.start,4) == 2;
	
	% Measure the variance treating these values as samples from a binomial
	% distribution.
	portionDiag = mean(isDiag);
	directionVariance = portionDiag * (1-portionDiag);
	
	% Measure how similar the sub-memos' directions are similar to the
	% parent memo's direction.
	directionConsistency = mean( isDiag == isDiagOrig );
	
end
function [parityVariance,parityConsistency,paritySplittingQuality] = measureParityMetrics(subMemos,parentMemo,orderArray)
	
	subParities = [ [subMemos.isOddH]; [subMemos.isOddW] ]; % 2xM
	origParity  = [ parentMemo.isOddH; parentMemo.isOddW ]; % 2x1
	
	% Generate a normalized histogram of these vertical/horizontal parties
	% in the domains (height=[even,odd])x(width=[even,odd])
	parityPMF = histcounts2( subParities(1,:),subParities(2,:), [-0.5,0.5,1.5],[-0.5,0.5,1.5]) / numel(subMemos);
	% Use this probability mass function to measure the centroid in each
	% dimension
	meanH = sum(sum(parityPMF .* [0;1]));
	meanW = sum(sum(parityPMF .* [0,1]));
	% Measure the literal variance for this distribution
	parityVariance = sum(sum(parityPMF .* (([0;1]-meanH).^2 + ([0,1]-meanW).^2)));
	
	% Measure how consistent the sub-memos' parities are wrt the parent
	% memo's parity.
	parityConsistency = mean(all(origParity == subParities,1));
	
	% Measure how consistent the splitting is compared to the parent. If
	% the parent was even in some dimension, it makes more sense to split
	% it into two pieces than 1 or 3 (although 1 is effectively a
	% non-split, so it feels like it shouldn't really count).
	sz = size(orderArray).';
	splitParity = mod(sz,2);
	isNonSplit = sz == 1;
	paritySplittingQuality = mean((splitParity == origParity) | isNonSplit);
	
end

function varargout = recordAndClose(f,~)
	persistent lastVec
	if isa(f,'matlab.ui.Figure')
		lastVec = arrayfun( @(c)c.Value, f.UserData);
		delete(f);
	else
		varargout{1} = lastVec;
	end
end