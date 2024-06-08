
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
			htcurve.MemoBase.orderArrayToList(...
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

%%

original = htcurve.MemoParity();
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
			
			subOrderings = cell(0,2);
			for splitInd = 1:numSplitCases
				tempOrderings = generateAllSplits( splittingCases(splitInd,:), original );
				
				% Simplify these solutions down as much as possible. We
				% don't need many representatives from each size class.
				tempOrderings = pruneSolutions(tempOrderings,original);
				
				% Accumulate these solutions
				subOrderings = [ subOrderings; tempOrderings ]; %#ok<AGROW>
			end
			
			
			
		end
	end
end

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
	list = htcurve.MemoParityList(s1*s2);
	
	% We'll gather valid suborderings as we find them
	subOrderings = cell(0,2);
	
	
	subSolns = splittingCase{4};
	for ssInd = 1:numel(subSolns)
		orderList = subSolns{ssInd};
		
		% Error checking and visualization
		orderArray = htcurve.MemoBase.orderListToArray( orderList ); % reconstruct, for error checking and visualization
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

function subOrderingBest = pruneSolutions(subOrderings, originalMemo)
	
	if isempty(subOrderings)
		subOrderingBest = subOrderings;
		return
	end
	
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
	[directionVariances,directionConsistencies] = cellfun(@(subMemos) measureDirectionUniformityConsistency(subMemos,originalMemo), subOrderings(:,2) );
	
	% Measure how uniformly the sub-solution parities are. Also measure the
	% portion of consistent parities with the parent parity.
	[parityVariances,parityConsistencies] = cellfun(@(subMemos) measureParityUniformityConsistency(subMemos,originalMemo), subOrderings(:,2) );
	
	% Summarize these qualities (weighting them differently). Here's the
	% priority: (#1 is most important)
	%   1. symmetric path                          MOST OBVIOUS "SYMMETRY"
	%   2. low direction variance                  INTER-CONSISTENCY
	%   3. low parity variance                     INTER-CONSISTENCY
	%   4. high direction consistency with parent  RECURSION CONSISTENCY
	%   5. high parity consistency with parent     RECURSION CONSISTENCY
	%   6. symmetric splitting parity              NICE TO HAVE
	%   7. solution guaranteed                     NICE TO HAVE (until we require it...)
	goodness = [isSymPath,-directionVariances,-parityVariances,directionConsistencies,parityConsistencies,isSymParity,isGuaranteed];
% Determine what the best is over the full set, elementwise. This will
% let us appreciate when certain qualities are impossible in this
% situation, so we can set our standards more appropriately.
maxGoodness = max(goodness,[],1);
	
% 	maxGoodness = [1,0,0,1,1,1,1];
	
% I THINK I NEED SOMETHING SMARTER THAN #2

	[~,order] = sortrows(goodness,'descend'); % put the best ones first
	
% 	close all
	kSet = find(all(goodness == goodness(order(1),:),2));
	for k = kSet.'
		figure;
		plot(subOrderings{k,2},subOrderings{k,1});
		title(sprintf('is guaranteed: %u',isGuaranteed(k)));
	end
	selection = order(1);
	
	if ~any(all(goodness==maxGoodness,2),1)
		disp('nothing perfect')
% 	else
% 		disp('something perfect')
	end
	
	% If we didn't pick out a guaranteed solution, we'll force it this time
	% while we grab a second pick
	if ~isGuaranteed(order(1))
		goodness = [isGuaranteed,goodness];
		[~,order] = sortrows(goodness,'descend'); % put the best ones first
		kSet = find(all(goodness == goodness(order(1),:),2));
		for k = kSet.'; figure; plot(subOrderings{k,2},subOrderings{k,1}); title(sprintf('is guaranteed: %u',isGuaranteed(k))); end
		selection(end+1) = order(1);
	end
	
	subOrderingBest = subOrderings(unique(selection),:);
	
end
function [directionVariance,directionConsistency] = measureDirectionUniformityConsistency(subMemos,parentMemo)
	starts = [subMemos.start];
	stops  = [subMemos.stop];
	isDiag = mod(stops-starts,4) == 2;
	isDiagOrig = mod(parentMemo.stop-parentMemo.start,4) == 2;
	portionDiag = mean(isDiag);
	directionVariance = portionDiag * (1-portionDiag);
	directionConsistency = mean( isDiag == isDiagOrig);
end
function [parityVariance,parityConsistency] = measureParityUniformityConsistency(subMemos,parentMemo)
	
	subParities = [ [subMemos.isOddH]; [subMemos.isOddW] ];
	origParity  = [ parentMemo.isOddH; parentMemo.isOddW ];
	
	parityPMF = histcounts2( subParities(1,:),subParities(2,:), [-0.5,0.5,1.5],[-0.5,0.5,1.5]) / numel(subMemos);
	mean1 = sum(sum(parityPMF .* [0;1]));
	mean2 = sum(sum(parityPMF .* [0,1]));
	parityVariance = sum(sum(parityPMF .* (([0;1]-mean1).^2 + ([0,1]-mean2).^2)));
	
	parityConsistency = mean(all(origParity == subParities,1));
	
end