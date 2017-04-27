function [ fres, correlation, calibration, fitRange, regressionEquations ] = SlopeStandardAdditionAnalysis(DATACELL, peakLocation, options)
%
% SlopeStandardAdditionAnalysis(DATACELL, peakLocation, options) is a function which tries
% to perform standard addition analysis using peak slope at the inflection points
% without any background correction. 
%
%
% DATACELL has four cell arrays (can be generated via prepareStructFromRawData):
%	Y: MATRIX with registered signals in each column
%	X: VECTOR with values to plot Y(:,i) against (the same number of rows as Y, and one
%	column) -- used only for plots.
%	CONC: VECTOR with values of concentration for each column of Y (so the same
%	number of columns and one row)
%	SENS: VECTOR with the INDEX (so 1:nrOfSens values) of sensitivity of each Y column (so the same
%	number of columns as Y and one row). There has to be minimum two
%	different sensitivities for each concentration.
%
% peakLocation is a scalar with rough estimate of the signal of interest
%	maximum. It is expressed as a field number i.e. of single curve cosists
%	of 100 data points peakLocation can be a number between 1 and 100.
%
% options is a structure, which may contain fields:
%	  .average - boolean, true means that columns with the same concentration
%				and sensitivity should be everaged
%	  .smooth - boolean, true means that data should be smoothed before further
%				processing
%	  .forceSamePoints - boolean, true means that inflection points should only
%				be calculated for the curve with highest concentration 
%				and used for all (as opposed to calculating the
%				inflection point for all curves independiently.

% Copyright (c) 2016, Filip Ciepiela <filip.ciepiela@agh.edu.pl> 
% All rights reserved.
%
%
% BSD-like license.
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
%
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

% Check input and set some default values (this may be tweaked)
%==============================================================
    close all;
    correlationTreshhold = 0.8;
    rowRemoveTresholdPercent = 0.3;
	slopeDiffRequired = 0.05;
	if ( ~isfield(options, 'average') )
		options.average = false;
	end
    
	if ( ~isfield(options,'smooth') )
		options.smooth = false;
	end

	if ( ~isfield(options, 'forceSamePoints') )
		options.forceSamePoints = true;
	end

	if ( size(DATACELL.SENS,2) ~= size(DATACELL.Y,2) || size(DATACELL.CONC,2) ~= size(DATACELL.Y,2) )
		disp('Number of columns in SENS CELL and Y has to be equal');
		return
	end
	
	dataY = DATACELL.Y;
	conc = DATACELL.CONC;
	sens = DATACELL.SENS;
	
	if ( options.average == true )
        % try to perform data averaging by sensitivity && concentration
        %==============================================================
		disp('Averaging curves with the same concentration and sensitivity');
		pos=1;
		for s=unique(sens)
			for c=unique(conc)
				sum = 0;
				num = 0;
				for i=1:size(conc,2)
					if ( s == sens(i) )
						if ( c == conc(i) )
							sum = sum+dataY(:,i);
							num = num+1;
						end
					end
				end
				newDataY(:,pos) = sum/num;
				newConc(1,pos)=c;
				newSens(1,pos)=s;
				pos=pos+1;
			end
		end
		dataY = newDataY;
		conc = newConc;
		sens = newSens;
	end
    
    % sort data by conc, so it is easier to follow the work flow
    %===========================================================
	[ concSort, sortOrder ] = sort(conc);
	sensSort = sens(1,sortOrder);
	dataYSort = zeros(size(dataY));
	dataYSort(:,:) = dataY(:,sortOrder);
    sensSort(:,:) = sens(1,sortOrder);

	if ( options.smooth )
        % smooth is required
        %===================
		for i=1:size(dataYSort,2)
			dataYSort(:,i) = smooth(dataYSort(:,i));
        end
	end

    % Very important step -- find the slope of each peak in inflection point. 
    % By default it reuses the point where the slope was taken from the peak
    % with the highest concentration -- if the location of the peak 
    % changes in subsequent data please use:
    % options.forceSamePoints = false;
    % Variables may end with R, L or AVG, which means they deal with
    % rigth, left or average slope in inflection point respectivly.
    % getSlopeInInflection is included in this file.
    %================================================================
    figure('Name','Plots');
    lineColors = hsv(length(unique(sens)));
	for i=size(dataYSort,2):-1:1
		if ( i == size(dataYSort,2) )
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, false, 0, lineColors(:,sensSort(i)));
		else
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, options.forceSamePoints, fitRange(:,i+1), lineColors(:,sensSort(i)));
		end
	end
	
	% It uses Normal Equation to find the optimal fit of calibration plot.
	%=====================================================================
	figure('Name', 'Fits');
	hold on;
	icons = [ 'o' '+' '*' 's' 'd' 'v' '^' '<' '>' 'p' 'h' ];
	if ( length(unique(sensSort)) > length(icons) ) 
		icons(1:length(unique(sensSort))) = 'o';
    end
	for s=unique(sensSort)
		list = (sensSort == s);
		% FIT using normal equation:
		X = [ ones(size(concSort(1,list),2),1) concSort(1,list)' ];
		normalFitX = pinv( X' * X );
		normalFitX = normalFitX * X';
		normalFitL(:,s) = normalFitX * slopeL(1,list)';
		normalFitR(:,s) = normalFitX * slopeR(1,list)';
		normalFitAVG(:,s) = normalFitX * slopeAVG(1,list)';
		concList(:,s) = concSort(1,list);
		plot(concList(:,s),slopeL(1,list),['b' icons(s) ]);
		plot(concList(:,s),slopeR(1,list),['r' icons(s) ]);
		tmp = corr( [(concList(:,s) .*normalFitL(2,s) + normalFitL(1,s)), slopeL(1,list)'] );
		correlation.L(s) = tmp(2,1);
		tmp = corr( [(concList(:,s) .*normalFitR(2,s) + normalFitR(1,s)), slopeR(1,list)'] );
		correlation.R(s) = tmp(2,1);
		tmp = corr( [(concList(:,s) .*normalFitAVG(2,s) + normalFitAVG(1,s)), slopeAVG(1,list)'] );
		correlation.AVG(s) = tmp(2,1);
	end
    
    calibration.L = slopeL(1,:);
    calibration.R = slopeR(1,:);
    calibration.AVG = slopeAVG(1,:);
	
    % Generate initial matrix of intercepts
    %======================================
	crossL = zeros(size(normalFitL,2),size(normalFitL,2),2);
	crossR = zeros(size(normalFitR,2),size(normalFitR,2),2);
	crossAVG = zeros(size(normalFitAVG,2),size(normalFitAVG,2),2);
    
	rpos = 1;
	fres = [];
    
    % We need to verify if slopes for each peak are found
    %====================================================
	avOK = true;
	rOK = true;
	lOK = true;
	if ( any(isnan(slopeR)))
		rOK=false;
		avOK=false;
	end
	if ( any(isnan(slopeL)))
		lOK=false;
		avOK=false;
    end
    
	for i=1:size(normalFitAVG,2)
		for ii=1:size(normalFitAVG,2)
			if ( i >= ii )
                crossAVG(i,ii,1) = NaN;
                crossR(i,ii,1) = NaN;
                crossL(i,ii,1) = NaN;
				continue;
            end
            % Note that it finds the intersect point for two lines at the
            % time, filling upper triangle of the matrix
            %============================================================
			if ( avOK )
				prop = normalFitAVG(2,i) / normalFitAVG(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
                    % Check if sensitivities are different enough, to small
                    % difference will result in serious problem with precission
                    % (and possible accuracy) of the result
                    %==========================================================
					crossAVG(i,ii,:) = pinv([normalFitAVG(2,i) -1; normalFitAVG(2,ii) -1]) * [-normalFitAVG(1,i);-normalFitAVG(1,ii)];
					plot([ crossAVG(i,ii,1) concSort(end) ], [ crossAVG(i,ii,1) concSort(end) ].*normalFitAVG(2,i) + normalFitAVG(1,i), 'g-');
                    plot([ crossAVG(i,ii,1) concSort(end) ], [ crossAVG(i,ii,1) concSort(end) ].*normalFitAVG(2,ii) + normalFitAVG(1,ii), 'g-');
					fresAVG(rpos) = -crossAVG(i,ii,1);
					plot(crossAVG(i,ii,1),crossAVG(i,ii,2),'gx','MarkerSize',20);
				else
					disp (['Sens1: ' num2str(normalFitAVG(2,i)) '; Sens2: ' num2str(normalFitAVG(2,ii)) '; Sens1/Sens2:' num2str(normalFitAVG(2,i)/normalFitAVG(2,ii)) ]);
					disp('Sensitivities are too similar for AVERAGE');
					avOK = false;
                end
                
            end

			
            % The same as above, but for L
            %=============================
			if ( lOK )
				prop = normalFitL(2,i) / normalFitL(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
					crossL(i,ii,:) = pinv([normalFitL(2,i) -1; normalFitL(2,ii) -1]) * [-normalFitL(1,i);-normalFitL(1,ii)];
					plot([ crossL(i,ii,1) concSort(end) ], [ crossL(i,ii,1) concSort(end) ].*normalFitL(2,i) + normalFitL(1,i), 'b-');
                    plot([ crossL(i,ii,1) concSort(end) ], [ crossL(i,ii,1) concSort(end) ].*normalFitL(2,ii) + normalFitL(1,ii), 'b-');
					fresL(rpos) = -crossL(i,ii,1);
					plot(crossL(i,ii,1),crossL(i,ii,2),'bx','MarkerSize',20);
				else
					disp (['Sens1: ' num2str(normalFitL(2,i)) '; Sens2: ' num2str(normalFitL(2,ii)) '; Sens1/Sens2:' num2str(normalFitL(2,i)/normalFitL(2,ii)) ]);
					disp('Sensitivities are too similar for LEFT');
					lOK = false;
                end
            end
			
            % The same as above but for R
            %============================
			if ( rOK )
				prop = normalFitR(2,i) / normalFitR(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
					crossR(i,ii,:) = pinv([normalFitR(2,i) -1; normalFitR(2,ii) -1]) * [-normalFitR(1,i);-normalFitR(1,ii)];
					plot([ crossR(i,ii,1) concSort(end) ], [ crossR(i,ii,1) concSort(end) ].*normalFitR(2,i) + normalFitR(1,i), 'r-');
                    plot([ crossR(i,ii,1) concSort(end) ], [ crossR(i,ii,1) concSort(end) ].*normalFitR(2,ii) + normalFitR(1,ii), 'r-');
					fresR(rpos) = -crossR(i,ii,1);
					plot(crossR(i,ii,1),crossR(i,ii,2),'rx','MarkerSize',20);
				else
					disp (['Sens1: ' num2str(normalFitR(2,i)) '; Sens2: ' num2str(normalFitR(2,ii)) '; Sens1/Sens2:' num2str(normalFitR(2,i)/normalFitR(2,ii)) ]);
					disp('Sensitivities are too similar for RIGHT');
					rOK = false;
                end
            end
			rpos = rpos+1;
		end
    end
    
    regressionEquations.AVG = normalFitAVG;
    regressionEquations.L = normalFitL;
    regressionEquations.R = normalFitR;
    
    % Here, is a little trick, to remove the intersection points
    % which are too far from average. It is done by the means of
    % Coefficient of Variance value, and can be tweaked in the setting at
    % the beggining of the file.
    % minimizeCV is included in this file.
    %====================================================================
    [ crossAVG, removedORDER ] = minimizeCV(squeeze(crossAVG(:,:,1)), 3, rowRemoveTresholdPercent);
    for i=1:length(removedORDER)
        correlation.AVG(removedORDER(i)) = [];
    end
    fresAVG = -crossAVG(logical(triu(ones(size(crossAVG)),1)));
    stdAVG = std(fresAVG);
    
    [ crossR, removedORDER ] = minimizeCV(squeeze(crossR(:,:,1)), 3, rowRemoveTresholdPercent);
    for i=1:length(removedORDER)
        correlation.R(removedORDER(i)) = [];
    end
    fresR = -crossR(logical(triu(ones(size(crossR)),1)));
    stdR = std(fresR);
    
     [ crossL, removedORDER ] = minimizeCV(squeeze(crossL(:,:,1)), 3, rowRemoveTresholdPercent);
    for i=1:length(removedORDER)
        correlation.L(removedORDER(i)) = [];
    end
    fresL = -crossL(logical(triu(ones(size(crossL)),1)));
    stdL = std(fresL);
    
    disp('Corelations Left:')
    for i=1:size(correlation.L,2)
        disp(correlation.L(i));
    end
    disp('===============')
    disp('Corelations Right:')
    for i=1:size(correlation.R,2)
        disp(correlation.R(i));
    end
    disp('===============')
    disp('Corelations AVG:')
    for i=1:size(correlation.AVG,2)
        disp(correlation.AVG(i));
    end
	
    % Here it selects, if it is possible to offer the final result, for which
    % set of data, the result is the best (left, right or average)
    %=====================================================================
	if ( lOK && min(correlation.L) > correlationTreshhold  ...
	&& ( isnan(stdAVG) || stdL <= stdAVG || min(correlation.AVG) <= correlationTreshhold ) ...
	&& ( isnan(stdR) || stdL <= stdR || min(correlation.R) <= correlationTreshhold ) )
		disp('Selecting left slope');
		fres = fresL;
	elseif ( rOK && min(correlation.R) > correlationTreshhold ...
	&& ( isnan(stdAVG) || stdR <= stdAVG || min(correlation.AVG) <= correlationTreshhold ) ...
	&& ( isnan(stdL) || stdR <= stdL || min(correlation.L) <= correlationTreshhold ) )
		disp('Selecting right slope');
		fres = fresR;
	elseif ( avOK && min(correlation.AVG) > correlationTreshhold )
		disp('Selecting average slope');
		fres = fresAVG;
	else 
		error('Could not select slope for calibration, please verify the data');
    end

	disp(sprintf('Partial results median: %0.6e ',median(fres)));
	disp(sprintf('The final result of standard addition: %0.6e ± %0.6e',mean(fres), (std(fres)/sqrt(numel(fres))*tinv(1-(.05/2),length(fres)-1))));

end

function [slopeL, slopeR, slopeAVGfitRange, fitRange] = getSlopeInInflection(signal, peak, forceFitRange, fitRange, lineColor )
% I find this as one of the most important steps, I has gone trou many
% iterations, so the code is a bit of mixture of different ideas.
% Some tweakalbe setting:
%======================================================================
	fitSize = 5; %How many points should be fitted to get slope%
	maxHit = 4;  %How many times the slope has to change to call it inflection point%
	verbose = true; %Draw some additional plots%
    
    % Prepare data structures:
    %=========================
    hitCnt = 0;
	sigLen = length(signal);
	prevNormalFit = [ NaN NaN ];
	finalNormalFit = [ NaN NaN ];
	%signal = smooth(signal,13,'sgolay',3);
	signal = smooth(signal,17,'sgolay',3);
	
	if ( forceFitRange ) 
        % If we have already the point range in which the slope has to be
        % found (i.e. this is not the plot with the highest concentration)
        %=================================================================
		[blackhole,pos] = max(isnan(fitRange));
		fitrangeL = fitRange(1:pos-1);
		fitrangeR = fitRange(pos+1:end);

		% Get linear fit on the left side
        %===========================
		X = [ ones(fitSize,1) (1:1:fitSize)' ];
		normalFitX = pinv( X' * X );
		normalFitX = normalFitX * X';
		yL = signal(fitrangeL,1);
		normalFitL = normalFitX * yL;
		
		% Get linear fit on the rigth side
        %===========================
		yR = signal(fitrangeR,1);
		normalFitR = normalFitX * yR;

        % Get slopes:
        %======================
		slopeL = normalFitL(2);
		slopeR = normalFitR(2);
		
	else
		% This is either the plot with the highest concentration or
        % option.forceSamePoints == false
        % So first, we need to find the inflection point, and then
        % compute its slope.
        %==========================================================
		X = [ ones(fitSize,1) (1:1:fitSize)' ];
		normalFitX = pinv( X' * X );
		normalFitX = normalFitX * X';
		fitrange = [];
% 		if ( verbose )
% 			figure;
% 		end
		for fitPos = peak:-1:fitSize
            % We start from the left side of the peak:
            %=========================================
			if ( peak < fitPos-fitSize )
				error('error');
			end
			fitrange = [ (fitPos-fitSize+1) :1: fitPos ];
			y = signal(fitrange,1);
			normalFit = normalFitX * y;

			if ( isnan(prevNormalFit(1)) )
				prevNormalFit = normalFit;
			elseif ( normalFit(2) < prevNormalFit(2) )
				if ( hitCnt == 0 ) 
					finalFitrange = fitrange;
					finalNormalFit = normalFit;
				end
				hitCnt = hitCnt+1;
				if ( hitCnt == maxHit ) 
					break;
				end
			elseif ( normalFit(2) >= prevNormalFit(2) )
				finalNormalFit = [ NaN NaN ];
				hitCnt = 0;
			end
			prevNormalFit = normalFit;

% 			if ( verbose ) 
% 				plot(X(:,2)+fitPos-5,y,'r*');
% 				hold on;
% 				plot(X(:,2)+fitPos-5,normalFit(1)+normalFit(2).*X(:,2));
% 			end
		end
		if ( isnan(finalNormalFit(2)) ) 
			% fit failed %
			disp('failed');
		else
			slopeL = finalNormalFit(2);
			fitrangeL = finalFitrange;
        end
		fitrange = [];
% 		if ( verbose )
% 			figure;
% 		end
		for fitPos = peak:1:(sigLen-fitSize)
            % And now rigth side:
            %====================
			fitrange= [ fitPos :1: (fitPos+fitSize-1) ];
			y = signal(fitrange,1);
			normalFit = normalFitX * y;
			if ( isnan(prevNormalFit(2)) )
				prevNormalFit = normalFit;
			elseif ( normalFit(2) > prevNormalFit(2) )
				if ( hitCnt == 0 ) 
					finalFitrange = fitrange;
					finalNormalFit = normalFit;
				end
				hitCnt = hitCnt+1;
				if ( hitCnt == maxHit ) 
					break;
				end
			elseif ( normalFit(2) <= prevNormalFit(2) )
				finalNormalFit = [ NaN NaN ];
				hitCnt = 0;
			end
			prevNormalFit = normalFit;

% 			if ( verbose ) 
% 				plot(X(:,2)+fitPos-1,y,'r*');
% 				hold on;
% 				plot(X(:,2)+fitPos-1,normalFit(1)+normalFit(2).*X(:,2));
% 			end
		end
		if ( isnan(finalNormalFit(2)) ) 
			% fit failed %
		else
			slopeR = finalNormalFit(2);
			fitrangeR = finalFitrange;
		end
		fitRange = [ fitrangeL NaN fitrangeR ];
	end %end catch
	
	% This part provide alternative method of getting the slope in the
	% infleciton, as any noise will make it much more difficult, here are
	% tested approaches. This overrites previous method, as is more
	% difficult to follow.
    %=====================================================================
	experimental = 1;

	if ( experimental == 1 )
		a=( signal(fitrangeR(end))-signal(fitrangeL(1)) )/( fitrangeR(end)-fitrangeL(1));
		b=signal(fitrangeL(1)) - a*fitrangeL(1);
		p = [ a b ];
		%p = polyfit( [ fitrangeL(1) fitrangeR(end) ]', signal([ fitrangeL(1) fitrangeR(end) ]),1), pause
		levelPeak = signal - polyval(p,[ 1 : numel(signal) ])';
		range2 = [ floor((fitrangeR(end)+fitrangeL(1))/2)-ceil(1.2*(fitrangeR(end)-fitrangeL(1))) ceil((fitrangeR(end)+fitrangeL(1))/2)+ceil(1.2*(fitrangeR(end)-fitrangeL(1))) ];
		%levelPeak = levelPeak - min(levelPeak(range2(1):range2(end)));
		f = fit( [range2(1):range2(end)]', levelPeak([range2(1):range2(end)]), 'smoothingspline' );
		%figure(99);plot(f,[1:numel(signal)], levelPeak);hold on; pause
		slopeL = f(fitrangeL(2)) - f(fitrangeL(1));
		slopeR = f(fitrangeR(end)) - f(fitrangeR(end-1));
	elseif ( experimental == 2 ) 
		d1 = sgsdf(signal,2,1,0,0);
		d2 = sgsdf(d1,2,1,0,0);
		plot(signal,'k');hold on;plot(d1,'b'); hold on; plot(d2,'r')
		vl=abs(d2(1));
		pl=0;
		vr=abs(d2(end));
		pr=0;
		fr = [ fitrangeL(1)-20 : fitrangeR(end)+20 ];
		for i= fr;
			if ( abs(d2(i)) < vl && d2(i)>d2(i-1) ) %signal is increasing 
				vl=abs(d2(i));
				pl=i;
			elseif ( abs(d2(i)) < vr && d2(i)<d2(i-1) ) %signal is decreasing
				vr = abs(d2(i));
				pr = i;
			end
        end
        
        %Set final values:
        %=================
		if ( pl ~= 0 ) 
			slopeL = d1(pl);
		else
			slopeL = NaN;
		end
		if ( pr ~= 0 ) 
			slopeR = d1(pr);
		else 
			slopeR = NaN;
		end
			
    end

    % Get average slope (rigth slope has negative slope (or should have):
    %==================================================================
	slopeAVGfitRange = (slopeL - slopeR)/2;
	
	if ( verbose )
		plot(signal,'Color', lineColor); hold on;
		plot(fitrangeL,signal(fitrangeL),'*b', 'MarkerSize', 2);
		plot(fitrangeR,signal(fitrangeR),'*r', 'MarkerSize', 2);
	end
	
end

function [ matrixNANDIAG, removedORDER ] = minimizeCV( matrixNANDIAG, minEntities, minCVchange)
    % We can try 3-sigma here, however I think, this needs to be more
    % tweakable:
    %================================================================
    removedORDER = [];
    while ( size(matrixNANDIAG,1) > minEntities )
        inrow = inrows(matrixNANDIAG);
        oldstdrows = std(inrow');
        oldmeancvtotal = std(inrow(:)) / abs(mean(mean(inrow)));
        [mval, mpos]=max(oldstdrows);
        newmatrixNANDIAG = matrixNANDIAG;
        newmatrixNANDIAG(:,mpos) = [];
        newmatrixNANDIAG(mpos,:) = [];
        newinrow = inrows(newmatrixNANDIAG);
        newstd = std(newinrow(:));
        newmeancvtotal = newstd / abs(mean(mean(newinrow)));
        if ( (oldmeancvtotal-newmeancvtotal) > minCVchange )
            matrixNANDIAG = newmatrixNANDIAG;
            removedORDER = [ removedORDER mpos ];
        else
            break;
        end
    end
end

function inrow = inrows(mat)
    inrow = zeros(size(mat,1), size(mat,1)-1);
    pos = ones(size(mat,1) , 1);
    for i=1:size(mat,1)
        for ii=1:size(mat,1)
            if ( ~isnan(mat(i,ii)) )
                inrow(i,pos(i)) = mat(i,ii);
                inrow(ii,pos(ii)) = mat(i,ii);
                pos(i) = pos(i)+1;
                pos(ii) = pos(ii)+1;
            end
        end
    end
end
                    
            
