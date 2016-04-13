function [ fres, correlation ] = standardAdditionSlope( DATACELL, peakLocation, options )
%%
% DATACELL has four cell arrays:
%	Y: MATRIX with registered signals in each column
%	X: VECTOR with values to plot Y(:,i) against (the same number of rows as Y, and one
%	column)
%	CONC: VECTOR with values of concentration for each column of Y (so the same
%	number of columns and one row)
%	SENS: VECTOR with the INDEX (so 1:nrOfSens values) of sensitivity of each Y column (so the same
%	number of columns as Y and one row)
%
% peakLocation is a scalar with rough estimate of the signal of interest
%	maximum. It is expressed as a field number i.e. of single curve cosists
%	of 100 data points peakLocation can be a number between 1 and 100.
%
% options is a structure, which may contain fields:
% .average - boolean, true means that columns with the same concentration
%			and sensitivity should be everaged
% .smooth - boolean, true means that data should be smoothed before further
%			processing
% .forceSamePoints - boolean, true means that inflection points should only
%			be calculated for the curve with highest concentration 
%			and used for all (as opposed to calculating the
%			inflection point for all curves independiently.
%

	try
		options.average;
	catch
		options.average = false;
	end

	try 
		options.smooth;
	catch
		options.smooth = false;
	end

	try
		options.forceSamePoints;
	catch
		options.forceSamePoints = false;
	end

	try
		options.slopeInfectionMethod;
	catch
		options.slopeInfectionMethod = 'clone';
	end

	if ( size(DATACELL.SENS,2) ~= size(DATACELL.Y,2) || size(DATACELL.CONC,2) ~= size(DATACELL.Y,2) )
		disp('Number of columns in SENS CELL and Y has to be equal');
		return
	end
	
	dataY = DATACELL.Y;
	conc = DATACELL.CONC;
	sens = DATACELL.SENS;
	
	if ( options.average == true )
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

	[ concSort, sortOrder ] = sort(conc);
	sensSort = sens(1,sortOrder);
	dataYSort = zeros(size(dataY));
	dataYSort(:,:) = dataY(:,sortOrder);

	if ( options.smooth ) 
		for i=1:size(dataYSort,2)
			dataYSort(:,i) = smooth(dataYSort(:,i));
		end
	end

	for i=size(dataYSort,2):-1:1
		if ( i == size(dataYSort,2) )
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, 'notClone');
		else
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, options.slopeInfectionMethod, fitRange(:,i+1));
		end
	end
	
	%
	% NORMAL FITS:
	%
	figure;
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
	
	crossL = zeros(size(normalFitL,2),size(normalFitL,2),2);
	crossR = zeros(size(normalFitR,2),size(normalFitR,2),2);
	crossAVG = zeros(size(normalFitAVG,2),size(normalFitAVG,2),2);
	rpos = 1;
	fres = [];
	for i=1:size(normalFitAVG,2)
		for ii=1:size(normalFitAVG,2)
			if ( i == ii )
				continue;
			end
			% MATRIX solution
			crossL(i,ii,:) = pinv([normalFitL(2,i) -1; normalFitL(2,ii) -1]) * [-normalFitL(1,i);-normalFitL(1,ii)];
			crossR(i,ii,:) = pinv([normalFitR(2,i) -1; normalFitR(2,ii) -1]) * [-normalFitR(1,i);-normalFitR(1,ii)];
			crossAVG(i,ii,:) = pinv([normalFitAVG(2,i) -1; normalFitAVG(2,ii) -1]) * [-normalFitAVG(1,i);-normalFitAVG(1,ii)];;
			plot([ crossR(i,ii,1) concSort(end) ], [ crossR(i,ii,1) concSort(end) ].*normalFitR(2,i) + normalFitR(1,i), 'r-');
			plot([ crossL(i,ii,1) concSort(end) ], [ crossL(i,ii,1) concSort(end) ].*normalFitL(2,i) + normalFitL(1,i), 'b-');
			plot([ crossAVG(i,ii,1) concSort(end) ], [ crossAVG(i,ii,1) concSort(end) ].*normalFitAVG(2,i) + normalFitAVG(1,i), 'g-');
			if ( i > ii ) 
				fresL(rpos) = -crossL(i,ii,1);
				fresR(rpos) = -crossR(i,ii,1);
				fresAVG(rpos) = -crossAVG(i,ii,1);
				plot(crossL(i,ii,1),crossL(i,ii,2),'bx')
				plot(crossR(i,ii,1),crossR(i,ii,2),'rx')
				plot(crossAVG(i,ii,1),crossAVG(i,ii,2),'gx');
				rpos = rpos+1;
			end
		end
	end
	
	stdL = std(fresL);
	stdR = std(fresR);
	stdAVG = std(fresAVG);
	
	if ( min(correlation.L) > 0.9  ...
	&& ( stdL <= stdAVG || min(correlation.AVG) <= 0.9 ) ...
	&& ( stdL <= stdR || min(correlation.R) <= 0.9 ) )
		disp('Selecting left slope');
		fres = fresL;
	elseif ( min(correlation.R) > 0.9 ...
	&& ( stdR <= stdAVG || min(correlation.AVG) <= 0.9 ) ...
	&& ( stdR <= stdL || min(correlation.L) <= 0.9 ) )
		disp('Selecting right slope');
		fres = fresR;
	elseif ( min(correlation.AVG) > 0.9 )
		disp('Selecting average slope');
		fres = fresAVG;
	else 
		error('Could not select slope for calibration, please verify the data');
	end
	
	disp(sprintf('mean: %0.6e', mean(fres)));
	disp(sprintf('median: %0.6e ',median(fres)));

end

function [slopeL, slopeR, slopeAVGfitRange, fitRange] = getSlopeInInflection(signal, peak, method, fitRange )
	fitSize = 5;
	maxHit = 4;
	hitCnt = 0;
	verbose = true;
	sigLen = length(signal);
	prevNormalFit = [ NaN NaN ];
	finalNormalFit = [ NaN NaN ];
	%signal = smooth(signal,13,'sgolay',3);
	%signal = smooth(signal,17,'sgolay',3);
	
	%
	% FITTING LEFT
	%
	% normal equation fit:
	
	if ( strcmp(method, 'clone') ) 
		fitRange;
		[blackhole,pos] = max(isnan(fitRange));
		fitrangeL = fitRange(1:pos-1);
		fitrangeR = fitRange(pos+1:end);

		% Left
		X = [ ones(fitSize,1) (1:1:fitSize)' ];
		normalFitX = pinv( X' * X );
		normalFitX = normalFitX * X';
		yL = signal(fitrangeL,1);
		normalFitL = normalFitX * yL;
		
		%Right
		yR = signal(fitrangeR,1);
		normalFitR = normalFitX * yR;

		slopeL = normalFitL(2);
		slopeR = normalFitR(2);
		
	else
		%
		% FIT LEFT
		%
		X = [ ones(fitSize,1) (1:1:fitSize)' ];
		normalFitX = pinv( X' * X );
		normalFitX = normalFitX * X';
		fitrange = [];
		
% 		if ( verbose )
% 			figure;
% 		end
		for fitPos = peak:-1:fitSize
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
			% nie udało się znaleźć dopasowania %
			disp('failed');
		else
			slopeL = finalNormalFit(2);
			fitrangeL = finalFitrange;
		end
		%
		% FITTING RIGHT
		%
		fitrange = [];
% 		if ( verbose )
% 			figure;
% 		end
		for fitPos = peak:1:(sigLen-fitSize)
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
			% nie udało się znaleźć dopasowania %
		else
			slopeR = finalNormalFit(2);
			fitrangeR = finalFitrange;
		end
		fitRange = [ fitrangeL NaN fitrangeR ];
	end %end catch

	slopeAVGfitRange = (slopeL + -1*slopeR)/2;
	
	if ( verbose )
		plot(signal); hold on;
		plot(fitrangeL,signal(fitrangeL),'*b');
		plot(fitrangeR,signal(fitrangeR),'*g');
	end
	
end
