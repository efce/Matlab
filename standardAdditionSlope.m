function [ fres, correlation ] = standardAdditionSlope( DATACELL, peakLocation, options )
%%
% DATACELL has four cell arrays:
%	Y: MATRIX with registered signals in each column
%	X: VECTOR with values to plot Y(:,i) against (the same number of rows as Y, and one
%	column)
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
% .average - boolean, true means that columns with the same concentration
%			and sensitivity should be everaged
% .smooth - boolean, true means that data should be smoothed before further
%			processing
% .forceSamePoints - boolean, true means that inflection points should only
%			be calculated for the curve with highest concentration 
%			and used for all (as opposed to calculating the
%			inflection point for all curves independiently.
%

	slopeDiffRequired = 0.01;
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
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, false);
		else
			[slopeL(1,i), slopeR(1,i), slopeAVG(1,i), fitRange(:,i)] = getSlopeInInflection(dataYSort(:,i), peakLocation, options.forceSamePoints, fitRange(:,i+1));
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
	avOK = true;
	rOK = true;
	lOK = true;
	for i=1:size(normalFitAVG,2)
		for ii=1:size(normalFitAVG,2)
			% Check if sensitivities are different enough
			if ( i >= ii )
				continue;
			end
			if ( avOK )
				prop = normalFitAVG(2,i) / normalFitAVG(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
					% MATRIX solution of intersection points for each
					% combination
					crossAVG(i,ii,:) = pinv([normalFitAVG(2,i) -1; normalFitAVG(2,ii) -1]) * [-normalFitAVG(1,i);-normalFitAVG(1,ii)];
					plot([ crossAVG(i,ii,1) concSort(end) ], [ crossAVG(i,ii,1) concSort(end) ].*normalFitAVG(2,i) + normalFitAVG(1,i), 'g-');
					fresAVG(rpos) = -crossAVG(i,ii,1);
					plot(crossAVG(i,ii,1),crossAVG(i,ii,2),'gx');
				else
					disp (['Sens1: ' num2str(normalFitAVG(2,i)) '; Sens2: ' num2str(normalFitAVG(2,ii)) '; Sens1/Sens2:' num2str(normalFitAVG(2,i)/normalFitAVG(2,ii)) ]);
					disp('Sensitivities are too similar for AVERAGE');
					avOK = false;
				end
			end
			
			if ( lOK )
				prop = normalFitL(2,i) / normalFitL(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
					% MATRIX solution of intersection points for each
					% combination
					crossL(i,ii,:) = pinv([normalFitL(2,i) -1; normalFitL(2,ii) -1]) * [-normalFitL(1,i);-normalFitL(1,ii)];
					plot([ crossL(i,ii,1) concSort(end) ], [ crossL(i,ii,1) concSort(end) ].*normalFitL(2,i) + normalFitL(1,i), 'b-');
					fresL(rpos) = -crossL(i,ii,1);
					plot(crossL(i,ii,1),crossL(i,ii,2),'bx');
				else
					disp (['Sens1: ' num2str(normalFitL(2,i)) '; Sens2: ' num2str(normalFitL(2,ii)) '; Sens1/Sens2:' num2str(normalFitL(2,i)/normalFitL(2,ii)) ]);
					disp('Sensitivities are too similar for LEFT');
					lOK = false;
				end
			end
			
			if ( rOK )
				prop = normalFitR(2,i) / normalFitR(2,ii);
				if ( prop > (1+slopeDiffRequired) || prop < (1-slopeDiffRequired) )
					% MATRIX solution of intersection points for each
					% combination
					crossR(i,ii,:) = pinv([normalFitR(2,i) -1; normalFitR(2,ii) -1]) * [-normalFitR(1,i);-normalFitR(1,ii)];
					plot([ crossR(i,ii,1) concSort(end) ], [ crossR(i,ii,1) concSort(end) ].*normalFitR(2,i) + normalFitR(1,i), 'r-');
					fresR(rpos) = -crossR(i,ii,1);
					plot(crossR(i,ii,1),crossR(i,ii,2),'rx');
				else
					disp (['Sens1: ' num2str(normalFitR(2,i)) '; Sens2: ' num2str(normalFitR(2,ii)) '; Sens1/Sens2:' num2str(normalFitR(2,i)/normalFitR(2,ii)) ]);
					disp('Sensitivities are too similar for RIGHT');
					rOK = false;
				end
			end
			
			rpos = rpos+1;
		end
	end
	toremove=[];
	if ( lOK )
		stdL = std(fresL)
		% Try Three-sigma to see if result are aligned;
		if ( length(fresL) >= 3 ) %only if there is three or more curves
			mfresL = mean(fresL);
			for si = 1:length(fresL)
				if ( abs(fresL(i)-mfresL) >= 3*stdL ) %three-sigma
					toremove = [ toremove i ];
				end
			end
		end
		fresL(toremove) = [];
		stdL=std(fresL);	
	else
		stdL = NaN;
	end
	if ( rOK )
		stdR = std(fresR)
		% Try Three-sigma to see if result are aligned;
		if ( length(fresR) >= 3 ) %only if there is three or more curves
			mfresR = mean(fresR);
			for si = 1:length(fresR)
				if ( abs(fresR(i)-mfresR) >= 3*stdR ) %three-sigma
					toremove = [ toremove i ];
				end
			end
		end
		fresR(toremove) = [];
		stdR=std(fresR)	
	else
		stdR = NaN;
	end
	if ( avOK )
		stdAVG = std(fresAVG)
		% Try Three-sigma to see if result are aligned;
		if ( length(fresAVG) >= 3 ) %only if there is three or more curves
			mfresAVG = mean(fresAVG);
			for si = 1:length(fresAVG)
				if ( abs(fresAVG(i)-mfresAVG) >= 3*stdAVG ) %three-sigma
					toremove = [ toremove i ];
				end
			end
		end
		fresAVG(toremove) = [];
		stdAVG=std(fresAVG)	
	else
		stdAVG = NaN;
	end
	
	if ( lOK && min(correlation.L) > 0.9  ...
	&& ( isnan(stdAVG) || stdL <= stdAVG || min(correlation.AVG) <= 0.9 ) ...
	&& ( isnan(stdR) || stdL <= stdR || min(correlation.R) <= 0.9 ) )
		disp('Selecting left slope');
		fres = fresL;
	elseif ( rOK && min(correlation.R) > 0.9 ...
	&& ( isnan(stdAVG) || stdR <= stdAVG || min(correlation.AVG) <= 0.9 ) ...
	&& ( isnan(stdL) || stdR <= stdL || min(correlation.L) <= 0.9 ) )
		disp('Selecting right slope');
		fres = fresR;
	elseif ( avOK && min(correlation.AVG) > 0.9 )
		disp('Selecting average slope');
		fres = fresAVG;
	else 
		error('Could not select slope for calibration, please verify the data');
	end
	
	disp(sprintf('mean: %0.6e', mean(fres)));
	disp(sprintf('median: %0.6e ',median(fres)));
	disp(sprintf('proponowany wynik: %0.6e ± %0.6e',mean(fres), (std(fres)/sqrt(numel(fres))*tinv(1-(.05/2),length(fres)-1))));

end

function [slopeL, slopeR, slopeAVGfitRange, fitRange] = getSlopeInInflection(signal, peak, forceFitRange, fitRange )
	fitSize = 5;
	maxHit = 4;
	hitCnt = 0;
	verbose = true;
	sigLen = length(signal);
	prevNormalFit = [ NaN NaN ];
	finalNormalFit = [ NaN NaN ];
	%signal = smooth(signal,13,'sgolay',3);
	signal = smooth(signal,17,'sgolay',3);
	
	%
	% FITTING LEFT
	%
	% normal equation fit:
	
	if ( forceFitRange ) 
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
	
	%
	% EXP
	%
	experimental = false;

	if ( experimental )
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
	end
	%
	%range2
	%

	slopeAVGfitRange = (slopeL + -1*slopeR)/2;
	
	if ( verbose )
		plot(signal); hold on;
		plot(fitrangeL,signal(fitrangeL),'*b');
		plot(fitrangeR,signal(fitrangeR),'*g');
	end
	
end
