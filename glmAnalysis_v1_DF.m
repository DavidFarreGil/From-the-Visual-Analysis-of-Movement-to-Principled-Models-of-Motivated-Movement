clear all;

clear;
clear; ignasiglob; cptf = 1; cptr = 0; beep off;
clearfigure = 1;  SaveFigure = PDFFIGURE; closewindows = 1;

if ( 1 )	% parameters
	%%% functions
	% ignasiaxisfun
	% ignasifig

	% ignasimean
	% ignasimean(x, stdtype, stdrange)
	%	x: data
	%	stdtype: 1=mean - 2=median
	%	stdrange: restriction to mean ? stdrange*std
	%		no restriction if stdrange = 0

	% ignasisem
	% ignasistd

	% ignasiplotcurve
	% ignasiplotline
	% ignasiplotmarkercurve
	% ignasiplotmarkercurveerrorbar
	% ignasiplotpoint
	% ignasiplotpointerrorbar

	% ignasiplot1
	%	plot of amplitude and duration as a function of period
	%	single subject, all subjects, drummers, nondrummers

	% ignasiplot2
	%	plot of percent
	%	only single subjects

	% ignasiplot3
	%	plot of value/std(value)
	%	only single subjects

	datadirectory = 'ignasidata';
    
    % list of boundaries
	minperiod = 0; maxperiod = 2.2;
	minamplitude = 0; maxamplitude = 1.2;
	minasynchrony = -.220; maxasynchrony = .080;
	minasynchronystd = -0; maxasynchronystd = .15;
	minpeakvelocity = 0; maxpeakvelocity = 5;
	mindwelltime = 0; maxdwelltime = 2;
	mindownswingtime = 0; maxdownswingtime = 1;
	minupswingtime = 0; maxupswingtime = 1;
	mindwelltimestd = 0; maxdwelltimestd = .2;
	minenergy = 0; maxenergy = 1500;
	minspatialvar = 0; maxspatialvar = 10;
	mintimevar = -.04; maxtimevar = .04;

    % list of xticks
    xtickamplitude = [0 0.5 1 1.5 2.0];
    ytickamplitude = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1];
	SubjectNumber = 1:12;
	SubjectNondrummerNumber = 1:8;
	SubjectDrummerNumber = 9:12;
	SubjectOrder = [ 1 2 3 5 6 7 8 9 10 11 12 13 ];
	SubjectId = [ 3 4 5 12 13 14 15 17 18 24 25 26 ];

	stdtype = STDNORMAL;
	stdrange = 2;
	plevel = 0;

	FrequencyList = [ -1 0.5 1.0 1.5 2.0 2.5 ];

	FREEFREQUENCY = 1;

	% plot
	AmplitudeSymbols = { DIAMOND, CIRCLE, MYSQUARE, PENTAGRAM, HEXAGRAM };
	FrequencyEdgeColor = { BLACK, RED, DARKYELLOW, GREEN,  BLUE, PURPLE } ;
	FrequencyFaceColor = { BLACK, RED, DARKYELLOW, GREEN,  BLUE, PURPLE } ;
	FreeFrequencyEdgeColor = BLACK;
	FreeFrequencyFaceColor = BLACK;
	AmplitudeLinestyles = { SOLID, SOLID, SOLID };

    FrequencyEdgeColorSubject = { ORANGE, DARKRED, CYAN, GREEN, RED, PURPLE, BLUE, WHITE, BLACK, RED, DARKYELLOW  };
	FrequencyFaceColorSubject = { ORANGE, DARKRED, CYAN, GREEN, RED, PURPLE, BLUE, WHITE, BLACK, RED,DARKYELLOW  };

    AmplitudeSubjectsSymbols = { DIAMOND, CIRCLE, MYSQUARE, DIAMOND, CIRCLE, MYSQUARE, DIAMOND, CIRCLE, MYSQUARE, DIAMOND };
	AmplitudeLineAllSubjectsstyles = { SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID };

    FrequencyFaceColor1 = { { BLACK, DARKRED, DARKYELLOW, DARKGREEN, DARKBLUE, DARKPURPLE }, ...
                    { DARKGRAY, RED, YELLOW, GREEN, BLUE, PURPLE }, ...
                    { MYGRAY, LIGHTRED, LIGHTYELLOW, LIGHTGREEN, LIGHTBLUE, LIGHTPURPLE } };

    FreeFrequencyEdgeColor1 = { BLACK, DARKGRAY, MYGRAY };

    BlackColor = { { BLACK, BLACK, BLACK, BLACK, BLACK, BLACK }, ...
                    { BLACK, BLACK, BLACK, BLACK, BLACK, BLACK }, ...
                    { BLACK, BLACK, BLACK, BLACK, BLACK, BLACK } };
     
    
	%DWELLTIMESYMBOL = LEFTTRIANGLE;
	%MOVEMENTTIMESYMBOL = RIGHTTRIANGLE;
	DWELLTIMESYMBOL = CIRCLE;
	MOVEMENTTIMESYMBOL = CIRCLE;
	CYCLETIMESYMBOL = PENTAGRAM;

	DWELLTIMECOLOR = MYGRAY;
	MOVEMENTTIMECOLOR = BLACK;
	CYCLETIMECOLOR = GREEN;
	RECONSTRUCTEDCYCLETIMECOLOR = ORANGE;
	DWELLMOVEMENTCOVARIANCECOLOR = PURPLE;

	colormap(jet);
	cmap = colormap;

% 	[myaxis, tx, ty] = axisfun(0, 1, 0, 1);
% 	%_______________________________________
% 	figpreps(0, 'Color', [], [], 0, 0, 0, 0, 6, 4, wfs, 0, myaxis, tx, ty, '', '');
% 
% 	n = length(cmap);
% 	for i=1:n
% 		plotline([0 i/n], [1 i/n], cmap(i,:), SOLID, 3*dlw);
% 	end
% 	%_______________________________________

	n = length(cmap);
	p = length(SubjectNumber);
	for i=1:p
		j = floor(i*(n/(p+1)));
		SubjectColors{i} = cmap(j,:);
	end

% 	[myaxis, tx, ty] = axisfun(0, 1, 0, 1);
% 	%_______________________________________
% 	figpreps(0, 'Color', [], [], 0, 0, 0, 0, 6, 4, wfs, 0, myaxis, tx, ty, '', '');
% 
% 	for i=1:p
% 		plotline([0 i/p], [1 i/p], SubjectColors{i}, SOLID, 3*dlw);
% 	end
% 	%_______________________________________

	markersize = 3;
	linewidth = 1;
	figuresizex = 5;
	figuresizey = 4;
	labelfontsize = 9;
	tickfontsize = 8;

end


taskID = 1;
taskName = 'precisionMovement';

nTrialsBlock = 100;

listCoefficients  = {'0','E','ST','#T','BCxC', ...
                     'ExST','Ex#T','ExBCxC','STx#T','STxBCxC', ...
                     '#TxBCxC','ExSTx#T','ExSTxBCxC','Ex#TxBCxC','STx#TxBCxC', ...
                     'ExSTx#TxBCxC'};

listCoefficients  = {'0','E','ST','#T', ...
                     'ExST','Ex#T','STx#T', ...
                     'ExSTx#T'};
% open connection
db = nmOpen('precisionMovement',{'ignasi','','com.mysql.jdbc.Driver','jdbc:mysql://localhost:3306/'},'STRUCT=1');

listSubjects = 2:12;

numS = numel(listSubjects)
nSubj = ceil(sqrt(numS));

% for each subject
nControl = 0:1;
nCoeffs = 8;

listVariables = {'PV1','MT','E'};

nBound = 0;

% level of control
% for each movement
XT = zeros(12*numel(listSubjects), 7*numel(listSubjects) + numel(listSubjects));
YT = zeros(12*numel(listSubjects), 3);

listTrials = {  [1 7 13 19], ... % E=0, ST
                [2 3 14 20], ... % E=0, LT
                [3 9 15 21], ... % E=1, ST
                [4 10 16 22], ... % E=1, LT
                [5 11 17 23], ... % E=2, ST
                [6 12 18 24]}; % E=2, LT

listE = [0 0 1 1 2 2];
listT = [1 -1 1 -1 1 -1];

listC = [-1 1];
% group analysis
for k=1:numel(listSubjects),

    dbss2 = nmSubset(db,sprintf('idSubject=%d',listSubjects(k)));

    X = [];
    Y = [];
    for p=1:6, 
        
        myt = listTrials{p};

        mys = sprintf('idSubject=%d and nChoice is NOT NULL and (trialTypeNumber=%d',listSubjects(k),myt(1));
        for q=2:4,
            mys = sprintf('%s or trialTypeNumber=%d',mys,myt(q));
        end
        mys = sprintf('%s)',mys);

        % choice
        for r=1:2,
            mys2 = sprintf('%s and nChoice=%d',mys,r);
            dbs = nmSubset(db,mys2);
            data = nmGet(dbs, {'tMvmtOnset','tMvmtOffset', 'MT=tEnterTarget-tMvmtOnset','PV=nPeakVelocity','tt=trialTypeNumber','E=nYourError','nB=nBlock','nT=nTrial','myE=nYourError'});

            X = [ X; listE(p) listT(p) listC(r)];
            Y = [ Y; nanmean(data.MT) nanmean(data.myE) nanmean(data.PV)];
        end

    end
    for p=1:size(Y,2),
        YZ(:,p) = (Y(:,p) - nanmean(Y(:,p)))/nanstd(Y(:,p));
    end
    
    %Xs 
    for p=1:size(X,2),
        XZ(:,p) = (X(:,p) - nanmean(X(:,p)))/nanstd(X(:,p));
    end
    X(:,4) = X(:,1).*X(:,2);
    X(:,5) = X(:,1).*X(:,3);
    X(:,6) = X(:,2).*X(:,3);
    X(:,7) = X(:,1).*X(:,2).*X(:,3);
    for p=1:size(X,2),
        XZ(:,p) = (X(:,p) - nanmean(X(:,p)))/nanstd(X(:,p));
    end

    YT((k-1)*12+1:k*12,:) = YZ;
    XT((k-1)*12+1:k*12, (k-1)*7+1:k*7 ) = XZ;
    XT((k-1)*12+1:k*12, 7*numel(listSubjects)+k) = ones(12,1);

    clear('X','Y');
end

% regression 
[b1,BINT,R,RINT,stats] = regress(YT(:,1),XT);
betas1(1,:) = b1(7*numel(listSubjects)+1:end);
betas1(2:8,:) = reshape(b1(1:7*numel(listSubjects)),7,numel(listSubjects));

[b2,BINT,R,RINT,stats] = regress(YT(:,2),XT);
betas2(1,:) = b2(7*numel(listSubjects)+1:end);
betas2(2:8,:) = reshape(b2(1:7*numel(listSubjects)),7,numel(listSubjects));

[b3,BINT,R,RINT,stats] = regress(YT(:,3),XT);
betas3(1,:) = b3(7*numel(listSubjects)+1:end);
betas3(2:8,:) = reshape(b3(1:7*numel(listSubjects)),7,numel(listSubjects));

for p=1:size(betas1,1),
    [h,pp1(p)] = ttest(betas1(p,:));
    [h,pp2(p)] = ttest(betas2(p,:));
    [h,pp3(p)] = ttest(betas3(p,:));
end


g1 = figure;
subplot(1,3,1),barweb(nanmean(betas1,2), nanstd(betas1,0,2)/sqrt(numel(listSubjects)), 1, [], [], {}, 'Beta Coeffs. (MT)', bone, [], {} );
subplot(1,3,2),barweb(nanmean(betas2,2), nanstd(betas2,0,2)/sqrt(numel(listSubjects)), 1, [], [], {}, 'Beta Coeffs. (Error)', bone, [], {} );
subplot(1,3,3),barweb(nanmean(betas3,2), nanstd(betas3,0,2)/sqrt(numel(listSubjects)), 1, [], [], {}, 'Beta Coeffs. (PV)', bone, [], {} );
set(g1,'color',[1 1 1]);
print(g1,'-dpdf','-r600',sprintf('Beta-CoeffsGroup2.pdf'));    
 
% %     MT = data.MT;
% %     PV1 = data.nPeakVelocity;
% %     E = data.E;
% 
% 
%     X1 = data.nMotivationLevel;
%     X2 = nB;
%     X3 = nTnB;
%     X4 = data.nMajorMinor;
%     X5 = data.nChoice;
%     X45 = X4.*X5; 
%     i45 = find(X45==4); X45(i45)=1;
% 
%     % Biomech=1 & Choice=1
%     i1x = find(X1==0 & X2==1 & X3<54 & X45==1);
%     i2x = find(X1==1 & X2==1 & X3<54 & X45==1);
%     i3x = find(X1==2 & X2==1 & X3<54 & X45==1);
%     i4x = find(X1==0 & X2==2 & X3<54 & X45==1);
%     i5x = find(X1==1 & X2==2 & X3<54 & X45==1);
%     i6x = find(X1==2 & X2==2 & X3<54 & X45==1);
%     i7x = find(X1==0 & X2==3 & X3<54 & X45==1);
%     i8x = find(X1==1 & X2==3 & X3<54 & X45==1);
%     i9x = find(X1==2 & X2==3 & X3<54 & X45==1);
% 
%     i10x = find(X1==0 & X2==1 & X3>54 & X45==1);
%     i11x = find(X1==1 & X2==1 & X3>54 & X45==1);
%     i12x = find(X1==2 & X2==1 & X3>54 & X45==1);
%     i13x = find(X1==0 & X2==2 & X3>54 & X45==1);
%     i14x = find(X1==1 & X2==2 & X3>54 & X45==1);
%     i15x = find(X1==2 & X2==2 & X3>54 & X45==1);
%     i16x = find(X1==0 & X2==3 & X3>54 & X45==1);
%     i17x = find(X1==1 & X2==3 & X3>54 & X45==1);
%     i18x = find(X1==2 & X2==3 & X3>54 & X45==1);
% 
%     i19x = find(X1==0 & X2==1 & X3<54 & X45==2);
%     i20x = find(X1==1 & X2==1 & X3<54 & X45==2);
%     i21x = find(X1==2 & X2==1 & X3<54 & X45==2);
%     i22x = find(X1==0 & X2==2 & X3<54 & X45==2);
%     i23x = find(X1==1 & X2==2 & X3<54 & X45==2);
%     i24x = find(X1==2 & X2==2 & X3<54 & X45==2);
%     i25x = find(X1==0 & X2==3 & X3<54 & X45==2);
%     i26x = find(X1==1 & X2==3 & X3<54 & X45==2);
%     i27x = find(X1==2 & X2==3 & X3<54 & X45==2);
% 
%     i28x = find(X1==0 & X2==1 & X3>54 & X45==2);
%     i29x = find(X1==1 & X2==1 & X3>54 & X45==2);
%     i30x = find(X1==2 & X2==1 & X3>54 & X45==2);
%     i31x = find(X1==0 & X2==2 & X3>54 & X45==2);
%     i32x = find(X1==1 & X2==2 & X3>54 & X45==2);
%     i33x = find(X1==2 & X2==2 & X3>54 & X45==2);
%     i34x = find(X1==0 & X2==3 & X3>54 & X45==2);
%     i35x = find(X1==1 & X2==3 & X3>54 & X45==2);
%     i36x = find(X1==2 & X2==3 & X3>54 & X45==2);
% 
%     % M, #B, BC, Ch
%     X = [ 0 1 27.5 1 ; ...
%           1 1 27.5 1 ; ...
%           2 1 27.5 1 ; ...
%           0 2 27.5 1 ; ...
%           1 2 27.5 1 ; ...
%           2 2 27.5 1 ; ...
%           0 3 27.5 1 ; ...
%           1 3 27.5 1 ; ...
%           2 3 27.5 1 ; ...
%           ...
%           0 1 81.5 1 ; ...
%           1 1 81.5 1 ; ...
%           2 1 81.5 1 ; ...
%           0 2 81.5 1 ; ...
%           1 2 81.5 1 ; ...
%           2 2 81.5 1 ; ...
%           0 3 81.5 1 ; ...
%           1 3 81.5 1 ; ...
%           2 3 81.5 1 ; ...
%           ...
%           0 1 27.5 2 ; ...
%           1 1 27.5 2 ; ...
%           2 1 27.5 2 ; ...
%           0 2 27.5 2 ; ...
%           1 2 27.5 2 ; ...
%           2 2 27.5 2 ; ...
%           0 3 27.5 2 ; ...
%           1 3 27.5 2 ; ...
%           2 3 27.5 2 ; ...
%           ...
%           0 1 81.5 2 ; ...
%           1 1 81.5 2 ; ...
%           2 1 81.5 2 ; ...
%           0 2 81.5 2 ; ...
%           1 2 81.5 2 ; ...
%           2 2 81.5 2 ; ...
%           0 3 81.5 2 ; ...
%           1 3 81.5 2 ; ...
%           2 3 81.5 2 ];
% 
%     X(:,5) = X(:,1).*X(:,2);
%     X(:,6) = X(:,1).*X(:,3);
%     X(:,7) = X(:,1).*X(:,4);
% 
%     nVars = size(X,2);
%     dimCases = size(X,1);
% 
%     % zscore
%     for m=1:nVars,
%        XZ(:,m) = (X(:,m)-nanmean(X(:,m)))/nanstd(X(:,m)); 
%     end
% 
%     % store global vars.
%     XT(dimCases*(k-1)+1:dimCases*k, nVars*(k-1)+1:nVars*k) = XZ;
%     XT(dimCases*(k-1)+1:dimCases*k, nVars*numel(listSubjects)+k) = ones(dimCases,1); 
% 
%     % list variables
%     % my var --- PV1
%     for t=1:dimCases,
%        listIndices{t} = eval(sprintf('i%dx',t));
%     end
% 
%     for m=1:numel(listVariables),
% 
%         Y(:,m) = evaluateVariable(eval(listVariables{m}), listIndices);
%         YZ(:,m) = (Y(:,m) - nanmean(Y(:,m)))/nanstd(Y(:,m));
%         YT(dimCases*(k-1)+1:dimCases*k,m) = YZ(:,m);
%         % eliminate nans
%         i1 = find(~isnan(Y(:,m)));        
%         [betasPV1(:,k,m,r), dev, statsPV1(k)] = glmfit(XZ(i1,:), YZ(i1,m));
%         pVM(:,k,m,r) = statsPV1(k).p;
%         betaVM(:,k,m,r) = betasPV1(:,k,m,r);
% 
%     end
% 
%     clear('Y','X','YZ','XZ');
% %        clear('Y','X','YZ','XZ','X1','X2','X3','X4','X5','X45');
% 
% end
% 
% % assess significance - for each variable
% XT = [XT ones(size(XT,1),1)];
% g1 = figure;
% 
% for m=1:numel(listVariables),
% 
%     for n=1:size(betasPV1,1),
%         [h,pvPV1(n,m,r)] = ttest(betasPV1(n,:,m,r));
%       %  [pvPV1(n,m,r)] = signrank(betasPV1(n,:,m,r));
%     end
% end
% 
% % with bootstraping
% for m=1:numel(listVariables),
%     [myp,XPT] = bootstrapVar(XT,YT(:,m),squeeze(betasPV1(:,:,m,:)),listSubjects,r,'myVar');
% end
% 
% 
% for m=1:numel(listVariables),
% % global GLM regression and statistical analysis
%     i2 = find(~isnan(YT(:,m))); 
%     [betasG, devG, R, rint, statsG] = regress(YT(i2,m), [XT(i2,:) ones(numel(i2),1)]);
% 
%     for n=1:numel(listSubjects),
%         betasG2(:,n,m,r) = [ betasG(nVars*numel(listSubjects)+n) betasG(nVars*(n-1)+1:nVars*n)' ]';
%     end
% 
%     %betasG2(:,:,m,r) = reshape(betasG,12,numel(listSubjects));
%     % assess significance
%     for n=1:size(betasG2,1),
% %            [h,pvPV2(n,m,r)] = ttest(betasG2(n,:,m,r));
%         [pvPV2(n,m,r),H,STATS(m)] = signrank(betasG2(n,:,m,r));
%     end
% 
%     % plot figures;
%     betasAve(:,m,r) = nanmean(betasG2(:,:,m,r),2);
%     betasStd(:,m,r) = nanstd(betasG2(:,:,m,r),0,2)/sqrt(numel(listSubjects));
% %        subplot(4,2,m), barweb(betasAve(:,m,r), betasStd(:,m,r), [], [], [], [], [], bone, [], listCoefficients);%, 1, 'axis');
%     subplot(4,2,m), barweb(betasAve(:,m,r), betasStd(:,m,r), [], [], [], [], [], bone, [], []);%, 1, 'axis');
%     title(listVariables{m},'fontname','Arial','fontsize',10);
% 
%     % calc stats        i1 = find(pvPV1(:,m,r)<0.1);
%  %   pp2(:,:,m,r) = bootstrapVar(XT,YT(:,m),betasG2(:,:,m,r),listSubjects,nControl(r),listVariables{m});
% 
% end
% set(0,'CurrentFigure',g1); hold on;
% set(g1,'color',[1 1 1]);
% print(g1,'-dpdf','-r600',sprintf('betasGeneral2-nC-%d.pdf',nControl(r)));    
% 
% 
% pvPV2(2:end,:,:)
% pvPV1(2:end,:,2)
% %    squeeze(pp2(:,end,:,2))
% 
% % plot effects
% % motivation 1-3 - for vaariables 4-7
% for p=1:3,
%     Y1(:,p,1:4) = YT(p:3:end,[1 2 5 7]);
% end
% gx = figure;
% for q=1:4,
%     subplot(1,4,q), boxplot(Y1(:,:,q)); ylim([-2 2]);
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('MotivationEffect-nC-%d.pdf',nControl(r)));    
% close(gx);
% 
% % figure 1 data
% figuresizex = 9; figuresizey = 9;
% labelfontsize = 8; tickfontsize = 8;
% 
% % % create figure
% % [myaxis, tx, ty] = ignasiaxisfun(minx, maxx, miny, maxy);
% % ignasifig(0, mytitle, mysubjecttitle, 0, 0, 0, 0, figuresizex, figuresizey, labelfontsize, tickfontsize, myaxis, tx, ty, xlabel, ylabel);
% 
% % plot raw data (force vs gain)
% if (0)
% 
%     PA1 = Y1(:,:,1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y1(:,:,2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y1(:,:,3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
%     % punch-through
%     maxy=0.5;  minx = 0; 
%     miny=-0.5; maxx = 3.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y1(:,:,4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     listVarsToPlot= {'E'}; %listVariables = {'E'};
%     plotVariableGroupS;
%     emptyfigure;
% 
% end
% 
% % block number 1-3 - for variables 4-7
% i1 = sort([1:9:size(YT,1) 2:9:size(YT,1) 3:9:size(YT,1)]);
% for p=1:3,
%     Y2(:,p,1:4) = YT(i1+3*(p-1),[1 2 5 7]);
% end
% gx = figure;
% for q=1:4
%     subplot(1,4,q), boxplot(Y2(:,:,q)); ylim([-2 3]);
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('BlockEffect-nC-%d.pdf',nControl(r)));    
% close(gx);
% 
% % plot raw data (force vs gain)
% if (0)
% 
%     PA1 = Y2(:,:,1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y2(:,:,2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y2(:,:,3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
%     % punch-through
%     maxy=0.5;  minx = 0; 
%     miny=-0.5; maxx = 3.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y2(:,:,4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     listVarsToPlot= {'E'};
%     plotVariableGroupS;
%     emptyfigure;
% 
% end
% 
% % motivationx#T 1-6 - for variables 4-7
% i0 = [1:3:18];
% i1 = {i0,i0+1,i0+2,i0+18,i0+19,i0+20};
% for p=1:numel(listSubjects),
%     i1{1} = [ i1{1} i0+36*(p-1) ];
%     i1{2} = [ i1{2} i0+36*(p-1)+1 ];
%     i1{3} = [ i1{3} i0+36*(p-1)+2 ];
%     i1{4} = [ i1{4} i0+36*(p-1)+3 ];
%     i1{5} = [ i1{5} i0+36*(p-1)+4 ];
%     i1{6} = [ i1{6} i0+36*(p-1)+5 ];
% end
% for p=1:6,
%     Y3(:,p,1:4) = YT(i1{p},[1 2 5 7]);
% end
% gx = figure;
% for q=1:4,
%     subplot(1,4,q), boxplot(Y3(:,[1 4 2 5 3 6],q));
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('MotivTrial-nC-%d.pdf',nControl(r)));    
% close(gx);
% 
% % plot again
% if (0)
% 
%     PA1 = Y3(:,[1 4 2 5 3 6],1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y3(:,[1 4 2 5 3 6],2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y3(:,[1 4 2 5 3 6],3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PA2 = Y3(:,[1 4 2 5 3 6],1); PA2MG = nanmean(PA2,1); PA2SG = nanstd(PA2,0,1)/sqrt(10);
% %         PV2 = Y3(:,[1 4 2 5 3 6],2); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PD2 = Y3(:,[1 4 2 5 3 6],3); PD2MG = nanmean(PD2,1); PD2SG = nanstd(PD2,0,1)/sqrt(10);
%     % punch-through
%     maxy=.5;  minx = 0; 
%     miny=-.5; maxx = 6.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y3(:,[1 4 2 5 3 6],4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     listVarsToPlot= {'E'}; 
%     plotVariableGroupS;
%     emptyfigure;
% 
% end
% 
% % motivationx#B 1-9 - for variables 4-7
% i0 = [ 1 10 19 28 ];
% i1 = {  i0,i0+1,i0+2, ...   % block 1  - Motiv 1-3
%         i0+3,i0+4,i0+5, ... % block 2  
%         i0+6,i0+7,i0+8};    % block 3
% for p=1:numel(listSubjects),
%     i1{1} = [ i1{1} i0+36*(p-1) ];
%     i1{2} = [ i1{2} i0+36*(p-1)+1 ];
%     i1{3} = [ i1{3} i0+36*(p-1)+2 ];
%     i1{4} = [ i1{4} i0+36*(p-1)+3 ];
%     i1{5} = [ i1{5} i0+36*(p-1)+4 ];
%     i1{6} = [ i1{6} i0+36*(p-1)+5 ];
%     i1{7} = [ i1{7} i0+36*(p-1)+6 ];
%     i1{8} = [ i1{8} i0+36*(p-1)+7 ];
%     i1{9} = [ i1{9} i0+36*(p-1)+8 ];
% end
% for p=1:9,
%     Y4(:,p,1:4) = YT(i1{p},[1 2 5 7]);
% end
% gx = figure;
% for q=1:4,
%     subplot(1,4,q), boxplot(Y4(:,1:9,q));
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('MotivBlock-nC-%d.pdf',nControl(r)));    
% close(gx);
% 
% if (0)
% 
%     PA1 = Y4(:,:,1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y4(:,:,2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y4(:,:,3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
%     % punch-through
%     maxy=1.5;  minx = 0; 
%     miny=-1.5; maxx = size(PA1,2)+.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y4(:,:,4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     listVarsToPlot= {'E'}; 
%     plotVariableGroupS;
% 
%     emptyfigure;
% 
% end
% 
% % Biomechanics
% clear('i1');
% i0 = [ 1:18 ];
% i1{1} = i0; i1{2} = i0+19;
% for p=2:numel(listSubjects),
%     i1{1} = [ i1{1} i0+36*(p-1) ];
%     i1{2} = [ i1{2} i0+36*(p-1)+18 ];
% end
% for p=1:2,
%     Y5(:,p,1:4) = YT(i1{p},[1 2 5 7]);
% end
% gx = figure;
% for q=1:4,
%     subplot(1,4,q), boxplot(Y5(:,1:2,q),'datalim',[-1 1]); %ylim([-2 2]);
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('Biomechanics-nC-%d.pdf',nControl(r)));    
% close(gx);
% 
% if (0)
% 
%     PA1 = Y5(:,:,1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y5(:,:,2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y5(:,:,3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% 
% %         PA2 = Y5(:,:,1); PA2MG = nanmean(PA2,1); PA2SG = nanstd(PA2,0,1)/sqrt(10);
% %         PV2 = Y5(:,:,2); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PD2 = Y5(:,:,3); PD2MG = nanmean(PD2,1); PD2SG = nanstd(PD2,0,1)/sqrt(10);
%     % punch-through
%     maxy=1;  minx = 0; 
%     miny=-1; maxx = size(PA1,2)+.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y5(:,:,4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     maxy=.4;  minx = 0; 
%     miny=-.4; maxx = size(PA1,2)+.5;
%     listVarsToPlot= {'E'}; 
%     plotVariableGroupS;
% 
%     emptyfigure;
% 
% end
% 
% 
% %     gx2 = figure;
% %     for q=1:4,
% %         subplot(1,4,q), hh1 = plot(1:2, nanmean(Y5(:,:,q),1),'o','markersize',5,'color',[0 0 1],'linewidth',1.5); 
% %         hold on; hh = errorbar(1:2, nanmean(Y5(:,:,q),1), nanstd(Y5(:,:,q),0,1)/sqrt(10)); set(hh,'color',[0 0 1],'linewidth',1.5);
% %         title(sprintf(listVariables{q+3})); box off;
% %         xlim([0.5 2.5]);
% %     end
% %     box off
% %     set(gx2,'color',[1 1 1]);
% %     print(gx2,'-dpdf','-r600',sprintf('Biomechanics2-nC-%d.pdf',nControl(r)));    
% 
% 
% % #Trial
% clear('i1');
% i0 = [ 1:9 19:27 ] ;
% i1{1} = i0; i1{2} = i0+9;
% for p=2:numel(listSubjects),
%     i1{1} = [ i1{1} i0+36*(p-1) ];
%     i1{2} = [ i1{2} i0+36*(p-1)+9 ];
% end
% for p=1:2,
%     Y6(:,p,1:4) = YT(i1{p},[1 2 5 7]);
% end
% gx = figure;
% for q=1:4,
%     subplot(1,4,q), boxplot(Y6(:,1:2,q));
%     hold on; title(sprintf(listVariables{q+3})); box off;
% end
% box off
% set(gx,'color',[1 1 1]);
% print(gx,'-dpdf','-r600',sprintf('NTrial-nC-%d.pdf',nControl(r))); 
% close(gx);
% 
% if (0)
% 
%     PA1 = Y6(:,:,1); PA1MG = nanmean(PA1,1); PA1SG = nanstd(PA1,0,1)/sqrt(10);
%     PV1 = Y6(:,:,2); PV1MG = nanmean(PV1,1); PV1SG = nanstd(PV1,0,1)/sqrt(10);
%     PV2 = Y6(:,:,3); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PA2 = Y5(:,:,1); PA2MG = nanmean(PA2,1); PA2SG = nanstd(PA2,0,1)/sqrt(10);
% %         PV2 = Y5(:,:,2); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PD2 = Y5(:,:,3); PD2MG = nanmean(PD2,1); PD2SG = nanstd(PD2,0,1)/sqrt(10);
% %         PA2 = Y6(:,:,1); PA2MG = nanmean(PA2,1); PA2SG = nanstd(PA2,0,1)/sqrt(10);
% %         PV2 = Y6(:,:,2); PV2MG = nanmean(PV2,1); PV2SG = nanstd(PV2,0,1)/sqrt(10);
% %         PD2 = Y6(:,:,3); PD2MG = nanmean(PD2,1); PD2SG = nanstd(PD2,0,1)/sqrt(10);
%     % punch-through
%     maxy=.35;  minx = 0; 
%     miny=-.35; maxx = size(PA1,2)+.5;
%     listVarsToPlot = {'PA1','PV1','PV2'}; plotVariableGroupS;
% 
%     E = Y6(:,:,4); EMG = nanmean(E,1); ESG = nanstd(E,0,1)/sqrt(10);
%     figuresizex = 5; figuresizey = 4;
%     maxy=.4;  minx = 0; 
%     miny=-.4; maxx = size(PA1,2)+.5;
%     listVarsToPlot= {'E'}; 
%     plotVariableGroupS;
% 
%     emptyfigure;
% 
% end
% 
% %     gx2 = figure;
% %     for q=1:4,
% %         subplot(1,4,q), hh1 = plot(1:2, nanmean(Y6(:,:,q),1),'o','markersize',5,'color',[0 0 1],'linewidth',1.5); 
% %         hold on; hh = errorbar(1:2, nanmean(Y6(:,:,q),1), nanstd(Y6(:,:,q),0,1)/sqrt(10)); set(hh,'color',[0 0 1],'linewidth',1.5);
% %         title(sprintf(listVariables{q+3})); box off;
% %         xlim([0.5 2.5]);
% %     end
% %     box off
% %     set(gx2,'color',[1 1 1]);
% %     print(gx2,'-dpdf','-r600',sprintf('NTriaal2-nC-%d.pdf',nControl(r)));    
% 
