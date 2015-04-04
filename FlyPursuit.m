
classdef FlyPursuit < handle
   % FlyPursuit - clustering-based fly tracker
   
   % TODO
   %  - keep track of current frame number and save tracks directly to object - DONE
   %  - export flybox video, DONE in external function
   %  - integrate with Pip's tracking code
   %  - USE REPELLENT FILTER: imagesc(fspecial('gaussian',30, 4) - fspecial('gaussian',30,6)), colorbar
   %  - Kalman Filter for assignment:
   %     - maybe use 'future values', e.g. predict label at time T from pos
   %     and vel at T-3:T:3
   properties
      vr, vFileName
      vw
      nFlies, samplesPerFly, initFrame
      w, h
      pos
      currentFrame, currentFrameOriginal, currentFrameIdx, NumberOfFrames
      
      medianFrame, meanFrame, stdFrame, madFrame
      foreGround, foreGroundOLD
      foreGroundThreshold
      colorChannels
      maxDist
      saveFlyBox, flyBoxWidthIdx, flyBoxHeightIdx
      temporalTrend, pixelNoiseLevel
      gmmStart, gmm, negLogLikelihood
      arena, arenaCrop, arenaX, arenaY, boundsX, boundsY
      seNew, seOld, se10, se2, se5, H
      mu, sigma, area, orientation, tracks, pathLabels
   end
   
   methods (Access='public')
      function obj = FlyPursuit(varargin)
         % get video reader
         obj.vFileName = varargin{1};
         try
            if isunix && ~ismac
               obj.vr = VideoReader3(obj.vFileName);
            else
               obj.vr = VideoReader3(obj.vFileName);
            end
         catch
            obj.vr = VideoReader3(obj.vFileName);
         end
         
         try
            [~ ,~, obj.colorChannels] = size(obj.getFrame(1));
         end
         
         obj.w = obj.vr.Width;
         obj.h = obj.vr.Height;
         obj.NumberOfFrames = obj.vr.NumberOfFrames;
         obj.currentFrameIdx = 0;
         % init arena to comprise full frame
         obj.arenaX = [1 1 obj.w obj.w 1];
         obj.arenaY = [1 obj.h obj.h 1 1];
         obj.boundsX = 1:obj.w;
         obj.boundsY = 1:obj.h;
         obj.arena = true(obj.w, obj.h);
         %
         obj.mu = zeros(obj.vr.NumberOfFrames,obj.nFlies,2);
         obj.sigma = zeros(obj.vr.NumberOfFrames,2,2,obj.nFlies);
         obj.tracks = zeros(obj.vr.NumberOfFrames, obj.nFlies,2);
         obj.pathLabels = zeros(obj.vr.NumberOfFrames, obj.nFlies);% track label for each centroid/sigma
         obj.pathLabels(1,:) = 1:obj.nFlies;% seed initial labels
         obj.orientation = zeros(obj.vr.NumberOfFrames,obj.nFlies,2);
         obj.area = zeros(obj.vr.NumberOfFrames,obj.nFlies);
         
         obj.foreGroundThreshold = 5;
         obj.maxDist = 20;
         %
         obj.seNew = strel('disk',10);
         obj.seOld = strel('disk',10);
         obj.se10 = strel('disk',10);
         obj.H = fspecial('gaussian',30,3);
         obj.se2 = strel('disk',2);
         obj.se5 = strel('disk',5);
      end
      
      function initTracker(obj, initFrame, gmmStart)
         % ARGS:
         %  initFrame
         %  gmmStart - struct with mu, Sigma, Pcomponents;
         obj.initFrame = initFrame;
         obj.gmmStart = gmmStart;
         %
         obj.mu = zeros(obj.vr.NumberOfFrames,obj.nFlies,2);
         obj.sigma = zeros(obj.vr.NumberOfFrames,2,2,obj.nFlies);
         obj.tracks = zeros(obj.vr.NumberOfFrames, obj.nFlies,2);
         obj.pathLabels = zeros(obj.vr.NumberOfFrames, obj.nFlies);% track label for each centroid/sigma
         obj.orientation = zeros(obj.vr.NumberOfFrames,obj.nFlies,2);
         obj.area = zeros(obj.vr.NumberOfFrames,obj.nFlies);
         obj.negLogLikelihood = zeros(obj.vr.NumberOfFrames,1);
         
         obj.foreGroundThreshold = 5;
         obj.maxDist = 20;
         obj.se2 = strel('disk',2);
         obj.se5 = strel('disk',5);
         % set initial track values
         obj.currentFrameIdx = max(obj.initFrame-1,1);
         obj.pathLabels(obj.currentFrameIdx,:) = 1:obj.nFlies;% seed initial labels
         obj.pathLabels(obj.currentFrameIdx+1,:) = 1:obj.nFlies;% seed initial labels
         obj.mu(obj.currentFrameIdx,:,:) = obj.gmmStart.mu;
         obj.sigma(obj.currentFrameIdx,:,:,1:obj.nFlies) = obj.gmmStart.Sigma;
         obj.area(obj.currentFrameIdx,:) = obj.gmmStart.PComponents;
         obj.pathLabels(obj.currentFrameIdx,:) = 1:obj.nFlies;% seed initial
         
      end
      
      
      %% __ TRACKER __
      function trackNextFrame(obj, plotResults)
         if nargin==1
            plotResults = false;
         end
         %% get current frame
         obj.getNextFrame(); %increments frame counter and gets new obj.currentFrame
         if plotResults
            subplot(211)
            imagesc(obj.currentFrame)
            set(gca,'DataAspectRatio',[1 1 1])
         end
         %% get new foreground
         % clean-up frame
         cleanFrame = obj.currentFrame;
         cleanFrame(~obj.arenaCrop) = 0;
         % calc forground
         [bw, obj.foreGround, clusterGroup, Nconn] = obj.getForeGroundAdaptive(cleanFrame, obj.gmmStart.mu, obj.foreGroundThreshold);
         % save foreGround as sparse matrix
         
         %% cluster
         obj.gmmStart = obj.clusterComponents(bw, obj.gmmStart, obj.foreGround, clusterGroup, Nconn);
         %% track parameters
         obj.mu(obj.currentFrameIdx,:,:) = obj.gmmStart.mu;
         obj.sigma(obj.currentFrameIdx,:,:,:) = obj.gmmStart.Sigma;
         obj.negLogLikelihood(obj.currentFrameIdx) = obj.gmmStart.negLogLikelihood;
         % get each fly's orientation and area from eig analysis of their xcov
         for fly = 1:obj.nFlies
            [eigVec, eigVal] = eig(obj.gmmStart.Sigma(:,:,fly));
            eigVal = diag(eigVal);
            obj.orientation(obj.currentFrameIdx,fly,:) = eigVec(:, argmax(eigVal));
            obj.area(obj.currentFrameIdx,fly) = pi*prod(sqrt(eigVal));
         end
         %% assign cluster centers to tracks/fly IDs
         posCurr = obj.mu(obj.currentFrameIdx,:,:);
         posPrev = obj.mu(obj.currentFrameIdx-1,:,:);
         
         newLabels = obj.assignCentroid2Track(posPrev, posCurr, obj.pathLabels(obj.currentFrameIdx-1,:),obj.maxDist);
         obj.pathLabels(obj.currentFrameIdx,:) = newLabels;% labels for each obj.centroid
         obj.tracks(obj.currentFrameIdx,obj.pathLabels(obj.currentFrameIdx,:),:) = squeeze(obj.mu(obj.currentFrameIdx,:,:));% actual obj.tracks
         if plotResults
            subplot(212)
            imagesc(obj.foreGround)
            hold on
            flyX = obj.gmmStart.mu(:,1)';
            flyY = obj.gmmStart.mu(:,2)';
            gscatter(flyX, flyY, obj.pathLabels(obj.currentFrameIdx,:), limit(jet(obj.nFlies)-0.2),'.',18,'off');%,'.','MarkerSize',16);
            plot([flyX;flyX + 50*obj.orientation(obj.currentFrameIdx,:,1)], [flyY;flyY + 50*obj.orientation(obj.currentFrameIdx,:,2)])
            title(obj.currentFrameIdx);
            set(gca,'DataAspectRatio',[1 1 1])
         end
      end
      
      
      function trackNextFrameSimple(obj, plotResults)
         if nargin==1
            plotResults = false;
         end
         %% get current frame
         obj.getNextFrame(); %increments frame counter and gets new obj.currentFrame
         if plotResults
            subplot(211)
            imagesc(obj.currentFrame)
            set(gca,'DataAspectRatio',[1 1 1])
         end
         %% get new foreground
         % clean-up frame
         cleanFrame = obj.currentFrame;
         cleanFrame(~obj.arenaCrop) = 0;
         % calc forground
         [obj.foreGround] = obj.getForeGround(cleanFrame, obj.foreGroundThreshold);
         % save foreGround as sparse matrix?
         
         %% cluster
         obj.gmmStart = obj.clusterFrame(obj.foreGround);
         %% track parameters
         obj.mu(obj.currentFrameIdx,:,:) = obj.gmmStart.mu;
         obj.sigma(obj.currentFrameIdx,:,:,:) = obj.gmmStart.Sigma;
         obj.negLogLikelihood(obj.currentFrameIdx) = obj.gmmStart.negLogLikelihood;
         % get each fly's orientation and area from eig analysis of their xcov
         for fly = 1:obj.nFlies
            [eigVec, eigVal] = eig(obj.gmmStart.Sigma(:,:,fly));
            eigVal = diag(eigVal);
            obj.orientation(obj.currentFrameIdx,fly,:) = eigVec(:, argmax(eigVal));
            obj.area(obj.currentFrameIdx,fly) = pi*prod(sqrt(eigVal));
         end
         %% assign cluster centers to tracks/fly IDs
         posCurr = obj.mu(obj.currentFrameIdx,:,:);
         posPrev = obj.mu(obj.currentFrameIdx-1,:,:);
         
         newLabels = obj.assignCentroid2Track(posPrev, posCurr, obj.pathLabels(obj.currentFrameIdx-1,:),obj.maxDist);
         obj.pathLabels(obj.currentFrameIdx,:) = newLabels;% labels for each obj.centroid
         obj.tracks(obj.currentFrameIdx,obj.pathLabels(obj.currentFrameIdx,:),:) = squeeze(obj.mu(obj.currentFrameIdx,:,:));% actual obj.tracks
         if plotResults
            subplot(212)
            imagesc(obj.foreGround)
            hold on
            flyX = obj.gmmStart.mu(:,1)';
            flyY = obj.gmmStart.mu(:,2)';
            gscatter(flyX, flyY, obj.pathLabels(obj.currentFrameIdx,:), limit(jet(obj.nFlies)-0.2),'.',18,'off');%,'.','MarkerSize',16);
            plot([flyX;flyX + 50*obj.orientation(obj.currentFrameIdx,:,1)], [flyY;flyY + 50*obj.orientation(obj.currentFrameIdx,:,2)])
            title(obj.currentFrameIdx);
            set(gca,'DataAspectRatio',[1 1 1])
         end
      end
      
      
      %% __ FUNCTIONS FOR GETTING FRAMES FROM VIDEO __
      function frames = getFrames(obj, varargin)
         % returns frame(s)
         % USAGE
         %  getFrames(idx); %get single frame
         %  getFrames(idx0,idx1)% get range of frames idx0:idx1
         %  getFrames([idx0:df:idx1])% get specified frames
         %
         % also - called by all other 'getFrame' functions - any change here will
         % be reflected in those, too.
         if length(varargin)==1
            framesToRead = varargin{1};
         elseif length(varargin)==2
            framesToRead = varargin{1}:varargin{2};
         end
         frames = zeros(obj.vr.Height, obj.vr.Width, length(framesToRead), 'uint8');
         for f = 1:length(framesToRead);
            tmpFrame = obj.vr.read(framesToRead(f));
            % make monochromatic
            if size(tmpFrame,3)>1
               tmpFrame = mean(tmpFrame,3);
            end
            frames(:,:,f) = tmpFrame;
         end
         obj.currentFrameOriginal = frames;
         frames = obj.cropFrameToArenaBounds(frames);
      end
      
      function frame = getFrame(obj, n)
         % returns single frame
         frame = obj.getFrames(n);
      end
      
      function varargout = getNextFrame(obj)
         % increments currentFrameIdx and loads new current frame
         % USAGE:
         %   obj.getNextFrame();% sets new obj.currentFrame
         %   currentFrame = obj.getNextFrame(); % also returns obj.currentFrame
         obj.currentFrameIdx = obj.currentFrameIdx + 1;
         obj.currentFrame = obj.getFrames(obj.currentFrameIdx);
         if nargout>0
            varargout{1} = obj.currentFrame;
         end
      end
      
      %% __ FUNCTIONS FOR INITIALIZING FLY POSITIONS __
      function drawArenaPoly(obj, frameNumber, radius)
         if nargin==1
            frameNumber = 1;
         end
         if nargin<=2
            radius = 1;
         end
         frame = obj.getFrames(frameNumber);
         
         [obj.arenaX, obj.arenaY] = ellipsePoints(obj.h./2*radius, obj.w/2*radius, 90,...
            obj.h./2, obj.w/2,12);
         obj.arena = poly2mask(obj.arenaX, obj.arenaY, obj.w, obj.h);
         obj.boundsY = limit(round(min(obj.arenaY):max(obj.arenaY)),1, obj.h);
         obj.boundsX = limit(round(min(obj.arenaX):max(obj.arenaX)),1, obj.w);
         obj.arenaCrop = obj.arena(obj.boundsY, obj.boundsX);% crop arena
         obj.w = length(obj.boundsX);
         obj.h = length(obj.boundsY);
         
         %          imagesc(frame);
         %          hold on
         %          plot(obj.arenaX, obj.arenaY,'k','LineWidth',3)
         %          drawnow
      end
      
      function drawArena(obj, frameNumber)
         %% draw arena and set initial positions
         if nargin==1
            frameNumber = 1;
         end
         frame = obj.getFrames(frameNumber);
         clf
         imagesc(frame);
         disp('please draw arena...')
         [obj.arena, obj.arenaX, obj.arenaY] = roipoly();
         obj.boundsY = limit(round(min(obj.arenaY):max(obj.arenaY)),1, obj.h);
         obj.boundsX = limit(round(min(obj.arenaX):max(obj.arenaX)),1, obj.w);
         obj.arenaCrop = obj.arena(obj.boundsY, obj.boundsX);% crop arena
         obj.w = length(obj.boundsX);
         obj.h = length(obj.boundsY);
         hold on
         plot(obj.arenaX, obj.arenaY)
         drawnow
      end
      
      function initFlies(obj,n)
         frame = obj.getFrames(n);
         imagesc(frame);
         hold on
         plot(obj.arenaX, obj.arenaY,'LineWidth',3)
         [~, x, y] = roipoly();
         obj.nFlies = max(1,size(x,1)-1);
         obj.pos = [x,y];
         obj.pos = unique(obj.pos, 'rows');
         % initialize gmm with positions and arbitrary guesses of fly size and orientation
         obj.gmmStart.Sigma = repmat(diag([10 10]),[1 1 obj.nFlies]);
         obj.gmmStart.mu = [x(1:obj.nFlies),y(1:obj.nFlies)];
         obj.gmmStart.PComponents = ones(obj.nFlies,1)/obj.nFlies;
      end
      
      %% __ FUNCTIONS FOR CLUSTERING FLIES __
      function gm = clusterFrame(obj, frame)
         % create samples from frame, cluster points, return structure MU, SIGMA, PCOMPONENTS
         [X, Y] = randp(single(frame), obj.nFlies*obj.samplesPerFly, 1);
         if isempty(obj.gmmStart)
            gmDist = gmdistribution.fit([X;Y]',obj.nFlies,'Replicates',500);
            [~, negLogLikelihood] = cluster(gmDist, [X;Y]');
         else
            try
               gmDist = gmdistribution.fit([X;Y]',obj.nFlies, 'Start',obj.gmmStart);
               [~, negLogLikelihood] = cluster(gmDist, [X;Y]');
            catch
               disp('   clustering failed - trying to recover...')
               gmDist= gmdistribution.fit([X;Y]',obj.nFlies, 'Replicates',500);
               [~, negLogLikelihood] = cluster(gmDist, [X;Y]');
            disp('     done.')
            end
            if obj.currentFrameIdx>obj.initFrame+100
               % try to recover from jumps
               historyFrameIdx = limit(obj.currentFrameIdx - 1 - (1:100),obj.initFrame+1,obj.NumberOfFrames);
               thres = mean(obj.negLogLikelihood(historyFrameIdx)) + 2.5*std(obj.negLogLikelihood(historyFrameIdx));
               if thres>0 && negLogLikelihood > thres
                  try
                     disp('   detected drop in cluster quality. trying to recover...')
                     gmDistR = gmdistribution.fit([X;Y]',obj.nFlies, 'Replicates',500);
                     [~, negLogLikelihoodR] = cluster(gmDist, [X;Y]');

                     if negLogLikelihoodR<negLogLikelihood
                        disp('      ll improved.')
                        gmDist = gmDistR;
                        negLogLikelihood = negLogLikelihoodR;
                     else
                        disp('      ll NOT improved.')
                     end
                  catch ME
                     disp(ME.getReport())
                  end
                  % check whether we improved
               end
            end
         end
         % extract fields from GM object
         gm.mu = gmDist.mu;
         gm.Sigma = gmDist.Sigma;
         gm.PComponents = gmDist.PComponents;
         gm.negLogLikelihood = negLogLikelihood;
      end
      
      function gmm = clusterConnComp(obj, frame)
         % 1 if number of conn comp==nFLies
         % take  centers as fly positions
         
         % 2. if n conn comp < nFlies
         % attempt to split lagest con comp
         
         % 3 if n conn cop>nFLies
         % use old fly positions to get rid of extraneous ones
      end
      
      function gmm = clusterComponents(obj, bw, gmm0, frame, clusterGroup, Nconn)
         % find connected components and cluster each separately - esp. useful for many flies, since
         % it divides one hard problem into a couple of easier ones
         gmm.mu = zeros(obj.nFlies,2);
         gmm.Sigma = zeros(2,2,obj.nFlies);
               
         cnt = 0;
         while any(gmm.mu(:)==0) && cnt<3
            for g = 1:Nconn % for each conn comp
               gFrame = frame;
               gFrame(bw~=g) = 0;% set pixels outside the group to zero
               
               gIdx = find(clusterGroup==g); % how many flies in the group?
               gNflies = length(gIdx);
               [X, Y] = randp(gFrame, gNflies*obj.samplesPerFly, 1);% get samples for centroids in the group
               if gNflies==0
                  %disp(['region ' int2str(g) ' contains no flies.'])
               elseif gNflies==1 % if there's only one fly we simply calc the mean and cov of the points
                  gmm.mu(gIdx,:) = mean([X;Y],2);
                  gmm.Sigma(:,:,gIdx) = cov([X;Y]');
                  gmm.negLogLikelihood(gIdx) = 0;
                  %                   S=ellipsefit(X,Y);
                  %                   mu(gIdx,:) = [S.Xc, S.Yc];
                  %                   cov([X;Y]')
                  %                   ss = rotateVec2D([S.A.*[1 0]'], S.Phi)
                  %                   ss = rotateVec2D([S.B.*[1 0]'], pi/2+S.Phi)
                  
                  %bwProps = regionprops(logical(bw),frame,{'WeightedCentroid'});
                  %wcent = reshape([bwProps.WeightedCentroid],2,[]);
               else % otherwise we cluster
                  % get initial conditions from previous positions for the subset
                  % of flies in the current region
                  
                  
                  % !!!
                  % IS THIS BLOCK SAME AS clusterFrame???
                  % can we call: gmmTmp = obj.clusterFrame(gFrame) instead??
                  gmmStart.mu = gmm0.mu(gIdx,:);
                  gmmStart.Sigma = gmm0.Sigma(:,:,gIdx);
                  gmmStart.PComponents = normalizeSum(gmm0.PComponents(gIdx));
                  % cluster
                  try
                     gmmTmp = gmdistribution.fit([X;Y]', gNflies, 'Start',gmmStart);
                     [~, negLogLikelihood] = cluster(gmmTmp, [X;Y]');
                  catch ME
                     disp(ME.getReport())
                     try % if an error occurs, try regularization
                        gmmTmp = gmdistribution.fit([X;Y]', gNflies, 'Start',gmmStart, 'Regularize',1);
                        [~, negLogLikelihood] = cluster(gmmTmp, [X;Y]');

                     catch ME2 % if that does not help, start from scratch by ignoring initial conditions (might be risky...)
                        disp(ME2.getReport())
                        gmmTmp = gmdistribution.fit([X;Y]', gNflies, 'Replicates',100);
                        [~, negLogLikelihood] = cluster(gmmTmp, [X;Y]');
                     end
                  end
                  % !!!
                  
                  % collect results over all connected components
                  gmm.mu(gIdx,:) = gmmTmp.mu;
                  gmm.Sigma(:,:,gIdx) = gmmTmp.Sigma;
                  gmm.negLogLikelihood(gIdx) = negLogLikelihood;

               end
            end
            if any(gmm.mu(:)==0)
               disp('centroids at 0 - doing it again.')
               cnt = cnt+1;
            end
         end
         gmm.PComponents = normalizeSum(gmm0.PComponents);
         gmm.negLogLikelihood = mean(gmm.negLogLikelihood);

      end
      
      function getBackGround(obj, nFrames)
         % USAGE: obj.getBackGround(nFrames)
         % estimate various statistics over a subset of frames (nFrames) including:
         %  - background frame (median and mean)
         %  - pixel noise (mad and std)
         %  - slow drift in overall luminance (temporalTrend)
         
         % framesToRead = randsample(obj.vr.NumberOfFrames, nFrames);
         framesToRead = round(linspace(obj.initFrame+10, obj.NumberOfFrames-10, nFrames));
         frames = obj.getFrames(framesToRead);
         frames = reshape(single(frames), [], size(frames,3));
         % background frame (median more robust)
         obj.medianFrame = reshape(median(frames,2), obj.h, obj.w);
         obj.meanFrame = reshape(mean(frames,2), obj.h, obj.w);
         % detect slow drift in luminance over time
         obj.temporalTrend = single(smooth(median(frames,1),100) - median(obj.medianFrame(:)));
         % remove trend from sample frames
         frames = bsxfun(@minus, frames, obj.temporalTrend');
         % calc detrended background frame (median more robust)
         obj.medianFrame = reshape(median(frames,2), obj.h, obj.w);
         % mad and std frames
         obj.stdFrame = reshape(std(single(frames),[],2), obj.h, obj.w);
         obj.madFrame = reshape(mad(frames',1), obj.h, obj.w);
         obj.madFrame = medfilt2(obj.madFrame,[2 2]);
         %obj.madFrame = limit(obj.madFrame,0,2);
         obj.madFrame(obj.madFrame(:)==0) = min(obj.madFrame(obj.madFrame>0));% avoid zero threshold
         % estimate pixel noise
         obj.pixelNoiseLevel = mean(obj.madFrame(:));%mean(std(frames));
         % interpolated temporalTrend for each frame
         obj.temporalTrend = interp1(framesToRead, obj.temporalTrend, obj.initFrame:obj.NumberOfFrames);
         
      end
      
      function frames = cropFrameToArenaBounds(obj, frames)
         frames = frames(obj.boundsY,obj.boundsX,:);
      end
      
      function frames = deleteFrameOutsideArena(obj, frames)
         if length(size(frames))==2
            frames(obj.arenaCrop) = 0;
         else
            for f = 1:size(frames,3)
               frames(obj.arenaCrop,f) = 0;
            end
         end
      end
      
      function frame = getForeGround(obj, frame, sense, frameNumber)
         % background subtraction
         % USAGE: frame = obj.getForeGround(frame, [sense=5], [frameNumber])
         % if frameNumber is provided, will correct for temporal drift in luminance
         if nargin==2 % set default sensitivity
            sense = 5;
         end
         if nargin>3 && frameNumber>0 %correct for temporal trend
            temporalTrend = obj.temporalTrend(frameNumber - obj.initFrame);
         else
            temporalTrend = 0;
         end
         frame = single(frame);
         frame = frame - temporalTrend;
         %% v3
         frame = abs(frame - obj.medianFrame);
         frame = frame>sense*obj.pixelNoiseLevel;
         if ~isempty(obj.arenaCrop)
            frame = frame.*obj.arenaCrop;               % delete everything outside arena
         end
         frame = medfilt2(frame,[3 3]);
         
         frame = imclose(imopen(frame,obj.se2),obj.se5);
         
         frame = single(frame);
         frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         frame = filter2(obj.H, frame);    % smooth - filter2 is faster that imfilter or conv2 and as fast as imgaussian - stick with filter2 to remove external dependency
         
         %% v1
         %          frame = (frame - obj.medianFrame);         % distance to background
         %          frame = limit(frame,0, inf);                % keep only pos deviations
         %          frame = imclose(frame, obj.se5);
         %          lum = mean(frame(:));
         %          frame = frame>lum*sense*(obj.madFrame);     % scale thres by mean luminance of the frame
         %
         %          frame = single(frame);
         %          frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         %          frame = medfilt2(frame,[3 3]);              % remove specks
         %          frame = imerode(frame, se);                 % erode...
         %          frame = medfilt2(frame,[3 3]);              % remove specks
         %          frame = filter2(obj.H, frame);    % smooth
         %% v2
         %          frame = abs(frame - obj.medianFrame);         % distance to background
         %          frame = imclose(imopen(log(frame)>sense*nanmedian(log(frame(:))),obj.se2),obj.se5);
         %          frame = single(frame);
         %          frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         %          frame = filter2(obj.H, frame);    % smooth
      end
      
      function [bw, frame, clusterGroup, Nconn] = getForeGroundAdaptive(obj, oriFrame, oldCentroids, sense,frameNumber)
         % background subtraction - treat everything that is not in the
         % neighbourhood of the fly as noise
         % if new foreground does not cover the position of the fly in the previous frame
         % sensitivity and neighbourhood are increased
         if nargin<4
            sense = 10;
         end
         if nargin<5
            frameNumber = 0;
         end
         senseOri = sense;
         lostFly = 1;
         cnt = 0;
         while ~isempty(lostFly) && cnt<50
            frame = obj.getForeGround(oriFrame, sense, frameNumber);
            % clean frame - remove specks that appear far away from old fly positions
            oldMask = false(size(frame));
            flyPosY = limit(round(oldCentroids(:,1)),1,obj.w);
            flyPosX = limit(round(oldCentroids(:,2)),1,obj.h);
            oldMask(sub2ind(size(oldMask), flyPosX, flyPosY)) = true;
            oldMask = imdilate(logical(oldMask), obj.seOld);% grow positions
            frame(~oldMask) = 0;% remove all pixels that are far away from the flies' previous positions
            % find isolated groups of flies
            mask = imdilate(logical(frame), obj.seNew);% get mask from current frame
            [bw, Nconn] = bwlabel(mask);% find conn comps in mask
            clusterGroup = diag(bw(...
               limit(round(oldCentroids(:,2)),1,max(obj.boundsY)),...
               limit(round(oldCentroids(:,1)),1,max(obj.boundsX))));% assign cluster center to conn comps
            lostFly = find(clusterGroup==0);% who did we loose?
            if ~isempty(lostFly)
               disp(['   we lost the following members of the crew: ' mat2str(lostFly') '.'])
               sense = 0.9*sense;
               disp(['   going back with lower threshold (new thres = ' num2str(sense) ').'])
               cnt = cnt+1;
            end
            %             subplot(221)
            %             cla
            %             imagesc(oriFrame)
            %             hold on
            %             %gscatter(oldCentroids(:,1), oldCentroids(:,2), obj.pathLabels(f-1,:),cols,[],[],'off')
            %             plot(oldCentroids(:,1), oldCentroids(:,2), '.k')
            %             hold off
            %             subplot(222)
            %             cla
            %             imagesc(mask-frame)
            %             hold on
            %             plot(oldCentroids(:,1), oldCentroids(:,2), '.r')
            %             hold off
            %             subplot(223)
            %             cla
            %             imagesc(bwlabel(mask))
            %             hold on
            %             plot(oldCentroids(:,1), oldCentroids(:,2), 'xk')
            %             drawnow
         end
         
         %% if we've lost all flies, perform simple background subtraction
         if length(lostFly)==obj.nFlies
            disp('we have lost all flies - starting from scratch')
            frame = obj.getForeGround(oriFrame,senseOri,frameNumber);
         end
      end
      
      function newLabels = assignCentroid2Track(obj, oldCentroid, newCentroid, oldLabels, maxDist)
         % assigns cluster centers to fly ids
         if nargin<5
            maxDist = 20;
         end
         if length(oldLabels)==1                               % if we have just a single fly, skip
            newLabels = oldLabels;
         else
            C1 = squeeze(oldCentroid);
            C2 = squeeze(newCentroid);
            D = pdist2(C1, C2);
            D(D>maxDist) = 1000*D(D>maxDist);                  % penalize jumps - works only partly
            assign = munkres(D);                               % calculates optimal assignment of centroids to paths
            [idx, ~] = find(assign);
            newLabels = oldLabels(idx);                        % labels for each obj.centroid
            if ~all(newLabels), newLabels = 1:obj.nFlies; end    % fallback in case everything goes wrong
         end
         if any(newLabels==20)
            disp(newLabels)
         end
      end
      
      function playTrack(obj, history, offset, plotFrame)
         if nargin==2
            offset = 1;
         end
         
         if nargin==3
            plotFrame = false;
         end
         for t = offset:size(obj.tracks,1)-history
            if plotFrame
               imagesc(obj.getFrame(t));
            end
            hold on
            plot(obj.tracks(t-(1:history),:,1), obj.tracks(t-(1:history),:,2))
            flyX = obj.tracks(t,:,1);
            flyY = obj.tracks(t,:,2);
            plot([flyX;flyX + 50*obj.orientation(t,:,1)], [flyY;flyY + 50*obj.orientation(t,:,2)])
            
            hold off
            title(t)
            set(gca,'XLim',[0 obj.w], 'YLim',[0 obj.h])
            set(gca,'DataAspectRatio',[1 1 1])
            drawnow
         end
      end
      
      function plotCluster(obj, gmm) %#ok<INUSL>
         hold on
         plot(gmm.mu(:,1), gmm.mu(:,2), '.')
         for nf = 1:size(gmm.mu,1)
            hee(nf) = error_ellipse(gmm.Sigma(:,:,nf), gmm.mu(nf,:));
         end
         hold off
         set(hee, 'LineWidth', 2, 'Color','r')
      end
      
      function flySpeed = getSpeed(obj)
         % get fly speeds
         flySpeed = sqrt(sum(diff(obj.tracks,1).^2,3));
      end
      
      function flyDist = getDist(obj)
         % get pairwise distances between flies
         flyCnt = 0;
         flyDist = zeros(size(obj.tracks,1), obj.nFlies^2/2 - obj.nFlies/2);
         for fly1 = 1:obj.nFlies
            for fly2 = fly1+1:obj.nFlies
               flyCnt = flyCnt+1;
               flyDist(:, flyCnt) = squeeze(sqrt(sum((obj.tracks(:,fly1,:) - obj.tracks(:,fly2,:)).^2,3)));
            end
         end
      end
      
      function flyAng = getAngle(obj)
         flyAng = zeros(size(obj.tracks,1), obj.nFlies^2/2 - obj.nFlies/2);
         flyCnt = 0;
         for fly1 = 1:obj.nFlies
            for fly2 = fly1+1:obj.nFlies
               flyCnt = flyCnt+1;
               flyAng(:, flyCnt) = 0;%squeeze(sqrt(sum((obj.tracks(:,fly1,:) - obj.tracks(:,fly2,:)).^2,3)));
            end
         end
      end
      
   end
end
