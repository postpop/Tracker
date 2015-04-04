classdef FlyPursuit2 < handle
   % FlyPursuit - clustering-based fly tracker
   
   % TODO
   %  - keep track of current frame number and save tracks directly to object
   %  - export flybox video
   %
   
   properties
      vr, vFileName
      nFlies, nSamp, initFrame
      w, h
      pos
      currentFrame, currentFrameIdx, NumberOfFrames
      medianFrame, meanFrame, stdFrame, madFrame
      foreGround, foreGroundOLD
      foreGroundThreshold
      maxDist
      temporalTrend, pixelNoiseLevel
      gmmStart, gmm
      arena, arenaCrop, arenaX, arenaY, boundsX, boundsY
      seNew, seOld, se10, se2, se5, H
      mu, sigma, area, orientation, tracks, pathLabels
      

   end
   
   methods (Access='public')
      function obj = FlyPursuit2(varargin)
         % get video reader
         obj.vFileName = varargin{1};
         try
            obj.vr = VideoReader(obj.vFileName);
         catch
            obj.vr = VideoReader(obj.vFileName);
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
         % set initial track values
         obj.currentFrameIdx = max(obj.initFrame-1,1);
         obj.mu(obj.currentFrameIdx,:,:) = obj.gmmStart.mu;
         obj.sigma(obj.currentFrameIdx,:,:,:) = obj.gmmStart.Sigma;
         obj.area(obj.currentFrameIdx,:) = obj.gmmStart.PComponents;
         obj.pathLabels(obj.currentFrameIdx,:) = 1:obj.nFlies;% seed initial 
         obj.foreGroundOLD = [];
      end
      
      
      %% __ TRACKER __
      function trackNextFrame(obj)
         %% get current frame
         obj.getNextFrame(); %increments frame counter and gets new obj.currentFrame
         %% check whether old frame exists
         if isempty(obj.foreGroundOLD)
            obj.foreGroundOLD = obj.getForeGround(obj.currentFrame);
            obj.foreGroundOLD(~obj.arenaCrop) = 0;
         else
            obj.foreGroundOLD = obj.foreGround;
         end
         %% get new forground
         oriFrame = obj.currentFrame;
         oriFrame(~obj.arenaCrop) = 0;
         [bw, obj.foreGround, clusterGroup, Nconn] = obj.getForeGroundAdaptive(oriFrame, obj.foreGrounOLD, obj.gmmStart.mu, obj.foreGroundThreshold);
         %%  cluster
         gmmNew = obj.clusterComponents(bw, obj.gmmStart, obj.foreGround, clusterGroup, Nconn);
         % save initial positions for next frame
         obj.mu(obj.currentFrameIdx,:,:) = gmmNew.mu;
         obj.sigma(obj.currentFrameIdx,:,:,:) = gmmNew.Sigma;
         % get each fly's orientation and area from eig analysis of their xcov
         for i = 1:obj.nFlies
            [eigVec, eigVal] = eig(gmmNew.Sigma(:,:,i));
            eigVal = diag(eigVal);
            obj.orientation(obj.currentFrameIdx,i,:) = eigVec(:, argmax(eigVal));
            obj.area(obj.currentFrameIdx,i) = pi*prod(sqrt(eigVal));
         end
         %% assign cluster centers to tracks
         newLabels = obj.assignCentroid2Track(obj.mu(max(f-1,1),:,:), obj.mu(f,:,:), obj.pathLabels(max(f-1,1),:), obj.maxDist);
         obj.pathLabels(obj.currentFrameIdx,:) = newLabels;% labels for each obj.centroid
         obj.tracks(obj.currentFrameIdx,obj.pathLabels(f,:),:) = squeeze(obj.mu(f,:,:));% actual obj.tracks
      end
         
      
      
      
      %% __ FUNCTIONS FOR GETTING FRAMES FROM VIDEO __
      function frames = getFrames(obj, varargin)
         % returns frame(s)
         % USAGE
         %  getFrames(idx); %get single frame
         %  getFrames(idx0,idx1)% get range of frames idx0:idx1
         %  getFrames([idx0:df:idx1])% get specified frames
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
         
         imagesc(frame);
         hold on
         plot(obj.arenaX, obj.arenaY,'k','LineWidth',3)
         drawnow
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
         % initialize gmm
         sigma = repmat(diag([10 10]),[1 1 obj.nFlies]);
         mu = [x(1:obj.nFlies),y(1:obj.nFlies)];
         p = ones(obj.nFlies,1)/obj.nFlies;
         obj.gmmStart = gmdistribution(mu, sigma, p);
      end
      
      function g = clusterFrame(obj, frame)
         % create samples from frame and custer points
         [X, Y] = randp(single(frame), obj.nSamp, 1);
         if isempty(obj.gmmStart)
            gmm = gmdistribution.fit([X;Y]',obj.nFlies,'Replicates',500);
         else
            g0.mu = obj.gmmStart.mu;
            g0.Sigma = obj.gmmStart.Sigma;
            g0.PComponents = obj.gmmStart.PComponents;
            try
               gmm = gmdistribution.fit([X;Y]',obj.nFlies, 'Start',g0);
            catch
               gmm = gmdistribution.fit([X;Y]',obj.nFlies, 'Replicates',500);
            end
         end
         g.mu = gmm.mu;
         g.Sigma =gmm.Sigma;
         g.PComponents =gmm.PComponents;
         
      end
      
      
      function gmm = clusterComponents(obj, bw, gmm0, frame, clusterGroup, Nconn)
         % find connected components and cluster each separately - esobj. useful for many flies, since
         % it divides one hard problem into a couple of easier ones
         mu = zeros(obj.nFlies,2);
         Sigma = zeros(2,2,obj.nFlies);
         cnt = 0;
         while any(mu(:)==0) && cnt<3
            for g = 1:Nconn % for each conn comp
               gFrame = frame;
               gFrame(bw~=g) = 0;% set pixels outside the group to zero
               
               gIdx = find(clusterGroup==g); % how many flies in the group?
               gNflies = length(gIdx);
               [X, Y] = randp(gFrame, gNflies*1000, 1);% get samples for centroids in the group
               if gNflies==0
                  disp(['region ' int2str(g) ' contains no flies.'])
               elseif gNflies==1 % if there's only one fly we simply calc the mean and cov of the points
                  mu(gIdx,:) = mean([X;Y],2);
                  Sigma(:,:,gIdx) = cov([X;Y]');
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
                  catch ME
                     disp(ME.getReport())
                     try % if an error occurs, try regularization
                        gmmTmp = gmdistribution.fit([X;Y]', gNflies, 'Start',gmmStart, 'Regularize',1);
                     catch ME2 % if that does not help, start from scratch by ignoring initial conditions (might be risky...)
                        disp(ME2.getReport())
                        gmmTmp = gmdistribution.fit([X;Y]', gNflies, 'Replicates',100);
                     end
                  end
                  % !!!
                  
                  % collect results over all connected components
                  mu(gIdx,:) = gmmTmobj.mu;
                  Sigma(:,:,gIdx) = gmmTmobj.Sigma;
               end
            end
            if any(mu(:)==0)
               disp('centroids at 0 - doing it again.')
               cnt = cnt+1;
            end
         end
         % superfluous?? return obj.gmmStart instead??
         gmm.mu = mu;
         gmm.Sigma = Sigma;
         gmm.PComponents = normalizeSum(gmm0.PComponents);
         % also save new positions in obj.gmmStart for next round
         obj.gmmStart = gmm;
      end
      
      function getBackGround(obj, nFrames)
         % USAGE: obj.getBackGround(nFrames)
         % estimate various statistics over a subset of frames (nFrames) including:
         %  - background frame (median and mean)
         %  - pixel noise (mad and std)
         %  - slow drift in overall luminance (temporalTrend)
         
         % framesToRead = randsample(obj.vr.NumberOfFrames, nFrames);
         framesToRead = linspace(obj.initFrame, obj.NumberOfFrames, nFrames);
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
         frame = medfilt2(frame,[3 3]);
         frame = imclose(imopen(frame,obj.se2),obj.se5);
         frame = single(frame);
         frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         frame = imgaussian(frame, 3, 30);    % smooth
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
         %          frame = imgaussian(frame, 3, 30);    % smooth
         %% v2
         %          frame = abs(frame - obj.medianFrame);         % distance to background
         %          frame = imclose(imopen(log(frame)>sense*nanmedian(log(frame(:))),obj.se2),obj.se5);
         %          frame = single(frame);
         %          frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         %          frame = imgaussian(frame, 3, 30);    % smooth
      end
      
      function [bw, frame, clusterGroup, Nconn] = getForeGroundAdaptive(obj, oriFrame, oldFrame, oldCentroids, sense,frameNumber)
         % background subtraction - treat everything that is not in the
         % neighbourhood of the fly as noise
         % if new foreground does not cover the position of the fly in the previous frame
         % sensitivity and neighbourhood are increased
         if nargin<5
            sense = 10;
         end
         if nargin<6
            frameNumber = 0;
         end
         senseOri = sense;
         lostFly = 1;
         cnt = 0;
         while ~isempty(lostFly) && cnt<50
            frame = obj.getForeGround(oriFrame, sense, frameNumber);
            % clean frame - remove specks that appear out of thin air...
            oldMask = imdilate(logical(oldFrame), obj.seOld);
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
         
         %% if we've lost all flies, perform simple
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
         if length(oldLabels)==1 % if we have just a single fly, skip
            newLabels = oldLabels;
         else
            C1 = squeeze(oldCentroid);
            C2 = squeeze(newCentroid);
            D = pdist2(C1, C2);
            D(D>maxDist) = 1000*D(D>maxDist);% penalize jumps - works only partly
            assign = munkres(D);
            [idx, ~] = find(assign);
            newLabels = oldLabels(idx);% labels for each obj.centroid
            if ~all(newLabels), newLabels = obj.nFlies; end
         end
      end
      
      function playTrack(obj, history, offset)
         if nargin==2
            offset = 1;
         end
         for t = offset:size(obj.tracks,1)-history
            plot(obj.tracks(t+(1:history),:,1), obj.tracks(t+(1:history),:,2))
            set(gca,'XLim',[0 obj.w], 'YLim',[0 obj.h])
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
      
   end
end
