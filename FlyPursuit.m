classdef FlyPursuit < handle
   
   properties
      vr, vFileName
      nFlies, nSamp, initFrame
      w, h
      pos
      currentFrame, previousFrame, frameIdx
      medianFrame, stdFrame, madFrame
      gmmStart, gmm
      arena, arenaCrop, arenaX, arenaY, boundsX, boundsY
      seNew, seOld, se10
      centroid, sigma, tracks, pathLabels
      
   end
   
   methods (Access='public')
      function obj = FlyPursuit(varargin)
         % get video reader
         obj.vFileName = varargin{1};
         obj.vr = VideoReader(obj.vFileName);
         obj.w = obj.vr.Width;
         obj.h = obj.vr.Height;
         % init arena to comprise full frame
         obj.arenaX = 1:obj.w;
         obj.arenaY = 1:obj.h;
         obj.boundsX = 1:obj.w;
         obj.boundsY = 1:obj.h;
         obj.arena = true(obj.w, obj.h);
         %
         obj.centroid = zeros(obj.vr.NumberOfFrames,obj.nFlies,2);
         obj.sigma = zeros(obj.vr.NumberOfFrames,2,2,obj.nFlies);
         obj.tracks = zeros(obj.vr.NumberOfFrames, obj.nFlies,2);
         obj.pathLabels = zeros(obj.vr.NumberOfFrames, obj.nFlies);% track label for each centroid/sigma
         obj.pathLabels(1,:) = 1:obj.nFlies;% seed initial labels
         %
         obj.seNew = strel('disk',10);
         obj.seOld = strel('disk',10);
         obj.se10 = strel('disk',10);
         
      end
      
      function frames = getFrames(obj, varargin)
         % getFrames(10), getFrames(10,20), getFrames([5:5:100])
         if length(varargin)==1
            framesToRead = varargin{1};
            frames = zeros(obj.vr.Height, obj.vr.Width, length(framesToRead), 'uint8');
            for f = 1:length(framesToRead);
               tmpFrame = obj.vr.read(framesToRead(f));
               if size(tmpFrame,3)>1
                  tmpFrame = mean(tmpFrame,3);
               end
               frames(:,:,f) = tmpFrame;
            end
         else
            frames = obj.vr.read([varargin{1} varargin{2}]);
            if length(size(frames))>4
               frames = squeeze(mean(frames,3));
            end
         end
         frames = obj.cropFrameToArenaBounds(frames);
      end
      
      function getFrame(obj, n)
         obj.currentFrame = obj.getFrames(n);
      end
      
      function getNextFrame(obj)
         disp('not implemented')
         %obj.previousFrame = obj.currentFrame;
         %obj.frameIdx = obj.frameIdx+1;
         % check for bounds!!
         %obj.currentFrame = obj.getFrames(obj.frameIdx);
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
         clf
         imagesc(frame);
         [~, x, y] = roipoly();
         obj.nFlies = size(x,1)-1;
         obj.pos = [x,y];
         % initialize gmm
         sigma = repmat(diag([10 10]),[1 1 obj.nFlies]);
         mu = [x(1:obj.nFlies),y(1:obj.nFlies)];
         p = ones(obj.nFlies,1)/obj.nFlies;
         obj.gmmStart = gmdistribution(mu, sigma, p);
      end
      
      function g = clusterFrame(obj, frame)
         % create samples from frame and custer points
         [X, Y] = randp(single(frame), obj.nSamp, 1);
         g0.mu = obj.gmmStart.mu;
         g0.Sigma = obj.gmmStart.Sigma;
         g0.PComponents = obj.gmmStart.PComponents;
         gmm = gmdistribution.fit([X;Y]',obj.nFlies, 'Start',g0);
         g.mu = gmm.mu;
         g.Sigma =gmm.Sigma;
         g.PComponents =gmm.PComponents;
         
      end
      
      
      function gmm = clusterComponents(obj, bw, gmm0, frame, clusterGroup, Nconn)
         % find connected components and cluster each separately - esp. uselfull for many flies, since
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
               %                [X, Y] = sampleGrid(gFrame);% get samples for centroids in the group
               %                X = X';
               %                Y = Y';
               
               hold on
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
                  % collect results over all connected components
                  mu(gIdx,:) = gmmTmp.mu;
                  Sigma(:,:,gIdx) = gmmTmp.Sigma;
                  
               end
            end
            if any(mu(:)==0)
               disp('centroids at 0 - doing it again.')
               cnt = cnt+1;
            end
         end
         gmm.mu = mu;
         gmm.Sigma = Sigma;
         
      end
      
      function getBackGround(obj, varargin)
         nFrames = varargin{1};
         framesToRead = randsample(obj.vr.NumberOfFrames, nFrames);
         frames = obj.getFrames(framesToRead);
         frames = reshape(single(frames), [], size(frames,3));
         obj.medianFrame = reshape(median(frames,2), obj.h, obj.w);
         obj.madFrame = reshape(mad(frames',1), obj.h, obj.w);
         obj.madFrame = medfilt2(obj.madFrame,[2 2]);
         obj.madFrame = limit(obj.madFrame,0,2);
         obj.madFrame(obj.madFrame(:)==0) = min(obj.madFrame(obj.madFrame>0));% avoid zero threshold
         obj.stdFrame = reshape(std(single(frames),[],2), obj.h, obj.w);
      end
      
      function frames = cropFrameToArenaBounds(obj, frames)
         frames = frames(obj.boundsY,obj.boundsX,:);
      end
      
      function frames = deleteFrameOutsideArena(obj, frames)
         for f = 1:size(frames,3)
            frames(obj.arena,f) = 0;
         end
      end
      
      function frame = getForeGround(obj, frame, sense)
         %frame = obj.getForeGround(frame, sense)
         if nargin==2
            sense = 5;
         end
         se = strel('disk',2);
         H = fspecial('gaussian',30,3);
         frame = single(frame);
         frame = -(frame - obj.medianFrame);         % distance to background
         frame = limit(frame,0, inf);                % keep only pos deviations
         lum = mean(frame(:));
         frame = frame>lum*sense*(obj.madFrame);     % scale thres by mean luminance of the frame
         
         frame = single(frame);
         frame(~obj.arenaCrop(:)) = 0;               % delete everything outside arena
         frame = medfilt2(frame,[3 3]);              % remove specks
         frame = imerode(frame, se);                 % erode...
         frame = imfilter(frame, H, 'replicate');    % smooth
         %%
      end
      
      function [bw, frame, clusterGroup, Nconn] = getForeGroundAdaptive(obj, oriFrame, oldFrame, oldCentroids, sense)
         if nargin<5
            sense = 10;
         end
         
         lostFly = 1;
         cnt = 0;
         while ~isempty(lostFly) && cnt<50
            frame = obj.getForeGround(oriFrame, sense);
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
            %             %gscatter(oldCentroids(:,1), oldCentroids(:,2), p.pathLabels(f-1,:),cols,[],[],'off')
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
      end
      
      function newLabels = assignCentroid2Track(obj, oldCentroid, newCentroid, oldLabels, maxDist)
         if nargin<5
            maxDist = 20;
         end
         C1 = squeeze(oldCentroid);
         C2 = squeeze(newCentroid);
         D = pdist2(C1, C2);
         D(D>maxDist) = 1000*D(D>maxDist);% penalize jumps - works only partly
         assign = munkres(D);
         [idx, ~] = find(assign);
         %%
         newLabels = oldLabels(idx);% labels for each p.centroid
         if ~all(newLabels), newLabels = obj.nFlies; end
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
      
      function plotCluster(obj, gmm)
         hold on
         plot(gmm.mu(:,1), gmm.mu(:,2), '.')
         for nf = 1:size(gmm.mu,1)
            hee(nf) = error_ellipse(gmm.Sigma(:,:,nf), gmm.mu(nf,:));
         end
         hold off
         set(hee, 'LineWidth', 2, 'Color','k')
      end
      
      
   end
end
