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
      
      function drawArena(obj)
         %% draw arena and set initial positions
         frame = obj.getFrames(1);
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
         [X, Y] = randp(single(frame), obj.nSamp, 1);
         
         g0.mu = obj.gmmStart.mu;
         g0.Sigma = obj.gmmStart.Sigma;
         g0.PComponents = obj.gmmStart.PComponents;
         gmm = gmdistribution.fit([X;Y]',obj.nFlies, 'Start',g0);
         g.mu = gmm.mu;
         g.Sigma =gmm.Sigma;
         g.PComponents =gmm.PComponents;
         
      end
      
      function gmm = clusterBW(obj)
         se = strel('disk',10);
         mask = imdilate(obj.currentFrame,se)>0;% get mask from current frame
         [bw, Nconn] = bwlabel(mask);% find conn comps in mask
         clusterGroup = diag(bw(round(fp.gmmStart.mu(:,2)), round(fp.gmmStart.mu(:,1))));% assign cluster center to conn comps
         
         clear mu
         for g = 1:Nconn % for each conn comp
            gFrame = obj.currentFrame;
            gFrame(bw~=g) = 0;% set pixels outside the group to zero
            gIdx = find(clusterGroup==g); % how many flies in the group?
            gNflies = length(gIdx);
            [X, Y] = randp(gFrame, gNflies*200, 1);% get samples for centroids in the group
            hold on
            if gNflies==1% if there's only one fly we simply calc the mean and cov of the points
               mu(gIdx,:) = mean([X;Y]');
               Sigma(:,:,gIdx) = cov([X;Y]');
            else% otherwise we have to cluster
               gmmStart = gmdistribution(obj.gmmStart.mu(gIdx,:), ...
                                         obj.gmmStart.Sigma(:,:,gIdx),...
                                         normalizeSum(obj.gmmStart.PComponents(gIdx)));
               gmm = gmdistribution.fit([X;Y]',gNflies, 'Start',gmmStart);
               mu(gIdx,:) = gmm.mu;
               Sigma(:,:,gIdx) = gmm.Sigma;
            end
         end
         for i = 1:p.obj.nFlies
            eigVal = eig(Sigma(:,:,i));
            area(i) = pi*prod(sqrt(diag(eigVal)));
         end
         gmm = gmdistribution(mu, sigma, normalizeSum(area));
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
         if nargin==2
            sense = 5;
         end
         se = strel('disk',2);
         H = fspecial('gaussian',30,3);
         frame = single(frame);
         frame = -(frame - obj.medianFrame);           % distance to background
         frame = limit(frame,0, inf);                           % keep only pos deviations
         %frame = bsxfun(@gt, frame, 5*(madFrame));             % thres
         lum = mean(frame(:));
         frame = frame>lum*sense*(obj.madFrame);             % scale thres by mean luminance of the frame
         frame = single(frame);
         frame(~obj.arenaCrop(:)) = 0;                                 % delete everything outside arena
         frame = medfilt2(frame,[3 3]);                  % remove specks
         frame = imerode(frame, se);                      % erode...
         frame = imfilter(frame, H, 'replicate');          % smooth
      end
      
      
      
      function playTrack(obj, history)
         for t = 1:size(obj.tracks,1)-history
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
         set(hee, 'LineWidth', 2)
      end
      
      
   end
end
