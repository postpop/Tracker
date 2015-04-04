classdef TrackerFrameFactory < handle
   % TrackerFrameFactory for 
   
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
      w, h
      currentFrame, currentFrameOriginal, currentFrameIdx, NumberOfFrames
      colorChannels
      arena, arenaCrop, arenaX, arenaY, boundsX, boundsY
   end
   
   methods (Access='public')
      function obj = FlyPursuit(varargin)
         % get video reader
         obj.vFileName = varargin{1};
         try
            if isunix && ~ismac
               obj.vr = VideoReader2(obj.vFileName);
            else
               obj.vr = VideoReader(obj.vFileName);
            end
         catch
            obj.vr = VideoReader(obj.vFileName);
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

      end
      
  
      
      %% __ FUNCTIONS FOR GETTING FRAMES FROM VIDEO __
      function frames = getFrames(obj, varargin)
         % returns frame(s)
         % USAGE
         %  getFrames(idx); %get single frame
         %  getFrames(idx0,idx1)% get range of frames idx0:idx1
         %  getFrames([idx0:df:idx1])% get specified frames
         % 
         % also - called by all other 'getFrame' functions - so any change here will
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
      
   end
end
