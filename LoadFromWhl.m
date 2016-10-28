function motion = LoadFromWhl(baseName, tos,lednum, varargin)
% Copied from StateEditor, by Andres Grosmark, 2012
% lednum specifies which LED you want to use... ie LED one uses whl columns
% 1&2...

colnums = [];
for a = 1:length(lednum);
    temp = lednum(a);
    colnums = [colnums temp*2-1 temp*2];
end
    
if length(varargin) >= 1
    whl = varargin{1};
else
    [whl,to,GoodRanges] = LoadFromWhlHelper1(baseName);
end
if size(whl,2) == 2 %if two LEDs then take the one with fewest NaN's, then fill in NaN's from other light if neccesary
    d = [abs(diff(whl(:, 1))) + abs(diff(whl(:, 2))), abs(diff(whl(:, 3))) + abs(diff(whl(:, 4)))];
    [n, mGood] = find(min(mean(isnan(d))));
    if mGood == 1
        mBad = 2;
        
    else
        mBad = 1;
    end
    motion1 = d(:, mGood);
    motion1(isnan(d(:, mGood))) = d(isnan(d(:, mGood)), mBad);
else
    motion1 = abs(diff(whl));
end
f = find(isnan(motion1));
motion = motion1;

for i = 1:length(f) %%%this looks 3 seconds before and after each NaN value to find a non-NaN estimate (mean of non-NaN neighbors)
    h = [f(i) - 39*3, f(i) + 39*3];
    if h(1) < 1
        h(1) = 1;
    else
        if h(2) > length(motion1)
            h(2) = length(motion1);
        end
    end
    c = motion1(h(1):h(2));
    motion(f(i)) = mean(c(~isnan(c)));
end

to2 = to(2:end) - (diff(to)/2);

motion2 = motion;
motion3 = motion2(:,colnums);
motion4 = sum(motion3,2);

motion = [];
% for i = 1:length(tos) %avearge over same bins as spectrogram
%     motion = [motion, mean(motion2(to2 > tos(i) & to2 <= (tos(i) + 1)))];
% end
for i = 1:length(tos) %avearge over same bins as spectrogram
    motion = [motion, mean(motion4(to2 > tos(i) & to2 <= (tos(i) + 1)))];
end


end

function [whl,t,GoodRanges] = LoadFromWhlHelper1(fbasename)
% USAGE
% [whl,t,GoodRanges,ep] = LoadPosition(fbasename)
%
% output:
%   whl: the position matrix
%   t: time vector
%   Good Ranges:

Fs = 1250/32;

whlt = dlmread([fbasename '.whl']);
[whl GoodRanges] = LoadFromWhlHelper2(whlt);

t = (1:size(whlt,1))'/Fs;
GoodRanges = GoodRanges/Fs;
end
function  [cWhl, GoodRanges_F] = LoadFromWhlHelper2(Whl, StretchLen, JumpSize, Gap)

% If the gap between the good strech is more than StrethcLen in terms of Whl row number,remove interporated values.
if nargin<2
    StretchLen = 30;
end

% If the Gap between the good strech is more than JumpSize, remove interporated values.
if nargin<3
    JumpSize = 30;
end

% if the distance between the two contimous rows are more than Gap centimeter, It's a big jump and do not use as an input for inpterp1.
if nargin<4,
    Gap = 30;
end

nWhl = size(Whl,1);

% interpolate missing values or large jumps.
% the value of whl(:,3:4) is also taken into account for the range of interpolation.
% A transision to and form (-1,-1) should be taken as a BigJump.

% I should use distance, not the one dimentinal projection of trajectory, by the way.

whltemp = Whl;
whltemp(find(whltemp)==-1) = -Gap;
dist_F = sqrt(diff(whltemp(:,1)).^2+diff(whltemp(:,2)).^2);
dist_R = sqrt(diff(whltemp(:,3)).^2+diff(whltemp(:,4)).^2);
BigJump_F = dist_F>Gap;
BigJump_R = dist_R>Gap;

Good_F = find(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F]));
Bad_F = find(~(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F])));
Good_R = find(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R]));
Bad_R = find(~(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R])));

whltemp(Bad_F,1:2) = -Gap;
whltemp(Bad_R,3:4) = -Gap;

WhlNaN = Whl;
WhlNaN(find(Whl==-1)) = NaN;

% Give -1 outside of the interpolation.

if length(Good_F)<2 || length(Good_R)<2;
    cWhl(:,1:2) = -ones(size(Whl,1),2);
else
    cWhl(:,1:2) = interp1(Good_F, Whl(Good_F,1:2), 1:nWhl, 'linear', -1);
    cWhl(:,3:4) = interp1(Good_R, Whl(Good_R,3:4), 1:nWhl, 'linear', -1);
end


% find missing stretches for Front LED
dGoodF = [-(whltemp(1,1)==-Gap) ; diff(whltemp(:,1)>-Gap)];
BadStartF = find(dGoodF<0);
BadEndF = find(dGoodF>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
    BadEndF = [BadEndF; nWhl];
end

if length(BadStartF)>length(BadEndF)
    BadEndF = [BadEndF; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartF>0)
    
    StartIndF = clip(BadStartF-1, 1, nWhl); % StartInd and EndInd give the
    EndIndF = clip(BadEndF+1, 1, nWhl);     % points you are interpolating between
    
    dist_F = sqrt((Whl(StartIndF,1)-Whl(EndIndF,1)).^2+(Whl(StartIndF,2)-Whl(EndIndF,2)).^2);
    ToChuckF = find(BadEndF-BadStartF>=StretchLen ...
        | dist_F > JumpSize);
    % chuck em
    
    for i=ToChuckF(:)'
        cWhl(BadStartF(i):BadEndF(i),1:2) = NaN;
    end
end

% find missing stretches for Rear LED
dGoodR = [-(whltemp(1,3)==-Gap) ; diff(whltemp(:,3)>-Gap)];
BadStartR = find(dGoodR<0);
BadEndR = find(dGoodR>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,3)==-1
    BadEndR = [BadEndR; nWhl];
end

if length(BadStartR)>length(BadEndR)
    BadEndR = [BadEndR; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartR>0)
    StartIndR = clip(BadStartR-1, 1, nWhl); % StartInd and EndInd give the
    EndIndR = clip(BadEndR+1, 1, nWhl);     % points you are interpolating between
    
    dist_R = sqrt((Whl(StartIndR,3)-Whl(EndIndR,3)).^2+(Whl(StartIndR,4)-Whl(EndIndR,4)).^2);
    ToChuckR = find(BadEndR-BadStartR>=StretchLen ...
        | dist_R > JumpSize);
    
    % chuck em
    for i=ToChuckR(:)'
        cWhl(BadStartR(i):BadEndR(i),3:4) = NaN;
    end
end


if 0 % OLD VERSION (BUG?)
    % % now find good ranges
    % dcGood = [-(Whl(1,1)==1) ; diff(cWhl(:,1)>-1)];
    % GoodStart = find(dcGood>0);
    % GoodEnd = find(dcGood<0)-1;
    % % if last point is good, need to finish GoodEnd
    % if cWhl(end,1)>-1
    %     GoodEnd = [GoodEnd; nWhl];
    % end
    % GoodRanges = [GoodStart, GoodEnd];
else
    dcGood_F = diff([0; cWhl(:,1)>-1; 0]);
    GoodStart_F = find(dcGood_F>0);
    GoodEnd_F = find(dcGood_F<0)-1;
    GoodRanges_F = [GoodStart_F, GoodEnd_F];
end


return

end