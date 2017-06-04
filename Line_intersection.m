% this attempt at linking of dikes will look at extending dikes to the Y
% axis or to the edges of the domain and then looking for the dikes at some
% close distance (or their midpoints are at some distance that is not too
% vast from the line created.  
%%
clc 
clear all

cd('/Users/Matthew/GitHub/Dike_linking');

%slope intercept script 3
% certain pre-processing of the dikes took place in Arcgis, removing
% vertices so that the lines were simplified to straight lines and only the
% end points were exported for use in this script.

%%


dikes = shaperead('simplified_dike_endpoints.shp');


%% now separate X's and Y's.

% i want to try and get out the x's and y's for the dike endpoints/start
% points

for i = 1:length(dikes)
    xtemp = [dikes(i).X]';
    xtemp(isnan(xtemp)) = []; %delete the NaN out of the xcoords
    
    ytemp = [dikes(i).Y]'; % delete the NaNs out of the ycoords
    ytemp(isnan(ytemp)) = [];
    
    x_coords{i} = xtemp;
    y_coords{i} = ytemp;
%     
%     x_coords{i} = [dikes(i).X]';
%     y_coords{i} = [dikes(i).Y]';
end


clear ytemp xtemp i

%% Perhaps  add a rotation matrix here



%%
% extract all of the start points, which are the odd numbered indices
a = 1;
for i = 1:2:(length(dikes))
    startX(a) = [x_coords{i}];
    startY(a) = [y_coords{i}];
    a = a + 1;
end

%extract all of the end points, which are the odd numbered indices
a = 1;
for i = 2:2:(length(dikes))
    
    
    endX(a) = [x_coords{i}];
    endY(a) = [y_coords{i}];
    a = a +1;
end

% %what if we normalized all variables right here? by normalizing the x and
% y locations of the points...?










%calculate slopes and y intercepts
a = 1;
for i = 1:((length(dikes)/2))
        
    slope(i) = ((endY(a)-startY(a)))./((endX(a)-startX(a))); % manuall checked appears right
    a = a + 1;
end

for i = 1:length(slope)
    Yintercept(i) = startY(i) - (slope(i).*startX(i));
%     Yintercept2(i) = endY(i) - (slope(i).*endX(i));
end


slope = slope';
Yintercept = Yintercept';



clear a i

%As of 5/17/17 I have confirmed that the slopes and y intercepts appear to
%be calculated correctly. One of the last things I'm going to check is
%whether or not the y intercept is the same regardless of whether the
%beginning for ending of the line is used, as it should be the same...
% 

% % % Below is a testing script for looking at the original mapped dikes
% and the lines from those lines to their y intercept to ensure correct
% linkage between all.

% xor = horzcat(startX', endX');
% yor = horzcat(startY', endY');
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% set(gca, 'Units','normalized');
% set(gcf, 'color', [1 1 1]);
% hold on;
% 
% 
% for i = 640
% 
% 
% line(xor(i,:), yor(i,:),'Color','r');
% plot(startX(i), startY(i),'*');
% plot(endX(i), endY(i),'r*');
% 
% line([startX(i) 0],[startY(i) Yintercept(i)])
% 
% % plot(0, Yintercept)
% 
% end
% 

%% Standardization 
% Try to normalize both the slope and y intercept to the same variance.
% The hope of the normalization step is to reduce the difference between
% the scales of Yintecept and slope

%This step eliminates those pesky infinite slope dikes 
% infs =  isinf(slope);
% slope(infs) = NaN;
% Yintercept(infs) = NaN;


% Call normc to normalized by the sum of the squares of each element in the
% column of slope and y intercept.  I opt for this method vs the
% normalization ased on max and min values because that transformation
% appears to change the relationship betwene the original data and the now
% transformed data, unlike the below functions.
normalizedSlope = normc(slope);
normalizedYint = normc(Yintercept);


%% Now look to find the minimization of slopes and Y intercepts so lines are 
%  as similar as possible 

% D will be the variable describing the distance between slope and y
% intercept, trying to find small values here
D = zeros(length(slope),length(slope));
h = waitbar(0, 'Please Be patient... wait .... keep waiting...');

for a = 1:length(slope)
    
    waitbar(a/numel(slope));
    for j = 1:length(slope)
        
%         D(a,j) = sqrt((slope(a) - slope(j))^2 + (Yintercept(a) - Yintercept(j))^2); % Euclidean distance
        
%         D(a,j) = abs((normalizedSlope(a) - normalizedSlope(j))) + abs((normalizedYint(a) - normalizedYint(j))); % Manhattan distance
       
    %Euclidian distance using normalized slope and intercept
        D(a,j) = sqrt((normalizedSlope(a) - normalizedSlope(j))^2 + (normalizedYint(a) - normalizedYint(j))^2); % Euclidean distance    end

    end
end
close(h)

%There are repeated values as you are comparing dikes with other dikes
%twice in the course of this analysis, so this command will eliminate the
%dikes above the diagonal.
D2 = triu(D);

clear h a j 

%%
clear matchDikes
% matchDikes = zeros(100*length(slope),7); % create dumby variable, speeds processing
a = 1;
h = waitbar(0, 'Please Be patient... wait .... keep waiting...');
for i= 1:length(slope)
    
    waitbar(i/numel(slope));
    
    for b = 1:length(slope)
        prob = D(i,b);
        
        %extract stored values to be put into match dikes array
%         s1 = slope(i);
%         s2 = slope(b);
%         y1 = Yintercept(i);
%         y2 = Yintercept(b);
        
        
        s1 = normalizedSlope(i);
        s2 = normalizedSlope(b);
        y1 = normalizedYint(i);
        y2 = normalizedYint(b);
        
        
        if D2(i,b) <= 1e-4 && D2(i,b) > 0 %0.00001 threshold for determining dike closeness
            
            matchDikes(a,:) = [i b prob s1 s2 y1 y2]; % will store the locations of matching dikes
            %for later in this variable
%             I and B are the indexes that refer back to which dikes are
%             matched
            % prob is the D or distance between the slopes and y intercepts
            % S1 and S2 are the slopes of the two dikes
%             Y1 and Y2 are the y intercepts for the two dikes
            a = a + 1;
        else
            D2(i,b) = NaN;
        end
    end
end

close(h)


% trim the empty cells off of matchDikes

clear i b a s1 s2 y1 y2 h

%%

%%  Ok so, now I have extracted the dikes that are possible matches
% but there are repetitions throughout the dataset.  Dikes that are matched
% with themselves and obviously have a very small D value.  I would like to
% eliminate these from the dataset.  Also there are going to repetitions of
% 1 vs 5 and 5 vs 1 are the same dike.  I need to consider how to eliminate
% these dikes.



%This loop shound remove all of the values where matchdikes(i,1) =
%matchdikes(1,2); This step may have been made irrelavent by removing the
%diagonal D values

for i = 1:length(matchDikes)
    if matchDikes(i,1) == matchDikes(i,2)
        matchDikes(i,:) = NaN;
    else
    end
end

    
%now I would like to delete all of the NaN values I just created.
matchDikes(~any(~isnan(matchDikes),2),:) = [];

%now clip off all those unnecesary zeroes
% matchDikes(matchDikes == 0) = [];

% And go ahead and order those matched dikes by their distance
prob = matchDikes(:,3);
[~, IX] = sort(prob);
matchDikesSorted = matchDikes(IX,:);

% matchDikesSorted = flipud(matchDikesSorted);
%% Ok so now I need some loop that figures out whether or not a dike has 
%  already been matched to another dike

% should be a for loop that goes through the list and compares if two dikes
% have been matched before...
matchA = matchDikesSorted(:,1);
matchB = matchDikesSorted(:,2);


% psuedo code, sort both lists and then look at the similarities (but the
% first one should be in the intact ordering... maybe that will help?

% for i = 1:length(matchedDikes)
%     for j = 1:length(matchedDikes)
%         
%         if(


% now I would like to run through all of the dike matches and calculate the
% hypothetical distances between the matched dikes and then sort the dikes
% in some way based off of their distance apart or eliminite some that seem
% a bit excessively far.

%%
phantomDikeLength = zeros(length(slope),1);

for i = 1:length(matchDikesSorted)
    dikeA = matchDikesSorted(i,1);
    dikeB = matchDikesSorted(i,2);
    
    phantomDikeLength(i) = sqrt((startX(dikeA) - startX(dikeB))^2 + (startY(dikeA) - startY(dikeB))^2);
end
    
% 
figure('units','normalized','outerposition',[0 0 1 1])
set(gca, 'Units','normalized');
set(gcf, 'color', [1 1 1]);
hold on;
histogram(phantomDikeLength)%,'normalization','probability');
% 
% histogram(dikeLength,'normalization','probability');

xlabel('Length (m)');
ylabel('Frequency');
vline(meanDikeLength);
vline(modeDikeLength);
vline(maxDikeLength)


%considering applying a threshold at the 4 - 5 km length scale for whether
%or not the connection is spurious, should help to ensure that more local
%dike matches occur.

%stores the index of the matched dikes that fulfill that requirement
thresholdLength = 30000; %4 km 

propDikeLengthIX = find(phantomDikeLength < thresholdLength & phantomDikeLength > 0);
propDikeLength = phantomDikeLength(propDikeLengthIX);

matchDikesShort = matchDikesSorted(propDikeLengthIX,:);

%% Ok, so now, I want to go to all of the dike pairings that are listed
% and draw lines between them, maybe pick a small subset to start with.

% startX and startY
% endX and endY
xor = horzcat(startX', endX');
yor = horzcat(startY', endY');

%plots the original dikes that I started with, working out from there.
figure('units','normalized','outerposition',[0 0 1 1])
set(gca, 'Units','normalized');
set(gcf, 'color', [1 1 1]);
hold on;
% xlim([4.5e5 4.9e5])
% ylim([4.98e6 5.04e6])

for i = 1:length(slope)
    line(xor(i,:), yor(i,:), 'LineWidth', 2);
end


% plot a certain number of the synthetic dikes
%will be the highest correlation dikes (proxy for probability first)
 for i = 100:120
%     dike1 = matchDikesShort(i,1); %identify the two dikes separately
%     dike2 = matchDikesShort(i,2);
    
%     dike1 = matchDikesShort(i,1); %identify the two dikes separately
%     dike2 = matchDikesShort(i,2);
    
%plot those dike matches with the longest distances
    dike1 =  matchDikesShortSorted(i,1); %identify the two dikes separately
    dike2 = matchDikesShortSorted(i,2);
    
    startX1 = startX(dike1); %find start xs
    startX2 = startX(dike2);
    
                            % and the end xs
    endX1 = endX(dike1);
    endX2 = endX(dike2);
    
    startY1 = startY(dike1); %find start locations
    startY2 = startY(dike2);
    
    endY1 = endY(dike1);
    endY2 = endY(dike2);
    
    
% %   then plot
%     line([endX1 startX2], [endY1 startY2], 'color', 'c');
%     line([startX1 endX2], [startY1 endY2], 'color','r');
    line([startX1 startX2], [startY1 startY2],'color','m');
    
    
end

xlabel('Easting'); 
ylabel('Northing');
title('Simplified dike dataset and synthetic dikes created through slope-intercept correlation');

%note on 2/13/17 This approach still has issues.  I think I am good on the
%method of comparing the dikes. I think the rub lies in how I'm calculating
%the slopes of the dikes maybe? as there appears to be great scatter in the
%true slope of the plotted dikes even if when I calculate the slope they
%appear similar.



%% Maybe try different ways of connecting disparate dikes?
% using a spline line to fit the end points, rather than forcing a straight
% line?