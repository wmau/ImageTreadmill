function bounds = sections_treadmill(x,y,direction,plotit)
%bounds = sections_treadmill(x,y,direction)
%
%   Imposes boundaries on the maze that correspond to arbitrary regions.
%   Also plots these boundaries for visualization. 
%
%   INPUTS
%       X&Y: Tracking data. 
%   
%       direction: String, specifying left, right, or alternation. 
%
%   OUTPUT
%       bounds: Structure array with fields that each have an x and y
%       field. The elements in the X and Y fields correspond to the bottom
%       left, right, top right, left corners of the boxes defining the
%       region. 
%           1. Base
%           2. Center
%           3. Choice
%           4. Approach for left and right
%           5. Left and right
%           6. Return for left and right 
%

%%
    xmax = max(x); xmin = min(x);
    ymax = max(y); ymin = min(y); 
    
    w = (xmax-xmin)/7;          %Width of arms.
    cx = 1.3;                   %Constants.
    cy = 1.8;
    
    switch direction
        case 'left'
            centerY =   [ymin+cy*w      ymin+cy*w       ymin+cy*1.5*w   ymin+cy*1.5*w];
        case 'right'
            centerY =   [ymin           ymin            ymin+w          ymin+w];
        case 'alternation'
            centerY =   [ymin+cy*w      ymin+cy*w       ymin+cy*1.5*w   ymin+cy*1.5*w];
    end
    
    choiceY =           [centerY(1)     centerY(1)      centerY(3)      centerY(3)]; 
  
    switch direction
        case 'left'
            approachY = [ymin           ymin            choiceY(1)      choiceY(1)];
        case 'right'
            approachY = [choiceY(3)     choiceY(3)      ymax            ymax];
    end
    
    %Center arm. 
    center.x =          [xmin+w         xmax-cx*w       xmax-cx*w       xmin+w];
    center.y =          centerY; 
    
    %Choice. 
    choice.x =          [xmin           center.x(1)     center.x(1)     xmin];
    choice.y =          choiceY; 
        
    %Base.
    base.x =            [center.x(2)    xmax            xmax            center.x(2)];
    base.y =            center.y;
    
    switch direction
    case 'left'
        %Left arm.
        left.x =        center.x;
        left.y =        [ymin           ymin            ymin+w          ymin+w];
        
        %Left return. 
        return_l.x =    base.x;
        return_l.y =    [ymin           ymin            center.y(1)     center.y(1)];

        %Approach.
        approach_l.x =  choice.x; 
        approach_l.y =  approachY;
    case 'right'
        %Right arm.
        right.x =       center.x; 
        right.y =       [ymax           ymax            ymax-w          ymax-w];
        
        %Right return.
        return_r.x =    base.x; 
        return_r.y =    [base.y(3)      base.y(3)       ymax            ymax];
        
        %Approach. 
        approach_r.x =  choice.x; 
        approach_r.y =  approachY; 
    case 'alternation'
        %Left arm.
        left.x =        center.x; 
        left.y =        [ymin           ymin            ymin+w          ymin+w];
        
        %Left return.
        return_l.x =    [center.x(2)    xmax            xmax            center.x(2)];
        return_l.y =    [ymin           ymin            center.y(1)     center.y(1)];
        
        %Right arm.
        right.x =       center.x; 
        right.y =       [ymax-w         ymax-w          ymax            ymax]; 
        
        %Right return.
        return_r.x =    return_l.x; 
        return_r.y =    [center.y(3)    center.y(3)     ymax            ymax]; 
        
        %Approach left.
        approach_l.x =  choice.x; 
        approach_l.y =  [ymin           ymin            choiceY(1)      choiceY(1)];
        
        %Approach right. 
        approach_r.x =  choice.x; 
        approach_r.y =  [choice.y(3)    choice.y(3)     ymax            ymax];
    end
    
    %Plot trajectory. 
    if plotit
    plot(x,y); 
    hold on;

        %Plot common sections. 
        plot(   [center.x center.x(1)],         [center.y center.y(1)],         'k-',...
                [choice.x choice.x(1)],         [choice.y choice.y(1)],         'k-',...               
                [base.x,base.x(1)],             [base.y,base.y(1)],             'k-');
    end
     
    %Plot sections specific to left/right.        
    switch direction
    case 'left'
        if plotit
        plot(   [left.x,left.x(1)],             [left.y,left.y(1)],             'k-',...
                [approach_l.x,approach_l.x(1)], [approach_l.y,approach_l.y(1)], 'k-',...
                [return_l.x,return_l.x(1)],     [return_l.y,return_l.y(1)],     'k-');
        end

        %Build struct.
        bounds.approach_l = approach_l;
        bounds.left = left; 
        bounds.return_l = return_l;
        
        %Zeros for right arm. 
        bounds.approach_r.x = zeros(1,4);       bounds.approach_r.y = zeros(1,4); 
        bounds.right.x = zeros(1,4);            bounds.right.y = zeros(1,4);
        bounds.return_r.x = zeros(1,4);         bounds.return_r.y = zeros(1,4);
    case 'right'
        if plotit
        plot(   [right.x,right.x(1)],           [right.y,right.y(1)],           'k-',...
                [approach_r.x,approach_r.x(1)], [approach_r.y,approach_r.y(1)], 'k-',...
                [return_r.x,return_r.x(1)],     [return_r.y,return_r.y(1)],     'k-');
        end
        
        bounds.approach_r = approach_r;
        bounds.right = right;
        bounds.return_r = return_r; 
        
        %Zeros for left arm. 
        bounds.approach_l.x = zeros(1,4);       bounds.approach_l.y = zeros(1,4); 
        bounds.left.x = zeros(1,4);            bounds.left.y = zeros(1,4);
        bounds.return_l.x = zeros(1,4);         bounds.return_l.y = zeros(1,4);
    case 'alternation'
        if plotit
        plot(   [approach_l.x,approach_l.x(1)], [approach_l.y approach_l.y(1)], 'k-',...
                [left.x,left.x(1)],             [left.y,left.y(1)],             'k-',...
                [return_l.x,return_l.x(1)],     [return_l.y,return_l.y(1)],     'k-',...
                [approach_r.x,approach_r.x(1)], [approach_r.y approach_r.y(1)], 'k-',...
                [right.x,right.x(1)],           [right.y,right.y(1)],           'k-',...
                [return_r.x,return_r.x(1)],     [return_r.y,return_r.y(1)],     'k-');
        end    
        %Build struct. 
        bounds.approach_l = approach_l;
        bounds.left = left; 
        bounds.return_l = return_l;
        
        bounds.approach_r = approach_r;
        bounds.right = right;
        bounds.return_r = return_r; 
    end
    
    %Build struct. 
    bounds.base = base;
    bounds.center = center;
    bounds.choice = choice;
    bounds.direction = direction; 
    
end