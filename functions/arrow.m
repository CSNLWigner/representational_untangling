function [] = arrow(P1,P2,h,w,lw,c,z)
% draw a vector with arrowhead using normalized coordinates
% (use for drawing vertical or horizontal arrows)
%  starting point: P1
%  ending point: P2
%  size of the arrowhead: h (height), w (weight)
%   inner height: z < h (deafault z = 0, just a triangle)
%  line width: lw, color: c
if nargin < 5, lw = 1; end
if nargin < 6, c = 'k'; end
if nargin < 7, z = 0; end
axes_position = [0 0 1 1];
axes('unit','normalized','position',axes_position,'visible','off','xlim',[0,1],'ylim',[0,1])
hold on
[phi,~] = cart2pol(P2(1)-P1(1),P2(2)-P1(2));
line([P1(1),P2(1)-cos(phi)*h],[P1(2),P2(2)-sin(phi)*h],'linewidth',lw,'color',c)
patch(P2(1)+[0 -cos(phi)*h+sin(phi)*w/2 -cos(phi)*(h-z) -cos(phi)*h-sin(phi)*w/2 0],P2(2)+[0 -sin(phi)*h-cos(phi)*w/2 -sin(phi)*(h-z) -sin(phi)*h+cos(phi)*w/2 0],c,'EdgeColor','none')