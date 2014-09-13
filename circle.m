function circle(R,x0,y0)
alpha=0:pi/50:2*pi;%½Ç¶È[0,2*pi]%R=2;%°ë¾¶
x=x0+R*cos(alpha);
y=y0+R*sin(alpha);
plot(x,y,'linewidth',2)
%fill(x,y,'b','edgealpha',0)
%^axis equal

function circle11(R,x0,y0)
alpha=0:pi/50:2*pi;%½Ç¶È[0,2*pi]%R=2;%°ë¾¶
x=x0+R*cos(alpha);
y=y0+R*sin(alpha);
plot(x,y,':','linewidth',2)
%fill(x,y,'b','edgealpha',0)
axis equal





% hold on;
% max=max(M);
% min=min(M);
% theata=0:0.01:6.3;
% x0=0.5*(max(1,1)+min(1,1));
% y0=0.5*(max(1,2)+min(1,2));
% a=x0+0.5*(max(1,1)-min(1,1));
% b=y0+0.5*(max(1,2)-min(1,2));
% x=a*cos(theata);
% y=b*sin(theata);
% plot(x,y)

