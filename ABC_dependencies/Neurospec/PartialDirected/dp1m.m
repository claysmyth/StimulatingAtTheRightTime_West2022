function dp1m(dat,rate,plots,plot_title,scale_flag,start_sec,sect_sec)

% function dp1m(dat,rate,plots,plot_title,scale_flag,start_sec,sect_sec)
%
% function to plot multiple sections of data.

if (nargin>5)
  start_pts=round(start_sec*rate)+1;
else	
  start_pts=1;
end
if (nargin>6)
  plot_pts=round(sect_sec*rate);
else
 plot_pts=round((length(dat)-start_pts+1)/plots);
end

stop_pts=start_pts+plot_pts;
end_pts=start_pts+(plot_pts*plots)-1;
if (end_pts>length(dat))
	end_pts=length(dat);
end	
dat_max=max(dat(start_pts:end_pts));
dat_min=min(dat(start_pts:end_pts));

for ind=1:plots;
  subplot(plots,1,ind)
  time=(start_pts-1:stop_pts-1)/rate;
  plot(time,dat(start_pts:stop_pts),'k-')
  if ((ind==1) & (nargin>3))
    title(plot_title);
  end	
  if (nargin<5) | (scale_flag==0)
    axis([-Inf,Inf,-Inf,Inf]);
  else	
    axis([-Inf,Inf,dat_min,dat_max]);
  end	
  start_pts=stop_pts;
  stop_pts=stop_pts+plot_pts;
  if(stop_pts>length(dat))
  	stop_pts=length(dat);
  end
end	
