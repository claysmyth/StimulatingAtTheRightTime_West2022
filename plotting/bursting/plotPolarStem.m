function p = plotPolarStem(X,cmap,mag,ls,mark,lw,ms) 
if nargin<6
    lw = 1.5;
end
if nargin<7
    ms = 150;
end
% X =circ_mean(aftdPhi{band,state},[],1);
if nargin<3
mag = 0.1*abs(mean(exp(1i.*X)));
end
phi = angle(mean(exp(1i.*X)));
phi_sem = circ_std(angle(exp(1i.*X)),[],[],2); %./sqrt(numel(X));

% polarplot([phi phi],[0 mag],'Color',cmap,'LineWidth',1.5); hold on % plots the main vector
theta = linspace(phi-phi_sem,phi+phi_sem,32);
p= polarplot(theta,repmat(mag,1,size(theta,2)),'Color',cmap,'LineWidth',lw,'LineStyle',ls); hold on
polarscatter(circ_mean(X'),mag,ms,mark,'filled','MarkerFaceColor',cmap)
