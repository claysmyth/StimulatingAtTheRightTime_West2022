function npdplotter_110717(NPD_data,NPD_sim,F,R,bestn,labelna)
if isempty(NPD_data)
    NPD_data = {zeros(size(NPD_sim{1}))};
end
if ~isfield(R.plot,'cmap')
    R.plot.cmap = [1 0 0];
end
if isempty(NPD_sim)
    NPD_sim= {zeros(size(NPD_data{1}))};
end
if nargin<5
    bestn = 1;
end
if isempty(labelna)
    labelna = 'NPD';
end


% Main Function Starts Here
NPD_data_n = NPD_data{1};
O = size(NPD_data_n,1);

for C = 1:O
    figure(C*10)
    clf
    for L = 1:length(NPD_sim)
        NPD_sim_n = NPD_sim{L};
        
        if L == bestn
            lwid = 2;
        else
            lwid = 0.5;
        end
        N = size(NPD_data_n,2); M = size(NPD_data_n,3);
        
        k = 0;
        for i = 1:N
            for j = 1:M
                k = k+1;
                subplot(N,M,k)
                
                if i == j
                    R2(i,j) = rsquare(squeeze(NPD_sim_n(C,i,j,1,:)),squeeze(NPD_data_n(C,i,j,1,:)));
                    ab(1) = plot(F,squeeze(NPD_sim_n(C,i,j,1,:)),'k--','linewidth',lwid); hold on
                    ab(1).Color = R.plot.cmap;
                    ab(2) = plot(F,squeeze(NPD_data_n(C,i,j,1,:)),'k','linewidth',2);
                    ab(2).Color = [0 0 0];
                    
                    plot(F,squeeze(abs(NPD_sim_n(C,i,j,1,:))),'r--','linewidth',lwid); hold on
                    plot(F,squeeze(imag(NPD_sim_n(C,i,j,1,:))),'b--','linewidth',lwid);
                    plot(F,squeeze(abs(NPD_data_n(C,i,j,1,:))),'r','linewidth',2);
                    plot(F,squeeze(imag(NPD_data_n(C,i,j,1,:))),'b','linewidth',2);
                    xlabel('Hz'); ylabel('Power'); %title(sprintf('Ch %1.f Pxx',i))
                else
                    if i<j
                        a(1) = plot(F,squeeze(NPD_sim_n(C,i,j,2,:)),'r--','linewidth',lwid); hold on
                        a(2) = plot(F,squeeze(NPD_data_n(C,i,j,2,:)),'r','linewidth',2);
                        
                        a(1).Color = R.plot.cmap;
                        a(2).Color = [0 0 0];
                    elseif i>=j
                        a(1) = plot(F,squeeze(NPD_sim_n(C,i,j,2,:)),'b--','linewidth',lwid); hold on
                        a(2) = plot(F,squeeze(NPD_data_n(C,i,j,2,:)),'b','linewidth',2);
                        a(1).Color = R.plot.cmap;
                        a(2).Color = [0 0 0];
                        
                    end
                    
                    xlabel('Hz'); ylabel(labelna); 
                end
                xlim([min(R.frqz) max(R.frqz)]);
                axis square
            end
        end
    end
    set(gcf,'Position',[380         235        1112         893])
end
