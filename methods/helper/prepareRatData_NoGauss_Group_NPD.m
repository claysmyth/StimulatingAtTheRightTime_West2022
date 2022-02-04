function R = prepareRatData_NoGauss_Group_NPD(R,plotop,nulldat)
if nargin<2
    plotop = 1;
end
if nargin<3
    nulldat = 0;
end
load([R.filepathn '\NPD_paper_RatNPD_150618.mat']);
NPDmat = fyA;
load([R.filepathn '\nsPow_paper_RatNPD_150618.mat']);
nsPowmat = fyA;
load([R.filepathn '\Frq_paper_RatNPD_150618.mat']);
meannpd_data = nan(1,6,6,4,311);
% condsel = [1 2];
R.condnames = {'OFF'};
R.Bcond = 1;
condsel = 2;
chsel = [1 3 2 4]'; % M1 STR GPe STN

for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                F_data = fxA(:,1);
                F_model = R.frqz;
                % Log transform of group average
                Pxy = abs(log10(mean(vertcat(nsPowmat{chsel(i),condsel(C),:}),1)));
                % Remove the powerline frequency
                Pxy(F_data>48 & F_data<52) = [];
                F_data(F_data>48 & F_data<52) = [];
                % Put in the same frequency space as the models
                Pxy = interp1(F_data,Pxy,F_model);
                % Regress out the 1/f background
                [xCalc yCalc b Rsq] = linregress(log10(F_model)',Pxy');
                if nulldat == 1
                    Pxy = 10.^yCalc;
                else
                    Pxy = Pxy-fliplr(yCalc');
                    clear ft gt
                    % Bring back to non-log space
                    Pxy = 10.^Pxy;
                end
                % Normalize and zero-base
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy-min(Pxy);
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                
                for k = 1:size(NPDmat,3)
                    F_data = fxA(:,1);
                    F_model = R.frqz;
                    Cxy = mean(horzcat(NPDmat{chsel(i),chsel(j),k,condsel(C),:})',1);
                    
                    % Remove the powerline frequency
                    Cxy(F_data>48 & F_data<52) = [];
                    F_data(F_data>48 & F_data<52) = [];
                    % Put in the same frequency space as the models
                    Cxy = interp1(F_data,Cxy,F_model);
                    if nulldat == 1
                        Cxy = zeros(size(Cxy));
                    else
                         clear ft gt
                    end
                    meannpd_data(C,i,j,k,:) = Cxy;
                end
                
            end
        end
    end
end

% Set data as working version
R.data.feat_emp = meannpd_data;
R.data.feat_xscale = R.frqz;

% Plot CSD
if strcmp('CSD',R.data.datatype)
    csdplotter_220517({meannpd_data},[],F_data,R)
elseif strcmp('NPD',R.data.datatype)
    npdplotter_110717({meannpd_data},[],R.frqz,R,[],[])
end