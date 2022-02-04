load('C:\David\Data\Tim\Data\JBon.mat')


dat_l=FTdata.trial{1}([8:10],:)';
rate=FTdata.fsample;

figure
dp1m(dat_l(:,1),rate,6,'STN\_L01',1)
figure
dp1m(dat_l(:,2),rate,6,'STN\_L12',1)
figure
dp1m(dat_l(:,3),rate,6,'STN\_L23',1)

[f12,t12,cl12]=sp2a2_R2(dat_l(:,1),dat_l(:,2),rate,11);
[f13,t13,cl13]=sp2a2_R2(dat_l(:,1),dat_l(:,3),rate,11);
[f23,t23,cl23]=sp2a2_R2(dat_l(:,2),dat_l(:,3),rate,11);
cl12.what=[FTdata.label{8},'-',FTdata.label{9}];
cl13.what=[FTdata.label{8},'-',FTdata.label{10}];
cl23.what=[FTdata.label{9},'-',FTdata.label{10}];

figure
psp2_R2(f12,t12,cl12,500,100,50,1)
figure
psp2_R2(f13,t13,cl13,500,100,50,1)
figure
psp2_R2(f23,t23,cl23,500,100,50,1)

return
