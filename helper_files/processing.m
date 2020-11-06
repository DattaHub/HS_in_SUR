clear
%mylist = [dir('Ysim5_*') dir('Ysim6_*') dir('Ysim3_*')];
%mylist = [dir('nograph_Ysim5_*') dir('nograph_Ysim6_*') dir('nograph_Ysim3_*')];
%mylist = [dir('nogamma_Ysim5_*') dir('nogamma_Ysim6_*') dir('nogamma_Ysim3_*')];

mylist = [dir('q300Ysim1_*') dir('q300Ysim2_*') dir('q300Ysim3_*') ...
          dir('q300Ysim4_*') dir('q300Ysim5_*') dir('q300Ysim6_*') ...
          dir('q300Ysim7_*') dir('q300Ysim8_*') dir('q300Ysim9_*') ...
          dir('q300Ysim10_*') ];

% $$$ mylist = [dir('nograph_q300Ysim1_*') dir('nograph_q300Ysim2_*') ...
% $$$           dir('nograph_q300Ysim3_*') dir('nograph_q300Ysim4_*') ...
% $$$           dir('nograph_q300Ysim5_*') dir('nograph_q300Ysim6_*') ...
% $$$           dir('nograph_q300Ysim7_*') dir('nograph_q300Ysim8_*') ...
% $$$           dir('nograph_q300Ysim9_*') dir('nograph_q300Ysim10_*')];
% $$$ 
mylist = [dir('nonograph_q300Ysim1_*') dir('nonograph_q300Ysim2_*') ...
          dir('nonograph_q300Ysim3_*') dir('nonograph_q300Ysim4_*') ...
          dir('nonograph_q300Ysim5_*') dir('nonograph_q300Ysim6_*') ...
         dir('nonograph_q300Ysim7_*') dir('nonograph_q300Ysim8_*') ...
          dir('nonograph_q300Ysim9_*') dir('nonograph_q300Ysim10_*')];

% $$$ mylist = [dir('nononograph_q300Ysim1_*') dir('nononograph_q300Ysim2_*') ...
% $$$           dir('nononograph_q300Ysim3_*') dir('nononograph_q300Ysim4_*') ...
% $$$           dir('nononograph_q300Ysim5_*') dir('nononograph_q300Ysim6_*') ...
% $$$          dir('nononograph_q300Ysim7_*') dir('nononograph_q300Ysim8_*') ...
% $$$           dir('nononograph_q300Ysim9_*') dir('nononograph_q300Ysim10_*')];


mylist = [dir('nogamma_q300Ysim1_*') dir('nogamma_q300Ysim2_*') ...
          dir('nogamma_q300Ysim3_*') dir('nogamma_q300Ysim4_*') ...
          dir('nogamma_q300Ysim5_*') dir('nogamma_q300Ysim6_*') ...
          dir('nogamma_q300Ysim7_*') dir('nogamma_q300Ysim8_*') ...
          dir('nogamma_q300Ysim9_*') dir('nogamma_q300Ysim10_*')];

threshold = (0:0.05:1);
TPR = zeros(length(threshold),1);
FPR = zeros(length(threshold),1);

zz=zeros(498, length(mylist));
adjzz=zeros(300,300,length(mylist));

for i = 1:length(mylist)
        load(mylist(1,i).name)
        zz(:,i) = z1_save;
        adjzz(:,:,i) = adj_save;
        clear(mylist(i).name)
end
meanz=mean(zz,2);

trueind=[30 40 57 62 161 239 269 322 335 399 457];
set(gca,'fontsize',20);
plot(meanz, 'linewidth',2);
hold on;
plot(trueind, meanz(trueind), '.', 'color', 'red', 'MarkerSize',30);
xlabel('Predictor','fontsize',20);
ylabel('Posterior Probability','fontsize',20);
hold off;

for k = 1:length(threshold )
    count = zeros(498,1); 
    z1_o = find(meanz>=threshold(k)); 
        for j =1:498
            if any(z1_o == j) 
                count(j) = count(j) +1;
            end
        end
    
    TPR(k) = ((count(30)>0) +(count(40)>0)+ (count(57)>0)+ (count(62)>0)+ ...
              (count(161)>0) + (count(239)>0)+ (count(269)>0) + ...
              (count(322)>0) + (count(335)>0) + (count(399)>0)+ (count(457)>0))/11;
    FPR (k) = (sum(count>0) - ((count(30)>0) +(count(40)>0)+ (count(57)>0)+ (count(62)>0)+ ...
              (count(161)>0) + (count(239)>0)+ (count(269)>0) + ...
              (count(322)>0) + (count(335)>0) + (count(399)>0)+ (count(457)>0)))/(498-11);
end

set(gca,'fontsize',20);
plot(FPR,TPR,'linewidth',2);
xlabel('FPR','fontsize',20);
ylabel('TPR','fontsize',20);
line('linewidth',2);
%save('ROC.mat','TPR','FPR')



meanadj = mean(adjzz,3);
%meanadj=meanadj.*1.8;
for i =1:300 meanadj(i,i)=1; end
imagesc(meanadj);
colormap(flipud(gray))
%colormap('default')
colorbar
set(gca,'fontsize',20);
%spy(meanadj>0.33,20);
xlabel('Responses','fontsize',20);
ylabel('Responses','fontsize',20);

threshold = (0:0.05:1);
TPR = zeros(length(threshold),1);
FPR = zeros(length(threshold),1);
           
for k = 1:length(threshold)
    count = zeros(300,300); 
    %[r,c] = find(meanadj>threshold(k)); 
    %count(r,c) = 1;
    count(find(meanadj>=threshold(k)))=1;
    
    TPR(k) = ( sum(sum(count(250:300,250:300))))./( 51*51);
    FPR (k) = ((sum(sum(count(:,:)>0))) - (sum(sum(count(250:300,250:300))))) ...
              ./(300*300 - ( 51*51));
end
TPR(length(TPR))=0;
FPR(length(FPR))=0;
TPR(1)=1;
FPR(1)=1;
set(gca,'fontsize',20);
plot(FPR,TPR,'--','linewidth',2);
%xlim([0 0.1]);
%ylim([0 0.1]);
xlabel('FPR','fontsize',20);
ylabel('TPR','fontsize',20);
line('linewidth',2);

% save('ROCadj.mat', 'TPR','FPR');
 