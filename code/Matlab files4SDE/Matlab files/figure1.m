
N=50000; M=700;
x1=zeros(N,1);x2=x1; xi1=x1; xi2=x1; eta=x1; 
for k=1:N
  [x1(k),x2(k),xi1(k),xi2(k),eta(k)]=SAV(M);
end

x0=-5:0.1:5;

figure(1)
h1=histogram(x1,x0,'Normalization','pdf');
hold on
h2=histogram(x2,x0,'Normalization','pdf');
con=3.5;
y=exp(-1/4*(x0.^2-1).^2)/con;
%y=x0.*exp(-2*x0)/con;
plot(x0,y,'r',linewidth=2)
legend('OAV method','TAV method','stationary measure')


[hd1,edges]= histcounts(x1,x0);
hd1=hd1/sum(hd1); hd1(hd1<1e-10)=1e-10;
bc = (edges(1:end-1) + edges(2:end)) / 2;
y1=exp(-1/4*(bc.^2-1).^2)/con;
y1(y1<1e-10)=1e-10; y1=y1/sum(y1);
KL1=sum(hd1.*log(hd1./y1))
%title(['One SAV: KL div=',num2str(KL1,'%.4f')])


hd2= histcounts(x2,x0);
hd2=hd2/sum(hd2); hd2(hd2<1e-10)=1e-10;
KL2=sum(hd2.*log(hd2./y1))
saveas(gcf,'distx','pdf')

% figure(2)
% subplot(3,1,1)
% h=histogram(xi1,'Normalization','pdf');
% title('\xi of OAV method')
% subplot(3,1,2)
% h=histogram(xi2,'Normalization','pdf');
% title('\xi of TAV method')
% subplot(3,1,3)
% h=histogram(eta,'Normalization','pdf');
% title('\eta of TAVs method')
%saveas(gcf,'distxieta','pdf')