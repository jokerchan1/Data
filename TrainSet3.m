A=hdfinfo('TrainSet3.hdf');
addpath('Meteorological data set1');
addpath('Meteorological data set2');
addpath('Precipitation data set1');
addpath('Precipitation data set2');
La=hdfread(A.Vgroup.Vgroup(1).SDS(1));%Latitude
Lo=hdfread(A.Vgroup.Vgroup(1).SDS(2));%Longtitude
Band19in14=hdfread(A.Vgroup.Vgroup(2).SDS(1));%Band 19 is on the 14th Band
Band20_36=hdfread(A.Vgroup.Vgroup(2).SDS(3));%Band 20-36
Band1_2=hdfread(A.Vgroup.Vgroup(2).SDS(5));%Band 1-2
Band3_7=hdfread(A.Vgroup.Vgroup(2).SDS(8));%Band 3-7
%%Latitude and longitude geographic location
LoB=121.5;
LoE=130;
LaB=55;
LaE=50;
temp=distance(La,Lo,LaB,LoB);
[xb,yb]=find(temp==min(min(temp)), 1 );
temp=distance(La,Lo,LaE,LoE);
[xe,ye]=find(temp==min(min(temp)), 1 );

if(xb>xe)
    temp=xb;
    xb=xe;
    xe=temp;
end
if(yb>ye)
    temp=yb;
    yb=ye;
    ye=temp;
end
xd=xe-xb;
yd=ye-yb;
%%Interception of band information
Band1=Band1_2(1,xb*5:xe*5-1,yb*5:ye*5-1);
Band1=squeeze(Band1);
Band2=Band1_2(2,xb*5:xe*5-1,yb*5:ye*5-1);
Band2=squeeze(Band2);
Band3=Band3_7(1,xb*5:xe*5-1,yb*5:ye*5-1);
Band3=squeeze(Band3);
Band4=Band3_7(2,xb*5:xe*5-1,yb*5:ye*5-1);
Band4=squeeze(Band4);
Band5=Band3_7(3,xb*5:xe*5-1,yb*5:ye*5-1);
Band5=squeeze(Band5);
Band26=Band19in14(15,xb*5:xe*5-1,yb*5:ye*5-1);
Band26=squeeze(Band26);
Band17=Band19in14(12,xb*5:xe*5-1,yb*5:ye*5-1);
Band17=squeeze(Band17);
Band18=Band19in14(13,xb*5:xe*5-1,yb*5:ye*5-1);
Band18=squeeze(Band18);
Band19=Band19in14(14,xb*5:xe*5-1,yb*5:ye*5-1);
Band19=squeeze(Band19);
Band7=Band3_7(5,xb*5:xe*5-1,yb*5:ye*5-1);
Band7=squeeze(Band7);
Band6=Band3_7(4,xb*5:xe*5-1,yb*5:ye*5-1);
Band6=squeeze(Band6);
Band31=Band20_36(11,xb*5:xe*5-1,yb*5:ye*5-1);
Band31=squeeze(Band31);
Band32=Band20_36(12,xb*5:xe*5-1,yb*5:ye*5-1);
Band32=squeeze(Band32);
P1=zeros(xd*5,yd*5);%Band1
P2=zeros(xd*5,yd*5);%Band2
P3=zeros(xd*5,yd*5);%Band3
P4=zeros(xd*5,yd*5);%Band4
P5=zeros(xd*5,yd*5);%Band5
P6=zeros(xd*5,yd*5);%Band6
P7=zeros(xd*5,yd*5); %%Band7
P26=zeros(xd*5,yd*5);%%Band26
P17=zeros(xd*5,yd*5);%%Band17
P18=zeros(xd*5,yd*5);%Band18
P19=zeros(xd*5,yd*5);%Band19
P31=zeros(xd*5,yd*5);%Band31
P32=zeros(xd*5,yd*5);%Band32

NDVI=zeros(xd*5,yd*5);


%Data preprocessing


%Calculation of reflectance
P26=double(Band26);
P1=double(Band1);
P2=double(Band2);
P3=double(Band3);
P4=double(Band4);
P5=double(Band5);
P6=double(Band6);
P7=double(Band7);
P17=double(Band17);
P18=double(Band18);
P19=double(Band19);
P31=double(Band31);
P32=double(Band32);
P17=(P17-A.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value(12)).*A.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value(12);
P18=(P18-A.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value(13)).*A.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value(13);
P19=(P19-A.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value(14))*A.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value(14);
P26=(P26-A.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value(15))*A.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value(15);
T31=(P31-A.Vgroup.Vgroup(2).SDS(3).Attributes(7).Value(11))*A.Vgroup.Vgroup(2).SDS(3).Attributes(6).Value(11);
T31=(6.6261e-34)*(3e8)./((11.03e-9)*(1.38e-23)*log(1+(2*(6.6261e-34)*(3e8)*(3e8)./(T31*(11.03e-9)^5))));
T31=T31./100;
T32=(P32-A.Vgroup.Vgroup(2).SDS(3).Attributes(7).Value(12))*A.Vgroup.Vgroup(2).SDS(3).Attributes(6).Value(12);
T32=(6.6261e-34)*(3e8)./((12.02e-9)*(1.38e-23)*log(1+(2*(6.6261e-34)*(3e8)*(3e8)./(T32*(12.02e-9)^5))));
T32=T32./100;
P1=P1.*A.Vgroup.Vgroup(2).SDS(5).Attributes(9).Value(1);
P2=P2.*A.Vgroup.Vgroup(2).SDS(5).Attributes(9).Value(2);
P3=P4.*A.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value(1);
P4=P4.*A.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value(2);
P5=P5.*A.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value(3);
P6=P6.*A.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value(4);
P7=P7.*A.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value(5);
ND7=(P2-P7)./(P2+P7);
NDVI=(P2-P1)./(P2+P1);%Vegetation coverage
FMC=P7./P2;


EVI=2.5*(P2-P1)./(P2+6*P1-7.5*P3+1);



% Calculate RG and FDDI
RG=zeros(xd*5,yd*5);
RG=(NDVI-min(min(NDVI)))./(max(max(NDVI))-min(min(NDVI)));


%Expand the latitude and longitude matrix
Latitude=zeros(xd,yd);
Lotitude=zeros(xd,yd);

Latitude(1:xd,1:yd)=La(xb:xe-1,yb:ye-1);
Lotitude(1:xd,1:yd)=Lo(xb:xe-1,yb:ye-1);

Latitude = imresize(Latitude, [xd*5 yd*5], 'bicubic');
Lotitude = imresize(Lotitude, [xd*5 yd*5], 'bicubic');


%%Fire Point
fire=zeros(xd*5,yd*5);
Fire=readtable('TrainSet3_fire.csv');
Fire=table2array(Fire);
la=0;
lo=0;
la=Fire(1,1);
lo=Fire(1,2);
for k=1:size(Fire,1)
    la=Fire(k,1);
    lo=Fire(k,2);
    temp=distance(Latitude,Lotitude,la,lo);
    [x,y]=find(temp==min(min(temp)), 1 );
    fire(x,y)=1;
end
fire(1,:)=0;
fire(:,1)=0;
fire(end,:)=0;
fire(:,end)=0;




temp2=find(fire==1&P1<0.3&NDVI>0.1&P26<0.02&RG>0.5&P2<0.2);
[row2, col2]=ind2sub([xd*5 yd*5],temp2);

%%Temperature

filelist=dir(fullfile('Meteorological data set1\*201704*.xls'));
n=length(filelist);
qiwen=zeros(n,1);
La_Xls=zeros(n,1);
Lo_Xls=zeros(n,1);
for i=1:n
    f=readtable(filelist(i).name,'ReadVariableNames',false);
    f=table2array(f);
    La_Xls(i,1)=str2double(cell2mat((f(2,12))));
    Lo_Xls(i,1)=str2double(cell2mat((f(2,13))));     
    qiwen(i,1)=(str2double(cell2mat((f(2,7))))-32)*5/9;
end
F=scatteredInterpolant(La_Xls,Lo_Xls,qiwen,'natural','nearest');
QiWen=F(Latitude,Lotitude);

%%Precipitation
filelist=dir(fullfile('Precipitation data set1\*201703*.xls'));
n=length(filelist);
jiangshui=zeros(n,7);
La_Xls=zeros(n,1);
Lo_Xls=zeros(n,1);
for i=1:n
    count=1;
    f=readtable(filelist(i).name,'ReadVariableNames',false);
    f=table2array(f);
    La_Xls(i,1)=f(31,2);
    Lo_Xls(i,1)=f(31,3); 
    for j=31:-1:25
        jiangshui(i,count)=f(j,5);
        count=count+1;
    end
        
end

F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,1),'natural','nearest');
JiangShui1=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,2),'natural','nearest');
JiangShui2=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,3),'natural','nearest');
JiangShui3=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,4),'natural','nearest');
JiangShui4=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,5),'natural','nearest');
JiangShui5=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,6),'natural','nearest');
JiangShui6=F(Latitude,Lotitude);
F=scatteredInterpolant(La_Xls,Lo_Xls,jiangshui(:,7),'natural','nearest');
JiangShui7=F(Latitude,Lotitude);

%%Relative humidity
filelist=dir(fullfile('Meteorological data set1\*201704*.xls'));

n=length(filelist);
La_Xls=zeros(n,1);
Lo_Xls=zeros(n,1);
rh=zeros(n,1);
for i=1:n
    f=readtable(filelist(i).name,'ReadVariableNames',false);
    f=table2array(f);
    a=7.63;
    b=241.9;
    La_Xls(i,1)=str2double(cell2mat((f(2,12))));
    Lo_Xls(i,1)=str2double(cell2mat((f(2,13))));   
    Td=(str2double(cell2mat((f(2,17))))-32)*5/9;
    T=(str2double(cell2mat((f(2,7))))-32)*5/9;
    cd=6.11*10^(a*Td/(b+Td));
    c=6.11*10^(a*T/(b+T));
    rh(i,1)=round(cd/c*100,2);
end

F=scatteredInterpolant(La_Xls,Lo_Xls,rh,'natural','nearest');
RH=F(Latitude,Lotitude);

%%Wind Speed
filelist=dir(fullfile('Meteorological data set1\*201704*.xls'));
n=length(filelist);
La_Xls=zeros(n,1);
Lo_Xls=zeros(n,1);
fengsu=zeros(n,1);
for i=1:n
    f=readtable(filelist(i).name,'ReadVariableNames',false);
    f=table2array(f);
    La_Xls(i,1)=str2double(cell2mat((f(2,12))));
    Lo_Xls(i,1)=str2double(cell2mat((f(2,13))));   
    fengsu(i,1)=str2double(cell2mat((f(2,5))));
end

F=scatteredInterpolant(La_Xls,Lo_Xls,fengsu,'natural','nearest');
FengSu=F(Latitude,Lotitude);

%%Altitude
[P,R]=geotiffread('Train3_Altitude.tif');
P=imresize(P,[1054 2494]);
LaB=R.LatitudeLimits(1);
LaE=R.LatitudeLimits(2);
LoB=R.LongitudeLimits(1);
LoE=R.LongitudeLimits(2);
La=LaB:R.CellExtentInLatitude*20:LaE-R.CellExtentInLatitude;
La=La';
Dem_Latitude=repelem(La,1,2494);
Lo=LoB:R.CellExtentInLongitude*20:LoE-R.CellExtentInLongitude;
Dem_Longitude=repelem(Lo,1054,1);

DEM = griddata(Dem_Latitude, Dem_Longitude, double(P), Latitude, Lotitude);

%%Slope
[S,R]=geotiffread('Train3_Slope.tif');
S=imresize(S,[1054,2494],'bicubic');
LaB=R.LatitudeLimits(1);
LaE=R.LatitudeLimits(2);
LoB=R.LongitudeLimits(1);
LoE=R.LongitudeLimits(2);
La=LaB:R.CellExtentInLatitude*20:LaE-R.CellExtentInLatitude;
La=La';
Slope_Latitude=repelem(La,1,2494);
Lo=LoB:R.CellExtentInLongitude*20:LoE-R.CellExtentInLongitude;
Slope_Longitude=repelem(Lo,1054,1);
Slope = griddata(Slope_Latitude, Slope_Longitude, double(S), Latitude, Lotitude);




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gaussian convolution
Test=imgaussfilt(fire,10);%6,10
Test=Test*10;%75,200
temp_Imag=find(Test>0);
[row3,col3]=ind2sub([xd*5 yd*5],temp_Imag);

%%%%%%%%%%%%%%%%%%%%%%%%

%%Training set data preprocessing
train=zeros(round(size(row3,1)*1.2),14);
for i=1:size(row3,1)
    train(i,1)=RH(row3(i),col3(i));
    train(i,2)=QiWen(row3(i),col3(i));
    train(i,3)=NDVI(row3(i),col3(i));
    train(i,4)=FengSu(row3(i),col3(i));
    train(i,5)=JiangShui1(row3(i),col3(i));
    train(i,6)=JiangShui2(row3(i),col3(i));
    train(i,7)=JiangShui3(row3(i),col3(i));
    train(i,8)=JiangShui4(row3(i),col3(i));
    train(i,9)=JiangShui5(row3(i),col3(i));
    train(i,10)=JiangShui6(row3(i),col3(i));
    train(i,11)=JiangShui7(row3(i),col3(i));    
    train(i,12)=Slope(row3(i),col3(i));
    train(i,13)=DEM(row3(i),col3(i));
    train(i,14)=Test(row3(i),col3(i));
end
temp1=find(Test==0&P1<0.3&NDVI>0.1&P26<0.02&RG>0.5&P2<0.2);
[row1,col1]=ind2sub([xd*5 yd*5],temp1);
% Non-fire data
A2=randperm(size(temp1,1),round(size(row3,1)*0.2));
[row4,col4]=ind2sub([xd*5 yd*5],A2);
for i=1:round(size(row3,1)*0.2)
    train(size(row3,1)+i,1)=RH(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,2)=QiWen(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,3)=NDVI(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,4)=FengSu(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,5)=JiangShui1(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,6)=JiangShui2(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,7)=JiangShui3(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,8)=JiangShui4(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,9)=JiangShui5(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,10)=JiangShui6(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,11)=JiangShui7(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,12)=Slope(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,13)=DEM(row1(A2(i)),col1(A2(i)));
    train(size(row3,1)+i,14)=0;
end