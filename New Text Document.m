clear
clc
data2Old= xlsread('052120191300.xlsx');
Data= xlsread('To_Dr_Li.csv');
distanceData= Data(360:447,6);
addpath('../Quaternions');
%time index (unit: sec)
k1 = 6738; % starting index
k2 = 14922; % ending index
data2=data2Old(k1:k2,:);
dataTime=data2(:,1)/1000;

dataTime = dataTime-dataTime(1);


dataAccX=data2(:,2);
dataAccY=data2(:,3);
dataAccZ=data2(:,4);
dataGyoX=data2(:,5);
dataGyoY=data2(:,6);
dataGyoZ=data2(:,7);


dataGyo=[];
aa = [50 50 50];
bb = [0.1 0.1 0.1];

for mm = 1:3
    data = data2(:,(mm+4));
    for i=2:(length(data)-1)
     if abs(data(i)-data(i-1))>aa(mm) & abs(data(i)-data(i+1))>aa(mm)
        data(i)=(data(i-1)+data(i+1))/2;
     end
    end
        
% meanBef=sum(data(1:50))/50;

meanBef=sum(data2Old(k1-100:k1-51,mm+4))/50;
data=data-meanBef;
%mechanical filter
for i=1:length(data)
    if abs(data(i))<=bb(mm) 
        data(i)=0;
    end
end

dataGyo=[dataGyo data];  %convert the unit of acceleration to g without moving filter 
       
end


dataAcc=[];
aa = [0.4 0.4 0.4];
bb = [0.02 0.02 0.02];

for kk = 1:3

    data = data2(:,(kk+1));
%eliminate the arbitrary peak

for i=2:(length(data)-1)
    if abs(data(i)-data(i-1))>aa(kk) & abs(data(i)-data(i+1))>aa(kk) & sign(data(i)-data(i-1))*sign(data(i)-data(i+1)) == 1
        data(i)=(data(i-1)+data(i+1))/2;
    end
end


% calibration
% meanBef=sum(data(1:50))/50;
meanBef=sum(data2Old(k1-100:k1-51,kk+1))/50;
data=data-meanBef;
%mechanical filter
for i=1:length(data)
    if abs(data(i))<=bb(kk) 
        data(i)=0;
    end
end

if kk==3
    data=data+meanBef;
    data = -data;
end
 dataAcc=[dataAcc data];  %convert the unit of acceleration to g without moving filter

end 
 
dataGyo = dataGyo * (pi/180);
quaternion = zeros(length(dataTime), 4);
tempQuaternion = [1 0 0 0];
quaternion(1,:) = tempQuaternion;
dataAccQuat = [smooth(dataAcc(:,1), 5) smooth(dataAcc(:,2), 5) -smooth(dataAcc(:,3), 5)];



for t = 2:length(dataTime)
  
tempQuaternion = estQuaternion(tempQuaternion, dataGyo(t,:), dataAccQuat(t,:),dataTime(t)-dataTime(t-1), 0);
% tempQuaternion = estQuaternionMahony(tempQuaternion, dataGyo(t,:), dataAccQuat(t,:),dataTime(t)-dataTime(t-1) , -0.5);

quaternion(t, :) = tempQuaternion;
end

euler = quatern2euler(quaternConj(quaternion)) * (180/pi); % euler angle

% quaternion = rotMat2quatern(euler2rotMat(pitch))

timeEncoder = 0:(length(distanceData)-1);
 
distance = interp1(timeEncoder,distanceData,dataTime,'liner')*0.3048;

% distance=0:0.753/length(data2):0.753;

%multiplying dv by the time difference:
for i=2:length(dataTime)
deltad=(distance(i)-distance(i-1));
dv(:,i)= quaternRotate([1 0 0]*deltad, quaternion(i,:));
end
pos=zeros(3,length(dataTime));
for i=2:length(dataTime)
pos(:,i)=pos(:,i-1)+dv(:,i);
end


figure
plot3(pos(1,:),pos(2,:),pos(3,:));
axis equal








