% Waveguiding video analysis code developed by Abhishek Kalpattu at the Fourkas Lab, UMD
clear all
close all

%% Loading data 
FF=0; % Load start timestamp
LF=12; % Load end timestamp

vidObj = VideoReader("DLinesHPol750.avi");
s = struct("cdata",zeros(vidObj.Height,vidObj.Width,3,"uint8"),colormap=[]);
vidObj.CurrentTime = FF;
k = 1;
while vidObj.CurrentTime <= LF
    s(k).cdata = readFrame(vidObj);
    k = k+1;
end

vidObj = VideoReader("Square patterning coupling big.avi");
s = struct("cdata",zeros(vidObj.Height,vidObj.Width,3,"uint8"),colormap=[]);
vidObj.CurrentTime = FF;
k = 1;
while vidObj.CurrentTime <= LF
    s(k).cdata = readFrame(vidObj);
    k = k+1;
end

%% Plot first frame 
figure(1)
for i=1:200
    VG=imagesc(s(i).cdata);
    title(i)
    pause(0.1)
end
%41
% Frame 18 - 59: First line, 60 - 88: Second line,  
figure(1)
for i=90:110
    imshow(s(i).cdata)
    title(i) 
    pause(1)
end
% 11 102
%% BG RGB sum 
figure(1);VG=imagesc(s(340).cdata);
[J,ROI1] = imcrop(VG);
[J2,ROI2] = imcrop(VG);
[J3,ROI3] = imcrop(VG);
[J4,ROI4] = imcrop(VG);
[J4,ROI5] = imcrop(VG);
S1=size(J(:,:,1));S2=size(J(:,:,2));S3=size(J(:,:,3));
EmIntBG=[sum(sum(J(:,:,1)))/(S1(1)*S1(2))+sum(sum(J(:,:,2)))/(S2(1)*S2(2))+sum(sum(J(:,:,3)))/(S3(1)*S3(2))];

% AuNP vid Comp 1: -37 px, Comp 2: 
ROI1=rect;
ROI2=rect;
ROI3=rect;
ROI4=rect;

% Centroid
ROI4Cx = (ROI4(1) + ROI4(1)+ROI4(3))/2
ROI4Cy = (ROI4(2) + ROI4(2)+ROI4(4))/2

ROI1Cx = (ROI1(1) + ROI1(1)+ROI1(3))/2
ROI1Cy = (ROI1(2) + ROI1(2)+ROI1(4))/2

V1A=[];
for i=1:29
    V1A=[V1A,(V1(2)+V1(4))-5.7*i];
end

DROI4=sqrt((ROI4Cy-(V1(2)+V1(4))).^2+(ROI4Cx-(V1(1)+5.34*(1-1))).^2);
for i=1:32
DROI4=[DROI4,sqrt((ROI4Cy-(V1(2)+V1(4))).^2+(ROI4Cx-(V1(1)+5.34*(i-1))).^2)]; %27,59,88,117,147
end
for i=1:29
DROI4=[DROI4,sqrt((ROI4Cy-(V1(2)+V1(4)-5.41*(i-1))).^2+(ROI4Cx-(V1(1)+V1(3))).^2)]; %27,59,88,117,147
end
for i=1:29
DROI4=[DROI4,sqrt((ROI4Cy-(V1(2))).^2+(ROI4Cx-(V1(1)+V1(3)-5.9*(i-1))).^2)]; %27,59,88,117,147
end
for i=1:30
DROI4=[DROI4,sqrt((ROI4Cy-(V1(2)+5.23*(i-1))).^2+(ROI4Cx-(V1(1))).^2)]; %27,59,88,117,147
end

% H1: 50-98, 3.48 px per frame, V1: 99-148, 3.46 px per frame, H2: 149-197,
% 3.41 px per frame, V2: 198-245, 3.60 px per frame
figure(1);VG=imagesc(s(1).cdata);
imh = findobj(figure(1), 'type', 'image');
YourData = get(imh(1),'CData');
Cr = imcrop(YourData, [ROI3]);
S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));

EmInt3=[];
for i=1:279;
    figure(1);VG=imagesc(s(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI3]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt3=[EmInt3,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end


% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
SS=s(11:102);
EmInt5=[];
for i=1:92
    figure(1);VG=imagesc(SS(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI5(1),ROI5(2)+1.835*(i-1),ROI5(3),ROI5(4)]);
%     figure(1);VP=imagesc(Cr);
%     title(i);pause(0.8)
%     title(sprintf('%0.2f ms', i*3*2.5));
%     axis off
%     subplot(4,4,i);
%     axis off;
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt5=[EmInt5,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end
for i=49:98;
    figure(1);VG=imagesc(SS(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1(1)-3.49*(47),ROI1(2)+3.46*(i-49),ROI1(3),ROI1(4)]);
    figure(1);VP=imagesc(Cr);
    title(i);pause(0.8)
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%     EmInt2=[EmInt2,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end
for i=99:147;
    figure(1);VG=imagesc(SS(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1(1)-3.49*(47)+3.41*(i-99),ROI1(2)+3.46*(98-49),ROI1(3),ROI1(4)]);
%     figure(1);VP=imagesc(Cr);
%     title(i);pause(0.8)
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%     EmInt2=[EmInt2,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end
for i=148:195;
    figure(1);VG=imagesc(SS(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1(1)-3.49*(47)+3.41*(147-99),ROI1(2)+3.46*(98-49)-3.60*(i-148),ROI1(3),ROI1(4)]);
%     figure(1);VP=imagesc(Cr);
%     title(i);pause(0.8)
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end

EmInt1=[];
for i=1:length(s);
    figure(1);VG=imagesc(s(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI2]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG2];
end

EmInt1=[];
for i=1:length(s);
    figure(1);VG=imagesc(s(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI3]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG3];
end


T=((LF-FF)/length(s))*(1:length(SS));
plot(T,EmInt1,'LineWidth',1.5);hold on
plot(T,EmInt1,'LineWidth',1.5);
plot(T,EmInt1,'LineWidth',1.5);
xlabel("Time (s)")
ylabel("Background subtracted RGB counts")
legend("ROI 1","ROI 2","ROI 3");
%% Multiple ROI runs
po=5;
ko=5;
RECT=[];
for i=1:po+1
    for j=1:ko+1
        RECT=[RECT;CRR(1)+(i-1)*round(CRR(3)/po),CRR(2)+(j-1)*round(CRR(4)/ko),round(CRR(3)/po),round(CRR(4)/ko)];
    end 
end
for p=17:1:19
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    V1=TotR;
    
    DR=acos(abs((ROCy-(V1(2)+V1(4)))/sqrt((ROCy-(V1(2)+V1(4))).^2+(ROCx-(V1(1))).^2)));
    for i=1:48
    DR=[DR,acos(abs((ROCy-(V1(2)+V1(4)))/sqrt((ROCy-(V1(2)+V1(4))).^2+((ROCx-(3.49*(i-1))-(V1(1)))).^2)))]; %27,59,88,117,147
    end
    for i=1:50
    DR=[DR,acos(abs(((ROCx-3.49*(47))-(V1(1)))/sqrt(((ROCy+3.46*(i-1))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    for i=1:49
    DR=[DR,acos(abs(((ROCy+3.46*(98-49))-(V1(2)+V1(4)))/sqrt(((ROCy+3.46*(98-49))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47)+3.41*(i-1))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    for i=1:48
    DR=[DR,acos(abs(((ROCx-3.49*(47)+3.41*(147-99))-(V1(1)))/sqrt(((ROCy+3.46*(98-49)-3.60*(i-1))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47)+3.41*(147-99))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    DR=(DR.*180)/pi;
    
    DRO=(abs(sqrt((ROCy-(V1(2)+V1(4))).^2+(ROCx-(V1(1))).^2)));
    for i=1:48
    DRO=[DRO,(abs(sqrt((ROCy-(V1(2)+V1(4))).^2+((ROCx-(3.49*(i-1))-(V1(1)))).^2)))]; %27,59,88,117,147
    end
    for i=1:50
    DRO=[DRO,(abs(sqrt(((ROCy+3.46*(i-1))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    for i=1:49
    DRO=[DRO,(abs(sqrt(((ROCy+3.46*(98-49))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47)+3.41*(i-1))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    for i=1:48
    DRO=[DRO,(abs(sqrt(((ROCy+3.46*(98-49)-3.60*(i-1))-(V1(2)+V1(4))).^2+((ROCx-3.49*(47)+3.41*(147-99))-(V1(1))).^2)))]; %27,59,88,117,147
    end
    hold on;scatter(DRO,EN(p,:),15,'o','filled')
end

EN=[];
for l=1:49
    ROI1=RECT(l,:);
    figure(1);VG=imagesc(s(300).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    EmInt1=[];
    for i=1:48;
        figure(1);VG=imagesc(SS(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1(1)-3.49*(i-1),ROI1(2),ROI1(3),ROI1(4)]);
        figure(1);VP=imagesc(Cr);
        title(i);pause(0.8)
    %     title(sprintf('%0.2f ms', i*3*2.5));
    %     axis off
    %     subplot(4,4,i);
    %     axis off;
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%         EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    for i=49:98;
        figure(1);VG=imagesc(SS(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1(1)-3.49*(47),ROI1(2)+3.46*(i-49),ROI1(3),ROI1(4)]);
    %     figure(1);VP=imagesc(Cr);
    %     title(i);pause(0.8)
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    for i=99:147;
        figure(1);VG=imagesc(SS(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1(1)-3.49*(47)+3.41*(i-99),ROI1(2)+3.46*(98-49),ROI1(3),ROI1(4)]);
    %     figure(1);VP=imagesc(Cr);
    %     title(i);pause(0.8)
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    for i=148:196;
        figure(1);VG=imagesc(SS(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1(1)-3.49*(47)+3.41*(147-99),ROI1(2)+3.46*(98-49)-3.60*(i-148),ROI1(3),ROI1(4)]);
    %     figure(1);VP=imagesc(Cr);
    %     title(i);pause(0.8)
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    EN=[EN;EmInt1];
end

EN=[];
for l=31:36
    ROI1=RECT(l,:);
    figure(1);VG=imagesc(s(300).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    EmInt1=[];
    for i=1:length(s);
    figure(1);VG=imagesc(s(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    EN=[EN;EmInt1]; % 4,5,4, 
end
TG=[];ER=[];
for p=1:6
ROI1=RECT(p+12,:);
ROCx = (ROI1(1) + ROI1(1)+ROI1(3))/2;
ROCy = (ROI1(2) + ROI1(2)+ROI1(4))/2;
TG=[TG,ROCx];
ER=[ER,max(EN3(p,:))];
end

EmInt1=[];
% for i=1:length(s);
%     figure(1);VG=imagesc(s(i).cdata);
%     imh = findobj(figure(1), 'type', 'image');
%     YourData = get(imh(1),'CData');
%     Cr = imcrop(YourData, [ROI1]);
%     S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%     EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
% end
%% 

EN=[];LArr=[];ANGarr=[];
for l=99:121
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    figure(1);VG=imagesc(s(1).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    EmInt1=[];LD=[];ANG=[];
    for i=1:length(s)-11
        LD=[LD,sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
        ANG=[ANG,acos(abs(DRF(1,i)-ROCx)/sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2))];
        figure(1);VG=imagesc(s(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
%         figure(1);VP=imagesc(Cr);
%         title(i);pause(0.8)
    %     title(sprintf('%0.2f ms', i*3*2.5));
    %     axis off
    %     subplot(4,4,i);
    %     axis off;
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    EN=[EN;EmInt1];
    LArr=[LArr;LD];
    ANGarr=[ANGarr;ANG];
end
plot(1:length(s)-11,EmInt1(1,:),'LineWidth',1.5)
plot(1:length(s)-11,EN3(1,:)+EN3(2,:)+EN3(3,:)+EN3(4,:)+EN3(5,:)+EN3(6,:)+EN3(7,:),'LineWidth',1.5)
scatter(LD,EmInt1(1,:),40,'o',"filled")
scatter(LArr(4,:),EN(4,:),40,'o',"filled")

% Parallel Processing
EN=[];LArr=[];
for i=1:length(s)-12
    EmInt1=[];LD=[];
    for l=1:30
        ROI1=RECT(l,:);p=l;
        ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
        ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
        LD=[LD;sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
        % BG Corr
        figure(1);VG=imagesc(s(1).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
        % ROI calc
        figure(1);VG=imagesc(s(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1;sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end 
    EN=[EN,EmInt1];
    LArr=[LArr,LD];
end
plot(1:length(s)-12,EN(1,:),'LineWidth',1.5)

MnD=[];
for i=1:30
    MnD=[MnD,min(LArr(i,:))];
end
MxA=[];
for i=1:15
    MxA=[MxA,max(EN_1(i,:))];
end

MP5=[];
for i=1:30
    MP5=[MP5,(EN(i,40))];
end
xx=[];yy=[];
for i=1:length(RECT)
    xx=[xx,(RECT(i,1) + RECT(i,1)+RECT(i,3))/2];
    yy=[yy,(RECT(i,2) + RECT(i,2)+RECT(i,4))/2];
end
DD5=[];
for i=1:length(MP)
    DD5=[DD5,sqrt((xx6(i)-xx6(1)).^2+(yy6(i)-yy6(1)).^2)*0.0581];
end

EnnA=[];LnA=[];ANGa=[];
pl=8;
MD=[EN(pl,:)];MX=[LArr(pl,:)];AN=[ANG(pl,:)];
for j=1:7
    MD=[MD;EN(8*j+pl,:)];
    MX=[MX;LArr(8*j+pl,:)];
    AN=[AN;ANG(8*j+pl,:)];
end
EnnA=[EnnA;MD];LnA=[LnA;MX];ANGa=[ANGa;AN];

% Rearrang LD
EN=[];LArr=[];ANGarr=[];
for l=1:6
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
%     figure(1);VG=imagesc(s(5).cdata);
%     imh = findobj(figure(1), 'type', 'image');
%     YourData = get(imh(1),'CData');
%     Cr = imcrop(YourData, [ROI1]);
%     S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%     EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    EmInt1=[];LD=[];
    for i=1:420
        LD=[LD,sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
%         ANG=[ANG,(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))];
%         figure(1);VG=imagesc(s(i).cdata);
%         imh = findobj(figure(1), 'type', 'image');
%         YourData = get(imh(1),'CData');
%         Cr = imcrop(YourData, [ROI1]);
% %         figure(1);VP=imagesc(Cr);
% %         title(i);pause(0.8)
%     %     title(sprintf('%0.2f ms', i*3*2.5));
%     %     axis off
%     %     subplot(4,4,i);
%     %     axis off;
%         S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
%         EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
%     EN=[EN;EmInt1];
    LArr=[LArr;LD];
%     ANGarr=[ANGarr;ANG];
end

DANG=acos(abs(DRF(2,:)-(TotR(2)+TotR(4)/2))./(sqrt((DRF(2,:)-(TotR(2)+TotR(4)/2)).^2+(DRF(1,:)-(TotR(1))).^2)));
ANG=[]
for l=1:30
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    ANG=[ANG;(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))-DANG];
end
ANG=(ANG)*(180/pi);

writerObj = VideoWriter('AuNP_CL.avi');
writerObj.FrameRate = 20;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(s)
    % convert the image to a frame
    frame = s(i).cdata ;    
    writeVideo(writerObj, frame);
end
close(writerObj)

% OLD Circle points
AX1=linspace(-TotR(4)/2,TotR(4)/2,201);
AX2=sqrt((TotR(4)/2)^2-AX1.^2);
DRF=[((TotR(1)-(TotR(4)/2))-AX2+(TotR(4)/2));(TotR(2)-AX1+(TotR(4)/2))];DRF1=[DRF(1,:);flip(DRF(2,:))];
DRF=[((TotR(1)-(TotR(4)/2))+AX2+(TotR(4)/2));(TotR(2)+AX1+(TotR(4)/2))];DRF2=[DRF(1,:);flip(DRF(2,:))];
DRF=[DRF1(1,:),DRF2(1,:);DRF1(2,:),DRF2(2,:)];

% NEW Circle points
x=0:0.015:pi;y=pi:0.015:2*pi;
G=(TotR(4)/2)*cos(x);J=(TotR(4)/2)*sin(x);
G1=(TotR(4)/2)*cos(y);J1=(TotR(4)/2)*sin(y);
A=[G+(TotR(2)+TotR(4)/2),G1+(TotR(2)+TotR(4)/2)];B=[J1+TotR(1),J+TotR(1)];
DRF=[B;A];
plot(B,A,'o')

x=0:0.0015:pi;y=pi:0.0015:2*pi;
G=(TotR(4)/2)*cos(x);J=(TotR(4)/2)*sin(x);
G1=(TotR(4)/2)*cos(y);J1=(TotR(4)/2)*sin(y);
A=[G+(TotR(2)+TotR(4)/2),G1+(TotR(2)+TotR(4)/2)];B=[J1+TotR(1),J+TotR(1)];
DP=[B;A];
%% GAUS plots

ELe=389;
% GXn=DRF(1,ELe)-100:0.5:DRF(1,ELe)+100;GX=GXn;
% GYn=DRF(2,ELe)-100:0.5:DRF(2,ELe)+100;GY=GYn;
GX=DRF(1,:);
GY=DRF(2,:);
Sig=33.7246;
GAUS=((1/(2*pi*Sig))*exp((-((GX-DRF(1,ELe)).^2+(GY-DRF(2,ELe)).^2))/(2*Sig))).^2;hold on;plot(ANGn(i,:),GAUS*((max(EN(i,:)))./max(GAUS)),'LineWidth',1.5,'color','k')
hold on;plot((1:length(DP(1,:))),GAUS,'LineWidth',1.5)

% DANG=acos(abs(DRF(2,:)-(TotR(2)+TotR(4)/2))./(sqrt((DRF(2,:)-(TotR(2)+TotR(4)/2)).^2+(DRF(1,:)-(TotR(1))).^2)));
x=0:0.015:pi;y=pi:0.015:2*pi;
DANG=[x,y];
ANG=[]
for l=1:72
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    if ROCx-TotR(1)>0 && ROCy-(TotR(2)+TotR(4)/2)>0
        QUAD=pi*2;
        ANG=[ANG;(QUAD-(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2))))-DANG];
    end
    if ROCx-TotR(1)>0 && ROCy-(TotR(2)+TotR(4)/2)<0
        QUAD=pi;
        ANG=[ANG;(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))-DANG+QUAD];
    end 
    if ROCx-TotR(1)<0 && ROCy-(TotR(2)+TotR(4)/2)>0
        QUAD=0;
        ANG=[ANG;(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))-DANG+QUAD];
    end
    if ROCx-TotR(1)<0 && ROCy-(TotR(2)+TotR(4)/2)<0
        QUAD=pi;
        ANG=[ANG;(QUAD-(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2))))-DANG];
    end
end
ANG=(ANG)*(180/pi);

% Global Coordinates
x=0:0.0015:pi;y=pi:0.0015:2*pi;
DANG1=[x,y];
ANGn=[]
for l=1:72
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    if ROCx-TotR(1)>0 && ROCy-(TotR(2)+TotR(4)/2)>0
        QUAD=pi*2;
        ANGn=[ANGn;(QUAD-(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2))))-DANG1];
    end
    if ROCx-TotR(1)>0 && ROCy-(TotR(2)+TotR(4)/2)<0
        QUAD=pi;
        ANGn=[ANGn;(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))-DANG1+QUAD];
    end 
    if ROCx-TotR(1)<0 && ROCy-(TotR(2)+TotR(4)/2)>0
        QUAD=0;
        ANGn=[ANGn;(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2)))-DANG1+QUAD];
    end
    if ROCx-TotR(1)<0 && ROCy-(TotR(2)+TotR(4)/2)<0
        QUAD=pi;
        ANGn=[ANGn;(QUAD-(acos(abs((TotR(2)+TotR(4)/2)-ROCy)/sqrt((TotR(1)-ROCx)^2+((TotR(2)+TotR(4)/2)-ROCy)^2))))-DANG1];
    end
end
ANGn=(ANGn)*(180/pi);
%% 
GX=1:1920;
GY=1:1200;
ZZ=[];
for i=1:length(GY)
    XX=[];
    for j=1:length(GX)
        XX=[XX,(1/(2*pi*Sig))*exp((-((GX(j)-DRF(1,ELe)).^2+(GY(i)-DRF(2,ELe)).^2))/(2*Sig))];
    end
    ZZ=[ZZ;XX];
end
% imagesc(GX,GY,ZZ);

ax1 = axes; 
im = imagesc(ax1,s(429).cdata); 
im.AlphaData = 0.5; % change this value to change the background image transparency 
axis square; 
hold all; 
ax2 = axes; 
im1 = imagesc(ax2,VG); 
im1.AlphaData = 0.5; % change this value to change the foreground image transparency 
hold all;
axis square
ax3 = axes; 
im1 = imagesc(ax3,ZZ); 
im1.AlphaData = 0.3; % change this value to change the foreground image transparency 
hold all;
axis square
linkaxes([ax1,ax2,ax3]) 
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = []; 
ax3.YTick = []; 
set([ax1,ax2,ax3],'Position',[.17 .11 .685 .815]); 

DRF=[linspace(TotR(1),TotR(1),500),linspace(TotR(1),TotR(1)+TotR(3),500),linspace(TotR(1)+TotR(3),TotR(1)+TotR(3),500),linspace(TotR(1)+TotR(3),TotR(1),500);linspace(TotR(2),TotR(2)+TotR(4),500),linspace(TotR(2)+TotR(4),TotR(2)+TotR(4),500),linspace(TotR(2)+TotR(4),TotR(2)+TotR(4),500),linspace(TotR(2),TotR(2),500)];
DANG=(abs(DRF(2,:)-(TotR(2)+TotR(4)/3))./(sqrt((DRF(2,:)-(TotR(2)+TotR(4)/3)).^2+(DRF(1,:)-TotR(1)).^2)));
%% Self-guidance NPs
dx=10;
SEN=[];SLArr=[];
for L=21:402
    l=L
    ROI1=[DRF(1,l)-(dx/2),DRF(2,l)-(dx/2),dx,dx];
    figure(1);VG=imagesc(s(300).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    EmInt1=[];LD=[];ANG=[];
    for i=1:length(s)-10
        LD=[LD,sqrt((DRF(1,i)-DRF(1,l))^2+(DRF(2,i)-DRF(2,l))^2)];
        ANG=[ANG,acos(abs(DRF(1,i)-ROCx)/sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2))];
        figure(1);VG=imagesc(s(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
%         figure(1);VP=imagesc(Cr);
%         title(i);pause(0.8)
    %     title(sprintf('%0.2f ms', i*3*2.5));
    %     axis off
    %     subplot(4,4,i);
    %     axis off;
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1,sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end
    SEN=[SEN;EmInt1];
    SLArr=[SLArr;LD];
end

for i=1:6
    hold on; plot(1:length(s)-11,EN(i,:),'LineWidth',1.5);
end
% ROIs: 25, 36, 44, 81
A = repmat(22.9637,[1 length(EN(21,:))])
[minValue,closestIndex] = min(abs(A-EN(21,:)));
closestValue = N(closestIndex) 
%[8.18182,8.84091,9.87879]

i=36
scatter(ANG(i,:),EN(i,:),30,'o','Filled','color','k');

ZZ=[];
for i=1:length(GY)
    XX=[];
    for j=1:length(GX)
        XX=[XX,(1/(2*pi*(1/sqrt(2))))*exp((-((GX(j)-0.4).^2+(GY(i)-4).^2))/(2*(1/sqrt(2))))];
    end
    ZZ=[ZZ;XX];
end
VG=(s(1).cdata);
for i=1:20
    VG=VG+(s(20+(i-1)*4).cdata-2.5);
end

% Mon Gaussian Fit
i=25
scatter(ANG(i,:),EN(i,:),30,'o','Filled','color','k');
hold on; plot(-100:0.01:260,a1*exp(-(((-100:0.01:260)-b1)/c1).^2)+d,'LineWidth',1.5,'color','r')
xlim([-100 260])
ylim([0 max(EN(i,:))+max(EN(i,:))*0.1])
xlabel("\phi_0 - \theta (degrees)")
ylabel("RGB counts")
legend(append('Data (ROI ',num2str(i),')'),'Single gaussian fit')

% Dual Gaussian fit
i=13
scatter(ANG(i,:),EN(i,:),30,'o','Filled','color','k');
hold on; plot(-100:0.01:260,a1*exp(-(((-100:0.01:260)-b1)/c1).^2)+a2*exp(-(((-100:0.01:260)-b2)/c2).^2)+d,'LineWidth',1.5,'color','r')
xlim([-100 260])
ylim([0 max(EN(i,:))+max(EN(i,:))*0.1])
xlabel("\phi_0 - \theta (degrees)")
ylabel("RGB counts")
legend(append('Data (ROI ',num2str(i),')'),'Double gaussian fit')

DA=[];IN=[];
for i=1:13
    Dist=[];Ints=[];
    for j=1:4
        Dist=[Dist,min(LArr(j+(i-1)*4,:))*0.0581];
        Ints=[Ints,max(EN(j+(i-1)*4,:))];
    end
    DA=[DA;Dist];
    IN=[IN;Ints];
end

for i=5:8
    hold on;plot(DA(i,:),IN(i,:),'*-','LineWidth',1.5);
end

% ROI at max: 2, 11, 18, 27, 33, 42, 43, 57, 66

% FLIP images clockwise 90 Degrees
AN=[];
for i=1:length(AD)
    AL=[];
    for j=1:length(AD)
        AL=[AL;AD(length(AD)-(i-1),j)-3];
    end
    AN=[AN,AL];
end

AX=[];
for i=1:length(AD)
    AL=[];
    for j=1:length(AD)
        if AN(j,i)<0
            AL=[AL;0];
        else 
            AL=[AL;AN(j,i)];
        end
    end
    AX=[AX,AL];
end

m=(ROI1Cy-ROI4Cy)./(ROI1Cx-ROI4Cx);
xx=linspace(ROI1Cx,ROI4Cx,30);
yy=m*xx+ROI4Cy-m*ROI4Cx;
xx2=linspace(ROI1Cx+0,ROI1Cx+(xx(30)-xx(1)+0),30);
yy2=m*xx2+(ROI1Cy+25)-m*(ROI1Cx+0);
xx3=linspace(ROI1Cx+36,ROI1Cx+(xx(30)-xx(1)+36),30);
yy3=m*xx3+(ROI1Cy+50)-m*(ROI1Cx+36);
xx4=linspace(ROI1Cx+60,ROI1Cx+(xx(30)-xx(1)+60),30);
yy4=m*xx4+(ROI1Cy+80)-m*(ROI1Cx+60);
xx5=linspace(ROI1Cx-80,ROI1Cx+(xx(30)-xx(1)-80),30);
yy5=m*xx5+(ROI1Cy-40)-m*(ROI1Cx-80);
xx6=linspace(ROI1Cx-100,ROI1Cx+(xx(30)-xx(1)-100),30);
yy6=m*xx6+(ROI1Cy-50)-m*(ROI1Cx-100);
RECT=[];
for i=1:length(xx)
    RECT=[RECT;[xx(i)-1,yy(i)-1,2,2]];
end
RECT1=[];
for i=1:length(xx2)
    RECT1=[RECT1;[xx2(i)-1,yy2(i)-1,2,2]];
end
RECT2=[];
for i=1:length(xx3)
    RECT2=[RECT2;[xx3(i)-1,yy3(i)-1,2,2]];
end
RECT3=[];
for i=1:length(xx4)
    RECT3=[RECT3;[xx4(i)-1,yy4(i)-1,2,2]];
end
RECT4=[];
for i=1:length(xx5)
    RECT4=[RECT4;[xx5(i)-1,yy5(i)-1,2,2]];
end
RECT5=[];
for i=1:length(xx6)
    RECT5=[RECT5;[xx6(i)-1,yy6(i)-1,2,2]];
end
RECT=[RECT;RECT1;RECT2;RECT3;RECT4;RECT5];
RECT=[RECT;RECT1];

qq=zeros(1200,1920);
ww=zeros(1200,1920);
ww1=ww;
ZH=cat(3,qq,ww,ww1);
for i=1:length(RECT)
    ZH(round(RECT(i,2)):round(RECT(i,2)+RECT(i,4)),round(RECT(i,1)):round(RECT(i,1)+RECT(i,3)),:)=1;
end
ZH=im2uint8(ZH);

qq=zeros(1200,1920);
ww=zeros(1200,1920);
ww1=ww;
ZH1=cat(3,qq,ww,ww1);
for i=1:length(RECT)
    ZH1(round(RECT(i,2)):round(RECT(i,2)+RECT(i,4)),round(RECT(i,1)):round(RECT(i,1)+RECT(i,3)),1)=100;
    ZH1(round(RECT(i,2)):round(RECT(i,2)+RECT(i,4)),round(RECT(i,1)):round(RECT(i,1)+RECT(i,3)),2)=100;
    ZH1(round(RECT(i,2)):round(RECT(i,2)+RECT(i,4)),round(RECT(i,1)):round(RECT(i,1)+RECT(i,3)),3)=10;
end
ZH1=im2uint16(ZH1);

SCL1=zeros(length(AD),length(AD));
SCL1(260:263,200:220)=200;
imagesc(AX+SCL1)
colormap 'hot'

%209,217,244,288,311
REn=[];RANG=[];Dist=[];RD=[];XD=[];u=60;
for i=1:70
    k=find(EN(u,:)==max(EN(u,:)))
    REn=[REn,max(EN(i,k))];
    RANG=[RANG,ANG(i,k)];
    Dist=[Dist,min(LArr(i,:))*0.0581];
    RD=[RD,LArr(i,k)*0.0581];
    p=i;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    p=u;
    ROCxx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCyy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    XD=[XD,sqrt((ROCx-ROCxx)^2+(ROCy-ROCyy)^2)];
end
hold on
scatter(RANG,REn+750,30,'o','Filled');xlabel("\phi_0 - \theta (degrees)")
ylabel("Waveguided emission (RGB counts)")

t = tiledlayout(1,1);
ax1 = axes(t);
scatter(1:20,REn,60,'o','Filled');xlabel('ROI number')
ylabel("Waveguided emission (RGB counts)")
ax2 = axes(t);
scatter(RANG,REn,60,'o','Filled')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax1.Box = 'off';
ax2.Box = 'off';

%ROIs: 2, 9, 15, 25, 31, 45

N=11; % The number of PreVid frames. These will not be analysed later because they don't contain fluorescence information
EN=[];LArr=[];

i=350
EmInt1=[];LD=[];
for l=1:length(RECT) % You can choose to analyse only a portion of the autogenerated ROIs 
    ROI1=RECT(l,:);p=l;
    ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
    ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
    LD=[LD;sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
    % BG Corr
    figure(1);VG=imagesc(s(1).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
    % ROI calc
    figure(1);VG=imagesc(s(i).cdata);
    imh = findobj(figure(1), 'type', 'image');
    YourData = get(imh(1),'CData');
    Cr = imcrop(YourData, [ROI1]);
    S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
    EmInt1=[EmInt1;sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
end 

RD=[2.18-1.21 3.39-2.18 4.59-3.39 5.57-4.59]; %750
RD=[2.42-1.21 3.39-2.42 mean([4.355,4.59])-3.39 5.57-mean([4.355,4.59]) 6.533-5.57]; %775
RD=[1.45-0.48394 2.42-1.45 3.39-2.42 mean([4.355,4.59])-3.39 5.57-mean([4.355,4.59]) mean([6.78,7.02])-5.57]; %825
RD=[1.94-0.72 3.14-1.94 4.36-3.14]; %850

RDA=[1.09 1.09 1.06 1.07 1.21];
Err=[0.07 0.07 0.05 0.06 0.01];

RD=[2.99-1.99 4.18-2.99 5.18-4.18 6.37-5.18 7.36-6.37]; %H Pol
RD=[1.98-0.79 2.97-1.98 3.96-2.97 4.96-3.96 mean([5.95,6.14])-4.96]; %V Pol
RD=[1.43-0.54 3.04-1.43 4.3-3.04 5.55-4.3 6.45-5.55 mean([7.7,7.88])-6.45]; %V Pol

RD=[mean([1.41,1.76])-0.35 mean([3.16,3.51])-mean([1.41,1.76]) 4.92-mean([3.16,3.51]) 6.68-4.92 8.08-6.68 13.0-11.59 14.05-13.0]; %Thick Video
P=[227,227,219;180,212,223;250,229,213;546,521,521;380,364,380;307,364,306;339,351,331;351,372,405;421,376,384;476,471,467;394,403,380;390,367,349;208,236,236;447,414,411]

P=[453,435,371;226,226,226;408,334,334;389,371,371;417,444,480;379,379,370;346,321,313];
P=[0.433
0.411
0.439
0.339
0.323
0.301
0.487
0.493
0.482
0.603
0.537
0.531
0.570
0.625
0.602];


Pavg=[];
for i=1:5
    Pavg=[Pavg;mean([P((i-1)*3+1) P(2+(i-1)*3) P(3+(i-1)*3)])];
end

Err=[];
for i=1:5
    Err=[Err;std([P((i-1)*3+1) P(2+(i-1)*3) P(3+(i-1)*3)])/sqrt(3)];
end

0.3226/0.138

0.1202/0.0564

0.1432/-0.1237

0.3328/0.2454

0.3289/-0.1995

0.2973/-0.2427

BG: 78.66/-242
%%%%

0.5317/-0.7766

0.285/-0.2288

0.7737/-1.202

0.6481/-0.7186

0.3928/-0.5746

0.239/0.3001
GH=[];
for i=1:120
    GH=[GH,min(LArr1(i,:))];
end

P=[];
for i=1:4
    P=[P;min(GH((i-1)*30+1:i*30))];
end

PJ=[45,57,70];
N=11; % The number of PreVid frames. These will not be analysed later because they don't contain fluorescence information
EN=[];LArr=[];
for o=1:3
    i=PJ(o);
    EmInt1=[];LD=[];
    for l=1:length(RECT(:,1))% You can choose to analyse only a portion of the autogenerated ROIs 
        ROI1=RECT(l,:);p=l;
        ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
        ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
        LD=[LD;sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
        % BG Corr
        figure(1);VG=imagesc(s(1).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmIntBG1=sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2));
        % ROI calc
        figure(1);VG=imagesc(s(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1;sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end 
    EN=[EN,EmInt1];
    LArr=[LArr,LD];
end
JK=[];
for i=1:50
    JK=[JK,max(EN(:,i))];
end

for i=11:20
    hold on;plot(1:65,EN(:,i),'o-','LineWidth',1.5,'MarkerFaceColor','b');
end


N=11; % The number of PreVid frames. These will not be analysed later because they don't contain fluorescence information
EN=[];LArr=[];
for i=1:70
    EmInt1=[];LD=[];
    for l=1:length(RECT(:,1))% You can choose to analyse only a portion of the autogenerated ROIs 
        ROI1=RECT(l,:);p=l;
        ROCx = (RECT(p,1) + RECT(p,1)+RECT(p,3))/2;
        ROCy = (RECT(p,2) + RECT(p,2)+RECT(p,4))/2;
%         LD=[LD;sqrt((DRF(1,i)-ROCx)^2+(DRF(2,i)-ROCy)^2)];
        % BG Corr
        EmIntBG1=2.5;
        % ROI calc
        figure(1);VG=imagesc(s(i).cdata);
        imh = findobj(figure(1), 'type', 'image');
        YourData = get(imh(1),'CData');
        Cr = imcrop(YourData, [ROI1]);
        S1=size(Cr(:,:,1));S2=size(Cr(:,:,2));S3=size(Cr(:,:,3));
        EmInt1=[EmInt1;sum(sum(Cr(:,:,1)))/(S1(1)*S1(2))+sum(sum(Cr(:,:,2)))/(S2(1)*S2(2))+sum(sum(Cr(:,:,3)))/(S3(1)*S3(2))-EmIntBG1];
    end 
    EN=[EN,EmInt1];
%     LArr=[LArr,LD];
end
