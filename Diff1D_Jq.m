
%%%%% Simple Jq Plot:
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
load('None_PP.mat');

%Setup limits:
[tmp tII] = min(abs(t/86400-60));
[tmp tFF] = min(abs(t/86400-120));
z = z_w(2:(end-1));
[tmp indtop] = min(abs(z+50));
[tmp indbot] = min(abs(z+150));

tvec = tII:(tFF-1);

kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq_nDIR_KPP = kt.*dTdz*Cp*rho0;
load('Kelvin_nTIW_PP.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq_nDIR_PP = kt.*dTdz*Cp*rho0;
load('TIW_PP.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq_nDIR_nTIW_KPP = kt.*dTdz*Cp*rho0;
load('Kelvin_vadvNUD.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq_nDIR_nTIW_PP = kt.*dTdz*Cp*rho0;

%Plotting:
%Setup figure:
figure;
set(gcf,'Position',[71 6 1845 999]);

linetype = {'.','+','s','x'};
dots = 720/4;
linesize = [20 20 7 20];
colors = {'k',[0    0.4471    0.7412]};%0.4941    0.1843    0.5569]};

Jqplot1 = subplot('Position',[0.2512    0.3423    0.1854    0.5495]);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
ylabel('Depth (m)','FontSize',25);
set(gca,'xtick',[-200:100:0]);
ylim([-150 0]);
xlim([0 270]);
box on;
grid on;
hold on;

Jqplot2 = subplot('Position',[0.45    0.3423    0.1854    0.5495]);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
set(gca,'ytick',[-150:50:0]);
set(gca,'yticklabel',[]);
set(gca,'xtick',[-200:100:0]);
ylim([-150 0]);
xlim([0 270]);
box on;
grid on;
hold on;

plot(mean(Jq_nDIR_KPP,2),z,'-','LineWidth',3,'color', ...
     colors{1});%,Jqplot1);
plot(mean(Jq_nDIR_PP,2),z,'-','LineWidth',3,'color', ...
     colors{2});%,Jqplot1);
plot(mean(Jq_nDIR_nTIW_KPP,2),z,'--','LineWidth',3,'color', ...
     colors{1});%,Jqplot1);
plot(mean(Jq_nDIR_nTIW_PP,2),z,'--','LineWidth',3,'color', ...
     colors{2});%,Jqplot1);
%Fill:
fcolors = {[0.7294    0.8314    0.9569],[ 0.9569  0.8314 0.7294]};
zF = z(1):1:z(end);
TIW1 = interp1(z',mean(Jq_nDIR_KPP,2)',zF,'linear');
nTIW1 = interp1(z',mean(Jq_nDIR_nTIW_KPP,2)',zF,'linear');
TIW2 = interp1(z',mean(Jq_nDIR_PP,2)',zF,'linear');
nTIW2 = interp1(z',mean(Jq_nDIR_nTIW_PP,2)',zF,'linear');

POS1 = (TIW1-nTIW1)<0;POS2 = (TIW2-nTIW2)<0;
%Separate blocks:
blocks1=cell(1,1);sign1=zeros(1,1);cnt = 1;iini = 1;
for zi=1:(length(zF)-1)
    if (POS1(zi+1)~=POS1(zi))
        sign1(cnt)=POS1(zi);
        blocks1{cnt} = zeros(size(zF));
        blocks1{cnt}(iini:zi) = 1;
        cnt = cnt+1;
        iini=zi+1;
    end
end
sign1(cnt)=POS1(zi);blocks1{cnt} = zeros(size(zF));
blocks1{cnt}(iini:zi) = 1;
blocks2=cell(1,1);sign2=zeros(1,1);cnt = 1;iini = 1;
for zi=1:(length(zF)-1)
    if (POS2(zi+1)~=POS2(zi))
        sign2(cnt)=POS2(zi);
        blocks2{cnt} = zeros(size(zF));
        blocks2{cnt}(iini:zi) = 1;
        cnt = cnt+1;
        iini=zi+1;
    end
end
sign2(cnt)=POS2(zi);blocks2{cnt} = zeros(size(zF));
blocks2{cnt}(iini:zi) = 1;

for ii=1:length(sign1)
X=[zF(logical(blocks1{ii})),fliplr(zF(logical(blocks1{ii})))];
Y=[nTIW1(logical(blocks1{ii})),fliplr(TIW1(logical(blocks1{ii})))];
h = fill(Y,X,fcolors{sign1(ii)+1},'Parent',Jqplot1,'EdgeColor','none');
uistack(h,'bottom');
end
for ii=1:length(sign2)
X=[zF(logical(blocks2{ii})),fliplr(zF(logical(blocks2{ii})))];
Y=[nTIW2(logical(blocks2{ii})),fliplr(TIW2(logical(blocks2{ii})))];
h = fill(Y,X,fcolors{sign2(ii)+1},'Parent',Jqplot1,'EdgeColor','none');
uistack(h,'bottom');
end

plot(mean(Jq_KPP,2),z,'-','LineWidth',3,'color', ...
     colors{1},'Parent',Jqplot2);
plot(mean(Jq_PP,2),z,'-','LineWidth',3,'color', ...
     colors{2},'Parent',Jqplot2);
plot(mean(Jq_nTIW_KPP,2),z,'--','LineWidth',3,'color', ...
     colors{1},'Parent',Jqplot2);
plot(mean(Jq_nTIW_PP,2),z,'--','LineWidth',3,'color', ...
     colors{2},'Parent',Jqplot2);

%Fill:
zF = z(1):1:z(end);
TIW1 = interp1(z',mean(Jq_KPP,2)',zF,'linear');
nTIW1 = interp1(z',mean(Jq_nTIW_KPP,2)',zF,'linear');
TIW2 = interp1(z',mean(Jq_PP,2)',zF,'linear');
nTIW2 = interp1(z',mean(Jq_nTIW_PP,2)',zF,'linear');

POS1 = (TIW1-nTIW1)>0;POS2 = (TIW2-nTIW2)>0;
%Separate blocks:
blocks1=cell(1,1);sign1=zeros(1,1);cnt = 1;iini = 1;
for zi=1:(length(zF)-1)
    if (POS1(zi+1)~=POS1(zi))
        sign1(cnt)=POS1(zi);
        blocks1{cnt} = zeros(size(zF));
        blocks1{cnt}(iini:zi) = 1;
        cnt = cnt+1;
        iini=zi+1;
    end
end
sign1(cnt)=POS1(zi);blocks1{cnt} = zeros(size(zF));
blocks1{cnt}(iini:zi) = 1;
blocks2=cell(1,1);sign2=zeros(1,1);cnt = 1;iini = 1;
for zi=1:(length(zF)-1)
    if (POS2(zi+1)~=POS2(zi))
        sign2(cnt)=POS2(zi);
        blocks2{cnt} = zeros(size(zF));
        blocks2{cnt}(iini:zi) = 1;
        cnt = cnt+1;
        iini=zi+1;
    end
end
sign2(cnt)=POS2(zi);blocks2{cnt} = zeros(size(zF));
blocks2{cnt}(iini:zi) = 1;

for ii=1:length(sign1)
X=[zF(logical(blocks1{ii})),fliplr(zF(logical(blocks1{ii})))];
Y=[nTIW1(logical(blocks1{ii})),fliplr(TIW1(logical(blocks1{ii})))];
h = fill(Y,X,fcolors{sign1(ii)+1},'Parent',Jqplot2,'EdgeColor','none');
uistack(h,'bottom');
end
for ii=1:length(sign2)
X=[zF(logical(blocks2{ii})),fliplr(zF(logical(blocks2{ii})))];
Y=[nTIW2(logical(blocks2{ii})),fliplr(TIW2(logical(blocks2{ii})))];
h = fill(Y,X,fcolors{sign2(ii)+1},'Parent',Jqplot2,'EdgeColor','none');
uistack(h,'bottom');
end


%Legend and text:
lg = legend('KPP TIW','PP TIW','KPP no TIW','PP no TIW');
set(lg,'Position',[0.2540    0.3550    0.1125    0.0975]);
text(260,-5,'No Diurnal','FontSize',25,'Parent',Jqplot1,'Backgroundcolor','w','HorizontalAlignment','right');
text(260,-5,'Diurnal','FontSize',25,'Parent',Jqplot2, ...
     'Backgroundcolor','w','HorizontalAlignment','right');

text(260,-145,'linear TIWs','FontSize',25,'Parent',Jqplot1,'Backgroundcolor','w','HorizontalAlignment','right');
text(260,-145,'linear TIWs','FontSize',25,'Parent',Jqplot2,'Backgroundcolor','w','HorizontalAlignment','right');
text(260,-145,'non-linear TIWs','FontSize',25,'Parent',Jqplot1,'Backgroundcolor','w','HorizontalAlignment','right');
text(260,-145,'non-linear TIWs','FontSize',25,'Parent',Jqplot2,'Backgroundcolor','w','HorizontalAlignment','right');
% $$$ lg2 = subplot('Position',[0.99 0.01 0.01 0.01]);
% $$$ hold on;
% $$$ plot(0,0,'-k','LineWidth',2,'visible','off');
% $$$ plot(0,0,'--k','LineWidth',2,'visible','off');
% $$$ set(lg2,'visible','off');
% $$$ lg2 = legend('TIWs','no TIW');
% $$$ set(lg2,'FontSize',25);
% $$$ set(lg2,'Position',[0.7594    0.3581    0.0908    0.0703]);


%%%%Table data:

z = z_w(2:(end-1));
dz = z_rho(2:end)-z_rho(1:(end-1));
zmin = -100;
zmax = -50;
zminp = min(z(z>zmin));
zmaxm = max(z(z<zmax));
zvec = find(z>=zminp & z<=zmaxm);

PARnames = {'_KPP','_PP'};
DIRnames = {'_nDIR',''};
TIWnames = {'_nTIW',''};
JqR = zeros(length(PARnames),length(DIRnames),length(TIWnames));

for PAR=1:2
    for DIR=1:2
        for TIW=1:2
            eval(['Jq = mean(Jq' DIRnames{DIR} TIWnames{TIW} ...
                  PARnames{PAR} ',2);']);
            
            Jqt = sum((Jq(zvec(2:end))+Jq(zvec(1:(end-1))))/2.*...
                      (z(zvec(2:end))-z(zvec(1:(end-1))))) + ...
                  (Jq(zvec(1))+Jq(zvec(1)-1))/2*...
                  (zminp-zmin) + ...
                  (Jq(zvec(end)+1)+Jq(zvec(end)))/2*...
                  (zmax-zmaxm);
            JqR(PAR,DIR,TIW) = Jqt/(zmax-zmin);
        end
    end
end

%Output table:
for PAR=1:2
    [PARnames{PAR}(2:end) ' & ' num2str(round(JqR(PAR,1,1))) ' ('  ...
     num2str(round(JqR(PAR,2,1))) ') & ' ...
     num2str(round(JqR(PAR,1,2))) ' ('  num2str(round(JqR(PAR,2,2))) ') & ' ...
     num2str(round((JqR(PAR,1,2)/JqR(PAR,1,1)-1)*100)) ' (' ...
     num2str(round((JqR(PAR,2,2)/JqR(PAR,2,1)-1)*100)) ') \\ ']
end

     
     
