%% Script to plot pole figure of bicrystal orientations
% written by Julian Mecklenburgh 2025
% Script require MTex 6.1.0 to be installed
clear all
%% load cif file
cs = loadCIF('forsterite');
% deifine boundary
bound=Miller(0,1,6,cs,'hkl');
bound_pole=round(Miller(vector3d(bound),cs,'uvw'));
% define misorientaion of boundary
orib=orientation.byAxisAngle(bound,180*degree,cs);
A=angle(orib);
% Define rotation to get boundary plane normal in centre of Stereogram
ori=orientation.byAxisAngle(Miller(1,0,0,cs,'uvw'),A/2,cs);
% Define rotatins of the compression directions
rot1=rotation.byAxisAngle(xvector,-20*degree)
rot2=rotation.byAxisAngle(yvector,20*degree)
% Define the two comrpession directions
cd1=rot1*zvector;
cd2=rot2*zvector;

% Define the crystal directions to plot
d100=Miller(1,0,0,cs,'uvw'); % For Crystal 1
d100r=orib*d100; % For crystal 2
% Rotate both into required position
d100a=ori*d100;
d100r=ori*d100r;
d010=Miller(0,1,0,cs,'uvw')% For Crystal 1
d010r=orib*d010;% For crystal 2
% Rotate both into required position
d010a=ori*d010;
d010r=ori*d010r;
% Define crystal plane to plot with all symmetric equivalents
p016=symmetrise(Miller(0,1,6,cs,'hkl'),'antipodal');% For Crystal 1
p016r=orib*p016;% For Crystal 2
% Rotate both into required position
p016a=ori*p016;
p016r=ori*p016r;
% Define the shear Directions
sd=Miller(0,6,-1,cs,'uvw');
sda=ori*sd

%% Plot figure
f=figure(2);
clf

h1=plot(d100a,'upper','antipodal','Marker','p','MarkerFaceColor',[0 0.4470 0.7410],'grid','label',...
    sprintf('[%i%i%i]',d100.u,d100.v,d100.w),'textAboveMarker',...
    'DisplayName','[100] Crystal 1');
hold on
h3=plot(d010a,'upper','antipodal','Marker','^','MarkerFaceColor',[0 0.4470 0.7410],'grid','label',...
    sprintf('[%i%i%i]',d010.u,d010.v,d010.w),...
    'DisplayName','[010] Crystal 1');
h5=plot(p016a,'upper','Marker','s','MarkerFaceColor',[0 0.4470 0.7410],'label',...
    ['\{',sprintf('%0.0f%0.0f%0.0f',p016.h(1),p016.k(1),p016.l(1)),'\}'], ...
    'DisplayName',['\{',sprintf('%0.0f%0.0f%0.0f',p016.h(1),p016.k(1),p016.l(1)),'\} Crystal 1']);
h2=plot(d100r,'upper','antipodal','Marker','p','MarkerFaceColor',[0.6350 0.0780 0.1840],'grid','label',...
    sprintf('[%i%i%i]',d100.u,d100.v,d100.w), ...
    'DisplayName','[100] Crystal 2');
h4=plot(d010r,'upper','antipodal','Marker','^','MarkerFaceColor',[0.6350 0.0780 0.1840],'grid','label',...
    sprintf('[%i%i%i]',d010.u,d010.v,d010.w),...
    'DisplayName','[010] Crystal 2');
h6=plot(p016r,'upper','Marker','s','MarkerFaceColor',[0.6350 0.0780 0.1840],'label',...
    ['\{',sprintf('%0.0f%0.0f%0.0f',p016.h(1),p016.k(1),p016.l(1)),'\}'],...
    'textAboveMarker','DisplayName',['\{',sprintf('%0.0f%0.0f%0.0f',p016.h(1),p016.k(1),p016.l(1)),'\} Crystal 2']);


h8=plot(cd1,'upper','Marker','+','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',2,...
    'label','CD1','DisplayName','Comp. Dir. 1');
h7=plot(sda,'upper','antipodal','Marker','x','MarkerEdgeColor',...
    [0.9290 0.6940 0.1250],'LineWidth',2,'label','[06\bar{1}]','textAboveMarker', ...
    'DisplayName','Shear Dir. 1');
h9=plot(cd2,'upper','Marker','+','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2,...
    'label','CD2','DisplayName','Comp. Dir. 2');
h10=plot(d100a,'upper','antipodal','Marker','x',...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2,'label',...
    'SD2','textBelowMarker',...
    'DisplayName','Shear Dir. 2');
hold off
legend
%f.Position=[2356 -240 900 900];