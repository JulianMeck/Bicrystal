%% load cif file
clear all
cs = loadCIF('forsterite');

bound=Miller(0,1,6,cs,'hkl');
ori = orientation.byAxisAngle(bound,180*degree,cs)
A=angle(ori)/degree;
bound_pole=round(Miller(vector3d(bound),cs,'uvw'));
d100=Miller(1,0,0,cs,'uvw');
d061=Miller(0,6,-1,cs,'uvw');
rot2=rotation.byAxisAngle(d100,20*degree);
v100=rot2*bound_pole;
v100r=ori*v100;
rot4=rotation.byAxisAngle(d061,20*degree);
v061=rot4*bound_pole;
v061r=ori*v061;

%% define miller of slip systems
d100r=ori*(d100);
m010=Miller(0,1,0,cs,'hkl');
m010r=ori*(m010);
d001=Miller(0,0,1,cs,'uvw');
d001r=ori*(d001);
m100=Miller(1,0,0,cs,'hkl');
m100r=ori*m100;
m001=Miller(0,0,1,cs,'hkl');
m001r=ori*m001;
m011=Miller(0,1,1,cs,'hkl');
m011=symmetrise(m011);
m011r=ori*m011;
m110=Miller(1,1,0,cs,'hkl');
m110r=ori*m110;
m130=Miller(1,3,0,cs,'hkl');
m130r=ori*m130;
m016=Miller(0,1,6,cs,'hkl');
m016r=ori*m016;

sstxt{1}='(010)[100]';
ss1=slipSystem(d100,m010);
ss1=ss1.symmetrise;
ss1r=ori*ss1;
Asstemp=angle(bound,ss1.n,'nosymmetry')/degree;
[Ass(1),ind1]=min(Asstemp);
Asstemp=angle(bound,ss1r.n,'nosymmetry')/degree;
[Assr(1),ind2]=min(Asstemp);
m_prime(1)=mPrime(ss1(ind1),ss1r(ind2));
SF_100(1)=SchmidFactor(ss1(ind1),v061);
SF_010(1)=SchmidFactor(ss1(ind1),v100);
SF_100r(1)=SchmidFactor(ss1(ind2),v061r);
SF_010r(1)=SchmidFactor(ss1(ind2),v100r);
% angle between Burger's vector and shear direction
phi1(1)=angle(ss1.b(ind1),d100,'nosymmetry');
phi2(1)=angle(ss1.b(ind1),d061,'nosymmetry');



sstxt{2}='(010)[001]';
ss2=slipSystem(d001,m010);
ss2=ss2.symmetrise;
ss2r=ori*ss2;
Asstemp=angle(bound,ss2.n,'nosymmetry')/degree;
[Ass(2),ind1]=min(Asstemp);
Asstemp=angle(bound,ss2r.n,'nosymmetry')/degree;
[Assr(2),ind2]=min(Asstemp);
m_prime(2)=mPrime(ss2(ind1),ss2r(ind2));
SF_100(2)=SchmidFactor(ss2(ind1),v061);
SF_010(2)=SchmidFactor(ss2(ind1),v100);
SF_100r(2)=SchmidFactor(ss2(ind2),v061r);
SF_010r(2)=SchmidFactor(ss2(ind2),v100r);
% angle between Brrger's vector and plane normal
phi1(2)=angle(ss2.b(ind1),d100,'nosymmetry');
phi2(2)=angle(ss2.b(ind1),d061,'nosymmetry');

sstxt{3}='(100)[001]';
ss3=slipSystem(d001,m100);
ss3=ss3.symmetrise;
ss3r=ori*ss3;
Asstemp=angle(bound,ss3.n,'nosymmetry')/degree;
[Ass(3),ind1]=min(Asstemp);
Asstemp=angle(bound,ss3r.n,'nosymmetry')/degree;
[Assr(3),ind2]=min(Asstemp);
m_prime(3)=mPrime(ss3(ind1),ss3r(ind2));
SF_100(3)=SchmidFactor(ss3(ind1),v061);
SF_010(3)=SchmidFactor(ss3(ind1),v100);
SF_100r(3)=SchmidFactor(ss3(ind2),v061r);
SF_010r(3)=SchmidFactor(ss3(ind2),v100r);
% angle between Brrger's vector and plane normal
phi1(3)=angle(ss3.b(ind1),d100);
phi2(3)=angle(ss3.b(ind1),d061);

sstxt{4}='(001)[100]';
ss4=slipSystem(d100,m001);
ss4=ss4.symmetrise;
ss4r=ori*ss4;
Asstemp=angle(bound,ss4.n,'nosymmetry')/degree;
[Ass(4),ind1]=min(Asstemp);
Asstemp=angle(bound,ss4r.n,'nosymmetry')/degree;
[Assr(4),ind2]=min(Asstemp);
m_prime(4)=mPrime(ss4(ind1),ss4r(ind2));
SF_100(4)=SchmidFactor(ss4(ind1),v061);
SF_010(4)=SchmidFactor(ss4(ind1),v100);
SF_100r(4)=SchmidFactor(ss4(ind2),v061r);
SF_010r(4)=SchmidFactor(ss4(ind2),v100r);
% angle between Brrger's vector and plane normal
phi1(4)=angle(ss4.b(ind1),d100);
phi2(4)=angle(ss4.b(ind1),d061);

sstxt{5}='{011}[100]';
ss5=slipSystem(d100,m011);
ss5=ss5.symmetrise;
ss5r=ori*ss5;
Asstemp=angle(bound,ss5.n,'nosymmetry')/degree;
[Ass(5),ind1]=min(Asstemp);
Asstemp=angle(bound,ss5r.n,'nosymmetry')/degree;
[Assr(5),ind2]=min(Asstemp);
m_prime5=mPrime(ss5,ss5r);
[m_prime(5),I]=max(m_prime5,[],'all');
SF=ss5.SchmidFactor(v061);
SF_100(5)=SF(I);
SF=ss5.SchmidFactor(v100);
SF_010(5)=SF(I);
SF=ss5r.SchmidFactor(v061r);
SF_100r(5)=SF(I);
SF=ss5r.SchmidFactor(v100r);
SF_010r(5)=SF(I);

% angle between Brrger's vector and plane normal
phi1(5)=angle(ss5.b(ind1),d100);
phi2(5)=angle(ss5.b(ind1),d061);


sstxt{6}='(110)[001]';
ss6=slipSystem(d001,m110);
ss6=ss6.symmetrise;
ss6r=ori*ss6;
Asstemp=angle(bound,ss6.n,'nosymmetry')/degree;
[Ass(6),ind1]=min(Asstemp);
Asstemp=angle(bound,ss6r.n,'nosymmetry')/degree;
[Assr(6),ind2]=min(Asstemp);
m_prime(6)=mPrime(ss6(ind1),ss6r(ind2));
SF_100(6)=SchmidFactor(ss6(ind1),v061);
SF_010(6)=SchmidFactor(ss6(ind1),v100);
SF_100r(6)=SchmidFactor(ss6(ind2),v061r);
SF_010r(6)=SchmidFactor(ss6(ind2),v100r);
% angle between Brrger's vector and plane normal
phi1(6)=angle(ss6.b(ind1),d100);
phi2(6)=angle(ss6.b(ind1),d061);

sstxt{7}='(130)[001]';
ss7=slipSystem(d001,m130);
ss7=ss7.symmetrise;
ss7r=ori*ss7;
Asstemp=angle(bound,ss7.n,'nosymmetry')/degree;
[Ass(7),ind1]=min(Asstemp);
Asstemp=angle(bound,ss7r.n,'nosymmetry')/degree;
[Assr(7),ind2]=min(Asstemp);
m_prime(7)=mPrime(ss7(ind1),ss7r(ind2));
SF_100(7)=SchmidFactor(ss7(ind1),v061);
SF_010(7)=SchmidFactor(ss7(ind1),v100);
SF_100r(7)=SchmidFactor(ss7(ind2),v061r);
SF_010r(7)=SchmidFactor(ss7(ind2),v100r);
% angle between Brrger's vector and plane normal
phi1(7)=angle(ss7.b(ind1),d100);
phi2(7)=angle(ss7.b(ind1),d061);





T016=table(sstxt',m_prime',Ass',phi1'/degree,phi2'/degree,SF_100',SF_100r',SF_010',SF_010r');
T016.Properties.VariableNames = {'Slip System','mPrime','Angle',...
    'phi_100','phi_061','SF [100]','SF [100] Rotated','SF [0 6 -1]','SF [0 6 -1] Rotated'};
T016
writetable(T016,'T016_new.xls')
m016=symmetrise(m016);
m016r=ori*m016;
%% Plot figure
f=figure(1);
clf


h1=plot(d100,'upper','antipodal','Marker','p','MarkerFaceColor',[0 0.4470 0.7410],'grid','label',...
    sprintf('[%i%i%i]',d100.u,d100.v,d100.w),'textAboveMarker',...
    'DisplayName','[100] Crystal 1');
hold on


h5=plot(m016,'upper','Marker','s','MarkerFaceColor',[0 0.4470 0.7410],'label',...
    ['\{',sprintf('%0.0f%0.0f%0.0f',m016.h(1),m016.k(1),m016.l(1)),'\}'], ...
    'DisplayName','\{014\} Crystal 1');
% h11=plot(m011,'upper','Marker','o','MarkerFaceColor',[0 0.4470 0.7410],'label',...
%     ['\{',sprintf('%0.0f%0.0f%0.0f',m011.h(1),m011.k(1),m011.l(1)),'\}'], ...
%     'DisplayName','\{011\} Crystal 1');
h2=plot(d100r,'upper','antipodal','Marker','p','MarkerFaceColor',[0.6350 0.0780 0.1840],'grid','label',...
    sprintf('[%i%i%i]',d100.u,d100.v,d100.w), ...
    'DisplayName','[100] Crystal 2');


h6=plot(m016r,'upper','Marker','s','MarkerFaceColor',[0.6350 0.0780 0.1840],'label',...
    ['\{',sprintf('%0.0f%0.0f%0.0f',m016.h(1),m016.k(1),m016.l(1)),'\}'],...
    'textAboveMarker','DisplayName','\{014\} Crystal 2');
% h12=plot(m011r,'upper','Marker','o','MarkerFaceColor',[0.6350 0.0780 0.1840],'label',...
%     ['\{',sprintf('%0.0f%0.0f%0.0f',m011.h(1),m011.k(1),m011.l(1)),'\}'],...
%     'textAboveMarker','DisplayName','\{011\} Crystal 2');


h8=plot(v100,'upper','Marker','+','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',2,...
    'label','CD1','DisplayName','Comp. Dir. 1');

h9=plot(v061,'upper','Marker','+','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2,...
    'label','CD2','DisplayName','Comp. Dir. 2');

hold off
legend
f.Position=[2356 -240 900 900];