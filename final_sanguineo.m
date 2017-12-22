model=fopen('E:\Guido\Tesis\Resultados\temp_P_para.txt','w');
c=0;
P=zeros(1,52);
param=zeros(1,35);
a=[1;1;1;1;1;1;2;1;1;1;1;0.04;1;10;1;1;1;1;10;0;10;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
am=a;
%y37=a1;y38=a2;y39=a3;y40=a4;y41=a5;y42=a6;y43=a7=20;y44=a8;y45=a9;y46=a10;y47=a11;y48=a12;y49=a13;y50=a14;y51=a15;
%y52=a16;y53=a17;y54=a18;y55=a19;y56=a20;y57=a21;y58=a22;y59=a23;y60=a24;y61=a25;y62=a26;y63=a27;y64=a28;y65=a29;y66=a30;y67=a31;y68=a32;y69=a33;y70=a34;y71=a35;


for j=0:0.01:0
for k=0:0.1:0
%Rate modifiers
ugt=a(1,1); 
uhk=a(2,1); 
ugpi=a(3,1); 
upfk=a(4,1); 
uald=a(5,1); 
utpi=a(6,1); 
ugapdh=1.8; 
ugpdh=1; 
ugpo=1; 
upgkg=a(10,1); 
upk=a(11,1); 
upyrsec=a(12,1); 
ugk=1; 
uatputc=a(14,1); 
upgm=a(15,1); 
ueno=a(16,1); 
udhap_t=1; 
ugly3p_t=1; 
ugly_sec=15; 
uppp=0; 
upga3_t=.5; 
upep_t=0; 
upepck=a(23,1); 
umdh=a(24,1); 
umal=a(25,1); 
ufum=a(26,1); 
ufred=a(27,1); 
uppdk=a(28,1); 
upyr_t=a(29,1); 
usuc_sec=a(30,1); 
uatputg=a(31,1); 
u13bpga_t=0; 
upgkc=a(33,1); 
ualadh=a(34,1); 
upyr_t2=a(35,1); 
    j,k
format long
tInitial=0;
tFinal=300; 
%yInitial=[1;1;1;1;1;1;1;1;1;1;1;0;1;2.36513927183603;1;1;1;0;0;0;0;1;1;0;1
%;1;1;1;0;0;0;0;1.24527023166839;0.28982498632157];
yInitial=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;j;k;ugt;uhk;ugpi;upfk;uald;utpi;ugapdh;ugpdh;ugpo;upgkg;upk;upyrsec;ugk;uatputc;upgm;ueno;udhap_t;ugly3p_t;ugly_sec;uppp;upga3_t;upep_t;upepck;umdh;umal;ufum;ufred;uppdk;upyr_t;usuc_sec;uatputg;u13bpga_t;upgkc;ualadh;upyr_t2;0;0;0;0;0;0];

tSpan=[tInitial tFinal];
fname='brucei_final_sanguineo_edo_definitivo';
options=odeset('AbsTol',0.0000000000000000001,'RelTol',0.000000000000000000001,'stats','on');
[t,y]=ode15s(fname,tSpan,yInitial);   
c=c+1;
% Index
P(c,1)=c;
P(c,2)=k;
P(c,3)=j;
param(c,1)=c;
param(c,2:35)=yInitial(37:70);
P(c,4:30)=(y((size(y,1)),1:27));
%Glucose intake
P(c,31)=((y((size(y,1)),28)-y((size(y,1)-4),28))/((t(size(t,1)))-t(size(t,1)-4)));
%Flux through Enolase
P(c,32)=((y((size(y,1)),29)-y((size(y,1)-4),29))/((t(size(t,1)))-t(size(t,1)-4)));
%Flux through PEP transporter
P(c,33)=((y((size(y,1)),30)-y((size(y,1)-4),30))/((t(size(t,1)))-t(size(t,1)-4)));
%Pyruvate secretion
P(c,34)=((y((size(y,1)),31)-y((size(y,1)-4),31))/((t(size(t,1)))-t(size(t,1)-4)));
%Succinate Secretion
P(c,35)=((y((size(y,1)),32)-y((size(y,1)-4),32))/((t(size(t,1)))-t(size(t,1)-4)));
%ATP utilization, cytosolic
P(c,36)=((y((size(y,1)),33)-y((size(y,1)-4),33))/((t(size(t,1)))-t(size(t,1)-4)));
%Flux through PPP
P(c,37)=((y((size(y,1)),34)-y((size(y,1)-4),34))/((t(size(t,1)))-t(size(t,1)-4)));
%Glycerol secretion
P(c,38)=((y((size(y,1)),35)-y((size(y,1)-4),35))/((t(size(t,1)))-t(size(t,1)-4)));
%Pyruvate proportion
P(c,39)=((P((c),34)/(P((c),34)+P((c),38)+P((c),35))));
%Glycerol proportion
P(c,40)=((P((c),38)/(P((c),34)+P((c),38)+P((c),35))));
%Succinate proportion
P(c,41)=((P((c),35)/(P((c),34)+P((c),38)+P((c),35))));
%Ratio Succinate/Pyruvate
P(c,42)=((P((c),35)/(P((c),34)+P((c),35))));
%Flux through PPDK
P(c,43)=((y((size(y,1)),73)-y((size(y,1)-4),73))/((t(size(t,1)))-t(size(t,1)-4)));
%Flux through ALADH
P(c,44)=((y((size(y,1)),74)-y((size(y,1)-4),74))/((t(size(t,1)))-t(size(t,1)-4)));
%Glycosomal Pyruvate secretion
P(c,45)=((y((size(y,1)),75)-y((size(y,1)-4),75))/((t(size(t,1)))-t(size(t,1)-4)));
%Pyruvate (measured as glucose consumption)
P(c,46)=P(c,34)/2/P(c,31)*100;
%Succinate (measured as glucose consumption)
P(c,47)=P(c,35)/2/P(c,31)*100;
%Flux through PPP (measured as glucose consumption)
P(c,48)=P(c,37)/P(c,31)*100;
%Proportion of PEP going inside the glycosome (measured as glucose consumption)
P(c,49)=P(c,33)/2/P(c,31)*100;
%Proportion of glucose metabolized to 2PGA
P(c,50)=P(c,32)/2/P(c,31)*100;
%Flux through PEPCK
P(c,51)=((y((size(y,1)),76)-y((size(y,1)-4),76))/((t(size(t,1)))-t(size(t,1)-4)));
%Flux Through GAPDH
P(c,52)=((y((size(y,1)),77)-y((size(y,1)-4),77))/((t(size(t,1)))-t(size(t,1)-4)));
P(c,1:51);

for i=1:1:(size(P,2)-30)
    if P(c,i)<(-0.1)
        display 'WARNING, NEGATIVE!!!!'
        P(c,i),i
        P(c,1:50)='*';
    end;
end;
for i=1:1:(size(P,2))
    if isnan(P(c,i))==1
        i
        display 'WARNING,\n Not a Number!!!!'
        P(c,1:50)='*';
    end;
end;
end;
a=am;
end;

parameters=yInitial(37:70)';
title='c\tJ\tK\tGlucose\tGlucose-6-Phosphate\tFructose-6-Phosphate\tFructose-1,6-Bisphosphate\tGlyceraldehyde-3-Phosphate\tDihydroxyacetone-3-Phosphate, glycosomal\tDihydroxyacetone-3-Phosphate, cytosolic\t1,3-Bisphophoglycerate\t3-Phosphoglycerate\t2-Phosphoglycerate\tPhosphoenol Pyruvate, cytosolic\tPyruvate, cytosolic\tNADH\tGlycerol-3-Phosphate, glycosomal\tGlycerol-3-Phosphate, cytosolic\tGlycerol\t3-Phosphoglycerate, cytosolic\tPhosphoenol Pyruvate, glycosomal\tOxalacetate\tPyruvate, glycosomal\tMalate\tFumarate\tSuccinate\t1,3-Bisphophoglycerate,cytosolic\tPg\tPc\tAlanine\tGlu_intake\tEno_flux\tPep_t_flux\tPyr_sec\tSuc_sec\tATP_ut\tPPP_flux\tGly_sec\tProp_pyr\tProp_gly\tProp_suc\tPerc. Sucinate\tFlux PPDK\tFlux ALADH\tPyr Secreted Flux\tP. Pyr (Glu)\tP. Suc (Glu)\tP. PPP (Glu)\tP. Glu in glycosome\tP. Glu metab. as 2PGA\tFlux PEPCK\n';
title2='Parameters\nc\tGT\tHK\tGPI\tPFK\tALD\tTPI\tGAPDH\tGPDH\tGPO\tPGKg\tPK\tPyrSec\tGK\tATPutC\tPGM\tENO\tDHAP_t\tGLY3P_t\tGly_sec\tPPP\tPGA3_t\tPEP_t\tPEPCK\tMDH\tMAL\tFUM\tFRED\tPPDK\tPyr_t\tSuc_sec\tATPutg\t13BPGA_t\tPGKc\tALADH\n';
save('E:\Guido\Tesis\Resultados\temp_para.txt','y','-ASCII','-tabs');
save('E:\Guido\Tesis\Resultados\tiempo_para.txt','t','-ASCII','-tabs');
fprintf(model,title);
fprintf(model,[repmat('%15.8e \t', 1, size(P, 2)) '\n'], P');
fprintf(model,'\n');
fprintf(model,title2);
fprintf(model,[repmat('%15.8e \t', 1, size(param, 2)) '\n'], param');
fclose(model);