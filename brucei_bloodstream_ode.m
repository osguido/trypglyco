function F=trypbruc(t,y)
format long

%Model of glycolysis of the bloodstream form of 'Trypanosoma brucei'   
%Variable Equivalences
    %Glucose;
    gluc_int=y(1);
    %Glucose-6-Phosphate;
    G6P=y(2);
    %Fructose-6-Phosphate;
    F6P=y(3);
    %Fructose-1,6-Bisphosphate;
    F16BP=y(4);   
    %Glyceraldehyde-3-Phosphate;
    GA3P=y(5);
    %Dihydroxyacetone-3-Phosphate, glycosomal;
    DHAPg=y(6);
    %Dihydroxyacetone-3-Phosphate, cytosolic;
    DHAPc=y(7);
    %1,3-Bisphophoglycerate;
    PGA13=y(8);
    %3-Phosphoglycerate;
    PGA3i=y(9);
    %2-Phosphoglycerate;
    PGA2=y(10);
    %Phosphoenol Pyruvate, cytosolic;
    PEPe=y(11);
    %Pyruvate, cytosolic;
    PYR=y(12);
    %NADH;
    NADH=y(13);
    %Glycerol-3-Phosphate, glycosomal;
    GLY3Pg=y(14);
    %Glycerol-3-Phosphate, cytosolic;
    GLY3Pc=y(15);
    %Glycerol;
    GLY=y(16);
    %3-Phosphoglycerate, cytosolic;
    PGA3e=y(17);
    %Phosphoenol Pyruvate, glycosomal;
    PEPi=y(18);
    %Oxalacetate;
    OAA=y(19);
    %Pyruvate, glycosomal;
    PYRg=y(20);
    %Malate;
    MAL=y(21);
    %Fumarate;
    FUM=y(22);
    %Succinate;
    SUC=y(23);
    %1,3-Bisphophoglycerate, cytosolic;
    PGA13c=y(24);
    %Conserved sum of adenosine moieties, glycosomal
    Pg=y(25);
    %Conserved sum of adenosine moieties, cytosolic
    Pc=y(26);
    %Alanine
    ALA=y(27);
for i=1:1:27
    if y(i)<1e-10
        y(i)=0;
    end
end

%Rate modifier constants

%y37=a1;y38=a2;y39=a3;y40=a4;y41=a5;y42=a6;y43=a7=20;y44=a8;y45=a9;y46=a10;y47=a11;y48=a12;y49=a13;y50=a14;y51=a15;
%y52=a16;y53=a17;y54=a18;y55=a19;y56=a20;y57=a21;y58=a22;y59=a23;y60=a24;y61=a25;y62=a26;y63=a27;y64=a28;y65=a29;y66=a30;
%y67=a31;y68=a32;y69=a33;y70=a34;y71=a35;
    %Glucose Transporter
    u1 = y(37);
    %Hexokinase
    u2 = y(38);
    %Glucose phosphate isomerase
    u3 = y(39);
    %Phosphofructokinase
    u4 = y(40);
    %Aldolase
    u5 = y(41);
    %Triose phosphate isomerase
    u6 = y(42);
    %Glyceraldehyde-3-Phosphate Dehydrogenase
    u7 = y(43);
    %Glicerol-3-Phosphate Dehydrogenase
    u13 = y(44);
    %Glicerol-3-Phospate Oxidase
    u15 = y(45);
    %Phospoglycerate Kinase, glycosomal
    u8 = y(46);
    %Pyruvate Kinase
    u11 = y(47);
    %Pyruvate Secretion
    u12 = y(48);
    %Glycerol Kinase
    u17 = y(49);
    %ATP Utilization, cytosolic
    u19 = y(50);
    %Phosphoglycerate mutase
    u9 = y(51);
    %Enolase
    u10 = y(52);
    %DHAP transport
    u18 = y(53);%
    %Glycerol-3-phosphate transport
    u20 = y(54);%
    %Glycerol secretion
    u21 = y(55);
    %Pentose Phosphate Sink
    u22 = y(56);
    %PGA3 Transport
    u23= y(57);
    %PEP Transport
    u24= y(58);
    %Phosphoenolpyruvate carboxikinase
    u25= y(59);
    %Malate Dehydrogenase
    u26 = y(60);
    %Malic Enzyme
    u27 = y(61);
    %Fumarase
    u28 = y(62);
    %Fumarate reductase
    u29 = y(63);
    %Phosphoenolpyruvate Dikinase
    u30 = y(64);
    %Pyruvate Transport
    u31 = y(65);
    %Succinate Secretion
    u32 = y(66);
    %ATP utilization (glicosomal)
    u33 = y(67);
    %13BPGA transport
    u35 = y(68);
    %PGK, cytosolic
    u36 = y(69);
    %Pyruvate secretion (glycosomal)
    u34 =y(71);
    %Alanine Dehydrogenase
    u37 =y(70);
    
    
%Compartments

    %cell volume

    vtot=5.7;

    %Glycosomal volume

    vgly=0.25;

    %Cytosolyc Volume

    vcyt=5.45;

    

%Equilibrium Constants

    ktim=0.045;

    kpgm=0.187;

    keno=6.7;

    keq=0.442;

    keq_atp=keq;

    

%Conserved sums & Algebraic equations

    %ATPg=(-b+sqrt(b^2-4*a*c))/(2*a);

    ratio=((vcyt/vgly)+1);

    sumt=120;

    sump=3.9;

    sumpc=3.9;

    a=1-4*keq_atp;

    b=sump-Pg*(1-4*keq_atp);% 3.9-0.768*y(9)

    c=-keq_atp*Pg^2;

    bc=sumpc-Pc*(1-4*keq_atp);

    cc=-keq_atp*Pc^2; 

    ATPg=(-b+sqrt(b^2-4*a*c))/(2*a);

    ADPg=Pg-2*ATPg;

    AMPg=3.9-ATPg-ADPg;

    AXP=ADPg+ATPg+AMPg;

    ATPc=(-bc+sqrt(bc^2-4*a*cc))/(2*a);

    AMPc=6.4391-0.5*Pc-0.6510*sqrt((sump+0.768*Pc)^2-1.3578*Pc^2);

    ADPc=-5.078125000+1.302083333*sqrt(bc^2-1.357824*Pc^2);

    NAD=6-NADH;

    

%Kinetic parameters

    %Glucose Transport

        %Bloodstream form

        vmax_gt=108.9;

        rv_gt=1;

        kg=1;

        

        %Procyclic form

        %vmax_gt=7;

        %rv_gt=0.35;

        %kg=0.059;

 

    %Hexokinase

        %Bloodstream form

        vmax_hk=1929.0;

        khk_atp=0.116;

        khk_glu=0.1;

        khk_adp=0.126;

        khk_g6p=12.0;

        %Procyclic form

    %   vmax_hk=35;

    %   khk_atp=0.116;

    %   khk_glu=.1;

    %   khk_adp=.126;

    %   khk_g6p=.126;

    

    %Phosphoglucose Isomerase

    %Bloodstream form

        vmax_gpi=1305;

        rvgpi=1;

        kgpi_g6p=0.4;

        kgpi_f6p=0.12;

    %Procyclic form

%         vmax_gpi=121;

%         rvgpi=1.857;

%         kgpi_g6p=0.263;

%         kgpi_f6p=0.142;    

%     



%Phosphofructokinase

    %Bloodstream form

        vmax_pfk=1708;

        kpfk_atp=2.6E-2;

        kpfk_f6p=0.82;

        k1i=15.8;

        k2i=10.7;  

    %Procyclic form

%     vmax_pfk=462;

%     kpfk_atp=0.026;

%     pkalpha=1.05;

%     pkbeta=0.32;

%     kpfk_amp=0.079;

%     kpfk_f6p=0.82;

%     nhpfk=1.41*(1-(pkbeta*AMPg)/(kpfk_amp+AMPg));



%Aldolase

    %Bloodstream form

        vmax_ald=560;

        rvald=1.19;

        keq=6.9E-2;

        kald_f16bp=9E-3*(1+ATPg/0.68+ADPg/1.51+AMPg/3.65);

        kald_dhap=1.5E-2*(1+ATPg/0.68+ADPg/1.51+AMPg/3.65);

        kald_ga3p=6.7E-2;

        ki_ga3p=9.8E-2;

    %Procyclic form

%         vmax_ald=11;

%         rvald=1.19;

%         kald_ga3p=0.067;

%         kald_f16bp=0.02;

%         kald_dhap=0.015;

%         ki_ga3p=0.098;





%Triosephosphate Isomerase

    %Bloodstream form

    vmax_tim=999.3;

    rvtim=5.7;

    ktim_dhap=1.2;

    ktim_ga3p=0.25;

    %Procyclic form

%         vmax_tim=5920;

%         rvtim=5.7;

%         ktim_dhap=1.3;

%         ktim_ga3p=0.3;



%Glyceraldehyde-3-phosphate dehydrogenase

        %Bloodstream form

            vmax_gapdh=720.9;

            rvgapdh=0.67;

            kgapdh_ga3p=0.15;

            kgapdh_nad=0.45;

            kgapdh_pga13=0.1;

            kgapdh_nadh=0.02;

        %Procyclic form

    %     vmax_gapdh=53;

    %     kgapdh_ga3p=0.15;

    %     kgapdh_nad=0.45;

    %     rvgapdh=0.67;

    %     kgapdh_nadh=0.02;

    %     kgapdh_pga13=0.1;

%Glycerol-3-Phosphate dehydrogenase

    %Bloodstream form

        vmax_gpdh=465;

        rvgpdh=0.28;

        kgpdh_dhap=0.1;

        kgpdh_nadh=0.01;

        kgpdh_gly3p=2;

        kgpdh_nad=0.4;

    %Procyclic form

%         vmax_gpdh=214;

%         rvgpdh=0.07;

%         kgpdh_dhap=0.1;

%         kgpdh_nadh=0.015;

%         kgpdh_gly3p=6.4;

%         kgpdh_nad=0.6;

  

%Glycerolphosphate Oxidase

    %Bloodstream form

        vmax_gpo=368;

        kgpo_gly3p=1.7;



    %Procyclic form

    %     vmax_gpo=10;

    %     kgpo_gly3p=.1;

    %     kgpo_dhap=.1;



%Phosphoglycerate Kinase (Glycosomal)

    %Bloodstream Form

        vmax_pgk=2862;

        kpgk_pga13=0.003;

        rvpgk=0.47;

        kpgk_atp=0.29;

        kpgk_pga3=1.62;

        kpgk_adp=0.1;

    %Procyclic form

    %     vmax_pgk=742;

    %     kpgk_pga13=0.003;

    %     rvpgk=0.47;

    %     kpgk_atp=0.29;

    %     kpgk_pga3=1.62;

    %     kpgk_adp=0.1;



    %Piruvate Kinase

        %Bloodstream form

            vmax_pk=1020;

            kpk_pep=0.34*(1+ATPc/0.57+ADPc/0.64);

            kpk_adp=0.114;

            n=2.5;        

        %Procyclic form

    %     vmax_pk=59;

    %     kpk_pep=(0.34*(1+ATPc/0.57+ADPc/0.64));

    %     kpk_adp=0.114;

    %     npk=2.5;

    %n=1.2;

    %Glycerol Kinase

        %Bloodstream form

            vmax_gk=200;

            rvgk=60.68;

            kgk_gly3p=3.83;

            kgk_adp=0.56;

            kgk_gly=0.44;

            kgk_atp=0.24;           

        %Procyclic form

    %     vmax_gk=16.92;

    %     kgk_gly3p=3.83;

    %     kgk_adp=0.56;

    %     rvgk=60.86;

    %     kgk_gly=0.44;

    %     kgk_atp=0.24;

    %ATP utilization

    katput=50;

    katput2=1;    

    %Pyruvate Secretion

    vmax_pyrsec=200;

    %kpyrsec_pyr=1.96;

    %External concentration of glucose

    gluc_ext=y(36);

    %External concentration of glycerol

    GLY_ext=0;

    

    %Phosphoglycerate Mutase

        %Bloodsstream form

            vmax_pgam=225;

            rvpgam=2.2;

            kpgam_pga3=0.27;

            kpgam_pga2=0.11;        

        %Procyclic form

%             vmax_pgam=486;

%             rvpgam=0.461;

%             kpgam_pga3=0.15;

%             kpgam_pga2=0.16;

    %Enolase

        %Bloodstream form

            vmax_eno=598;

            rveno=0.66;

            keno_pga2=0.054;

            keno_pep=0.24;

        %Procyclic form

%             vmax_eno=486;

%             rveno=0.68;

%             keno_pga2=0.054;

%             keno_pep=0.244;

    %DHAP transport across glycosomal membrane

    vmax_tradhap=1000;

    rvtradhap=vmax_tradhap;

    ktradhap=0.1;

    %GLY3P transport across glycosomal membrane

    vmax_tragly3p=1000;

    rvtragly3p=vmax_tragly3p;

    ktragly3p=0.1;

    %GLY transport

    vmax_glyt=200;

    kglyt=.196;

    %Pentose Phosphate Pathway

    vmaxppp=0.5;

    kppp_fbp=1;

    %PGA3 Transport

    vmax_pga3t=1000;

    kpga3t=1;

    rvpga3t=vmax_pga3t;

    %PEP Transport

    vmax_pept=1000;

    kpept=0.1;

    rvpept=vmax_pga3t;

    %PEPCK

    vmf_pepck=353;

    vmr_pepck=3318;

    kpepck_pep=0.035;

    kpepck_co2=2.77;

    kpepck_atp=0.027;

    kpepck_oaa=0.041;

    kpepck_adp=0.017;

    npepck=2.06;

    co2=25;



    %Phosphoenol piruvate dikinase

    vmf_ppdk=15;

    vmr_ppi=1.986;

    kppdk_pep=0.027;

    kppdk_amp=0.009;

    kppdk_pyr=0.1;

    kppdk_atp=0.19;

    kppdk_ppi=0.089;

    ppi=0.02;



    %Malate dehydrogenase

    vmf_mdh=1057;

    vmr_mdh=144;

    kmdh_oaa=0.029;

    kmdh_mal=.252;

    kmdh_nadh=0.047;

    kmdh_nad=0.099;



    %Fumarase

    vmr_fase=100;

    vmf_fase=100;

    kfase_mal=0.15;

    kfase_fum=0.63;



    %Fumarato Reductase

    vmf_fumrd=41;

    kfumrd_fum=.41;

    kfumrd_nadh=0.02;

    vmr_fumrd=.174;

    kfumrd_suc=.047;

    kfumrd_nad=.099;



    %PYR transport

    vmax_pyrt=1000;

    kpyrt=.1;

    rvpyrt=vmax_pyrt;

    

    %SUC transport

    vmax_suct=200;

    ksuct=1.96;

    

    %Adenilate Kinase, cytosolic

    kakg=10000;



    %PGK, cytosolic

    vmax_pgkc=742;

    kpgkc_adp=0.1;

    kpgkc_pga13=0.003;

    vmr_pgkc=349;

    kpgkc_pga3=1.62;

    kpgkc_atp=0.29;

    %Malic Enzyme

    vmax_mal=49;

    vmr_mal=0.061;

    kmal_pyr=4.4;

    kmal_mal=0.2;

    kmal_co2=0.1;

    %Alanine Dehydrogenase (Streptomyces fradiae)
    vmax_aladh=13.4;
    vmaxr_aladh=0.95;
    kaladh_pyr=0.23;
    kaladh_nadh=0.05;
    kaladh_ala=10;
    kaladh_nad=0.18;
    kaladh_nh3=11.6;
    NH3=0.5;
    
    %Alanine Secretion
    vmax_alasec=50;

%Rate Equations:

    %Glucose Transport

    %Bloodstream form

    vgt = vmax_gt*(gluc_ext-gluc_int)/((kg+gluc_ext+gluc_int+.75*gluc_ext*gluc_int/kg))*u1;

    %Procyclic form

    %vgt =

    %((vmax_gt*gluc_ext-rv_gt*gluc_int)/kg)/((1+gluc_ext/kg+gluc_int/kg))*u1;

    %Hexokinase

    vhk = vmax_hk*(ATPg/khk_atp*gluc_int/khk_glu)/((1+ATPg/khk_atp+ADPg/khk_adp)*(1+gluc_int/khk_glu+G6P/khk_g6p))*u2;

    %Phosphoglucose Isomerase

    vpgi = (vmax_gpi*(G6P/kgpi_g6p-rvgpi*F6P/kgpi_f6p))/(1+G6P/kgpi_g6p+F6P/kgpi_f6p)*u3;

    %Phosphofructokinase

    %Bloodstream form

    vpfk = (vmax_pfk*(k1i/(k1i+F16BP))*(F6P*ATPg/(kpfk_f6p*kpfk_atp)))/((1+F6P/kpfk_f6p+F16BP/k2i)*(1+ATPg/kpfk_atp))*u4;

    %Procyclic form

%     vpfk = (vmax_pfk*(((F6P/kpfk_f6p))^nhpfk)*ATPg)/(kpfk_atp*(1+(sqrt((F6P/kpfk_f6p)^2))^nhpfk)*(1+ATPg/kpfk_atp))*u4;

  

    %Aldolase

        %Bloodstream form

        vald = (vmax_ald*(F16BP-GA3P*DHAPg/keq))/(kald_f16bp+F16BP+kald_dhap*GA3P/(rvald*keq)+kald_ga3p*DHAPg/(rvald*keq)+GA3P*DHAPg/(rvald*keq)+F16BP*GA3P/ki_ga3p)*u5;

        %Procyclic form

        %vald = (vmax_ald*(F16BP/kald_f16bp-GA3P*DHAPg*rvald/(kald_dhap*kald_ga3p)))/(1+F16BP/kald_f16bp+GA3P/kald_ga3p+DHAPg/kald_dhap+F16BP*GA3P/(kald_f16bp*ki_ga3p)+GA3P*DHAPg/(kald_dhap*kald_ga3p))*u5;

    %Triosephosphate Isomerase

    vtpi = vmax_tim*(DHAPg/ktim_dhap-rvtim*GA3P/ktim_ga3p)/(1+DHAPg/ktim_dhap+GA3P/ktim_ga3p)*u6;

    %Glyceraldehyde-3-phosphate dehydrogenase

    vgapdh = vmax_gapdh*(GA3P*NAD/(kgapdh_ga3p*kgapdh_nad)-rvgapdh*PGA13*NADH/(kgapdh_nadh*kgapdh_pga13))/((1+GA3P/kgapdh_ga3p+PGA13/kgapdh_pga13)*(1+NAD/kgapdh_nad+NADH/kgapdh_nadh))*u7;

    %Glycerol-3-Phosphate dehydrogenase

    vgpdh = vmax_gpdh*(DHAPg/kgpdh_dhap*NADH/kgpdh_nadh-rvgpdh*GLY3Pg/kgpdh_gly3p*NAD/kgpdh_nad)/((1+DHAPg/kgpdh_dhap+GLY3Pg/kgpdh_gly3p)*(1+NADH/kgpdh_nadh+NAD/kgpdh_nad))*u13;

    %Glycerolphosphate Oxidase

    vgo = (vmax_gpo*GLY3Pc/kgpo_gly3p)*u15;

    %vgo = (vmax_gpo*GLY3Pc/kgpo_gly3p-rvgpo*DHAPc/kgpo_dhap)/(1+GLY3Pc/kgpo_gly3p+DHAPc/kgpo_dhap)*u15;

    %Phosphoglycerate Kinase (Glycosomal)

    vpgk= vmax_pgk*(PGA13*ADPg/(kpgk_pga13*kpgk_adp)-rvpgk*PGA3i*ATPg/(kpgk_atp*kpgk_pga3))/((1+PGA13/kpgk_pga13+PGA3i/kpgk_pga3)*(1+ADPg/kpgk_adp+ATPg/kpgk_atp))*u8;

    %Piruvate Kinase

        %Bloodstream form

        vpyk = (vmax_pk*((PEPe/kpk_pep)^n)*ADPc/kpk_adp)/(((1+(PEPe/kpk_pep)^n))*(1+ADPc/kpk_adp))*u11;

        %Procyclic form

        %vpyk = (vmax_pk*((PEPe/kpk_pep))^npk)*ADPc/kpk_adp)/(((1+(PEPe/kpk_pep)^npk))*(1+ADPc/kpk_adp))*u11;

    %Glycerol Kinase

    vgk = vmax_gk*(GLY3Pg/kgk_gly3p*ADPg/kgk_adp-rvgk*GLY/kgk_gly*ATPg/kgk_atp)/((1+GLY3Pg/kgk_gly3p+GLY/kgk_gly)*(1+ADPg/kgk_adp+ATPg/kgk_atp))*u17;

    %ATP utilization

    vatput = (katput*ATPc/ADPg)*u19;

    vatput2 = (katput2*ATPg/ADPg)*u33;

    %Pyruvate Secretion

    vpyrsec = vmax_pyrsec*PYR*u12;

    %Phosphoglycerate Mutase

    vpgam = vmax_pgam*(PGA3e/kpgam_pga3-rvpgam*PGA2/kpgam_pga2)/(1+PGA3e/kpgam_pga3+PGA2/kpgam_pga2)*u9;

    %Enolase

    veno = vmax_eno*(PGA2/keno_pga2-rveno*PEPe/keno_pep)/(1+PGA2/keno_pga2+PEPe/keno_pep)*u10;

    %DHAP transport

    vdhapt = (vmax_tradhap*DHAPg/ktradhap-rvtradhap*DHAPc/ktradhap)/(1+DHAPg/ktradhap+DHAPc/ktradhap)*u18;

    %GLY3P transport

    vgly3pt = (vmax_tragly3p*GLY3Pg/ktragly3p-rvtragly3p*GLY3Pc/ktragly3p)/(1+GLY3Pc/ktragly3p+GLY3Pg/ktragly3p)*u20;

    %Glycerol transport

    vglyt=(vmax_glyt*GLY/kglyt-vmax_glyt/5*GLY_ext/kglyt)/(1+GLY/kglyt+GLY_ext/kglyt)*u21;

    %Pentose Phosphate sink

    %vppp=vmaxppp*F16BP/kppp_fbp/(1+F16BP/kppp_fbp)*u22;

    vppp=vgt*u22*G6P;

    %PGA3 Transport

    vpga3t=(vmax_pga3t*PGA3i/kpga3t-rvpga3t*PGA3e/kpga3t)/(1+PGA3e/kpga3t+PGA3i/kpga3t)*u23;

    %PEP Transport

    vpept=(vmax_pept*PEPe/kpept-rvpept*PEPi/kpept)/(1+PEPe/kpept+PEPi/kpept)*u24;

    %Phosphoenol Pyruvate Carboxikinase

    vpepck=((vmf_pepck*co2*ADPg*PEPi)/(kpepck_co2*kpepck_adp*kpepck_pep)-(vmr_pepck*ATPg*((OAA/kpepck_oaa))^npepck))/((1+PEPi/kpepck_pep+((OAA/kpepck_oaa))^npepck)*(1+ATPg/kpepck_atp+ADPg/kpepck_adp)*(1+co2/kpepck_co2))*u25;

    %Malate Dehydrogenase

    vmdh=(vmf_mdh*OAA*NADH/(kmdh_oaa*kmdh_nadh)-vmr_mdh*MAL*NAD/(kmdh_mal*kmdh_nad))/((1+OAA/kmdh_oaa+MAL/kmdh_mal)*(1+NADH/kmdh_nadh+NAD/kmdh_nad))*u26;

    %Malic Enzyme

    vmalic=(vmax_mal*MAL/kmal_mal-vmr_mal*PYRg*co2/(kmal_pyr*kmal_co2))/((1+MAL/kmal_mal+PYRg/kmal_pyr)*(1+co2/kmal_co2))*u27;

    %Fumarase

    vfase=(vmf_fase*MAL/kfase_mal-vmr_fase*FUM/kfase_fum)/(1+MAL/kfase_mal+FUM/kfase_fum)*u28;

    %Fumarate reductase

    vfumrd=(vmf_fumrd*FUM*NADH/(kfumrd_fum*kfumrd_nadh)-vmr_fumrd*SUC*NAD/(kfumrd_suc*kfumrd_nad))/((1+FUM/kfumrd_fum+SUC/kfumrd_suc)*(1+NADH/kfumrd_nadh+NAD/kfumrd_nad))*u29;

    %Phosphoenolpyruvate Dikinase

    vppdk=(vmf_ppdk*PEPi*AMPg*ppi/(kppdk_pep*kppdk_amp*kppdk_ppi)-vmr_ppi*PYRg*ATPg/(kppdk_pyr*kppdk_atp))/((1+PEPi/kppdk_pep+PYRg/kppdk_pyr)*(1+AMPg/kppdk_amp+ATPg/kppdk_atp)*(1+ppi/kppdk_ppi))*u30;

    %Pyruvate transport

    vpyrt=(vmax_pyrt*PYR/kpyrt-rvpyrt*PYRg/kpyrt)/(1+PYRg/kpyrt+PYR/kpyrt)*u31;

    vpyrt2=(0*PYRg/ksuct)/(1+PYRg/ksuct)*u31;

    %Succinate transport

    vsucsec=vmax_suct*SUC/ksuct/(1+SUC/ksuct)*u32;

    %13BPGA transport

    vpga13t=(vmax_pga3t*PGA13/kpga3t-rvpga3t*PGA13c/kpga3t)/(1+PGA13/kpga3t+PGA13c/kpga3t)*u35;

    %Cytosolic PGK

    vpgkc=(vmax_pgkc*PGA13c*ADPc/(kpgkc_adp*kpgkc_pga13)-vmr_pgkc*ATPc*PGA3e/(kpgkc_pga3*kpgkc_atp))/((1+PGA3e/kpgkc_pga3+ADPc/kpgkc_adp)*(1+PGA13c/kpgkc_pga13+ATPc/kpgkc_atp))*u36;

    %Alanine Dehydrogenase
    valadh=((vmax_aladh*PYRg*NADH*NH3/(kaladh_pyr*kaladh_nadh*kaladh_nh3)-vmaxr_aladh*ALA*NAD/(kaladh_ala*kaladh_nad)))/((1+PYRg/kaladh_pyr+ALA/kaladh_ala)*(1+NADH/kaladh_nadh+NAD/kaladh_nad)*(1+NH3/kaladh_nh3))*u37;
    %Alanine Secretion
    valasec = vmax_alasec*ALA;

%Diferential Equations    

%Glucose;

    dglu=(vgt-vhk)/vtot;

%Glucose-6-Phosphate;

    dg6p = (vhk-vpgi-vppp)/vgly;

%Fructose-6-Phosphate;

    df6p = (vpgi-vpfk)/vgly;

%Fructose-1,6-Bisphosphate;

    df16bp=(vpfk-vald)/vgly;

%Glyceraldehyde-3-Phosphate;

    dga3p=(vald-vgapdh+vtpi)/vgly;

%Dihydroxyacetone-3-Phosphate, glycosomal;

    ddhapg=(vald-vgpdh-vdhapt-vtpi)/vgly;

%Dihydroxyacetone-3-Phosphate, cytosolic;

    ddhapc=(vdhapt+vgo)/vcyt;

%1,3-Bisphophoglycerate;

    dbpga=(vgapdh-vpgk-vpga13t)/vgly;

%3-Phosphoglycerate;

    dpga3i=(vpgk-vpga3t)/vgly;

%2-Phosphoglycerate;

    dpga2=(vpgam-veno)/vcyt;

%Phosphoenol Pyruvate, cytosolic;

    dpepe=(veno-vpyk-vpept)/vcyt;

%Pyruvate, cytosolic;

    dpyr=(vpyk-vpyrsec+vpyrt2)/vcyt;

%NADH;

    dnadh=(vgapdh-vgpdh-vmdh-vfumrd-valadh)/vgly;

%Glycerol-3-Phosphate, glycosomal;

    dgly3pg=(vgpdh-vgk-vgly3pt)/vgly;

%Glycerol-3-Phosphate, cytosolic;

    dgly3pc=(-vgo+vgly3pt)/vcyt;

%Glycerol;

    dgly=(vgk-vglyt)/vgly;;

%3-Phosphoglycerate, cytosolic;

    dpga3e=(vpga3t+vpgkc-vpgam)/vcyt;

%Phosphoenol Pyruvate, glycosomal;

    dpepi=(vpept-vppdk-vpepck)/vgly;

%Oxalacetate;

    doaa=(vpepck-vmdh)/vgly;

%Pyruvate, glycosomal;

    dpyrg=(vppdk-vpyrt2+vmalic-valadh)/vgly;

%Malate;

    dmal=(vmdh-vfase-vmalic)/vgly;

%Fumarate;

    dfum=(vfase-vfumrd)/vgly;

%Succinate;

    dsuc=(vfumrd-vsucsec)/vgly;

%1,3-Bisphophoglycerate, cytosolic;

    dbpgac=(vpga13t-vpgkc)/vcyt;

%d(Pg)/dt

    dpg=(-vhk-vpfk+vpgk+vgk+vpepck+2*vppdk-vatput2)/vgly;

%d(Pc)/dt

    dpc=(vpyk-vatput-vpgkc)/vcyt;

%Alanine
    dala=(valadh-valasec)/vgly;

F = [dglu,dg6p,df6p,df16bp,dga3p,ddhapg,ddhapc,dbpga,dpga3i,dpga2,dpepe,dpyr,dnadh,dgly3pg,dgly3pc,dgly,dpga3e,dpepi,doaa,dpyrg,dmal,dfum,dsuc,dbpgac,dpg,dpc,dala,vgt,veno,vpept,vpyrsec,vsucsec,vatput2,vppp,vglyt,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,vppdk,valadh,vpyrt2,vpepck,vgapdh]';
%t,y';
