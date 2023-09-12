clc
clf
clear
format shortg
n_reponses=9;
data=load('data.txt');
%transformation FWHM
% data(:,20)=log(data(:,20))

%---------------------------------------------------------------------------------------------------------------
disp('Répétitions au centre');
X1 = data([1:3,20],:)

% Etape 1 variabilité
disp('Moyennes');
y0=mean(X1(:,17:17+n_reponses-1))
disp('Ecarts-type');
sigma=std(X1(:,17:17+n_reponses-1));
S0=2*std(X1(:,17:17+n_reponses-1))

disp('Matrice screening H16 + réponses');
X2=data(4:19,:)

disp('Matrice de calcul');
EXP=[ones(16,1) X2(:,(2:16))] % ajout colonne de 1 pour constante
H16=hadamard(16)

% pour les w réponses
for w=1:n_reponses %nombre de réponses
   for v=1:16 %nombre de paramètres
        coeffR(v,w)=sum(EXP(:,v).*X2(:,16+w))/16;
   end  
end
disp('Matrice de coefficients');
coeffR
   
for i=1:n_reponses;
figure(i)
hold on
bar ((coeffR((2:end),i)),'y')
plot([0 16], [S0(i) S0(i)],'b')
plot([0 16], [-S0(i) -S0(i)],'b')

if i==1;
title('Raman shift (cm-1)');
end
if i==2;
title('FWHM Raman shift E2(h) (cm-1)');
end
if i==3;
title('1/Curvature (m-1)');
end
if i==4;
title('FWHM (0002) (arcsec)');
end
if i==5;
title('Thickness (nm)');
end
if i==6;
title('Roughness RMS (nm)');
end
if i==7;
title('Roughness RA (nm)');
end
if i==8;
title('Island size (m, from curvature)');
end
if i==9;
title('Island size (m, from Raman)');
end

xlabel('Parameters');
ylabel('Parameter weigth');

text(1,0,'Pellet age','FontSize',24,'Rotation',-90*sign(coeffR(2,i)))
text(2,0,'Cooling down gas','FontSize',24,'Rotation',-90*sign(coeffR(3,i)))
text(3,0,'Mounting delay','FontSize',24,'Rotation',-90*sign(coeffR(4,i)))
text(4,0,'Pellet size','FontSize',24,'Rotation',-90*sign(coeffR(5,i)))
text(5,0,'Sequential or simultaneous','FontSize',24,'Rotation',-90*sign(coeffR(6,i)))
text(6,0,'Ammonia flow rate','FontSize',24,'Rotation',-90*sign(coeffR(7,i)))
text(7,0,'Purge time','FontSize',24,'Rotation',-90*sign(coeffR(8,i)))
text(8,0,'Pre-chlorination time','FontSize',24,'Rotation',-90*sign(coeffR(9,i)))
text(9,0,'Etching at 1100°C','FontSize',24,'Rotation',-90*sign(coeffR(10,i)))
text(10,0,'Pellet cleaning (hydrogen)','FontSize',24,'Rotation',-90*sign(coeffR(11,i)))
text(11,0,'Hydrogen flow rate','FontSize',24,'Rotation',-90*sign(coeffR(12,i)))
text(12,0,'Furnace power','FontSize',24,'Rotation',-90*sign(coeffR(13,i)))
text(13,0,'Cooling down rate','FontSize',24,'Rotation',-90*sign(coeffR(14,i)))
text(14,0,'Temperature','FontSize',24,'Rotation',-90*sign(coeffR(15,i)))
text(15,0,'Pressure','FontSize',24,'Rotation',-90*sign(coeffR(16,i)))

text(16,S0(i),'2\sigma','FontSize',24);
text(16,-S0(i),'-2\sigma','FontSize',24);
text(16,0,num2str(coeffR(1,i),4),'FontSize',24);

set(gca,'FontSize',24)
%cas particulier de la FWHM
if i==4;
X1_alt = data(1:3,:);
        
% Etape 1 variabilité
disp('Moyennes');
y0_alt=mean(X1_alt(:,17:17+n_reponses-1))
disp('Ecarts-type');
S0_alt=2*std(X1_alt(:,17:17+n_reponses-1))
text(16,S0_alt(i),'2\sigma','FontSize',24);
text(16,-S0_alt(i),'-2\sigma','FontSize',24);
plot([0 16], [S0_alt(i) S0_alt(i)],'--b')
plot([0 16], [-S0_alt(i) -S0_alt(i)],'--b')
end

hold off
saveas(gcf,['Response_',num2str(i),'.png']);
end;

for i=1:n_reponses;
figure(i+10)
plot(data(:,1),data(:,i+16),'bd');

if i==1;
ylabel('Raman shift cm-1');
end
if i==2;
ylabel('FWHM Raman shift E2(h) cm-1');
end
if i==3;
ylabel('1/Curvature m-1');
end
if i==4;
ylabel('FWHM (0002) arcsec');
end
if i==5;
ylabel('Thickness nm');
end
if i==6;
ylabel('Roughness nm RMS');
end
if i==7;
ylabel('Roughness nm RA');
end
if i==8;
ylabel('Island size m (Curvature)');
end
if i==9;
ylabel('Island size m (Raman)');
end
xlabel('Experiment order');

set(gca,'FontSize',16)
saveas(gcf,['Ordre_',num2str(i),'.png']);
end;

%---------------------------------------------------------------------------------------------------------------
disp('Expérience -1');
X3 =data(21,:)

for i=1:n_reponses
disp(['Expérience -1, calcul de la réponse ' num2str(i)]);
 R= sum(coeffR((1:end),i)'.*[1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]);
 disp('Modèle - Experience - modèle-expérience - 2*écart-type');
 [R  X3(1,i+16) R-X3(1,i+16) S0(i)]
 if abs(R-X3(1,i+16))<S0(i)
disp(['Absence d''intéractions pour la réponse ', num2str(i)]);
disp(['Répétions au centre, calcul de la réponse ' num2str(i)]);
 R= sum(coeffR((1:end),i)'.*[1 1 1 1 -1 -1 0 0 0 0 0 0 0 0 0 0]);
 disp('Modèle - répet au centre - Modèle-répet au centre - 2*écart-type');
 [R  y0(1,i) R-y0(1,i) S0(i)]
 if abs(R-y0(1,i))<S0(i)
disp(['Absence d''effets quadratiques pour la réponse ', num2str(i)]);
else
disp(['Effets quadratiques suspectés pour la réponse ', num2str(i)]);    
 end
else
disp(['Intéractions suspectées pour la réponse ', num2str(i)]);    
end
end;


% simplification du modèle
for i=1:n_reponses;
    for j=1:1:21
    
    rep_predite(j)=sum(coeffR((1:end),i)'.*[1,data(j,2:16)]);
    
    end
figure(i+20)
hold on
plot(data(:,i+16),rep_predite,'bd');

if i==1;
title('Raman shift cm-1');
end
if i==2;
title('FWHM Raman shift E2(h) cm-1');
end
if i==3;
title('1/Curvature m-1');
end
if i==4;
title('FWHM (0002) arcsec');
end
if i==5;
title('Thickness nm');
end
if i==6;
title('Roughness nm RMS');
end
if i==7;
title('Roughness nm RA');
end
if i==8;
title('Island size m (Curvature)');
end
if i==9;
title('Island size m (Raman)');
end
ylabel('Predicted');
xlabel('Experimental');
set(gca,'FontSize',16)
hold off
saveas(gcf,['Simulated_',num2str(i),'.png']);
end;
% graphe des effets
hadam_matrix=X2(:,2:end);
for j=1:1:n_reponses;
for i =1:1:15;
    mat_effect=sortrows(hadam_matrix,i);
    moins(i,j)=mean(mat_effect(1:8,15+j));
    plus(i,j)=mean(mat_effect(9:16,15+j));
end;
end;

for i=1:1:n_reponses;
y=[moins(:,i), plus(:,i)];    
% y=[moins(:,j), ones(16,1)*(y0(j)), plus(:,j)];
figure(i+30)
bar(y)
if i==1;
title('Raman shift cm-1');
end
if i==2;
title('FWHM Raman shift E2(h) cm-1');
end
if i==3;
title('1/Curvature m-1');
end
if i==4;
title('FWHM (0002) arcsec');
end
if i==5;
title('Thickness nm');
end
if i==6;
title('Roughness nm RMS');
end
if i==7;
title('Roughness nm RA');
end
if i==8;
title('Island size m Curvature)');
end
if i==9;
title('Island size m (Raman)');
end
text(1,0,'Pellet age','FontSize',16,'Rotation',90)
text(2,0,'Cooling down gas','FontSize',16,'Rotation',90)
text(3,0,'Mounting delay','FontSize',16,'Rotation',90)
text(4,0,'Pellet size','FontSize',16,'Rotation',90)
text(5,0,'Sequential or simultaneous','FontSize',16,'Rotation',90)
text(6,0,'Ammonia flow rate','FontSize',16,'Rotation',90)
text(7,0,'Purge time','FontSize',16,'Rotation',90)
text(8,0,'Pre-chlorination time','FontSize',16,'Rotation',90)
text(9,0,'Etching at 1100°C','FontSize',16,'Rotation',90)
text(10,0,'Pellet cleaning under hydrogen','FontSize',16,'Rotation',90)
text(11,0,'Hydrogen flow rate','FontSize',16,'Rotation',90)
text(12,0,'Furnace power','FontSize',16,'Rotation',90)
text(13,0,'Cooling down rate','FontSize',16,'Rotation',90)
text(14,0,'Temperature','FontSize',16,'Rotation',90)
text(15,0,'Pressure','FontSize',16,'Rotation',90)
set(gca,'FontSize',16)
saveas(gcf,['Effect_',num2str(i),'.png']);
end;

%calcul de la contrainte selon stoney
Y_AlN=340*1E9/(1-0.21); %Pa
Y_saphir=431*1E9/(1-0.285); %Pa
ts=430e-6;
tf=data(:,21)*1E-9;
un_sur_r=data(:,19);
raman_shift=data(:,17);
stress_stoney=1E-9*Y_saphir.*ts.^2.*un_sur_r./(6.*tf);%GPa
stress_raman=(raman_shift-657.67)/-4.04;%GPa
t_griffith_stoney=2.*Y_AlN.*0.4./(pi.*(stress_stoney.*1E9).^2)*1e9; %nm
t_griffith_raman=2.*Y_AlN.*0.4./(pi.*(stress_raman.*1E9).^2)*1e9; %nm
casse=tf*1E9>t_griffith_stoney;
casse_maeva=[1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 0 1 0 1 1 0]';
figure(41)
scatter(stress_stoney,stress_raman,[],casse);
colormap jet
%title('Stoney vs Raman');
xlabel('Stress from curvature (GPa)');
ylabel('Stress from raman (GPa)');
set(gca,'FontSize',16);
saveas(gcf,'Stoney_vs_raman.png');

temperature=1250+data(:,15)*50;
figure(42);
plot(temperature, stress_stoney, 'bd');
%title('Stoney vs Temperature');
xlabel('Temperature °C');
ylabel('Stress Stoney (GPa)');
set(gca,'FontSize',16);
saveas(gcf,'Stoney_vs_Temperature.png');

figure(43);
plot(temperature, raman_shift, 'bd');
%title('Raman Shift vs Temperature');
xlabel('Temperature °C');
ylabel('Raman shift (cm-1)');
set(gca,'FontSize',16);
saveas(gcf,'Raman_vs_Temperature.png');

figure(44);
scatter(temperature, un_sur_r, [], casse_maeva);
colormap jet
%title('1/r vs Temperature');
xlabel('Temperature °C');
ylabel('1/r (m-1)');
set(gca,'FontSize',16);
saveas(gcf,'1_r_vs_Temperature.png');

figure(45);
scatter(raman_shift, un_sur_r, [], casse_maeva);
colormap jet
%title('Raman shift vs 1/r');
xlabel('Raman shift (cm-1)');
ylabel('1/r (m-1)');
set(gca,'FontSize',16);
saveas(gcf,'1_r_vs_Raman.png');

% fissuration=[casse casse_maeva tf*1E9-t_griffith_stoney];

Temp=[1200 1250 1300];
raman_stress_T=[mean(stress_raman([6 7 8 9 12 13 18 19 21])),mean(stress_raman([1 2 3 20])),mean(stress_raman([4 5 10 11 14 15 16 17]))];
stoney_stress_T=[mean(stress_stoney([6 7 8 9 12 13 18 19 21])),mean(stress_stoney([1 2 3 20])),mean(stress_stoney([4 5 10 11 14 15 16 17]))];

ilots_curvature=data(:,24);
ilots_raman=data(:,25);
ilots_curvature_T=[mean(ilots_curvature([6 7 8 9 12 13 18 19 21])),mean(ilots_curvature([1 2 3 20])),mean(ilots_curvature([4 5 10 11 14 15 16 17]))];
A=mean(ilots_curvature([6 7 8 9 12 13 18 19 21]));
B=mean(ilots_curvature([1 2 3 20]));
C=mean(ilots_curvature([4 5 10 11 14 15 16 17]));
disp('ilots courbure');
Temp;
[A B C]
A=std(ilots_curvature([6 7 8 9 12 13 18 19 21]));
B=std(ilots_curvature([1 2 3 20]));
C=std(ilots_curvature([4 5 10 11 14 15 16 17]));
[A B C]
A=mean(ilots_raman([6 7 8 9 12 13 18 19 21]));
B=mean(ilots_raman([1 2 3 20]));
C=mean(ilots_raman([4 5 10 11 14 15 16 17]));
disp('ilots Raman');
Temp;
[A B C]
A=std(ilots_raman([6 7 8 9 12 13 18 19 21]));
B=std(ilots_raman([1 2 3 20]));
C=std(ilots_raman([4 5 10 11 14 15 16 17]));
[A B C]

FWHM=data(:,20);

figure(46)
scatter(ilots_curvature, FWHM, [], casse_maeva);
colormap jet
%title('FWHM vs island size (curvature)');
xlabel('Island size (m)');
ylabel('FWHM (0002) arcsec');
set(gca,'FontSize',16);
saveas(gcf,'FWHM_vs_island_curv.png');

figure(47)
scatter(ilots_raman, FWHM, [], casse_maeva);
colormap jet
%title('FWHM vs island size (Raman)');
xlabel('Island size (m)');
ylabel('FWHM (0002) arcsec');
set(gca,'FontSize',16);
saveas(gcf,'FWHM_vs_island_raman.png');

figure(49)
hold on
%H16
%errorbar(Temp, [3.5698e-08   3.4705e-08    3.822e-08],[2.5516e-09   2.2967e-09   3.6596e-09], 'bd');
errorbar(Temp, [3.1303e-08   3.6131e-08   3.4496e-08],[1.0249e-08   2.5055e-09   6.9551e-09], 'kd');
% %H8
% errorbar([1150 1200 1250]+4, [3.5698e-08   3.4705e-08    3.822e-08],[2.5516e-09   2.2967e-09   3.6596e-09], 'gd');
% errorbar([1150 1200 1250]+6, [2.6998e-08    2.462e-08    3.497e-08],[8.3381e-09   3.3734e-09   9.3363e-09], 'kd');

%balaji et Wu
errorbar([650 850 1100 1300]+8, [17.9e-9 18.9e-9 25e-9 35e-9], [5e-9 5e-9 0 0], 'kd');

colormap jet
%title('Island size (curvature) vs T');
xlabel('T(°C)');
ylabel('Island size (m)');
set(gca,'FontSize',16);
saveas(gcf,'T_vs_island_curv.png');
hold off

FWHM_raman=data(:,18);
figure(50)
plot(FWHM_raman, FWHM, 'bd');
%title('FWHM Raman vs FWHM (0002)');
xlabel('FWHM Raman E2(h) cm-1');
ylabel('FWHM (0002) arcsec');
set(gca,'FontSize',16);
saveas(gcf,'FWHM_vs_FWHM.png');

figure(51)
scatter(ilots_raman, ilots_curvature, [], casse_maeva);
colormap jet
%title('Island size (Raman) vs Island size (curvature)');
xlabel('Island size (Raman, m)');
ylabel('Island size (curvature, m)');
set(gca,'FontSize',16);
saveas(gcf,'ilots_vs_ilots.png');

figure(52)
scatter(stress_raman, FWHM, [], casse_maeva);
colormap jet
%title('Stress (Raman) vs FWHM (0002)');
xlabel('Stress (Raman) GPa');
ylabel('FWHM (0002) arsec');
set(gca,'FontSize',16);
saveas(gcf,'raman_vs_FWHM.png');

figure(53)
scatter(stress_stoney, FWHM, [], casse_maeva);
colormap jet
%title('Stress (Curvature) vs FWHM (0002)');
xlabel('Stress (curvature) GPa');
ylabel('FWHM (0002) arsec');
set(gca,'FontSize',16);
saveas(gcf,'stoney_vs_FWHM.png');

figure(54)
scatter(tf, FWHM, [], casse_maeva);
colormap jet
%title('Thickness vs FWHM (0002)');
xlabel('Thickness (m)');
ylabel('FWHM (0002) arsec');
set(gca,'FontSize',16);
saveas(gcf,'Thickness_vs_FWHM.png');

figure(60);
hold on
error_ilots_Raman=std(ilots_raman([1:3,20]));
error_ilots_Curvature=std(ilots_curvature([1:3,20]));
errorbarxy(ilots_raman, ilots_curvature,error_ilots_Raman*ones(21), error_ilots_Curvature*ones(21),'k');
scatter(ilots_raman, ilots_curvature, [], casse_maeva, 'filled', 'MarkerEdgeColor','k');
s.LineWidth = 10;
s.MarkerEdgeColor = 'k';
colormap gray
hold off

%title('Island size (Raman) vs Island size (curvature)');
xlabel('Island size Raman shift (m)');
ylabel('Island size curvature (m)');
set(gca,'FontSize',16);
axis([1e-8 7e-8 1e-8 7e-8]);
saveas(gcf,'ilots_vs_ilots_error.png');


figure(61);
hold on
error_stress_Raman=std(stress_raman([1:3,20]));
error_FWHM=std(FWHM([1:3,20]));
errorbarxy(stress_raman, FWHM,error_stress_Raman*ones(21), error_FWHM*ones(21),'k');
scatter(stress_raman, FWHM, [], casse_maeva, 'filled', 'MarkerEdgeColor','k');
s.LineWidth = 10;
s.MarkerEdgeColor = 'k';
colormap gray
hold off
axis([-1 3.5 0 6000]);
%title('Stress (Raman) vs FWHM (0002)');
xlabel('Stress from Raman shift (GPa)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'raman_vs_FWHM_error.png');

figure(62);
hold on
error_stress_stoney=std(stress_stoney([1:3,20]));
error_FWHM=std(FWHM([1:3,20]));
errorbarxy(stress_stoney, FWHM,error_stress_stoney*ones(21), error_FWHM*ones(21),'k');
scatter(stress_stoney, FWHM, [], casse_maeva, 'filled', 'MarkerEdgeColor','k');
s.LineWidth = 10;
s.MarkerEdgeColor = 'k';
colormap gray
hold off
axis([-1 1 0 6000]);
%title('Stress (Curvature) vs FWHM (0002)');
xlabel('Stress from Curvature (GPa)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'raman_vs_Stoney_error.png');

figure(63);
hold on
error_tf=std(tf([1:3,20]));
error_FWHM=std(FWHM([1:3,20]));
errorbarxy(tf, FWHM,error_tf*ones(21), error_FWHM*ones(21),'k');
scatter(tf, FWHM, [], casse_maeva, 'filled', 'MarkerEdgeColor','k');
s.LineWidth = 10;
s.MarkerEdgeColor = 'k';
colormap gray
hold off
axis([0 2.5e-6 0 6000]);
%title('Thickness vs FWHM (0002)');
xlabel('Thickness (m)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'Thickness_vs_FWHM_error.png');

figure(64)
hold on
error_ilots_Raman=std(ilots_raman([1:3,20]));
error_FWHM=std(FWHM([1:3,20]));
errorbarxy(ilots_raman, FWHM,error_ilots_Raman*ones(21), error_FWHM*ones(21),'k');
scatter(ilots_raman, FWHM, [], casse_maeva, 'filled', 'MarkerEdgeColor','k');
s.LineWidth = 10;
s.MarkerEdgeColor = 'k';
colormap gray
hold off
axis([1e-8 5e-8 0 6000]);
%title('FWHM vs island size (Raman)');
xlabel('Island size Raman shift (m)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'FWHM_vs_island_raman_error.png');

time=load('time.txt');
figure(65)
plot(time,FWHM,'kd')
%title('FWHM vs time');
xlabel('time (min)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'FWHM vs time.png');

time=load('time.txt');
figure(66)
plot((tf*1e6)./(time/60),FWHM,'kd')
%title('GR vs time');
xlabel('Growth rate (µm/h)');
ylabel('FWHM 0002 (arcsec)');
set(gca,'FontSize',16);
saveas(gcf,'GR vs time.png');

%création du tableau résumé analysé avec regstats__________________________
input_reg=data([4:19,21],2:16);
outputs_reg=data([4:19,21],17:25);

whichstats = {'tstat','beta','fstat'};
for i = 1:1:n_reponses;
out=regstats(outputs_reg(:,i),input_reg,'linear',whichstats);
reg_p_value_table(:,i)=out.tstat.pval;
reg_coeff_table(:,i)=out.beta;
end;
figure(67)
imagesc(reg_p_value_table)
colorbar
colormap jet

%création du tableau résumé analysé avec la règle des sigmas_______________

for i = 1:1:n_reponses;
    for j=2:1:16
    if abs(coeffR(j,i))<=2*sigma(i); influence(j-1,i)=0;end;
    if abs(coeffR(j,i))>2*sigma(i) && abs(coeffR(j,i))<=3*sigma(i); influence(j-1,i)=1*sign(coeffR(j,i));end;
    if abs(coeffR(j,i))>3*sigma(i); influence(j-1,i)=2*sign(coeffR(j,i));end;
    end;
end;
figure(68)
imagesc(influence)
colorbar
colormap jet
inverted=influence';
figure (69)
imagesc(inverted)
colorbar
colormap jet