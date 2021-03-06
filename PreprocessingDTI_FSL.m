function NomFichierNiiX = PreprocessingDTI_FSL(Indir,Outdir,NomFichierNii,B0dir)
% function NomFichierNiiX = PreprocessingDTI_FSL(Indir,Outdir,NomFichierNii,B0dir, idx1, idx2)
% FSLcommand='/usr/bin/fsl5.0-';
FSLcommand='/usr/share/fsl/5.0/bin/';
%Outdir=[Outdir '/'];
% Nombvec= [Indir strrep(NomFichierNii,'_d.nii', '.bvec')];
Nombvec= TrouverNomFichier(Indir,'bvec');
if isempty(Nombvec)
    warning('Le fichier bvec est introuvable, vérification requise');
    Nombvec = input('Renseigner le nom du fichier bvec dans le dossier DTI_orig','s');
end
Nombvec=fullfile(Indir,Nombvec);

Nombvecold= TrouverNomFichier(Indir,'oldvec');
if ~isempty(Nombvecold)
    Nombvecold=fullfile(Indir, Nombvecold);
else
    Nombvecold=Nombvec;
end
% Nombval= [Indir strrep(NomFichierNii, '_d.nii', '.bval')];
Nombval= TrouverNomFichier(Indir,'bval');
if isempty(Nombval)
    warning('Le fichier bval est introuvable, vérification requise');
    Nombval = input('Renseigner le nom du fichier bval dans le dossier DTI_orig','s');
end
Nombval=fullfile(Indir,Nombval);

bvecf = importdata(Nombvec,' ');
bvecfold= importdata(Nombvecold,' ');
bvalf=importdata(Nombval,' ');
disp(bvecf);
NomFichierNiiX = ['X' NomFichierNii];

copyfile(Nombval,fullfile(Outdir,strrep(NomFichierNiiX,'.nii', '.bval')));
copyfile(Nombvec,fullfile(Outdir,strrep(NomFichierNiiX,'.nii', '.bvec')));

%Modification des fichier bvec et bval en enlevant le valeur
%correspondant au volume supplémentaire philips

Vol = 1:size(bvecf,2);
idx = find(bvecf(1,:)==0 & bvecf(2,:)==0 & bvecf(3,:)==0 & bvalf(1,:)~=0);

if ~isempty(idx)
    Nombval=fullfile(Outdir,strrep(NomFichierNiiX,'.nii', '.bval'));
    Nombvec=fullfile(Outdir,strrep(NomFichierNiiX,'.nii', '.bvec'));
    
    bvecf(:,idx(1))=[]; % élimine le volume supplémentaire
    Vol(idx(1))=[];
    
    bvalf(idx(1))=[];
    
    disp('Extraction du volume supplémentaire...');
    %     system([FSLcommand 'fslroi ' Indir NomFichierNii ' ' Outdir NomFichierNiiX ' 0 ' num2str(idx(2)-1)]); %Exclue le volume supplémentaire des Nii
    
    
    for j= 1 : size(Vol,2)
        matlabbatch{1}.spm.util.cat.vols{j,1} = [fullfile(Indir,NomFichierNii)  ',' num2str(Vol(j)) ];
    end
    
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.name = fullfile(Outdir,NomFichierNiiX);
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    
    
    disp('Extraction du volume supplémentaire OK');
    fmt= '%d'; % Définit le format pour l'écriture dans le fichier bvec et bval
    for j= 1 :(size(bvecf,2)-1)
        fmt=[fmt ' %6.4E'];
    end
    fmt=[fmt '\n'];
    
    fid=fopen(Nombvec,'wt');
    fprintf(fid,fmt,bvecf');
    fclose(fid);
    
    fid2=fopen(Nombval,'wt');
    fprintf(fid2,'%d ',bvalf');
    fclose(fid2);
    clear fid fid2
    
    disp('Modification des bvec et bval OK');
else
    copyfile(fullfile(Indir, NomFichierNii),fullfile(Outdir,NomFichierNiiX));
end

% %On modifie XinFile pour avoir les mêmes dimensions que b0_P
% h1=spm_vol([Outdir NomFichierNiiX]);
% Xin=spm_read_vols(h1);
% 
% if (idx1 ~= 0) || (idx2 ~= 0) 
%     Xin=Xin(:,:,idx2+1:(size(Xin,3)-idx1),:);
%     [h1.dim]=deal([size(Xin,1) size(Xin,2) size(Xin,3)]);
%     spm_write_vol(h1, Xin); %Error! Can only handle a maximum of 3 dimensions.
% end


%----------------------------------------------
% plot des direction sur sphère 3D
if exist(fullfile(Outdir,'QC.fig'),'file')
    dirsphere=open(fullfile(Outdir,'QC.fig'));
else
    disp('Attention, abscence de figure de contrôle qualité dans le dossier DTI_Pretrt, création d''une nouvelle figure');
    dirsphere=figure;
end
subplot(3,2,2)
for i = 1 :size(bvecfold,2)
    X= bvecfold(1,i) ;
    Y= bvecfold(2,i) ;
    Z= bvecfold(3,i) ;
    
    plot3(X,Y,Z,'r-*')
    hold on
    grid on
    axis equal
end

for i = 1 :size(bvecf,2)
    X= bvecf(1,i) ;
    Y= bvecf(2,i) ;
    Z= bvecf(3,i) ;
    
    plot3(X,Y,Z,'g-*')
    hold on
    grid on
    axis equal
    
end

[x1,y1,z1] = sphere();
r = 0.99;
surf( r*x1, r*y1, r*z1 ) % sphere with radius 5 centred at (0,0,0)
colormap('gray')
title('Direction gradients de diffusion');
saveas(dirsphere,fullfile(Outdir,'QC.fig'));

%--------------------------------------------------
fid=fopen(fullfile(Outdir,'index.txt'),'w');
ind=ones(1,size(bvecf,2));
fprintf(fid,'%d ',ind);
fclose(fid);
clear ind fid
if exist(fullfile(Outdir,'acqparamsDTI.txt'),'file')>0 && exist(fullfile(Outdir,'acqparamsB0.txt'),'file')>0
    acqDTI = importdata(fullfile(Outdir,'acqparamsDTI.txt'),' ');
    acqB0 =importdata(fullfile(Outdir,'acqparamsB0.txt'),' ');
    acq = [acqDTI ; acqB0];
    fid = fopen(fullfile(Outdir,'acqparams.txt'),'w');
    fprintf(fid,'%d %d %d %7.5f\n',acq');
    fclose(fid);
    clear fid;
else
    warning('les fichiers acqparamsDTI et B0 n''ont pas pu être concaténés, Vérification requise');
    
end
cd (Outdir) % Change le répertoire courant pour

% NomDossierB0map = TrouverNomDossier(PathSujet,'B0map');
% if isempty(NomDossierB0map) && ~isempty(NomDossierB0);
if ~isempty(B0dir);
    disp('Estimation de la carte de champs ...')
    v1 = spm_vol('b0_P.nii');
    if mod(v1.dim(3),2)==0
        status=system([FSLcommand 'topup --imain=b0_PA.nii --datain=acqparams.txt --config=/home/veronica/Donnees/DTIPark/b02b0.cnf --out=topup_results']);
    else
        status=system([FSLcommand 'topup --imain=b0_PA.nii --datain=acqparams.txt --config=/home/veronica/Donnees/DTIPark/b02b02.cnf --out=topup_results']);% b02b02.cnf sans subsampling dans le répertoire contenant tous les centres
    end
    if status==0
        disp('Estimation de la carte de champs  OK')
    else
        warning('L''estimation de la carte de champs par topup a échoué, Vérification requise (dbcont pour continuer)');
        dbstop in PreprocessingDTI_FSL.m
    end
    disp('Correction des distorsions sur les b0 des images de diffusions ...')
    status=system([FSLcommand 'applytopup --imain=b0_P.nii,rb0_A.nii --topup=topup_results --datain=acqparams.txt --inindex=1,2 --out=b0_corr.nii']); % Correction de la b0 pour servir de cible de coregistration pour les images anatomiques
    if status ==0
        disp('Correction des distorsions sur les b0 des images de diffusions  OK')
    else
        warning('La correction des distorsions sur la B0 a échoué, Vérification requise (dbcont pour continuer)');   
    end
else
    disp('Pas de Bo valide, La correction n''inclura pas la correction d''artefact de suceptibilité')
    
end
