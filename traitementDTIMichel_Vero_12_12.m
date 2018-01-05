clear all;
close all;

%% Path et variables d'environnement
clc
dicomdict('set','/media/veronica/DATAPART1/Code/Code_DTI/dictionnaire_dicom.txt'); %charge le dictionnaire contenant également les tags privés

Path = '/media/veronica/DATAPART1/Documents/MATLAB/';
Path2 = '/media/veronica/DATAPART1/Donnees/DTIPark/Park/';
cd(Path2);
addpath([Path 'spm12/']);
addpath([Path 'DWIDenoisingPackage/']);

FSLcommand='/usr/lib/fsl/5.0/';
setenv('FSLOUTPUTTYPE','NIFTI');
s=1; % Compteur de nouveaux sujets

Doss= dir([ Path2 '*']);
Doss = Doss(arrayfun(@(x) ~strcmp(x.name(1),'.'),Doss));% pour supprimer les . et .. du résultat du dir

%% Traitement de tous les sujets dans le répertoire renseigné
for i = 1 : size(Doss,1)
    if Doss(i,1).isdir==1
        PathSujetS=[Path2 Doss(i,1).name '/'];
        % -----------CAS 1a: images DICOM----------------
        if exist([PathSujetS 'DICOM'],'dir')
            a=dir([PathSujetS 'DICOM']);
            a = a(arrayfun(@(x) ~strcmp(x.name(1),'.'),a));
            s=1;
            % Ne pas traiter le pseudo-sujet DTI_64dir
            while  s< size(a,1)+1
                pname=a(s,1).name;
                hdr2=dicominfo(fullfile(PathSujetS,'DICOM',pname)); %extrait les informations du fichier DICOM
                if isfield(hdr2,'ProtocolName') && ~isempty(strfind(hdr2.ProtocolName,'DTI_64dir'))
                    s=size(a,1)+1;
                else
                    s= s+1;
                end
                hdr2.ProtocolName
            end
            
            % Extraction de paramètres wfsHz, EchoSpacing, TRtime
            nomacqparam='acqparamsDTI.txt';
            if isfield(hdr2,'MagneticFieldStrength') && isfield(hdr2,'WaterFatShift') && isfield(hdr2,'EchoTrainLength') && isfield(hdr2,'EPIFactor')
                wfsHz=3.4*42.57* double(hdr2.MagneticFieldStrength);
                EchoSpacing = ((1000*double(hdr2.WaterFatShift))/(wfsHz*(double(hdr2.EchoTrainLength)+1)));%/2.2;% EchoSpacing en ms 
                if size(EchoSpacing,2)>1
                    disp('Possible passage par le PACS, verification par python  et entrer valeur EchoSpacing ou alerter le centre pour réimportation des images');
                end
                save([PathSujetS 'EchoSpacing.mat'],'EchoSpacing');
                
                TRtime = (EchoSpacing * double(hdr2.EPIFactor))/1000; %Total readout time en secondes
                if TRtime<0.01
                    disp('TRtime anormalement faible, vérification requise, dbcont pour continuer')
                end
                
                fid=fopen(fullfile(PathSujetS,nomacqparam),'w');
                fprintf(fid,'%d %d %d %7.5f',0,1,0,TRtime); % Attention valable pour une direction d'encodage de phase P vers A
                fclose(fid);
            else
                disp('Des champs importants sont manquants, vérifier manuellement, dbcont pour continuer');
            end
            
            % Conversion de fichier DICOM en nifti
            status = ConversionNifti([PathSujetS 'DICOM/'],PathSujetS,'DTI');
            disp(status);
            
            
            Nomdti= TrouverNomFichier2(PathSujetS,'DTI_64dir','.nii');
            %inFile=MainDWIDenoising2(Path,Nomdti,indir);
            inFile= TrouverNomFichier([PathSujetS '/' NomDossier],'.nii');
            indir=PathSujetS;
            outdir =[PathSujetS 'DTI/'];
            mkdir(outdir)
            if exist([ PathSujetS 'acqparamsDTI.txt'],'file')
                copyfile([ PathSujetS 'acqparamsDTI.txt'],[outdir 'acqparamsDTI.txt']);
            end
            
        % -------------CAS 1b: Images nifti ou autre-----------------    
        else
            NomDossier = TrouverNomDossier(PathSujetS,'DTI_64dir');
            if isempty(NomDossier)
                nomseq=input('Nom spécifique de la séquence  pour la DTI: ','s');
                NomDossier = TrouverNomDossier(PathSujetS,nomseq);
            end
            Nomdti = TrouverNomFichier([PathSujetS '/' NomDossier],'.nii');
            %Nomdti = MainDWIDenoising2(Path,Nomdti,[PathSujetS '/' NomDossier]);
            %inFile = TrouverNomFichier([PathSujetS NomDossier],'_d.nii');
            inFile = TrouverNomFichier([PathSujetS NomDossier],'.nii');
            indir = [PathSujetS NomDossier '/'];
            outdir =[PathSujetS 'DTI/'];
            mkdir(outdir)
            if exist([ Path2 'acqparamsDTI.txt'],'file')
                copyfile([ Path2 'acqparamsDTI.txt'],[outdir 'acqparamsDTI.txt']);
            end
        end
        
        %
        if exist([indir strrep(inFile,'_d.nii', '.bval')],'file')==0  %find and replace string 
            movefile([indir strrep(inFile,'_d.nii', '_bvals')],[indir strrep(inFile,'_d.nii', '.bval')]);
            movefile([indir strrep(inFile,'_d.nii', '_bvecs')],[indir strrep(inFile,'_d.nii', '.bvec')]);
        end
        
        nombval=TrouverNomFichier(indir,'bval');
        bvaltest= importdata([indir nombval]);
        posb0= find(bvaltest==0);
        % fslroi: extract region of interest from an image 
        system([FSLcommand 'fslroi ' strrep(indir,' ','\ ') inFile ' ' strrep(indir,' ','\ ') 'b0_P.nii '  num2str(posb0(1)-1) ' 1']);
        copyfile([ indir 'b0_P.nii'],[outdir 'b0_P.nii']);
        clear bvaltest
       
        % ----------------Cas 2a: Images DICOM-------------------
        if exist([PathSujetS 'DICOM'],'dir')
            B0file = TrouverNomFichier2(PathSujetS,'PrepAPA','.nii');
            nombval=TrouverNomFichier(PathSujetS ,'bval');
            if ~isempty(nombval)
                bvaltest= importdata([PathSujetS nombval]);
            else
                bvaltest= 0;
            end
            posb0A= find(bvaltest== 0);
            
            if ~isempty(posb0A) % Si le B0 n'est pas valide car b=/= 0
                system([FSLcommand 'fslroi ' strrep(PathSujetS,' ','\ ') '/' B0file ' ' strrep(PathSujetS,' ','\ ') '/b0_A.nii ' num2str(posb0A(1)-1) ' 1']);
                clear bvaltest
                
                nomacqparam='acqparamsB0.txt';
                pname=a(17,1).name;
                hdr2=dicominfo(fullfile(PathSujetS,'DICOM',pname)); %extrait les informations du fichier DICOM
                if isfield(hdr2,'MagneticFieldStrength') && isfield(hdr2,'WaterFatShift') && isfield(hdr2,'EchoTrainLength') && isfield(hdr2,'EPIFactor')
                    wfsHz=3.4*42.57* double(hdr2.MagneticFieldStrength);
                    
                    EchoSpacing = ((1000*double(hdr2.WaterFatShift))/(wfsHz*(double(hdr2.EchoTrainLength)+1)));%/2.2;% EchoSpacing en ms
                    
                    if size(EchoSpacing,2)>1
                        disp('Possible passage par le PACS, verification par python  et entrer valeur EchoSpacing ou alerter le centre pour réimportation des images');
                        dbstop in GUIImport_CQData.m
                    end
                    
                    TRtime = (EchoSpacing * double(hdr2.EPIFactor))/1000; %Total readout time en secondes
                    save([PathSujetS 'EchoSpacing.mat'],'EchoSpacing');
                    if TRtime<0.01
                        disp('TRtime anormalement faible, vérification requise, dbcont pour continuer')
                        dbstop in GUIImport_CQData.m
                    end
                    
                    fid=fopen(fullfile(PathSujetS,nomacqparam),'w');
                    fprintf(fid,'%d %d %d %7.5f',0,-1,0,TRtime); % Attention valable pour une direction d'encodage de phase P vers A
                    fclose(fid);
                else
                    disp('Des champs importants sont manquants, vérifier manuellement, dbcont pour continuer');
                    dbstop in GUIImport_CQData.m
                end
                copyfile([ PathSujetS 'acqparamsB0.txt'],[outdir 'acqparamsB0.txt']);
                copyfile([ PathSujetS 'b0_A.nii'],[outdir 'b0_A.nii']);
                
                clear matlabbatch
                spm_jobman('initcfg');
                
                % Coregistration du b0_A sur le b0_P pour assurer une coregistration parfaite pour l'estimation de la carte de champs par topup
                matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[outdir 'b0_P.nii,1']};
                matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[outdir 'b0_A.nii,1']};
                matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                
                spm('defaults', 'FMRI');
                spm_jobman('run', matlabbatch);
                clear matlabbatch
                
                spm_check_registration([outdir 'b0_P.nii'],[outdir 'rb0_A.nii'])
                disp('Vérification de la coregistration, pressez une touche pour continuer');
                pause
                
                
                %-------------------------------------------
                system([FSLcommand 'fslmerge -t ' strrep(outdir,' ','\ ') 'b0_PA.nii ' strrep(outdir,' ','\ ') 'b0_P.nii ' strrep(outdir,' ','\ ') 'rb0_A.nii']); %merge les 2 b0 pour topup
                NomDossierB0=1;
            else
                NomDossierB0=[];
            end
        % -------------------Cas 2b: Images nifti ou autres----------------
        else    
            NomDossierB0 = TrouverNomDossier(PathSujetS,'PrepAPA');
            if isempty(NomDossierB0)
                nomseq=input('Nom spécifique de la séquence pour le B0 : ','s');
                NomDossierB0 = TrouverNomDossier(PathSujetS,nomseq);
            end
            if ~isempty(NomDossierB0)
                B0file = TrouverNomFichier([PathSujetS NomDossierB0],'_d.nii');
                
                if isempty(B0file)
                    B0file = TrouverNomFichier([PathSujetS NomDossierB0],'.nii');
                    %   bvaltest= importdata([PathSujetS '/' NomDossierB0 '/' strrep(B0file,'.nii', '_bvals')]);
                    % else
                    %   bvaltest= importdata([PathSujetS '/' NomDossierB0 '/' strrep(B0file,'_d.nii', '_bvals')]);
                end
                nombval=TrouverNomFichier([PathSujetS NomDossierB0],'bval');
                if ~isempty(nombval)
                    bvaltest= importdata([PathSujetS NomDossierB0 '/' nombval]);
                else
                    bvaltest= 0;
                end
                posb0A= find(bvaltest== 0);
                
                if ~isempty(posb0A) % Si le B0 n'est pas valide car b=/= 0
                    system([FSLcommand 'fslroi ' strrep(PathSujetS,' ','\ ') '/' NomDossierB0 '/' B0file ' ' strrep(PathSujetS,' ','\ ') '/' NomDossierB0 '/b0_A.nii ' num2str(posb0A(1)-1) ' 1']);
                    clear bvaltest
                    
                    copyfile([ Path2 'acqparamsB0.txt'],[outdir 'acqparamsB0.txt']);
                    copyfile([ PathSujetS NomDossierB0 '/b0_A.nii'],[outdir 'b0_A.nii']);
                    clear matlabbatch
                    spm_jobman('initcfg');
                    
                    
                    % Coregistration du b0_A sur le b0_P pour assurer une coregistration parfaite pour l'estimation de la carte de champs par topup
                    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[outdir 'b0_P.nii,1']};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[outdir 'b0_A.nii,1']};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                    
                    spm('defaults', 'FMRI');
                    spm_jobman('run', matlabbatch);
                    clear matlabbatch
                    
                    spm_check_registration([outdir 'b0_P.nii'],[outdir 'rb0_A.nii'])
                    disp('Vérification de la coregistration, pressez une touche pour continuer');
                    %pause
                    
                    
                    %-------------------------------------------
                    system([FSLcommand 'fslmerge -t ' strrep(outdir,' ','\ ') 'b0_PA.nii ' strrep(outdir,' ','\ ') 'b0_P.nii ' strrep(outdir,' ','\ ') 'rb0_A.nii']); %merge les 2 b0 pour topup
                else
                    NomDossierB0=[];
                end
            end
        end
        XinFile = PreprocessingDTI_FSL(indir,outdir,inFile,NomDossierB0); % Prépare les données DTI
        
     
        %% Rééchantillonnage de la b0_P ou b0_corr.nii en 1x1x1
        C3Dcommand='/home/veronica/Downloads/Programs/c3d/bin/';
        if exist([PathSujetS 'DTI/b0_corr.nii'],'file')
            system([C3Dcommand 'c3d ' strrep(PathSujetS,' ','\ ') 'DTI/b0_corr.nii -resample-mm 1.0x1.0x1.0mm -o ' strrep(PathSujetS,' ','\ ') 'DTI/rb0_corr.nii']);
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ PathSujetS 'DTI/rb0_corr.nii,1']};
        else
            system([C3Dcommand 'c3d ' strrep(PathSujetS,' ','\ ') 'DTI/b0_P.nii -resample-mm 1.0x1.0x1.0mm -o ' strrep(PathSujetS,' ','\ ') 'DTI/rb0_P.nii']);
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ PathSujetS 'DTI/rb0_P.nii,1']};
        end
        if exist([PathSujetS 'DICOM'],'dir')
            NomAnat = TrouverNomFichier2(PathSujetS,'Anat','.nii');
            PathAnat = [PathSujetS NomAnat];
        else
            NomDossierAnat = TrouverNomDossier(PathSujetS,'Anat');
            NomAnat = TrouverNomFichier([PathSujetS NomDossierAnat],'.nii');
            PathAnat =[PathSujetS NomDossierAnat '/' NomAnat];
            %Coregistration entre FLAIR T1 et T2* sur la DTI corrigée
        end
        disp('Corgistration des FLAIR et T1 sur la DTI...');
        spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[PathAnat ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
        
        cd(outdir)
        %         system([FSLcommand 'bet b0_corr.nii mask2.nii -f 0.1 -g 0 -m']);
        ROBEXcommand='/home/veronica/Downloads/Programs/ROBEXv12/ROBEX/';
        status=system([ROBEXcommand 'ROBEX -i b0_corr.nii -r ' ROBEXcommand 'ref_vols -m atlas_mask.nii']);
        %status=system('/usr/local/robex-ubuntu/ROBEX -i b0_corr.nii -r /usr/local/robex-ubuntu/ -m mask1.nii');
        
        %     system([FSLcommand 'eddy_correct ' strrep(indir,' ','\ ') Nomdti ' ' strrep(out_ec,' ','\ ') ' 1 - spline']);
        FSLcommand2= '/usr/share/fsl/5.0/bin/';
        sinFile= [outdir XinFile];
        
        %         system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' strrep(outdir,' ','\ ') 'mask2_mask.nii --index=' strrep(outdir,' ','\ ') 'index.txt --acqp=' strrep(outdir,' ','\ ') 'acqparams.txt --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0  --slm=linear --topup=' strrep(outdir,' ','\ ') 'topup_results --flm=quadratic --out=' strrep(outdir,' ','\ ') 'DTI_esc.nii']);
        status=system([FSLcommand2 'eddy_openmp --imain=' sinFile ' --mask=' fullfile(strrep(outdir,' ','\ '), 'mask1.nii') ' --index=' fullfile(strrep(outdir,' ','\ '),'index.txt') ' --acqp=' fullfile(strrep(outdir,' ','\ '),'acqparams.txt') ' --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0 --resamp=jac --interp=spline --dont_peas --fep --repol --flm=quadratic --slm=linear --topup=' fullfile(strrep(outdir,' ','\ '),'topup_results') ' --out=' fullfile(strrep(outdir,' ','\ '),'DTI_esc.nii')]);
        
        
        %     system([FSLcommand 'eddy_correct ' Nomdti ' ec_data.nii 1 - spline']);
        %         Anat=TrouverNomFichier(PathDoss,'anatHR.nii');
        %
        Nombvec=TrouverNomFichier(outdir,'_bvecs');
        Nombval=TrouverNomFichier(outdir,'bval');
        %     system([FSLcommand 'bet ' Anat ' brain.nii -f 0.1 -g 0 -m']);
        %         system([FSLcommand 'bet b0.nii braindti.nii -f 0.1 -g 0 -m']);
        
        system([FSLcommand 'dtifit --data=DTI_esc.nii --out=dti --mask=mask1.nii --bvecs=' Nombvec ' --bvals=' Nombval ' --save_tensor'])
        
        nomFichierMD = TrouverNomFichier(outdir,'MD');
        nomFichierFA = TrouverNomFichier(outdir,'FA');
        ValArtefactDTI(outdir,nomFichierMD,nomFichierFA)
        
        %         %Coregistration sur T1
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[outdir Anat ',1']};
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[outdir 'ec_data.nii,1']};
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
        %             [outdir 'dti_MD.nii']
        %             [outdir 'dti_FA.nii']
        %             [outdir 'dti_L1.nii']
        %             };
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        %         spm('defaults', 'FMRI');
        %         spm_jobman('run', matlabbatch);
        %
        %         clear matlabbatch
        %
        cd(Path2)
%alias matlab='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/matlab -desktop'
        
     end
end