clear all;
close all;
dicomdict('set','/media/veronica/DATAPART1/Code/Code_DTI/dictionnaire_dicom.txt'); %charge le dictionnaire contenant également les tags privés

Path = '/media/veronica/DATAPART1/Documents/MATLAB/';
Path2 = '/media/veronica/DATAPART1/Donnees/DTIPark/Park/';
cd(Path2);
addpath([Path 'spm12/']);
addpath([Path 'DWIDenoisingPackage/']);

FSLcommand='/usr/lib/fsl/5.0/';
setenv('FSLOUTPUTTYPE','NIFTI');
C3Dcommand='/home/veronica/Downloads/Programs/c3d/bin/';
s=1; % Compteur de nouveaux sujets

Doss= dir([ Path2 '*']);
Doss = Doss(arrayfun(@(x) ~strcmp(x.name(1),'.'),Doss));% pour supprimer les . et .. du résultat du dir

for i = 1 : size(Doss,1)
    if Doss(i,1).isdir==1
       PathSujetS=[Path2 Doss(i,1).name '/'];
       
       if exist([PathSujetS 'DICOM'],'dir')
           disp('DICOM');
           %disp('Converting DICOM to nifti');
           %disp('Organizing files');
           %Prepare_dicom4fsl(PathSujetS);
       end
       
       NomDossier = TrouverNomDossier(PathSujetS,'DTI_64dir');
           
       if ~isempty(NomDossier)
           disp('=============');
           disp(Doss(i,1).name);
           disp('=============');
           indir= [PathSujetS NomDossier '/'];
           outdir= [PathSujetS 'DTI/'];
           inFile= TrouverNomFichier(indir, '.nii');
           NomDossierB0 = TrouverNomDossier(PathSujetS,'PrepAPA');
           Extract_b0s_FSL(PathSujetS, Path2, FSLcommand);
           
           if exist([outdir 'b0_PA.nii'], 'file')
               disp('Extract_b0s_FSL a bien fini');
               XinFile = PreprocessingDTI_FSL(indir,outdir,inFile,NomDossierB0); % Prépare les données DTI
               %XinFile = ['X' inFile];
               cd(outdir);
               if exist([outdir 'b0_corr.nii'],'file')
                   disp('PreprocessingDTI_FSL a bien fini');
                   system([C3Dcommand 'c3d ' strrep(PathSujetS,' ','\ ') 'DTI/b0_corr.nii -resample-mm 1.0x1.0x1.0mm -o ' strrep(PathSujetS,' ','\ ') 'DTI/rb0_corr.nii']);
               else
                   disp('Le fichier b0_corr.nii na pas été géneré par PreprocessingDTI_FSL');
                   disp('Le script continue avec b0_P');
                   system([C3Dcommand 'c3d ' strrep(PathSujetS,' ','\ ') 'DTI/b0_P.nii -resample-mm 1.0x1.0x1.0mm -o ' strrep(PathSujetS,' ','\ ') 'DTI/rb0_P.nii']);
               end
               if exist([PathSujetS 'DTI/rb0_corr.nii'],'file')
                   system([FSLcommand 'bet rb0_corr.nii b0_ss.nii -f 0.1 -m -Z']);
               elseif exist([PathSujetS 'DTI/rb0_P.nii'],'file') 
                   system([FSLcommand 'bet rb0_P.nii b0_ss.nii -f 0.1 -m -Z']); 
               else
                   disp('Aucun résultat de c3d disponible');
               end
               
               %     system([FSLcommand 'eddy_correct ' strrep(indir,' ','\ ') Nomdti ' ' strrep(out_ec,' ','\ ') ' 1 - spline']);
               FSLcommand2= '/usr/share/fsl/5.0/bin/';
               sinFile= [outdir XinFile];

               %         system([FSLcommand2 'eddy --imain=' sinFile ' --mask=' strrep(outdir,' ','\ ') 'mask2_mask.nii --index=' strrep(outdir,' ','\ ') 'index.txt --acqp=' strrep(outdir,' ','\ ') 'acqparams.txt --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0  --slm=linear --topup=' strrep(outdir,' ','\ ') 'topup_results --flm=quadratic --out=' strrep(outdir,' ','\ ') 'DTI_esc.nii']);
               status=system([FSLcommand2 'eddy_openmp --imain=' sinFile ' --mask=' [outdir 'b0_ss_mask.nii'] ' --index=' [outdir 'index.txt'] ' --acqp=' [outdir 'acqparams.txt'] ' --bvecs=' strrep(sinFile,'.nii', '.bvec') ' --bvals=' strrep(sinFile,'.nii', '.bval')  ' --fwhm=0 --resamp=jac --interp=spline --dont_peas --fep --repol --flm=quadratic --slm=linear --topup=' [outdir 'topup_results'] ' --out=' [outdir 'DTI_esc.nii']]);
               Nombvec=TrouverNomFichier(outdir,'bvec');
               Nombval=TrouverNomFichier(outdir,'bval');
               system([FSLcommand 'dtifit --data=DTI_esc.nii --out=dti --mask=b0_ss_mask.nii --bvecs=' Nombvec ' --bvals=' Nombval ' --save_tensor'])
               nomFichierMD = TrouverNomFichier(outdir,'MD');
               nomFichierFA = TrouverNomFichier(outdir,'FA');
               if ~isempty(nomFichierMD)
                   ValArtefactDTI(outdir,nomFichierMD,nomFichierFA)
               end
           end
           cd(Path2)
       end
    end
end