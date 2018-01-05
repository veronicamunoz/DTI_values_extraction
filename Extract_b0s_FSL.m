function Extract_b0s_FSL(PathSujetS, Path2, FSLcommand)

% ============= EXTRACTION b0_P =============
NomDossier = TrouverNomDossier(PathSujetS,'DTI_64dir');
if isempty(NomDossier)
    nomseq=input('Nom spécifique de la séquence  pour la DTI: ','s');
    NomDossier = TrouverNomDossier(PathSujetS,nomseq);
end
%Nomdti = TrouverNomFichier([PathSujetS '/' NomDossier],'.nii');
%Nomdti = MainDWIDenoising2(Path,Nomdti,[PathSujetS '/' NomDossier]);
indir = [PathSujetS NomDossier '/'];
inFile = TrouverNomFichier(indir,'.nii');
outdir =[PathSujetS 'DTI/'];
mkdir(outdir);

if exist([ Path2 'acqparamsDTI.txt'],'file')
    copyfile([ Path2 'acqparamsDTI.txt'],[outdir 'acqparamsDTI.txt']);
end

if exist([indir strrep(inFile,'.nii', '.bval')],'file')==0  %find and replace string 
    movefile([indir strrep(inFile,'.nii', '_bvals')],[indir strrep(inFile,'.nii', '.bval')]);
    movefile([indir strrep(inFile,'.nii', '_bvecs')],[indir strrep(inFile,'.nii', '.bvec')]);
end

nombval=TrouverNomFichier(indir,'bval');
bvaltest= importdata([indir nombval]);
posb0= find(bvaltest==0);

% fslroi: extract region of interest from DTI_64dir file to b0_P file
system([FSLcommand 'fslroi ' strrep(indir,' ','\ ') inFile ' ' strrep(indir,' ','\ ') 'b0_P.nii '  num2str(posb0(1)-1) ' 1']);
copyfile([ indir 'b0_P.nii'],[outdir 'b0_P.nii']);
clear bvaltest

% ================ EXTRACTION DE b0_A ====================
% ================= CALCUL DE b0_PA ======================
NomDossierB0 = TrouverNomDossier(PathSujetS,'PrepAPA');
if isempty(NomDossierB0)
    nomseq=input('Nom spécifique de la séquence pour le B0 : ','s');
    NomDossierB0 = TrouverNomDossier(PathSujetS,nomseq);
end

if ~isempty(NomDossierB0)
    B0file = TrouverNomFichier([PathSujetS NomDossierB0],'.nii');
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
        
        % Coregistration du b0_A sur le b0_P
        % Pour assurer une coregistration parfaite pour l'estimation de la carte de champs par topup
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

        %Recoupage de b0_P par rapport à rb0_A
        h1=spm_vol([outdir 'b0_P.nii']);
        h2=spm_vol([outdir 'rb0_A.nii']);
        b0_P=spm_read_vols(h1);
        rb0_A=spm_read_vols(h2);

        idx1= input('Coupes supérieures à rogner : ');
        
        if idx1 ~= 0 
        b0_P=b0_P(:,:,1:(size(b0_P,3)-idx1));
        rb0_A=rb0_A(:,:,1:(size(rb0_A,3)-idx1));

            if size(b0_P)==size(rb0_A)
                h1.dim=size(b0_P);
                h2.dim=size(rb0_A);
                spm_write_vol(h1, b0_P);
                spm_write_vol(h2, rb0_A);
            else
                disp('ERREUR de dimensions entre b0_P et rb0_A');
            end
        end
        
        idx2= input('Coupes inférieures à rogner : ');
        if idx2 ~= 0
        b0_P=b0_P(:,:,idx2:size(b0_P,3));
        rb0_A=rb0_A(:,:,idx2:size(rb0_A,3));
        b0_P(isnan(rb0_A))=0;
        rb0_A(isnan(rb0_A))=0; 

            if size(b0_P)==size(rb0_A)
                h1.dim=size(b0_P);
                h2.dim=size(rb0_A);
                spm_write_vol(h1, b0_P);
                spm_write_vol(h2, rb0_A);
            else
                disp('ERREUR de dimensions entre b0_P et rb0_A');
            end
        end
        
        %system([FSLcommand 'flirt -in ' [outdir 'b0_A.nii'] ' -ref ' [outdir 'b0_P.nii'] '-interp spline -out ' [outdir 'rb0_A.nii']]);
        clear matlabbatch
        spm_check_registration([outdir 'b0_P.nii'],[outdir 'rb0_A.nii'])
        disp('Vérification de la coregistration, pressez une touche pour continuer');
        pause    
        
        %fslmerge:merge b0_P and rb0_A into b0_PA
        system([FSLcommand 'fslmerge -t ' strrep(outdir,' ','\ ') 'b0_PA.nii ' strrep(outdir,' ','\ ') 'b0_P.nii ' strrep(outdir,' ','\ ') 'rb0_A.nii']); %merge les 2 b0 pour topup
    else
        disp('posb0A vide');
    end
end

