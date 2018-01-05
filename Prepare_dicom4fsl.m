function Prepare_dicom4fsl(PathSujetS)
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
    wfsHz=3.4*42.57*double(hdr2.MagneticFieldStrength);
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

indir=[PathSujetS strep(NomD, '.nii', '')]
% 
%             indir=PathSujetS;
%             Nomdti= TrouverNomFichier2(PathSujetS,'DTI_64dir','.nii');
%             inFile=MainDWIDenoising2(Path,Nomdti,indir);
%             outdir =[PathSujetS '/DTI/'];
%             mkdir(outdir)
%             if exist([ PathSujetS 'acqparamsDTI.txt'],'file')
%                 copyfile([ PathSujetS 'acqparamsDTI.txt'],[outdir 'acqparamsDTI.txt']);
%             end
%         if exist([indir strrep(inFile,'_d.nii', '.bval')],'file')==0
%             movefile([indir strrep(inFile,'_d.nii', '_bvals')],[indir strrep(inFile,'_d.nii', '.bval')]);
%             movefile([indir strrep(inFile,'_d.nii', '_bvecs')],[indir strrep(inFile,'_d.nii', '.bvec')]);
%         end
%         
%         nombval=TrouverNomFichier(indir,'bval');
%         bvaltest= importdata([indir nombval]);
%         posb0= find(bvaltest== 0);
%         system([FSLcommand 'fslroi ' strrep(indir,' ','\ ') inFile ' ' strrep(indir,' ','\ ') 'b0_P.nii '  num2str(posb0(1)-1) ' 1']);
%         copyfile([ indir 'b0_P.nii'],[outdir 'b0_P.nii']);
%         clear bvaltest
    
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

    % Merge de b0_P et rb0_A 
    system([FSLcommand 'fslmerge -t ' strrep(outdir,' ','\ ') 'b0_PA.nii ' strrep(outdir,' ','\ ') 'b0_P.nii ' strrep(outdir,' ','\ ') 'rb0_A.nii']); %merge les 2 b0 pour topup
    NomDossierB0=1;
else
    NomDossierB0=[];
end
disp('DICOM2');