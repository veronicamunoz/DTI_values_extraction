function Prepare_dicom4fsl(PathSujetS)
a=dir([PathSujetS 'DICOM']);
a = a(arrayfun(@(x) ~strcmp(x.name(1),'.'),a));
s=1;

%Lis soit un fichier DTI_64dir soit plusieurs fichiers avec quelques
%directions chacun
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
status = system( ['mcverter -o ' PathSujetS ' -f fsl -d -n -F ''+PatientName,-PatientId,-SeriesDate,-SeriesTime,-StudyId,-StudyDescription,-SeriesNumber,+SequenceName,-SeriesDescription,+ProtocolName'' -v ' PathSujetS 'DICOM/']);
% status = system(['mrconvert ' PathSujetS 'DICOM/ ext.nii']);
disp(status);

cd(PathSujetS);
a=dir(PathSujetS);
a = a(arrayfun(@(x) ~strcmp(x.name(1),'.'),a));
for i = 1:size(a,1)
    if ~isdir(a(i).name) && isempty(strfind(a(i).name,'DICOM'))
        tempname=a(i).name(1:24);
        if exist(tempname,'dir')==0
            mkdir(tempname);
        end
        movefile(fullfile(PathSujetS, a(i).name), fullfile(PathSujetS,tempname,a(i).name));
    end
end

disp('DICOM2');