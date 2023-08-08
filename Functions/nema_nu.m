function [NU, CV, SUV, slice_flags, positions] = nema_nu(image, dicom_headers)
% Function [NU, CV, SUV, slice_flags, positions] = nema_nu(image, dicom_headers)
% Processes a 3D image from a 20cm cylindrical phantom  and estimates SUV
% error and non-uniformity according to NEMA NU-2001 standards. 
% Copyright 2023 Mauro Namías - mnamias@gmail.com
% image: pixels
% dicom_headers: DICOM headers
% NU: array of maximum slice non-uniformity 
% CV: array of slice coefficient of variation
% SUV: array of slice SUV error



[a,b,c] = size(image);
[N,X] = hist(image(:),100);
slice_flags = [];

[peaks,pos] = findpeaks(N,'NPEAKS',2,'sortstr','descend');

mean_phantom = X(pos(1));

%parámetros del fantoma;

t1 = dicom_headers(1).RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
A0 = dicom_headers(1).RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;

prompt = {'Enter phantom initial activity [Bq]:'};
dlg_title = 'Phantom activity';
num_lines = 1;
def = {A0};

prompt = {'Enter phantom initial activity time of measurement:'};
dlg_title = 'Phantom activity time';
num_lines = 1;
def = {t1};


t2 = dicom_headers(1).SeriesTime;
%t2 = dicom_headers(1).AcquisitionTime;

hs = str2num(t2(1:2))-str2num(t1(1:2));
mins = str2num(t2(3:4))-str2num(t1(3:4));
dt = hs*60+mins;


fd = exp(-log(2)*dt/110);

A1 = A0*fd;


prompt = {'Enter phantom volume [ml]:'};
dlg_title = 'Phantom volume';
num_lines = 1;
def = {'5640'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
vol = str2num(cell2mat(answer));

%vol = 5640;

CA = A1/vol; % concentración de actividad esperada en Bq/ml

modo = dicom_headers(1).ReconstructionMethod;
%modo = modo(1:2);
fecha = dicom_headers(1).AcquisitionDate;

xlswrite(filename,{fecha},3,'C7')
progressbar()
for i = 1:c

% mascara = medfilt2(image(:,:,i),[9;9]) > mean_phantom*0.8;
% se = strel('disk',1);
% mascara = imerode(mascara,se);

pix_spacing = dicom_headers(1).PixelSpacing;
blk_size = round(10/pix_spacing(1));

%fun = @medias_blk;
fun = @medias_block;
%medias = blkproc(image(:,:,i),[blk_size blk_size],fun);
medias = blockproc(image(:,:,i),[blk_size blk_size],fun);

%fun = @min_blk;
fun = @min_block;
%minis = blkproc(image(:,:,i),[blk_size blk_size],fun);
minis = blockproc(image(:,:,i),[blk_size blk_size],fun);

mascara_small = medfilt2(minis) > mean_phantom*0.8;
se = strel('disk',1);
mascara_small = imerode(mascara_small,se);

medias_mask = medias(mascara_small > 0);
media_slice(i) = mean(medias_mask(:));

NRois(i) = sum(mascara_small(:)>0);

if(NRois(i) > 100)
    slice_flags(i) = 1;
    NU1 = 100*(max(medias_mask(:)) - media_slice(i))/media_slice(i);
    NU2 = 100*(media_slice(i) - min(medias_mask(:)))/media_slice(i);
    NU(i) = max(NU1,NU2);
    CV(i) = std(medias_mask(:))/media_slice(i)*100;
    SUV(i) = (media_slice(i)-CA)/CA*100;
else
    slice_flags(i) = 0;
NU(i) = 0;
CV(i) = 0;
SUV(i) = 0;
end 
  
progressbar(i/c)
end

positions = {dicom_headers.SliceLocation};
positions = cell2mat(positions);
%positions = fliplr(positions);


%% Busca la porción central del fantoma
rango = abs((positions(end)-positions(1)));

indices = find(slice_flags==1);
indices = indices(4:end-3);

comienzo = indices(1);
fin = indices(end);
centro = (positions(comienzo)+positions(fin))/2;

posicion_inicio = positions(comienzo)
posicion_fin = positions(fin)



if(rango > 170)
disp('2 camillas')
else
disp('1 camilla')
end


%% SUV
figure;
plot( positions(indices), SUV(indices), 'k -');
%ylim([-10 10]);

SUVs_ok = SUV(indices);
media_SUV = mean(SUVs_ok(:));
max_SUV = max(SUVs_ok(:));
std_SUV = std(SUVs_ok(:));
min_SUV = min(SUVs_ok(:));

hold on;
plot(positions(indices), media_SUV,'k --');

posi = positions(indices);

title( ['Error SUV ' modo '  ' fecha(7:8) '/' fecha(5:6) '/' fecha(1:4)  ] )
xlabel('Posición Axial [mm]');
ylabel('Error SUV [%]');

f = '%6.2g';
text(posicion_inicio+10, 9, ['SUV: ' num2str(media_SUV,f) ' +- ' num2str(std_SUV,f) '  [' num2str(min_SUV,f) '-' num2str(max_SUV,f) ']  %'])
disp(['SUV: ' num2str(media_SUV,f) ' +- ' num2str(std_SUV,f) '  [' num2str(min_SUV,f) '-' num2str(max_SUV,f) ']  %']);

%%


%% NU
figure
plot( positions(indices), NU(indices), 'k -');
%ylim([0 25]);

NUs_ok = NU(indices);
media_NU = mean(NUs_ok(:));
max_NU = max(NUs_ok(:));
std_NU = std(NUs_ok(:));
min_NU = min(NUs_ok(:));

hold on;
plot(positions(indices),media_NU,'k --');


title( ['Uniformidad Tomográfica ' modo '  ' fecha(7:8) '/' fecha(5:6) '/' fecha(1:4)  ] )
xlabel('Posición Axial [mm]');
ylabel('No Uniformidad [%]');

f = '%6.2g';
text(posicion_inicio+10, 23, ['NU: ' num2str(media_NU,f) ' +- ' num2str(std_NU,f) '  [' num2str(min_NU,f) '-' num2str(max_NU,f) ']  %']);
%text(mean(positions(:))/2,min_NU*1.2, ['NU: ' num2str(media_NU,f) ' +- ' num2str(std_NU,f) '  [' num2str(min_NU,f) '-' num2str(max_NU,f) ']  %'])
disp(['NU: ' num2str(media_NU,f) ' +- ' num2str(std_NU,f) '  [' num2str(min_NU,f) '-' num2str(max_NU,f) ']  %']);

%% CV
figure
plot( positions(indices), CV(indices), 'k -');
%ylim([0 5]);

CVs_ok = CV(indices);
media_CV = mean(CVs_ok(:));
max_CV = max(CVs_ok(:));
std_CV = std(CVs_ok(:));
min_CV = min(CVs_ok(:));

hold on;
plot(positions(indices),media_CV,'k --');


title( ['Coeficientes de Variación ' modo '  ' fecha(7:8) '/' fecha(5:6) '/' fecha(1:4)  ] )
xlabel('Posición Axial [mm]');
ylabel('CV [%]');

text(posicion_inicio+10, 1, ['CV: ' num2str(media_CV,f) ' +- ' num2str(std_CV,f) '  [' num2str(min_CV,f) '-' num2str(max_CV,f) ']  %']);
disp(['CV: ' num2str(media_CV,f) ' +- ' num2str(std_CV,f) '  [' num2str(min_CV,f) '-' num2str(max_CV,f) ']  %'])



disp(['Actividad presente al momento del scan [uCi]: ' num2str(A1/37000) ]);





