%%
%             B: [1x53940 double]
%            Xv: [53940x6 double]
%            Yv: [53940x6 double]
%            Av: [53940x5 double]
%      xy_names: {'CoM'  'head'  'head-mid'  'mid'  'tail-mid'  'tail'}
%     ang_names: {'CoM'  'head-mid'  'tail-mid'  'body-bend'  'tail-bend'}
%     
   
load D130304_jove.mat;

%% frame time conversion
FLIM = round(interp1([10 550], [129 8230], [294.5 296.5])); %[4397, 4427]
aviobj = VideoReader(D.S.vidFname);
MOV = read(aviobj, FLIM);
MOV = MOV(:,:,1,:);
implay(MOV);

%%
IMa = MOV(:,:,:,1); IMa = IMa(:,:,1);
IMb = MOV(:,:,:,end); IMb = IMb(:,:,1);
intrng = stretchlim(IMa, [0 .85]);
IMa = imadjust(IMa, intrng);
IMb = imadjust(IMb, intrng);
figure; imshow(IMa);
figure; imshow(IMb);

%%
IDX = find(D.T>=294.5 & D.T<=296.55); %sample rate 100 Hz
XV = D.Xv(IDX(1:20:end),:);
YV = D.Yv(IDX(1:20:end),:);
nt = size(XV, 1);

%

figure; hold on;
axis square;
plot(75*cos(0:.01:2*pi), 75*sin(0:.01:2*pi), 'k'); %draw aquarium

% time point color
mrColor = jet(nt);

for i=1:nt
    x = XV(i, 2:end);
    y = YV(i, 2:end);
    n = numel(x);
    ni = linspace(1, n, n*8); %8x interp
    xi = interp1(1:n, x, ni, 'cubic');
    yi = interp1(1:n, y, ni, 'cubic');
%     plot(xi, yi, 'color', .25*[1 1 1]);    
    plot(xi, yi, 'color', mrColor(i,:));
end

XV = D.Xv(IDX(1:end),:);
YV = D.Yv(IDX(1:end),:);

x = XV(:, 2);
y = YV(:, 2);
plot(x, y, 'k:');    

x = XV(:, 6);
y = YV(:, 6);
plot(x, y, 'w:');   

% axis([-50 -10 25 50]);

%% superimpose first and last images
IMc = IMa/2 + IMb/2;
intrng = stretchlim(IMc, [.2 .85]);
IMc = imadjust(IMc, intrng);
figure; imshow(IMc);

%% eod rate
R = D.R(IDX);
figure; plot(R);
% [~, RANK] = sort(R);
% RANK = RANK / max(RANK);
R = (R-min(R))/(max(R)-min(R));
figure; plot(R);

%% plot
dt = 2;
IDX = find(D.T>=294.5 & D.T<=296.5); %sample rate 100 Hz
XV = D.Xv(IDX(1:dt:end),:);
YV = D.Yv(IDX(1:dt:end),:);
nt = size(XV, 1);

R = D.R(IDX(1:dt:end));
R = (R-min(R))/(max(R)-min(R));

figure; hold on;
axis square;
plot(75*cos(0:.01:2*pi), 75*sin(0:.01:2*pi), 'k'); %draw aquarium

% time point color
mrColor = jet(255);
RIDX = ceil(R * 255);
RIDX(RIDX<=0)=1;

for i=1:nt
    x = XV(i, 2:end);
    y = YV(i, 2:end);
    n = numel(x);
    ni = linspace(1, n, n*8); %8x interp
    xi = interp1(1:n, x, ni, 'cubic');
    yi = interp1(1:n, y, ni, 'cubic');
%     plot(xi, yi, 'color', mrColor(RIDX(i),:));
    plot(xi(1), yi(1), '.', 'color', mrColor(RIDX(i),:));
    mrColor(RIDX(i),:)
end