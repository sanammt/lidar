clear all
close all
% Load the two files extracted from FDTD run
% The farfield |E^2| data
% Angles of the farfield (in farfield plot we are looking at |E^2| vs Angle
farfieldE2File='f2d_full_26.txt';
farfieldAnlgle='fangle_26.txt';
%
farfieldE2Data=load(farfieldE2File);
farfieldAngleData=load(farfieldAnlgle);

% Get sizi of the data and check for a mismath in datasize
% Datasize repesents number of lambda/points in the FOV in grating coupler
% scanning
sizeE2Data=size(farfieldE2Data)
sizeAngleData=size(farfieldAngleData)
if (sizeE2Data(2) ~=sizeAngleData(2))
    disp('WARNNIG, size mismath in data')
end
%
tol=0.002;  % the tol is used to calculate 3db FWHM
lamda=linspace(1335,1360,sizeE2Data(2))

for i=1:1:sizeE2Data(2)
% for i=1:1:1
    y{i}=farfieldE2Data(:,i); % farfield |E2| 
    maxy{i} = max(y{i});      % peak value for farfield |E2|  
    f{i}=find(y{i}==maxy{i}); % index of farfield |E2|  peak
    y1{i}=y{i}/maxy{i};       % Normalized farfield |E2|  
    y_fwhm=1/sqrt(2);         % -3db power for a normalized |E2|
    theta_peakPower{i}=farfieldAngleData(f{i},1);   % radiation angle at each lambda
    f_fwhm{i}=find((y1{i}<y_fwhm+tol) & ((y1{i}>y_fwhm-tol))); % index of farfield |E2|  for -3db loss (fwhm) it must be only 2 indecies
    while (size(f_fwhm{i}) < 2)
        tol=tol+0.0005;
        f_fwhm{i}=find((y1{i}<y_fwhm+tol) & ((y1{i}>y_fwhm-tol)));
    end
    while (size(f_fwhm{i}) > 2)
        tol=tol-0.00005;
        f_fwhm{i}=find((y1{i}<y_fwhm+tol) & ((y1{i}>y_fwhm-tol)));
    end
    if (size(f_fwhm{i}) ~= 2)
        disp('WARNNIG')
    end
    tol=0.002;
    theta1=farfieldAngleData(f_fwhm{i}(1),1);
    theta2=farfieldAngleData(f_fwhm{i}(2),1);
    fwhm{i}=abs(theta2-theta1);
    
    disp('----------')
end
fov=abs(theta_peakPower{sizeE2Data(2)} - theta_peakPower{1})
plot(farfieldAngleData(:,1),y1{1})
hold on
for i=2:1:sizeE2Data(2)
    plot(farfieldAngleData(:,1),y1{i})
end
hold on



x=linspace(-100,100,1000);
y=0+0*x+11.3827*x.^2+0.0005*x.^3-0.0011*x.^4;
fwhm



