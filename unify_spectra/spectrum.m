clear all

% Doesnt work on linux
% opengl software

% Prepare the layout (size, position, etc)
varble='V';

plane180=[ 7]'; % yp = 4.7
plane180=[15]'; % yp ~ 15
plane550=[10]';
nplanes=numel(plane180);
iplane=1;

case0=5;
caseF=8;
ncase=caseF-case0+1;
NVert_plot=3;
NHorz_plot=4;

nribs550=[0 216 128 72];

% icase180=[1 3 5 8];

pos{1}=[.130 .700 .17 .24];
pos{2}=[.315 .700 .17 .24];
pos{3}=[.500 .700 .17 .24];
pos{4}=[.685 .700 .17 .24];

contours0=(0:.044:.22);

maxdark =0.15;
maxclear=1.00;
mapBW=(maxclear-maxdark)/(numel(contours0)-1);
mapBW=(maxclear:-mapBW:maxdark);
mapBW=[mapBW' mapBW' mapBW'];

figure(2),clf
set(gcf,'position',[400 200 1120 270*NVert_plot])
set(gcf,'PaperPositionMode','auto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jpanel=1:3
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Loading data of the 180 case
  % load utauRe
  % jribn = (icase180(jpanel));
  % utau  =  utau(jribn);
  % Retau = Retau(jribn);
  % jribn = (icase180(jpanel));
  utau  = 0.0617; % approx 
  Retau = 182; % approx

  switch jpanel
    case 1,
      file = ['spU016x032' '.dat'];
      grid = 'U';
    case 2,
      file = ['spV016x032' '.dat'];
      grid = 'V';
    case 3,
      file = ['spW016x032' '.dat'];
      grid = 'W'; 
  end

  disp(' ');
  disp(['READING ' file]);

  fid = fopen(file,'r');
          fread(fid,5,'int32');
  alp   = fread(fid,1,'float64');
  bet   = fread(fid,1,'float64');
          fread(fid,2,'int32');
  nband = fread(fid,1,'int32');
  mband = (nband+1)/2;
          fread(fid,91,'int32');
  N     = fread(fid,[4 nband+2],'int32');
  
  
  nyu = N(4,nband+1)-N(4,1)+2;
  nyv = N(3,nband+1)-N(3,1)+2;
  

          fread(fid,2,'int32');
  yu    = fread(fid,nyu,'float64')';
          fread(fid,nyu+2,'float64')';
          
          
  yv    = fread(fid,nyv,'float64')';
          fread(fid,nyv+2,'float64')';
          
  istat = fread(fid,1,'int32');
  nx  = N(1,mband+1)/2;
  nz  = N(2,mband+1)/2+1; %The plus 1 is due to a rouge +1 in spectra.f90...
  
  if(grid=='U' | grid=='W')
    ny = nyu;
    y=yu;
    
    % jpl: index of the plane
    jpl = -N(4,1)+1+plane180(1:nplanes);
          fread(fid,2,'int32');
    % spT: Spectrum
    spT   = fread(fid,[(nx+1)*nz*ny],'float64');
    fclose(fid); 
    
    
    % Divide by the number of samples
    spT = spT/istat;
    spT = reshape(spT,nx+1,nz,ny);
    % Sum statistics of the two halves of the channel
    spT = .5*(spT(:,:,1:end/2)+spT(:,:,end:-1:end/2+1));
    % Extract one plane only
    spT = squeeze(spT(:,:,jpl));

    Vrms2 = (sum(spT(:))-spT(1,1))/utau^2
    cntrlines = Vrms2*contours0;
    
    
    
elseif(grid=='V')
    ny = nyv;
    y=yv;
    
    % jpl: index of the plane
    jpl = -N(3,1)+1+plane180(1:nplanes);
          fread(fid,2,'int32');
    % spT: Spectrum
    spT   = fread(fid,[(nx+1)*nz*ny],'float64');
    fclose(fid); 
    
    % Divide by the number of samples
    spT = spT/istat;
    spT = reshape(spT,nx+1,nz,ny);
    % Sum statistics of the two halves of the channel
    spT = .5*(spT(:,:,1:ceil(end/2))+spT(:,:,end:-1:ceil(end/2)));
    % Extract one plane only
    spT = squeeze(spT(:,:,jpl));

    Vrms2 = (sum(spT(:))-spT(1,1))/utau^2
    cntrlines = Vrms2*contours0;


elseif(grid=='P')
    ny = nyu-2;
    y=yu;
    
    % jpl: index of the plane
    jpl = -N(4,1)+1+plane180(1:nplanes);
          fread(fid,2,'int32');
    % spT: Spectrum
    spT   = fread(fid,[(nx+1)*nz*ny],'float64');
    fclose(fid); 
    
    spT = spT/istat;
    spT = reshape(spT,nx+1,nz,ny);
    % Sum statistics of the two halves of the channel
    spT = .5*(spT(:,:,1:end/2)+spT(:,:,end:-1:end/2+1));
    % Extract one plane only
    spT = squeeze(spT(:,:,jpl));

    Vrms2 = (sum(spT(:))-spT(1,1))/utau^2
    cntrlines = Vrms2*contours0;
else
    break
end
  
  
  % Premultiplied spectrum
  kx  = [0:nx]';
  kz  = [0:nz-1];
  spT = spT.*(kx*kz);

  y = 1+y(1:ceil(end/2));

  spT180 = spT(2:end,2:end)/utau^2;
  % Rz180  = 2*pi/bet./nribn(jribn)*Retau;

  ypl = y(jpl)*Retau;

  disp(' ');
  disp(['PLOTTING ' varble ' SPECTRA']);
  disp(['AT y' num2str(ypl)]);
  disp(' ');
  Lx = 2*pi/alp./[1:nx  ]*Retau;
  Lz = 2*pi/bet./[1:nz-1]*Retau;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(1,4,jpanel),hold on
  %contourf(Lx,Lz,spT550',cntrlines)
  contourf(Lx,Lz,spT180',cntrlines)
  %contour(Lx,Lz,spT180',cntrlines,'k','LineWidth',2)
  shading flat
  set(gca,'clim',[cntrlines(1) cntrlines(end)])
  colormap(mapBW)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if jpanel ~= 1
  %   plot([20 40],Rz550*[1 1],'k','LineWidth',3)
  %   fill([20 30 30 20],Rz180*[1.05 1.05 0.95 0.95],'w','LineWidth',1)
  % end
  %pause

  subplot(1,4,jpanel)
  set(gca,'xsca','log','ysca','log')
  axis([20 Lx(1) 10 Lz(1)])
  set(gca,'Fontn','Times','FontSize',19,'linew',1.2)
  set(gca,'XTick',[100 1000])
  set(gca,'XMinorTick','on')
  set(gca,'YTick',[10 100 1000])
  set(gca,'YMinorTick','on')
  set(gca,'TickLe',[.03 0])
  xlabel('\lambda_x^+','position',[350 5.5 0]);
  if jpanel == 1
    ylabel('   \lambda_z^+','Ver','t','Hor','l');
  else
    set(gca,'YTickLabel',[])
  end
  box on
  set(gca,'Position',pos{jpanel})
  set(gca,'layer','top') 
  
  if(grid=='U')
      cntrlinesU = cntrlines;
      spTU = spT180;
      spTLx = Lx;
      spTLz = Lz;
      save('SpTU','spTLx','spTLz','spTU','cntrlinesU');
  elseif(grid=='W')
      cntrlinesW = cntrlines;
      spTW = spT180;
      spTLx = Lx;
      spTLz = Lz;
      save('SpTW','spTLx','spTLz','spTW','cntrlinesW');
  end
  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print -depsc2 V2pof
