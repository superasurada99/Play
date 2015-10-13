% Experiment:
    % ---------------
    % load processed experimental run:
    matName1 = frameName(['processedExp_' matName1Series],rpmAll(iRun));
    matName2 = '';
    matName = fullfile(dirNameSeries,matName1,matName2);
    eval(['load ' matName])
    
    % convert to metric:
    xMesh = xMesh/1000;
    yMesh = yMesh/1000;
    uMesh = uMesh/1000;
    vMesh = vMesh/1000;
    xCentre = xCentre/1000;
    yCentre = yCentre/1000;
    XGrid = XGrid/1000;
    YGrid = YGrid/1000;
    PsiGrid = PsiGrid/10^6;
    xyBot = xyBot/1000;
    xyTop = xyTop/1000;
    xyFlow = xyFlow/1000;
    
    % obtain velocity field from streamfunction:
    dx = xMesh(1,2)-xMesh(1,1);
    dy = yMesh(2,1)-yMesh(1,1);
    uMesh = ( 0.5*(PsiGrid(2:end,1:end-1)+PsiGrid(2:end,2:end)) ...
        - 0.5*(PsiGrid(1:end-1,1:end-1)+PsiGrid(1:end-1,2:end)) )/dy;
    vMesh = -( 0.5*(PsiGrid(1:end-1,2:end)+PsiGrid(2:end,2:end)) ...
        - 0.5*(PsiGrid(1:end-1,1:end-1)+PsiGrid(2:end,1:end-1)) )/dx;
    PsiRefExp = omegaExp/2*R^2; % use actual experimental value
    
    % tilt to inclined coordinate system:
    xMeshTilt = xMesh*cos(alpha) - yMesh*sin(alpha);
    yMeshTilt = xMesh*sin(alpha) + yMesh*cos(alpha);
    xyTopTilt = [xyTop(:,1)*cos(alpha)-xyTop(:,2)*sin(alpha),...
        xyTop(:,1)*sin(alpha)+xyTop(:,2)*cos(alpha)];
               
    % obtain hydrostatic pressure and effective viscosity:
    ySurfMeshTilt = interp1(xyTopTilt(:,1),xyTopTilt(:,2),xMeshTilt);
    pMesh = rho*g*cosAlpha*(ySurfMeshTilt-yMeshTilt);
    pMesh(pMesh<0) = nan;
    dynViscosityMesh = chi*D*sqrt(rho*pMesh); % effective viscosity
    
    % obtain shear rate:
    [dudxMesh,dudyMesh] = gradient(uMesh,dx,dy);
    [dvdxMesh,dvdyMesh] = gradient(vMesh,dx,dy);
    % Need to use invariant dissipation function ! (see Aris, 1962, p. 118):
    PhiMesh = -4*( dudxMesh.*dvdyMesh - (1/2*(dudyMesh+dvdxMesh)).^2 );
    gammaDotMesh = sqrt(PhiMesh);
    %gammaDotMesh = dudyMesh + dvdxMesh; % Wrong: not rotation invariant!
    dissipMesh = dynViscosityMesh.*gammaDotMesh.^2;
    
    % obtain divergence of kinetic energy flux:
    fluxKEx = 0.5*rho*(uMesh.^2+vMesh.^2).*uMesh;
    fluxKEy = 0.5*rho*(uMesh.^2+vMesh.^2).*vMesh;
    [dfluxKExdx,dfluxKExdy] = gradient(fluxKEx,dx,dy);
    [dfluxKEydx,dfluxKEydy] = gradient(fluxKEy,dx,dy);
    divKEfluxMesh = dfluxKExdx+dfluxKEydy;
    
    % form closed contours encircling drum atmosphere and rigidly rotating zone:
    thetaBoundary = 2*pi*(0:0.01:1)';
    xyBoundary = ones(size(thetaBoundary))*[xCentre yCentre] + R*[cos(thetaBoundary), sin(thetaBoundary)];
    [~,iUp] =min( (xyBoundary(:,1)-xyTop(1,1)).^2 + (xyBoundary(:,2)-xyTop(1,2)).^2 );
    [~,iDown] =min( (xyBoundary(:,1)-xyTop(end,1)).^2 + (xyBoundary(:,2)-xyTop(end,2)).^2 );
    xyAtm = [xyTop; xyBoundary(iDown:end,:); xyBoundary(1:iUp,:)];
    xyRigid = [xyBot(1,:); xyBoundary(iUp:iDown,:); flipud(xyBot)];
    
    % obtain distance functions:
    [m,n] = size(uMesh);
    xyMesh = [reshape(xMesh,m*n,1) reshape(yMesh,m*n,1)];
    DAtm = dpoly(xyMesh,xyAtm); % from Per-Olof Persson's distMesh toolbox
    DAtmMesh = reshape(DAtm,m,n);
    DRigid = dpoly(xyMesh,xyRigid);
    DRigidMesh = reshape(DRigid,m,n);
    
    % set atmosphere values to nan and rigid body rotation values to 0:
    dissipMesh(DAtmMesh<0) = nan;
    dissipMesh(DRigidMesh<0) = 0;
    divKEfluxMesh(DAtmMesh<0) = nan;
    divKEfluxMesh(DRigidMesh<0) = 0;
    
    figure
%    subplot(2,2,2*(iRun-1)+2)
    switch mapOption
        case 1
            mypcolor(xMesh*1000,yMesh*1000,dissipMesh)
        case 2
            mypcolor(xMesh*1000,yMesh*1000,divKEfluxMesh) 
        case 3
            mypcolor(xMesh*1000,yMesh*1000,divKEfluxMesh+dissipMesh) 
    end
    hold on
    contour(XGrid*1000,YGrid*1000,PsiGrid,PsiRefExp*[-1:0.125:-0.125],'k-','lineWidth',lineWidth)
    plot(xyFlow(:,1)*1000,xyFlow(:,2)*1000,'k-','lineWidth',2*lineWidth)
    plot(xyRigid(:,1)*1000,xyRigid(:,2)*1000,'k-','lineWidth',2*lineWidth)
    shading flat
    axis equal
    axis(xyAxis*1000)
    caxis(caxScaleAll(iRun)*[-1 1])
    axis off
    