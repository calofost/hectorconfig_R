#############################################################################
### Purpose: Configure the position and orientation of HECTOR probes after
### greedy field and target selection. This version uses optimisation to make
### an initial guess of angles for non-overlapping probes.
### Date: 15 December 2018
### Author: Caroline Foster
#############################################################################

# Dependencies are imported by the hectorconfig package. This file is
# evaluated in a private engine environment by .hector_new_engine(), rather
# than attaching packages to the user's search path.

#Probe shape definition as global variables (currently has a head/tip, a body/probe and a cable zone of avoidance (all should be in mm):
#plate_scale"/mm scaling
#plate_scale <<- 15.008
#! Added by Sam to use the same constants file as Ayoan's code uses
# The factory supplies hector_constants before evaluating this file. Keeping
# the check here makes accidental direct sourcing fail clearly instead of
# silently using a machine-specific path.
if (!exists("hector_constants", inherits = FALSE)) {
  stop("HECTOR constants were not supplied to the configuration engine.", call. = FALSE)
}

# Small Cartesian helpers are kept local to the engine. The original script
# obtained these from several attached packages, which made package namespace
# loading fragile. All distances used by this configuration code are in the
# plate's Cartesian millimetre coordinates (lonlat=FALSE).
cosd <- function(x) cos(x * pi / 180)
sind <- function(x) sin(x * pi / 180)
asind <- function(x) asin(x) * 180 / pi
atand <- function(x) atan(x) * 180 / pi
atan2d <- function(y, x) atan2(y, x) * 180 / pi

Rotate <- function(x, y = NULL, mx = NULL, my = NULL, theta = pi / 3,
                   asp = 1) {
  if (is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  }
  if (is.null(mx)) mx <- 0
  if (is.null(my)) my <- 0
  x0 <- x - mx
  y0 <- y - my
  list(
    x = x0 * cos(theta) - y0 * sin(theta) + mx,
    y = x0 * sin(theta) + y0 * cos(theta) + my
  )
}

.as_xy <- function(points) {
  if (is.data.frame(points) && all(c("x", "y") %in% names(points))) {
    return(as.matrix(points[, c("x", "y"), drop = FALSE]))
  }
  if (is.null(dim(points))) {
    if (length(points) != 2) stop("Cartesian points must have two coordinates.", call. = FALSE)
    return(matrix(as.numeric(points), nrow = 1, ncol = 2))
  }
  points <- as.matrix(points)
  if (ncol(points) == 2) return(points)
  if (nrow(points) == 2) return(t(points))
  stop("Cartesian points must have two coordinates.", call. = FALSE)
}

pointDistance <- function(p1, p2 = NULL, lonlat = FALSE) {
  if (isTRUE(lonlat)) {
    stop("The HECTOR engine only supports Cartesian distances.", call. = FALSE)
  }
  a <- .as_xy(p1)
  if (is.null(p2)) {
    return(sqrt(outer(a[, 1], a[, 1], "-")^2 +
                outer(a[, 2], a[, 2], "-")^2))
  }
  b <- .as_xy(p2)
  if (nrow(a) == 1) {
    return(sqrt(rowSums((b - matrix(a[1, ], nrow(b), 2, byrow = TRUE))^2)))
  }
  if (nrow(b) == 1) {
    return(sqrt(rowSums((a - matrix(b[1, ], nrow(a), 2, byrow = TRUE))^2)))
  }
  if (nrow(a) == nrow(b)) return(sqrt(rowSums((a - b)^2)))
  sqrt(outer(a[, 1], b[, 1], "-")^2 +
       outer(a[, 2], b[, 2], "-")^2)
}

polyarea <- function(x, y) {
  x_next <- c(x[-1], x[1])
  y_next <- c(y[-1], y[1])
  sum(x * y_next - x_next * y) / 2
}

spDists <- function(x, y = NULL, diagonal = FALSE, longlat = FALSE) {
  if (is.null(y)) return(pointDistance(x, lonlat = longlat))
  pointDistance(x, y, lonlat = longlat)
}

spDistsN1 <- function(pts, pt, longlat = FALSE) {
  pointDistance(pts, pt, lonlat = longlat)
}


#Field dimensions and parameters:
fov <- hector_constants$HECTOR_plate_radius*2  #diameter of the field in mm
skybuffer <- 8 #clear-edge buffer
wceg <- 2.* 15  #Width of the cable exit gap (this is the 2*n where n=15mm initially in Julia's notes)
ceg_pos <- cbind(x=fov/2.*cosd(c(-90,30,150)),y=fov/2.*sind(c(-90,30,150)), deg=c(-90,30,150))
E0 <- sqrt((wceg/2.)**2-(9)**2) #E is the line joining the solid ferule end to the centre of the exit point. 
#E0=SQRT((n^2)-81). If E<E0, then need to point to the next least contaminated exit. This step prevents a ferule from 
#blocking an exit such that other cables will not get past.
fbr <- 45 #Fibre bend radius

#Galaxy and standard star bundle dimensions:
tip_w <- 14.5/sqrt(2)  
tip_l <- 14.5/sqrt(2) 
# Testing- make this 18.5 from 14.5
probe_w <- 18.5
probe_l <- 45 
cable_w <- 8.5 
cable_l <- 45 #value changed from 38 on 22 August 2022
excl_radius <- sqrt(2.*((tip_w/2)^2)) * 1.05 #Added by Sam after an issue where two probes were ever so slightly too close in Feb2022 commissioning
dngalprobes <- 19
nstdprobes <- 2
ngprobesmin <- 3
ngprobesmax <- 6
ngalprobes <- dngalprobes
nprobes <- dngalprobes #Start with the number of galaxy targets. There will be an additional 2 probes for standard stars.

#Guide bundle dimensions:
#Julia: The guide bundles should be the same 14mm width but shorter than the hexabundle ferules by 20mm ~50mm total plus cable bit.
gtip_w <- 14.5/sqrt(2)
gtip_l <- 14.5/sqrt(2)
gprobe_w <- 14.5
gprobe_l <- 45
gcable_w <- 8.5
gcable_l <- 45 #value changed from 38 on 22 August 2022
gexcl_radius <- sqrt(2.*((gtip_w/2)^2))

#Setting timeout time (in cpu seconds, I think) for how long each attempt lasts before moving on. 
#The code will try for at most 15 times that before starting to swap a probe.
#timeout<<-180 #180 seconds x 15 is about half an hour and is about right. Less than about 90 seconds led to problems in highly conflicted field..
timeout <- 600 #Increased to account for longer conflict calculation time due to rounded heads so code still has enough time time try enough angles variations.
timeout_behaviour = 'error'

#angle step going around the head to approximate circle as a polygon. The smaller angles are more precice but lead to longer configuration time.
delta_poly <- 5 #5 degrees leads to a 0.1% radius buffer around the head and configuration seems to take a bit longer than a square head.

#Compute the (horizontal/RA/x-axis) distance to the field edge.
dist2fieldedge <- function(x,y){
  dist=sqrt((fov/2.)**2-y**2)-abs(x)
  return(dist)
}

#Compute range of allowed angles for given probe:
computeAnglesRange <- function(x, y, ceg, delta_ang=20,guide=FALSE){
  
  ang=as.numeric(first_guess_angs(pos=cbind(x,y),cegs=ceg))*180./pi
  
  # From Julia: " fbr is the minimum fibre bend radius, fbr=45mm. If E<2fbr then the maximum configuration 
  #angle that should be set for that ferule is +/- arcsin(E/2fbr).
  #Implemented this limit on the fibre bend radius in the next 2 lines. This will be applied regardless of the allowed delta_ang
  E=pointDistance(p1=ceg_pos[ceg,c('x','y')], p2=c(x,y), lonlat=FALSE) - (probe_l+tip_l/2.) ##E is the line joining the solid ferule end to the centre of the exit point.
  if (E < 2*fbr){
    new_delta_ang=asind(E/(2*fbr))
    if (new_delta_ang < delta_ang) delta_ang=new_delta_ang #Only apply correction if the angle is more restrictive.
  }
  
  minang=ang-delta_ang
  maxang=ang+delta_ang
  
  r0=fov/2
  if (guide){
    r1=gcable_l+gprobe_l+gtip_l
  }else{
    r1=cable_l+probe_l+tip_l
  }
  d=sqrt(x^2 + y^2)
  
  #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
  #Runs between -90 and 270.
  azAng=(pi+atan2(y,x))/pi*180
  if (azAng > 270) azAng=azAng-360
  if (azAng < -90) azAng=azAng+360
  
  #If the probe can touch the edge of the field:
  if (d+r1 >= r0){
    a=(r0^2-r1^2+d^2)/(2*d)
    h=sqrt(r0^2-a^2)
    x1=a*x/d+h*(y)/d      
    y1=a*y/d-h*(x)/d       
    x2=a*x/d-h*(y)/d      
    y2=a*y/d+h*(x)/d
    
    ang1=atan2(y1-y,x1-x)/pi*180-180
    ang2=atan2(y2-y,x2-x)/pi*180-180
    
    ang1=adjust_angle(theta=ang1,fgang=ang)
    ang2=adjust_angle(theta=ang2,fgang=ang)
    
    elminang=min(ang1,ang2) #edge limit minimum angle
    elmaxang=max(ang1,ang2) #edge limit maximum angle
    
    if (guide) {
      dcegth=(gprobe_l+gcable_l+gtip_l/2.)
    } else {
      dcegth=(probe_l+cable_l+tip_l/2.)
    }
    
    #computing the extremes of the relevant cable exit gap.
    dx=as.numeric(wceg/2.* sind(ceg_pos[ceg,'deg']))
    dy=as.numeric(-1.*wceg/2.* cosd(ceg_pos[ceg,'deg']))
    p2=c(ceg_pos[ceg,'x']-dx, ceg_pos[ceg,'y']-dy)
    p3=c(ceg_pos[ceg,'x']+dx, ceg_pos[ceg,'y']+dy)
    
    #Computing distance between target (p1), and p2 and p3.
    d12=pointDistance(p1=c(x,y),p2=p2, lonlat=FALSE)
    d13=pointDistance(p1=c(x,y),p2=p3, lonlat=FALSE)
    
    #We have 2 cases:
    
    #Case 1: probe goes through the cable exit gap: With change to E0 from 22 Aug 2022, case 1 should never occur.
    #allow for maximal angle change across the gap, unless it is larger than delta_ang.
    #In this case, the angle range should include the initial ang guess.
    if ((d12 < dcegth) | (d13 < dcegth)) {

      #Computing angles range within the cable exit gap
      psi2=asind(cable_w/2./d12)
      psi3=asind(cable_w/2./d13)

      ang1=atan2d((y-p2[2]),(x-p2[1]))-psi2
      ang2=atan2d((y-p3[2]),(x-p3[1]))+psi3

      ang1=adjust_angle(theta=ang1,fgang=ang)
      ang2=adjust_angle(theta=ang2,fgang=ang)

      minang=min(ang1,ang2)
      maxang=max(ang1,ang2)
      
      # temp=min(ang1,ang2)
      # ang2=max(ang1,ang2)
      # ang1=temp
      # 
      # if (ang1 > elminang) {
      #   minang=ang1
      # }
      # if (ang2 < elmaxang) {
      #   maxang=ang2
      # }
    #Case 2: the probe does not go through the exit gap but may hit the side of the FoV.
    #Pick the angles range that contains the first guess ang..
    } else{
      if ((ang < elminang) & (ang < elmaxang)) {
        if (maxang > elminang) maxang=elminang
      } else if ((ang > elminang) & (ang > elmaxang)){
        if (minang < elmaxang) minang=elmaxang
      } else if ((ang > minang) & (ang < maxang)){
        if (minang < elminang) {
          minang = elminang
        } else if (maxang > elmaxang) {
          maxang = elmaxang
        } 
      }
    }
  }
  if (minang > maxang) print(paste('Houston, we have a problem... minang=', as.character(minang),' and maxang=',as.character(maxang)))
  return(c(minang,maxang))
}

#Wanna pick a 360 degree range that includes that first-guess angle (fgang).
adjust_angle <- function(theta,fgang){
  if (theta < fgang-180){
    theta=theta+360
  } else if (theta > fgang + 180){
    theta=theta-360
  }
  return(theta)
}

#############################################################################
#test code: delete when done with finalising the above.
#for (probe in TPs){
#angRange=computeAnglesRange(pos[probe,'x'],pos[probe,'y'], ceg=cegs[probe], ang=angs[probe], delta_ang=90)/180*pi
#ang_ranges[probe,]=angRange
#}
#draw_pos(x=pos[,'x'],y=pos[,'y'])
#defineprobe(xs=pos[,'x'],ys=pos[,'y'], angs=angs, interactive=TRUE)

#defineprobe(xs=pos[,'x'],ys=pos[,'y'], angs=ang_ranges[,1], interactive=TRUE, col='green')
#defineprobe(xs=pos[,'x'],ys=pos[,'y'], angs=ang_ranges[,2], interactive=TRUE, col='red')

#############################################################################
### Randomly allocate nprobes positions on 2 degree field of view for initial 
### simulation and code testing.Will need to test on real data (GAMA G23?) later
### to get clustering right..
#############################################################################
random_allocation <- function (nrandprobes=ngalprobes+6, nstdprobes=20, nguideprobes=20) {
  #nprobes is the number of positions. HECTOR only has nprobes, but also adding 6 spare ones.
  #fov is the field of view in degrees (global variable).
  #excl_radius is the exclusion radius in arcsec for each object (hexabundle head size, global variable).
  found=FALSE
  while (found==FALSE){
    r <- sqrt(runif(n=nrandprobes/2, min=0, max=fov/2.-excl_radius*2))
    r <- c(r,runif(n=nrandprobes/2, min=0, max=fov/4.-excl_radius*2))
    t <- 2*pi*runif(n=nrandprobes)
    x <- r*cos(t) ; y <- r*sin(t)
    test=as.data.frame(spDists(cbind(x,y)))
    if (length(test[test>0 & test<excl_radius*2.])==0) {found=TRUE}
  }
  
  sx=rep(NaN,n=nstdprobes)
  sy=rep(NaN,n=nstdprobes)
  found=FALSE
  i=1
  while (found==FALSE) {
    r <- sqrt(runif(n=1, min=0,max=fov/2.- excl_radius*2))
    t <- 2*pi*runif(n=1)
    sx[i] <- r*cos(t) ; sy[i] <- r*sin(t)
    test=as.data.frame(spDistsN1(pts=cbind(x,y),pt=c(sx[i],sy[i])))
    if (length(test[test<excl_radius*2.])==0) {
      i=i+1
      if (i>nstdprobes) {found=TRUE}
    }
  }
  
  gx=rep(NaN,n=nguideprobes)
  gy=rep(NaN,n=nguideprobes)
  found=FALSE
  i=1
  while (found==FALSE) {
    r <- sqrt(runif(n=1, min=0,max=fov/2.- gexcl_radius*2))
    t <- 2*pi*runif(n=1)
    gx[i] <- r*cos(t) ; gy[i] <- r*sin(t)
    test=as.data.frame(spDistsN1(pts=cbind(x,y),pt=c(gx[i],gy[i])))
    if (length(test[test<gexcl_radius*2.])==0) {
      i=i+1
      if (i>nguideprobes) {found=TRUE}
    }
  }
  
  return(list(pos=cbind(x,y), stdpos=cbind(x=sx,y=sy), guidepos=cbind(x=gx,y=gy)))
}

#Draw buttons around each position to be allocated.
draw_pos <- function(x,y, label_probe=TRUE, overwrite=TRUE){
  if (overwrite){ #set overwrite to true to start the plots from scratch.
    magplot(0,0,pch='x',col="red",xlim=c(-fov/2.,fov/2.),ylim=c(-fov/2.,fov/2.), asp=1, xlab='X (mm)', ylab='Y (mm)')
    draw.circle(0,0,radius=fov/2.,border='red')
    draw.circle(0,0,radius=fov/2.+skybuffer,border='black')
    #Draw the cable exit gaps in cyan
    for (cceg in c(1:length(ceg_pos[,1]))){
      dx=wceg/2.* sind(ceg_pos[cceg,'deg'])
      dy=-1.*wceg/2.* cosd(ceg_pos[cceg,'deg'])
      lines(c(ceg_pos[cceg,'x']-dx,ceg_pos[cceg,'x']+dx),c(ceg_pos[cceg,'y']-dy,ceg_pos[cceg,'y']+dy), lwd=3, col='cyan')
    }
    text(as.character(c(1:3)), x=ceg_pos[,'x'],y=ceg_pos[,'y'])
  }
  for (i in c(1:length(x))){
    if (i<=ngalprobes){
      draw.circle(x[i],y[i],radius=excl_radius, col='lightblue')
    } else{
      draw.circle(x[i],y[i],radius=excl_radius, col='lightblue', lty=3)
    }
    if (label_probe){
      text(x[i],y[i],as.character(i), cex=0.5)
    }
  }
}

#Add buttons around potential standard stars:.
draw_stds <- function(xs,ys,label_probe=TRUE){
  for (i in c(1:length(xs))){
    draw.circle(xs[i],ys[i],radius=excl_radius, col='lightpink', lty=3)
    if (label_probe) text(xs[i],ys[i],as.character(i), cex=0.5)
  }
}

#Add buttons around potential standard stars:.
draw_guides <- function(xg,yg){
  for (i in c(1:length(xg))){
    draw.circle(xg[i],yg[i],radius=gexcl_radius, col='yellow', lty=3)
  }
}

#Draw all pos:
refresh_plot <- function(){
  if (isTRUE(visualise)) {
    try(dev.flush(), silent=TRUE)
    try(utils::flush.console(), silent=TRUE)
    if (is.finite(plot_pause) && plot_pause > 0) {
      Sys.sleep(plot_pause)
    }
  }
  invisible(NULL)
}

draw_all <- function(x,y,xs,ys,xg,yg, interactive=TRUE){
  if (interactive){
    draw_pos(x,y)
    draw_guides(xg,yg)
    draw_stds(xs,ys)
    refresh_plot()
  }
}

#This function computes the probe polygon shape defined by a series of data points. This function needs to remain flexible enough 
#to allow for changes in the probe shape as the HECTOR design evolves.
defineprobe <- function(xs,ys,angs,interactive=FALSE, col='black'){
  probes=c()
  for (i in c(1:length(xs))) {
    x=xs[i]
    y=ys[i]
    rbuff=2-cosd(delta_poly/2)#buffer to account for the approximation of small angles difference in radius.
    #probe includes the whole exclusion area, so including the exclusion radius around the probe head.
    probe=cbind(c(x+rbuff*excl_radius*sind(seq(-1.*atand(tip_l/tip_w)+delta_poly,180+atand(tip_l/tip_w)-delta_poly,delta_poly)), 
                  x-tip_l/2.,
                  x-tip_l/2.,
                  x-tip_l/2.-probe_l, 
                  x-tip_l/2.-probe_l, 
                  x-tip_l/2.-probe_l-cable_l, 
                  x-tip_l/2.-probe_l-cable_l, 
                  x-tip_l/2.-probe_l, 
                  x-tip_l/2.-probe_l, 
                  x-tip_l/2., 
                  x-tip_l/2.),
                c(y-rbuff*excl_radius*cosd(seq(-1.*atand(tip_l/tip_w)+delta_poly,180+atand(tip_l/tip_w)-delta_poly,delta_poly)), 
                  y+tip_w/2., 
                  y+probe_w/2., 
                  y+probe_w/2., 
                  y+cable_w/2.,
                  y+cable_w/2.,
                  y-cable_w/2.,
                  y-cable_w/2.,
                  y-probe_w/2., 
                  y-probe_w/2., 
                  y-tip_w/2.))
    #footprint is just for display purposes when plotting, it is not used outside this function.
    lpoly=length(probe[,1])
    footprint=cbind(c(x+tip_w/2., 
                      probe[(lpoly-9):lpoly,1], 
                      x+tip_w/2.),
                    c(y+tip_w/2., 
                      probe[(lpoly-9):lpoly,2],
                      y-tip_w/2.))
    
    #Rotating probe.
    probetemp=cbind(Rotate(probe,theta=angs[i], mx=x, my=y)$x, Rotate(probe,theta=angs[i], mx=x, my=y)$y)
    probe=probetemp
    colnames(probe)=c(paste('X',as.character(i), sep=''),paste('Y',as.character(i),sep=''))
    
    footprintemp=cbind(Rotate(footprint,theta=angs[i], mx=x, my=y)$x, Rotate(footprint,theta=angs[i], mx=x, my=y)$y)
    footprint=footprintemp
    colnames(footprint)=c(paste('X',as.character(i), sep=''),paste('Y',as.character(i),sep=''))
    
    if (interactive) {polygon(footprint, border=col)}#Drawing the actual shape of the probe, not the exclusion area.
    assign(paste('probe',as.character(i), sep=''), probe)
    probes=cbind(probes,probe)
  }
  if (interactive) {refresh_plot()}
  return(as.data.frame(probes))
}

#This function computes the guide probe polygon shape defined by a series of data points.
defineguideprobe <- function(gxs,gys,gangs,interactive=FALSE){
  gprobes=c()
  for (i in c(1:length(gxs))) {
    x=gxs[i]
    y=gys[i]
    rbuff=2-cosd(delta_poly/2)#buffer to account for the approximation of small angles difference in radius.
    #gprobe includes the whole exclusion area, so including the exclusion radius around the probe head.
    gprobe=cbind(c(x+rbuff*gexcl_radius*sind(seq(-1.*atand(gtip_l/gtip_w)+delta_poly,180+atand(gtip_l/gtip_w)-delta_poly,delta_poly)), 
                   x-(gtip_l-gtip_w/2.),
                   x-(gtip_l-gtip_w/2.),
                   x-(gtip_l-gtip_w/2.)-gprobe_l, 
                   x-(gtip_l-gtip_w/2.)-gprobe_l, 
                   x-(gtip_l-gtip_w/2.)-gprobe_l-gcable_l, 
                   x-(gtip_l-gtip_w/2.)-gprobe_l-gcable_l, 
                   x-(gtip_l-gtip_w/2.)-gprobe_l, 
                   x-(gtip_l-gtip_w/2.)-gprobe_l, 
                   x-(gtip_l-gtip_w/2.), 
                   x-(gtip_l-gtip_w/2.)),
                 c(y-rbuff*gexcl_radius*cosd(seq(-1.*atand(gtip_l/gtip_w)+delta_poly,180+atand(gtip_l/gtip_w)-delta_poly, delta_poly)), 
                   y+gtip_w/2., 
                   y+gprobe_w/2., 
                   y+gprobe_w/2., 
                   y+gcable_w/2.,
                   y+gcable_w/2.,
                   y-gcable_w/2.,
                   y-gcable_w/2.,
                   y-gprobe_w/2., 
                   y-gprobe_w/2., 
                   y-gtip_w/2.))
    
    #footprint is just for display purposes when plotting, it is not used outside this function.
    lpoly=length(gprobe[,1])
    footprint=cbind(c(x+gtip_w/2., 
                      gprobe[(lpoly-9):lpoly,1], 
                      x+gtip_w/2.),
                    c(y+gtip_w/2., 
                      gprobe[(lpoly-9):lpoly,2],
                      y-gtip_w/2.))
    
    #Rotating probe.
    gprobetemp=cbind(Rotate(gprobe,theta=gangs[i], mx=x, my=y)$x, Rotate(gprobe,theta=gangs[i], mx=x, my=y)$y)
    gprobe=gprobetemp
    colnames(gprobe)=c(paste('X',as.character(i), sep=''),paste('Y',as.character(i),sep=''))
    footprintemp=cbind(Rotate(footprint,theta=gangs[i], mx=x, my=y)$x, Rotate(footprint,theta=gangs[i], mx=x, my=y)$y)
    footprint=footprintemp
    colnames(footprint)=c(paste('X',as.character(i), sep=''),paste('Y',as.character(i),sep=''))
    
    if (interactive) {polygon(footprint, border='darkgrey')} #Drawing the actual shape of the probe, not the exclusion area.
    assign(paste('gprobe',as.character(i), sep=''), gprobe)
    gprobes=cbind(gprobes,gprobe)
  }
  if (interactive) {refresh_plot()}
  return(as.data.frame(gprobes))
}

# Check the complete exclusion polygon, not only the probe centre, against
# the circular HECTOR field.  The circle is convex, so checking every polygon
# vertex is sufficient to ensure that all straight polygon edges remain inside
# it.  A tiny tolerance avoids rejecting a point solely because of round-off at
# the plate boundary.
probe_fov_status <- function(polygons, radius=fov/2., tolerance=1e-7){
  xcols=grep('^X[0-9]+$', names(polygons), value=TRUE)
  if (length(xcols) == 0) {
    return(list(ok=TRUE, outside=integer(0)))
  }

  outside=vapply(xcols, function(xcol) {
    ycol=sub('^X', 'Y', xcol)
    if (!(ycol %in% names(polygons))) return(TRUE)
    x=polygons[[xcol]]
    y=polygons[[ycol]]
    any(!is.finite(x) | !is.finite(y)) ||
      any((x^2 + y^2) > (radius + tolerance)^2)
  }, logical(1))

  list(ok=!any(outside), outside=which(outside))
}

configuration_fov_status <- function(probes, gprobes=NULL){
  target_status=probe_fov_status(probes)
  guide_status=if (is.null(gprobes)) {
    list(ok=TRUE, outside=integer(0))
  } else {
    probe_fov_status(gprobes)
  }
  list(
    ok=target_status$ok && guide_status$ok,
    targets=target_status,
    guides=guide_status
  )
}

#Check overlapping probes (This is a bit slow, could it be optimised?):
find_probe_conflicts <- function(probes, pos, angs){
  probe_conflicts=data.frame()
  if (dim(pos)[1] == 1) {
    return(probe_conflicts)
  }
  lpoly=length(probes[,1]) #length of polygon defining probe.
  for (i in c(1:(length(pos[,'x'])-1))) {
    p1 <- data.frame(PID=rep(i, lpoly), POS=(1:lpoly), X=probes[,paste('X',as.character(i), sep='')], Y=probes[,paste('Y',as.character(i), sep='')])
    for (j in c((i+1):length(pos[,'x']))) {
      p2 <- data.frame(PID=rep(j, lpoly), POS=(1:lpoly), X=probes[,paste('X',as.character(j), sep='')], Y=probes[,paste('Y',as.character(j), sep='')])
      if (!is.null(joinPolys(p1,p2))) {
        junction=joinPolys(p1,p2)
        #Compute the percentage area of a probe that overlaps.100% means they overlap perfectly. 0% means no overlap.
        areapoly= round(polyarea(junction$X,junction$Y)/polyarea(p1$X,p1$Y) * 100, digits=1)
        probe_conflicts=rbind(probe_conflicts, c(i,j,areapoly))
      }
    }
  }
  nconflicts=dim(probe_conflicts)[1]
  if (nconflicts>0) {
    colnames(probe_conflicts)=c('Probe_1st','Probe_2nd','PercentArea')
  }
  return(as.data.frame(probe_conflicts))
}

#Looking for guide/target/standard conflicts with the target/standard probes fixed.
find_guide_conflicts <- function(probes, pos, angs, gprobes, gpos, gangs, verbose=TRUE){
  probe_conflicts=data.frame()
  lpoly=length(probes[,1]) #length of polygon defining probe.
  glpoly=length(gprobes[,1])
  if (verbose) print("Identifying clashes with target and/or standard probes")
  #First identify clashes with target/standard probes:
  for (i in c(1:(length(pos[,'x'])))) {
    p1 <- data.frame(PID=rep(i, lpoly), POS=(1:lpoly), X=probes[,paste('X',as.character(i), sep='')], Y=probes[,paste('Y',as.character(i), sep='')])
    for (j in 1:length(gpos[,'x'])) {
      p2 <- data.frame(PID=rep(j, glpoly), POS=(1:glpoly), X=gprobes[,paste('X',as.character(j), sep='')], Y=gprobes[,paste('Y',as.character(j), sep='')])
      if (!is.null(joinPolys(p1,p2))) {
        junction=joinPolys(p1,p2)
        #Compute the percentage area of a probe that overlaps.100% means they overlap perfectly. 0% means no overlap.
        areapoly= round(polyarea(junction$X,junction$Y)/polyarea(p1$X,p1$Y) * 100, digits=1)
        probe_conflicts=rbind(
          probe_conflicts,
          data.frame(Probe_1st=i, Probe_2nd=j, PercentArea=areapoly,
                     ConflictType='target-guide')
        )
      }
    }
  }
  if (verbose) print('Identifying conflict with other guide probes.')
  #Now identify clashes between guide probes.
  if (length(gpos[,'x']) > 1) { #Only if there is more than one guide probe.
    for (i in c(1:(length(gpos[,'x'])))) {
      p1 <- data.frame(PID=rep(i, lpoly), POS=(1:glpoly), X=gprobes[,paste('X',as.character(i), sep='')], Y=gprobes[,paste('Y',as.character(i), sep='')])
      for (j in c(1:length(gpos[,'x']))) { 
        if (j!=i){
          p2 <- data.frame(PID=rep(j, lpoly), POS=(1:glpoly), X=gprobes[,paste('X',as.character(j), sep='')], Y=gprobes[,paste('Y',as.character(j), sep='')])
          if (!is.null(joinPolys(p1,p2))) {
            junction=joinPolys(p1,p2)
            #Compute the percentage area of a probe that overlaps.100% means they overlap perfectly. 0% means no overlap.
            areapoly=round(polyarea(junction$X,junction$Y)/polyarea(p1$X,p1$Y) * 100, digits=1)
            probe_conflicts=rbind(
              probe_conflicts,
              data.frame(Probe_1st=i, Probe_2nd=j, PercentArea=areapoly,
                         ConflictType='guide-guide')
            )
          }
        }
      }
    }
  }
  nconflicts=dim(probe_conflicts)[1]
  if (nconflicts>0) {
    colnames(probe_conflicts)=c('Probe_1st','Probe_2nd','PercentArea',
                                'ConflictType')
  }
  return(as.data.frame(probe_conflicts))
}

#Determine the appropriate cable exit gap for a position by picking
#the exit from the direction with the fewest positions in front. This is 
#done by dividing the field into three using the position and the central points
#between the cable exit gaps and counting the number of targets within.
choose_cegs <- function(pos) {
  
  npos=length(pos[,1])
  cegs=rep(0,npos)
  
  #Visualising what's going on (may be deleted when done with editing this bit of code):
  #x=pos[1,'x']
  #y=pos[1,'y']
  #draw_pos(x=pos_master$pos[c(1:19),'x'],y=pos_master$pos[c(1:19),'y'])
  
  for (j in c(1:npos)){
    #Define area or "pizza slice" for each cable exit gap.
    pizza1 = data.frame(PID=rep(1, 122), POS=1:122, X=c(pos[j,'x'],fov/2.*cosd(c(-30:-150))), Y=c(pos[j,'y'],fov/2.*sind(c(-30:-150))))
    pizza2 = data.frame(PID=rep(2, 122), POS=1:122, X=c(pos[j,'x'],fov/2.*cosd(c(-30:90))), Y=c(pos[j,'y'],fov/2.*sind(c(-30:90))))
    pizza3 = data.frame(PID=rep(3, 122), POS=1:122, X=c(pos[j,'x'],fov/2.*cosd(c(90:210))), Y=c(pos[j,'y'],fov/2.*sind(c(90:210))))
    nceg=c(-1,-1,-1) #Starting at -1 to avoid counting the probe itself.
    for (i in c(1:npos)){
      cx=pos[i,'x']
      cy=pos[i,'y']
      #Approximating each position as a tiny square for this purpose only.
      polpos= data.frame(PID=rep(1, 4), POS=1:4, X=c(cx+excl_radius,cx, cx-excl_radius, cx), Y=c(cy,cy+excl_radius,cy,cy-excl_radius))
      if (!is.null(joinPolys(pizza1,polpos))) {
        nceg[1]=nceg[1]+1
      }
      if(!is.null(joinPolys(pizza2,polpos))){
        nceg[2]=nceg[2]+1
      } 
      if (!is.null(joinPolys(pizza3,polpos))){
        nceg[3]=nceg[3]+1
      }
    }
    #Need to account for when there are two identically good options.
    if (length(which(min(nceg)==nceg))==1){
      cegs[j]=which(min(nceg)==nceg)
    } else { #If two exits are equally good, pick the closest.
      cegs[j]=which(min(nceg)==nceg)[which.min(pointDistance(p1=ceg_pos[which(min(nceg)==nceg),c('x','y')],p2=c(cx,cy), lonlat=FALSE))]
    }
    #If position is too close to the exit, move them to next closest exit.
    E=pointDistance(p1=ceg_pos[cegs[j],c('x','y')], p2=pos[j,c('x','y')], lonlat=FALSE) - (probe_l+cable_l+tip_l/2.) ##E is the line joining the solid ferule + cable length end to the centre of the exit point. 
    #Was increased on 24 August 2022 to avoid ferule being too close to exit gap.
    if (E < E0){ #E0 is the threshold to avoid blocking a cable exit gap. It is set above.
      cegs[j]=order(pointDistance(p1=ceg_pos[,c('x','y')],p2=pos[j,c('x','y')], lonlat=FALSE))[2]
    }
  }
  return(cegs)
  
}

#Initial guess of the angles. Should point to the appropriate exit.
first_guess_angs <- function(pos, cegs, guide=FALSE){ 
  
  nprobes=length(pos[,'x'])
  
  angs=rep(0,nprobes)
  
  for (probe in c(1:nprobes)){
    angs[probe]=atan2(ceg_pos[cegs[probe],'y'] - pos[probe,'y'], ceg_pos[cegs[probe],'x'] - pos[probe,'x']) + pi
    if (angs[probe] > 2*pi) angs[probe]=angs[probe]-2.*pi
  }
  
  #Determine probes near the edge of the field that cannot have the initial guess angle.
  if (guide){
    mindist=sqrt((gcable_l+gprobe_l+gtip_l-gtip_w/2.)^2 + (gcable_w/2.)^2)
  }else{
    mindist=sqrt((cable_l+probe_l+tip_l-tip_w/2.)^2 + (cable_w/2.)^2)
  }
  
  #Computing the allowed range... this has changed in v2.0
  #for (probe in c(1:nprobes)){
    #dist=(fov/2.)-sqrt(pos[probe,'y']**2+pos[probe,'x']**2)
    #if (dist < mindist){
    #  arange=computeAnglesRange(x=pos[probe,'x'],y=pos[probe,'y'], ceg=cegs[probe], guide=guide)
    #  if (angs[probe]/pi*180 < min(arange)) angs[probe]=min(arange/180*pi)
    #  else if ((angs[probe]/pi*180 > max(arange))) angs[probe]=max(arange/180*pi)
    #  # else do nothing, you're in a cable gap exit, so just keep pointing to the middle.
    #}
  #}
  
  return(angs)
}

MinimizeConflicts <- function(pos=pos, angs=angs, cegs=cegs, delta_ang=20, probe_conflicts=conflicts){
  
  #Find all conflicted probes:
  CPs=unique(c(probe_conflicts[,1],probe_conflicts[,2]))
  #Identify all probes within reach of the conflicted probes as a list of probes to tweak.
  TPs= CPs #List of tweakable probes:
  for (probe in CPs) {
    dists=sqrt( (pos[,'x']- pos[probe,'x'])^2 + (pos[,'y']- pos[probe,'y'])^2)
    newTPs=which(dists < (probe_l/3))
    TPs=c(TPs,newTPs)
  }
  TPs=unique(TPs)
  
  ang_ranges=cbind(angs,angs)
  for (probe in TPs){
    angRange=computeAnglesRange(x=pos[probe,'x'],y=pos[probe,'y'], ceg=cegs[probe], delta_ang=delta_ang)/180*pi
    ang_ranges[probe,]=angRange
  }
  
  draw_all(x=pos[,'x'],y=pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  
  nit=length(TPs)*500
  angs_tries=c()
  for (probe in c(1:nprobes)){
    angs_tries= rbind(angs_tries,runif(n=nit,min=ang_ranges[probe,1],max=ang_ranges[probe,2]))
  }
  
  angs_tries=as.data.frame(angs_tries)
  #View(angs_tries)
  angs_tries[nprobes+1,]=1000
  areamin=1000
  for (i in c(1:nit)){
    angs_tries[nprobes+1,i]=compute_total_conflict_area(angs=angs_tries[c(1:nprobes),i],pos=pos)
    if (angs_tries[nprobes+1,i] < areamin){
      areamin=angs_tries[nprobes+1,i]
      angsmin=angs_tries[c(1:nprobes),i]
      #print(areamin)
      angsmin[angsmin>3*pi/2]=angsmin[angsmin>3*pi/2]-2*pi
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angsmin, interactive=FALSE)
      if (areamin<=0) {return(angsmin)}
    }
  }
  angsmin[angsmin>3*pi/2]=angsmin[angsmin>3*pi/2]-2*pi
  angsmin <<- angsmin
  
  draw_all(x=pos[,'x'],y=pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angsmin, interactive=visualise)
  
  #Identify remaining conflicts.
  conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angsmin)
  nconflicts=dim(conflicts)[1]
  
  #Try expanding the angles of probes still conflicted.
  if (nconflicts>0){
    #Recurrently run but first fix probes that are now no longer conflicted...
    angsmin <<- MinimizeConflicts(pos=pos, cegs=cegs, delta_ang=delta_ang, angs=angsmin, probe_conflicts=conflicts)
  }
  
  draw_all(x=pos[,'x'],y=pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angsmin, interactive=visualise)
  
  angs=angsmin
  return(angs)
}

wrapMinimizeConflicts <- function(pos=pos, angs=angs, cegs=cegs, init_delta_ang=20, probe_conflicts=conflicts){
  print('Wrapping around MinimizeConflicts')
  conflicts=probe_conflicts
  nconflicts=dim(conflicts)[1]
  
  #Look for a solution.
  i=1 
  while ((nconflicts > 0) & (i<4)){ #We will wrap around this loop at most 3 times.
    #Start by allowing probes to move by the initial delta_ang.
    
    delta_ang=init_delta_ang
    print(paste('Starting (again) with delta_ang =', as.character(delta_ang)))
    
    newangs=angs 
    try(withTimeout(assign('newangs',MinimizeConflicts(pos=pos, angs=angs, cegs=cegs, delta_ang=delta_ang, probe_conflicts=conflicts)), timeout=timeout, onTimeout=timeout_behaviour), silent=TRUE)
    if (FALSE %in% (angs==newangs)){
      angs=newangs
    } else {
      angs=angsmin
    }
    draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
    probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
    #Check if could configure:
    conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
    nconflicts=dim(conflicts)[1]

    #Trying swapping the cable exit gap first as this often resolves conflicts and can speed up the process..
    if (nconflicts > 0) {
      #Trying to swap some cable exit gaps first
      nceg=try_next_ceg(pos=pos, angs=angsmin, cegs=cegs, probe_conflicts=conflicts)
      cegs=nceg$cegs
      angs=nceg$angs
      newangs=angs 
      
      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
      if (nconflicts >0){
        try(withTimeout(assign('newangs',MinimizeConflicts(pos=pos, angs=angs, cegs=cegs, delta_ang=delta_ang, probe_conflicts=conflicts)), timeout=timeout, onTimeout=timeout_behaviour))
      }
      if (FALSE %in% (angs==newangs)){
        angs=newangs
      } else {
        angs=angsmin
      }
      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      #Check if could configure:
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
    }
    
    if (nconflicts > 0) {
      print('Could not find a solution with delta_ang=20, trying delta_ang=45')
      delta_ang=45
      newangs=angs
      try(withTimeout(assign('newangs',MinimizeConflicts(pos=pos, angs=angs, cegs=cegs, delta_ang=delta_ang, probe_conflicts=conflicts)), timeout=timeout, onTimeout=timeout_behaviour))
      if (FALSE %in% (angs==newangs)){
        angs=newangs
      } else {
        angs=angsmin #If it timed out, use global angsmin.
      }
      
      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      #Check if could configure:
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
    }
    
    if (nconflicts > 0) {
      print('Could not find a solution with delta_ang=45, trying delta_ang=90')
      delta_ang=90
      newangs=angs
      try(withTimeout(assign('newangs',MinimizeConflicts(pos=pos, angs=angs, cegs=cegs, delta_ang=delta_ang, probe_conflicts=conflicts)), timeout=timeout, onTimeout=timeout_behaviour))
      if (FALSE %in% (angs==newangs)){
        angs=newangs
      } else {
        angs=angsmin
      }

      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      #Check if could configure:
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
    }
    
    if (nconflicts > 0) {
      print('Could not find a solution with delta_ang=90, trying to swap some cable exit gaps and delta_ang=20')
      nceg=try_next_ceg(pos=pos, angs=angsmin, cegs=cegs, probe_conflicts=conflicts)
      cegs=nceg$cegs
      angs=nceg$angs
      
      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
    }
    
    i=i+1
  }
  return(angs)
}

MinimizeGuideConflicts <- function(pos=pos, angs=angs, cegs=cegs, gpos=gpos, gangs=gangs, gcegs=gcegs, delta_ang=20, guide_conflicts){
  
  # Find conflicted guide probes.  For target-guide conflicts only the guide
  # angle is free here.  For guide-guide conflicts both guide IDs must be
  # included, so both probes can move in the stochastic search.
  if ('ConflictType' %in% names(guide_conflicts)) {
    guide_guide=guide_conflicts$ConflictType == 'guide-guide'
    CGPs=unique(c(
      guide_conflicts$Probe_2nd[!guide_guide],
      guide_conflicts$Probe_1st[guide_guide],
      guide_conflicts$Probe_2nd[guide_guide]
    ))
  } else {
    # Backward-compatible fallback for conflict tables created by older
    # versions of the engine.
    CGPs=unique(guide_conflicts[,2])
  }
  CGPs=as.integer(CGPs)
  CGPs=CGPs[is.finite(CGPs) & CGPs >= 1 & CGPs <= nrow(gpos)]
  TPs=integer(0)
  if ('ConflictType' %in% names(guide_conflicts)) {
    target_guide=guide_conflicts$ConflictType == 'target-guide'
    TPs=unique(as.integer(guide_conflicts$Probe_1st[target_guide]))
    TPs=TPs[is.finite(TPs) & TPs >= 1 & TPs <= nrow(pos)]
  }
  print('Minimizing conflicts with guide probe(s):')
  print(as.character(CGPs))
  if (length(TPs) > 0) {
    print('Also considering target/standard probe(s):')
    print(as.character(TPs))
  }

  # Return the representation of a 180-degree alternative that lies within
  # the physically allowed angular range.  Angles are deliberately not
  # wrapped to [-pi, pi], because computeAnglesRange can use an equivalent
  # representation outside that interval.
  flipped_angle <- function(angle, x, y, ceg, guide=FALSE) {
    angle_range=computeAnglesRange(
      x=x, y=y, ceg=ceg, delta_ang=180, guide=guide
    ) / 180*pi
    candidates=angle + pi + ((-2:2) * 2*pi)
    valid=candidates[
      is.finite(candidates) &
      all(is.finite(angle_range)) &
      candidates >= min(angle_range) &
      candidates <= max(angle_range)
    ]
    if (length(valid) > 0) valid[1] else NA_real_
  }

  # Search locally around a proposed configuration.  In particular, this is
  # used after a 180-degree flip: the exact flip is often still marginally in
  # conflict, while a few degrees of additional rotation resolves it.
  search_guide_neighbourhood <- function(base_gangs, deltas=c(delta_ang,45,90)) {
    probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                       interactive=FALSE)
    for (try_delta in unique(deltas)) {
      gang_ranges=cbind(base_gangs,base_gangs)
      valid_ranges=TRUE
      for (i in CGPs){
        allowed_range=computeAnglesRange(
          x=gpos[i,'x'], y=gpos[i,'y'], ceg=gcegs[i],
          delta_ang=180, guide=TRUE
        ) / 180*pi
        local_range=c(base_gangs[i] - try_delta/180*pi,
                      base_gangs[i] + try_delta/180*pi)
        if (!all(is.finite(allowed_range))) {
          valid_ranges=FALSE
          break
        }
        gang_ranges[i,1]=max(min(allowed_range), local_range[1])
        gang_ranges[i,2]=min(max(allowed_range), local_range[2])
        if (gang_ranges[i,1] > gang_ranges[i,2]) {
          valid_ranges=FALSE
          break
        }
      }
      if (!valid_ranges) next

      nit=max(100, length(CGPs)*50)
      gangs_tries=matrix(NA_real_, nrow=length(base_gangs), ncol=nit)
      gangs_tries[,1]=base_gangs
      for (probe in seq_along(base_gangs)) {
        if (probe %in% CGPs && nit > 1) {
          if (gang_ranges[probe,1] == gang_ranges[probe,2]) {
            gangs_tries[probe,2:nit]=gang_ranges[probe,1]
          } else {
            gangs_tries[probe,2:nit]=runif(
              n=nit-1, min=gang_ranges[probe,1], max=gang_ranges[probe,2]
            )
          }
        } else {
          gangs_tries[probe,]=base_gangs[probe]
        }
      }

      for (i in seq_len(nit)) {
        tempgangs=gangs_tries[,i]
        gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],
                                 gangs=tempgangs, interactive=FALSE)
        fov_status=configuration_fov_status(probes, gprobes)
        if (!fov_status$ok) next
        tempconflicts=find_guide_conflicts(
          probes=probes, pos=pos, angs=angs, gprobes=gprobes,
          gpos=gpos, gangs=tempgangs, verbose=FALSE
        )
        if (nrow(tempconflicts) == 0) return(tempgangs)
      }
    }
    return(NULL)
  }

  # If a guide still clashes with a target after all guide-only searches,
  # allow the implicated target and guide(s) to move together.  Target-target
  # conflicts are checked for every trial, so this cannot silently damage the
  # already configured target field.
  search_joint_neighbourhood <- function(base_angs, base_gangs,
                                         deltas=c(10,delta_ang,45),
                                         samples_per_probe=50) {
    if (length(TPs) == 0) return(NULL)

    for (try_delta in unique(deltas)) {
      ang_ranges=cbind(base_angs,base_angs)
      gang_ranges=cbind(base_gangs,base_gangs)
      valid_ranges=TRUE

      for (i in TPs){
        allowed_range=computeAnglesRange(
          x=pos[i,'x'], y=pos[i,'y'], ceg=cegs[i],
          delta_ang=180, guide=FALSE
        ) / 180*pi
        local_range=c(base_angs[i] - try_delta/180*pi,
                      base_angs[i] + try_delta/180*pi)
        if (!all(is.finite(allowed_range))) {
          valid_ranges=FALSE
          break
        }
        ang_ranges[i,1]=max(min(allowed_range), local_range[1])
        ang_ranges[i,2]=min(max(allowed_range), local_range[2])
        if (ang_ranges[i,1] > ang_ranges[i,2]) {
          valid_ranges=FALSE
          break
        }
      }
      if (!valid_ranges) next

      for (i in CGPs){
        allowed_range=computeAnglesRange(
          x=gpos[i,'x'], y=gpos[i,'y'], ceg=gcegs[i],
          delta_ang=180, guide=TRUE
        ) / 180*pi
        local_range=c(base_gangs[i] - try_delta/180*pi,
                      base_gangs[i] + try_delta/180*pi)
        if (!all(is.finite(allowed_range))) {
          valid_ranges=FALSE
          break
        }
        gang_ranges[i,1]=max(min(allowed_range), local_range[1])
        gang_ranges[i,2]=min(max(allowed_range), local_range[2])
        if (gang_ranges[i,1] > gang_ranges[i,2]) {
          valid_ranges=FALSE
          break
        }
      }
      if (!valid_ranges) next

      # This is deliberately a small search: normally only one target and
      # one guide are being varied, so a few dozen draws per probe are enough
      # to explore each local neighbourhood without delaying the field run.
      nit=max(100, (length(TPs)+length(CGPs))*samples_per_probe)
      angs_tries=matrix(NA_real_, nrow=length(base_angs), ncol=nit)
      gangs_tries=matrix(NA_real_, nrow=length(base_gangs), ncol=nit)
      angs_tries[,1]=base_angs
      gangs_tries[,1]=base_gangs

      for (probe in seq_along(base_angs)) {
        if (probe %in% TPs && nit > 1) {
          if (ang_ranges[probe,1] == ang_ranges[probe,2]) {
            angs_tries[probe,2:nit]=ang_ranges[probe,1]
          } else {
            angs_tries[probe,2:nit]=runif(
              n=nit-1, min=ang_ranges[probe,1], max=ang_ranges[probe,2]
            )
          }
        } else {
          angs_tries[probe,]=base_angs[probe]
        }
      }
      for (probe in seq_along(base_gangs)) {
        if (probe %in% CGPs && nit > 1) {
          if (gang_ranges[probe,1] == gang_ranges[probe,2]) {
            gangs_tries[probe,2:nit]=gang_ranges[probe,1]
          } else {
            gangs_tries[probe,2:nit]=runif(
              n=nit-1, min=gang_ranges[probe,1], max=gang_ranges[probe,2]
            )
          }
        } else {
          gangs_tries[probe,]=base_gangs[probe]
        }
      }

      for (i in seq_len(nit)) {
        tempangs=angs_tries[,i]
        tempgangs=gangs_tries[,i]
        temp_probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=tempangs,
                                interactive=FALSE)
        gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],
                                 gangs=tempgangs, interactive=FALSE)
        fov_status=configuration_fov_status(temp_probes, gprobes)
        if (!fov_status$ok) next
        target_conflicts=find_probe_conflicts(
          probes=temp_probes, pos=pos, angs=tempangs
        )
        if (nrow(target_conflicts) > 0) next

        tempconflicts=find_guide_conflicts(
          probes=temp_probes, pos=pos, angs=tempangs,
          gprobes=gprobes, gpos=gpos, gangs=tempgangs, verbose=FALSE
        )
        if (nrow(tempconflicts) == 0) {
          return(list(angs=tempangs, gangs=tempgangs))
        }
      }
    }
    return(NULL)
  }

  # Try 180-degree alternatives before the broader guide-only search.  For
  # the usual one-target/one-guide clash this means exactly three starts:
  # flip the target, flip the guide, or flip both.  If more probes are listed
  # in the conflict table, single and pair flips are sufficient here; trying
  # every higher-order combination is disproportionately expensive.
  search_joint_flip_neighbourhood <- function(
      deltas=c(10,delta_ang,45), samples_per_probe=50) {
    if (length(TPs) == 0) return(NULL)
    movable=c(TPs,CGPs)
    nmovable=length(movable)
    if (nmovable == 0) return(NULL)

    if (nmovable <= 4) {
      mask_ids=seq_len(2^nmovable)-1L
    } else {
      mask_ids=integer(0)
      for (i in seq_len(nmovable)) {
        mask_ids=c(mask_ids, 2^(i-1))
      }
      if (nmovable > 1) {
        for (i in seq_len(nmovable-1)) {
          for (j in (i+1):nmovable) {
            mask_ids=c(mask_ids, 2^(i-1)+2^(j-1))
          }
        }
      }
    }

    for (mask_id in mask_ids) {
      tempangs=angs
      tempgangs=gangs
      flip_bits=as.logical(intToBits(as.integer(mask_id)))[seq_along(movable)]
      valid_flip=TRUE

      for (j in seq_along(TPs)) {
        if (flip_bits[j]) {
          new_angle=flipped_angle(
            angle=angs[TPs[j]], x=pos[TPs[j],'x'], y=pos[TPs[j],'y'],
            ceg=cegs[TPs[j]], guide=FALSE
          )
          if (!is.finite(new_angle)) {
            valid_flip=FALSE
            break
          }
          tempangs[TPs[j]]=new_angle
        }
      }
      if (!valid_flip) next

      for (j in seq_along(CGPs)) {
        bit=j+length(TPs)
        if (flip_bits[bit]) {
          new_angle=flipped_angle(
            angle=gangs[CGPs[j]], x=gpos[CGPs[j],'x'],
            y=gpos[CGPs[j],'y'], ceg=gcegs[CGPs[j]], guide=TRUE
          )
          if (!is.finite(new_angle)) {
            valid_flip=FALSE
            break
          }
          tempgangs[CGPs[j]]=new_angle
        }
      }
      if (!valid_flip) next

      joint_solution=search_joint_neighbourhood(
        tempangs, tempgangs, deltas=deltas,
        samples_per_probe=samples_per_probe
      )
      if (!is.null(joint_solution)) return(joint_solution)
    }
    return(NULL)
  }

  # Try the joint target/guide search from the current configuration first.
  # This avoids spending the full guide-only search budget when the target
  # itself needs to move to make room.
  if (length(TPs) > 0) {
    joint_solution=search_joint_neighbourhood(angs,gangs)
    if (!is.null(joint_solution)) {
      print('Resolved guide-target conflict with an early joint target/guide angular search.')
      return(joint_solution)
    }
  }

  joint_solution=search_joint_flip_neighbourhood()
  if (!is.null(joint_solution)) {
    print('Resolved guide-target conflict after an early 180-degree flip and local angular search.')
    return(joint_solution)
  }

  # Try the discrete 180-degree alternatives before stochastic optimisation.
  # These are often the useful solutions when two guide cables are crossing:
  # flipping either probe, or both, can resolve the conflict immediately.
  if (length(CGPs) <= 10) {
    nflip_combinations=2^length(CGPs)
    for (mask_id in (seq_len(nflip_combinations)-1L)) {
      tempgangs=gangs
      flip_bits=as.logical(intToBits(as.integer(mask_id)))[seq_along(CGPs)]
      if (any(flip_bits)) {
        for (j in which(flip_bits)) {
          angle_range=computeAnglesRange(
            x=gpos[CGPs[j],'x'], y=gpos[CGPs[j],'y'],
            ceg=gcegs[CGPs[j]], delta_ang=180, guide=TRUE
          ) / 180*pi
          # Keep the angle representation used by computeAnglesRange.  The
          # same physical direction can be written with several 2*pi-shifted
          # values, and the range may legitimately extend beyond [-pi, pi].
          # Choose the equivalent 180-degree flip that lies in that range.
          flip_candidates=tempgangs[CGPs[j]] + pi + ((-2:2) * 2*pi)
          valid_flip=flip_candidates[
            is.finite(flip_candidates) &
            all(is.finite(angle_range)) &
            flip_candidates >= min(angle_range) &
            flip_candidates <= max(angle_range)
          ]
          if (length(valid_flip) > 0) {
            tempgangs[CGPs[j]]=valid_flip[1]
          } else {
            flip_bits[j]=FALSE
          }
        }
      }
      if (any(flip_bits)) {
        probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                           interactive=FALSE)
        gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],
                                 gangs=tempgangs, interactive=FALSE)
        fov_status=configuration_fov_status(probes, gprobes)
        tempconflicts=find_guide_conflicts(
          probes=probes, pos=pos, angs=angs, gprobes=gprobes,
          gpos=gpos, gangs=tempgangs, verbose=FALSE
        )
        if (nrow(tempconflicts) == 0 && fov_status$ok) {
          print('Resolved guide conflict with a 180-degree probe flip.')
          return(tempgangs)
        }
        local_solution=search_guide_neighbourhood(tempgangs)
        if (!is.null(local_solution)) {
          print('Resolved guide conflict after a 180-degree probe flip and local angular search.')
          return(local_solution)
        }
      }
    }
  }
  
  # The old routine sampled only 50 angles in a narrow range. Try the
  # requested range first, then wider ranges, with enough samples that a
  # narrow conflict-free interval is not missed merely by chance.
  delta_angles=unique(c(delta_ang,45,90))
  for (try_delta in delta_angles) {
    gang_ranges=cbind(gangs,gangs)
    for (i in CGPs){
      gang_ranges[i,] = computeAnglesRange(x=gpos[i,'x'], y=gpos[i,'y'],
                                           ceg=gcegs[i],
                                           delta_ang=try_delta,
                                           guide=TRUE) / 180*pi
    }

    draw_all(x=pos[,'x'],y=pos[,'y'],xs=pos_master$stdpos[,'x'],
             ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],
             yg=pos_master$guidepos[,'y'], interactive=visualise)
    probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                       interactive=visualise)
    gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],
                             gangs=gangs, interactive=visualise)

    nit=max(100, length(CGPs)*50)
    gangs_tries=c()
    for (probe in c(1:length(gangs))){
      if (probe %in% CGPs){
        gangs_tries = rbind(gangs_tries,
                            runif(n=nit,min=gang_ranges[probe,1],
                                  max=gang_ranges[probe,2]))
      }else{
        gangs_tries = rbind(gangs_tries,rep(gangs[probe],times=nit))
      }
    }
    gangs_tries=as.data.frame(gangs_tries)

    for (i in c(1:nit)){
      tempgangs=gangs_tries[,i]

      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                         interactive=FALSE)
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],
                               gangs=tempgangs, interactive=FALSE)
      fov_status=configuration_fov_status(probes, gprobes)
      if (!fov_status$ok) next

      tempconflicts=find_guide_conflicts(probes=probes, pos=pos,
                                         angs=angs,gprobes=gprobes,
                                         gpos=gpos, gangs=tempgangs,
                                         verbose=FALSE)
      nconflicts=dim(tempconflicts)[1]
      if (nconflicts <=0){
        return(tempgangs)
      }
    }
  }

  # Last chance for genuinely difficult fields.  This larger search is only
  # reached after the fast current-state and flip-start searches have failed.
  # It is still restricted to the implicated target/guide probes.
  if (length(TPs) > 0) {
    print('Trying extended joint target/guide search for this hard conflict.')
    joint_solution=search_joint_neighbourhood(
      angs, gangs, deltas=c(5,10,20,45,90), samples_per_probe=250
    )
    if (is.null(joint_solution)) {
      joint_solution=search_joint_flip_neighbourhood(
        deltas=c(5,10,20,45,90), samples_per_probe=250
      )
    }
    if (!is.null(joint_solution)) {
      print('Resolved guide-target conflict with the extended joint search.')
      return(joint_solution)
    }
  }

  print('Could not configure this guide probe. Will try the next one in the list.')
  
  return(gangs)
}

compute_total_conflict_area <- function(angs,pos=pos){
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=FALSE)
  lpoly=length(probes[,1]) #length of polygon defining probe.
  sprobes=data.frame()
  for (i in c(1:(length(pos[,'x'])))) {
    sprobes <- rbind(sprobes,cbind(PID=rep(i, lpoly), POS=1:lpoly, X=probes[,paste('X',as.character(i), sep='')], Y=probes[,paste('Y',as.character(i), sep='')]))
  }
  sprobes=as.PolySet(sprobes)
  junction=joinPolys(sprobes,sprobes,operation='INT')
  nconflicts=(length(unique(junction[,'PID']))-nprobes)/2
  
  area=c()
  for (i in unique(junction[,'PID'])){
    area=c(area,polyarea(junction[junction[,'PID']==i,'X'],junction[junction[,'PID']==i,'Y']))
  }
  total_area=sum(area)-nprobes*Mode(area)
  
  return(total_area)
}

#Removing standard stars that are conflicting with a probe within the exclusion radius.
cleanconflictedstds <- function(angs,pos=pos, stdpos=pos_master$stdpos){
  
  draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=stdpos[,'x'],ys=stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  lpoly=length(probes[,1]) #length of polygon defining probe.
  
  #Setting up into sp objects that package rgeos will understand.
  sprobes=data.frame()
  for (i in c(1:ngalprobes)) {
    sprobes <- rbind(sprobes,cbind(PID=rep(i, lpoly), POS=1:lpoly, X=probes[,paste('X',as.character(i), sep='')], Y=probes[,paste('Y',as.character(i), sep='')]))
  }
  sprobes=as.PolySet(sprobes)
  for (j in c(1:ngalprobes)) {
    assign(paste('cprobe',as.character(j),sep=''),Polygon(coords=as.matrix(sprobes[sprobes[,'PID']==j,c('X','Y')]), hole=FALSE))
  }
  cprobes=Polygons(srl=mget(paste('cprobe',as.character(c(1:ngalprobes)),sep='')), ID='meh')
  pols=SpatialPolygons(Srl=list(cprobes)) #sp object for the probes as spatial polygons.
  spts = SpatialPoints(coords=as.matrix(pos_master$stdpos[,c('x','y')])) #sp objects for the standard stars as simple points.
  #Following line replaced Sept 2026 by CFoster as used deprecated rgeos.
  #mindists=apply(gDistance(spts, pols, byid=TRUE),2,min) #The minimum distance to the nearest probe for each standard. 
  d = st_distance(st_as_sf(spts), st_as_sf(pols))   # rows = standards, cols = polygons
  mindists = apply(d, 1, min)                       # min distance to nearest probe, per standard
  
  clstdpos=stdpos[order(mindists, decreasing=TRUE),] #Ordering standards according to their distance from the closest probe.
  mindists=mindists[order(mindists, decreasing=TRUE)]
  clstdpos=clstdpos[mindists > excl_radius,] #Removing standards that overlap with nearest probe.
  
  #Reordering the clean standard stars by distance from the centre (closest at the top):
  dist2cen=sqrt(clstdpos[,'x']^2+clstdpos[,'y']^2)
  #removing standard probes that are very close to the edge of the field.
  clstdpos=clstdpos[which(dist2cen < (fov/2. - 1.5* excl_radius)),]
  clstdpos=clstdpos[order(dist2cen, decreasing=FALSE),]
  print(dist2cen[order(dist2cen, decreasing=FALSE)])
  
  return(clstdpos)
}

#Removing guide stars that are conflicting with a probe within the exclusion radius.
cleanconflictedguides <- function(angs,pos=pos, guidepos=pos_master$guidepos){
  
  draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  lpoly=length(probes[,1]) #length of polygon defining probe.
  
  #Setting up into sp objects that package rgeos will understand.
  sprobes=data.frame()
  for (i in c(1:(ngalprobes+nstdprobes))) {
    sprobes <- rbind(sprobes,cbind(PID=rep(i, lpoly), POS=1:lpoly, X=probes[,paste('X',as.character(i), sep='')], Y=probes[,paste('Y',as.character(i), sep='')]))
  }
  sprobes=as.PolySet(sprobes)
  for (j in c(1:(ngalprobes+nstdprobes))) {
    assign(paste('cprobe',as.character(j),sep=''),Polygon(coords=as.matrix(sprobes[sprobes[,'PID']==j,c('X','Y')]), hole=FALSE))
  }
  cprobes=Polygons(srl=mget(paste('cprobe',as.character(c(1:(ngalprobes+nstdprobes))),sep='')), ID='meh')
  pols=SpatialPolygons(Srl=list(cprobes)) #sp object for the probes as spatial polygons.
  spts = SpatialPoints(coords=as.matrix(pos_master$guidepos[,c('x','y')])) #sp objects for the guide stars as simple points.
  #Following line replaced Sept 2026 by CFoster as used deprecated rgeos.
  #mindists=apply(gDistance(spts, pols, byid=TRUE),2,min) #The minimum distance to the nearest probe for each guide.
  d = st_distance(st_as_sf(spts), st_as_sf(pols))   # rows = guides, cols = polygons
  mindists = apply(d, 1, min)
  
  clguidepos=guidepos[order(mindists, decreasing=TRUE),] #Ordering guides according to their distance from the closest probe.
  mindists=mindists[order(mindists, decreasing=TRUE)]
  clguidepos=clguidepos[mindists > gexcl_radius,] #Removing guides that overlap with nearest probe.
  
  #Reordering guides by distance from the centre:
  dist2cen=sqrt(clguidepos[,'x']^2+clguidepos[,'y']^2)
  #removing guide probes that are very close to the edge of the field.
  clguidepos=clguidepos[which(dist2cen < (fov/2. - 1.5* excl_radius)),]
  clguidepos=clguidepos[order(dist2cen, decreasing=TRUE),]
  
  #Dividing into quadrants:
  first_q=clguidepos[clguidepos[,'x']<0 & clguidepos[,'y']<0,]
  #print(first_q)
  second_q=clguidepos[clguidepos[,'x']<0 & clguidepos[,'y']>=0,]
  #print(second_q)
  third_q=clguidepos[clguidepos[,'x']>=0 & clguidepos[,'y']>=0,]
  #print(third_q)
  fourth_q=clguidepos[clguidepos[,'x']>=0 & clguidepos[,'y']<0,]
  #print(fourth_q)
  
  oclguidepos=data.frame(x=double(),y=double()) #ordered and cleaned guide positions.
  #Cycling through by quadrant to ensure guides from all corners of the plate get prioritised:
  gui=1
  print(gui)
  
  
  remaining_first_q = nrow(as.matrix(first_q))
  remaining_second_q = nrow(as.matrix(second_q))
  remaining_third_q = nrow(as.matrix(third_q))
  remaining_fourth_q = nrow(as.matrix(fourth_q))
                            
  while ((length(oclguidepos[,1]) < length(clguidepos[,1])) & gui < length(clguidepos[,1])){
    
    if (remaining_first_q > 0) {
      if (nrow(as.matrix(first_q))>1){
        assign('oclguidepos',rbind(oclguidepos, first_q[gui,]))
        remaining_first_q = remaining_first_q -1
    } else{
      assign('oclguidepos',rbind(oclguidepos, first_q[1, , drop=FALSE]))
      remaining_first_q = remaining_first_q -1
    }
      
    }
    if (remaining_second_q > 0){
      
    if (nrow(as.matrix(second_q))>1 & remaining_second_q > 0){
      assign('oclguidepos',rbind(oclguidepos, second_q[gui,]))
      remaining_second_q = remaining_second_q -1
    } else {
      assign('oclguidepos',rbind(oclguidepos, second_q[1, , drop=FALSE]))
      remaining_second_q = remaining_second_q - 1
    }
    }
    
    if (remaining_third_q > 0){
    if (nrow(as.matrix(third_q))>1){
      assign('oclguidepos',rbind(oclguidepos, third_q[gui,]))
      remaining_third_q = remaining_third_q -1
    } else {
      assign('oclguidepos',rbind(oclguidepos, third_q[1, , drop=FALSE]))
      remaining_third_q = remaining_third_q -1
    }
    }
    if (remaining_fourth_q > 0){
    if (nrow(as.matrix(fourth_q))>1){
      assign('oclguidepos',rbind(oclguidepos, fourth_q[gui,]))
      remaining_fourth_q = remaining_fourth_q -1
    } else{
      assign('oclguidepos',rbind(oclguidepos, fourth_q[1, , drop=FALSE]))
      remaining_fourth_q = remaining_fourth_q -1
    }
    }
    gui=gui+1
    #print(c('oclguidepos is length', length(oclguidepos[,1])))
  }
  #print(c("OCLguidepos has ", nrow(oclguidepos), "rows"))
  colnames(oclguidepos)=c('x','y')
  
  return(oclguidepos)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

probe_density <- function(pos=pos) {
  nneighbours=rep(0,ngalprobes)
  #Counting the number of neighbours within the total probe length:
  for (probe in c(1:ngalprobes)){
    dists=sqrt( (pos[,'x']- pos[probe,'x'])^2 + (pos[,'y']- pos[probe,'y'])^2)
    nneighbours[probe]=length(which(dists < 3.*excl_radius))
  }
  return(nneighbours)
}

try_next_ceg <- function(pos, cegs, angs, probe_conflicts=conflicts){
  
  print('Trying to swap the cable exit gap for some probes to see if that helps.')
  #Identify clusters of conflicted probes.
  nconflicts=length(probe_conflicts[,1])
  for (i in c(1:nconflicts)) {
    if (i==1){
      CP_clusters=list(as.numeric(probe_conflicts[i,c(1:2)]))
    }else{
      found=FALSE
      ncl=length(CP_clusters)
      for (cl in c(1:ncl)){
        cCPs=as.numeric(probe_conflicts[i,c(1:2)])
        #If any of the conflicted probes in the conflict is in existing cluster, add them to the cluster.
        if ((cCPs[1] %in% CP_clusters[[cl]]) | (cCPs[2] %in% CP_clusters[[cl]])){
          CP_clusters[[cl]]=unique(c(as.numeric(probe_conflicts[i,c(1:2)]), CP_clusters[[cl]]))
          found=TRUE
        }
      }
      if (!found) { # else, create a new cluster.
        CP_clusters[[ncl+1]]=as.numeric(probe_conflicts[i,c(1:2)])
      }
    }
  }
  
  #Identify clusters of conflicted probes pointing to the same exit and change the ceg for the furthest probe
  # if allowed, otherwise try second furthest probe. If probes pointing at different exit gaps, try swap one.
  cegs_new=cegs
  angs_new=angs
  
  ncl=length(CP_clusters)
  for (cl in c(1:ncl)){
    cCPs=CP_clusters[[cl]]
    if (length(cCPs)-length(unique(cegs[cCPs])) > 0) { #if more than 2 in the cluster point to the same exit gap
      cceg=Mode(cegs[cCPs])
      cCPs=cCPs[cegs[cCPs]==cceg] #focusing on probes in common cable exit gap.
      
      d=pointDistance(p1=pos[cCPs,c('x','y')],p2=ceg_pos[cegs[cCPs][1],c('x','y')], lonlat=FALSE)
      
      #identify probe furthest from the common cable exit gap to swap to second closest cable exit gap:
      toswap=cCPs[which(d==max(d))][1]
    } else { #If all conflicted probes in cluster pointing to different exit gaps, move the probe with the most limited ang_ranges
      angRanges=c()
      for (i in cCPs){
        angRanges=c(angRanges, abs(diff(computeAnglesRange(x=pos[i,'x'],y=pos[i,'y'],ceg=cegs[i]))))
      }
      toswap=cCPs[which(angRanges==min(angRanges))][1]
    }
    d2cegs=pointDistance(p1=ceg_pos[,c('x','y')], p2=pos[toswap,c('x','y')], lonlat=FALSE)
    if ((cegs_new[toswap] == order(d2cegs)[2]) & (d2cegs[order(d2cegs)[1]] > (probe_l+tip_l+cable_l))) {
      cegs_new[toswap]=order(d2cegs)[1]
    } else{
      cegs_new[toswap]=order(d2cegs)[2]
    }
    #angs_new[toswap]=first_guess_angs(pos=pos,cegs=cegs_new)[toswap]
    #Better to reset all other angles.
    angs_new=first_guess_angs(pos=pos,cegs=cegs_new)
  }
  angs_new[angs_new>3*pi/2]=angs_new[angs_new>3*pi/2]-2*pi
  return(list(cegs=cegs_new, angs=angs_new))
}

try_swaps <- function(pos_master=pos_master){
  if (length(pos_master$pos[,1])<=dngalprobes){
    print('WARNING: There is no spare target to swap in this field.')
    #* SAM- added in pos= and angs= here, as otherwise we can't read the output list and fail with a really obscure error
    #* I have no idea why temp$pos gives NA rather than an error if temp has no attribute pos. REALLY HARD TO DEBUG!!!
    return(list(pos=pos_master$pos, angs=first_guess_angs(pos=pos_master$pos, cegs=choose_cegs(pos_master$pos))))
  }
  #Computing the probe density.
  nspare=length(pos_master$pos[,1])-ngalprobes
  pos=pos_master$pos[1:ngalprobes,]
  pd=probe_density(pos=pos)
  if (nspare<=3) {
    pos2swap=order(pd,decreasing=TRUE)[1:3]
  } else{
    pos2swap=order(pd,decreasing=TRUE)[1:nspare]
  }
  for (p in pos2swap) {
    print(paste('Swapping probe ', as.character(p)))
    for (i in c(1:nspare)){
      print(paste('with spare ',as.character(i)))
      cpos=pos
      cpos[p,]=pos_master$pos[ngalprobes+i,]
      cegs=choose_cegs(pos=cpos)
      angs=first_guess_angs(pos=cpos, cegs=cegs)
      probes=defineprobe(x=cpos[,'x'],y=cpos[,'y'],angs=angs, interactive=FALSE)
      
      print('Identifying initial conflicts')
      conflicts=find_probe_conflicts(probes=probes, pos=cpos, angs=angs)
      nconflicts=dim(conflicts)[1]
      if (nconflicts <1) {
        print('Bam! No initial conflicts, this one is done!')
        return(list(pos=cpos,angs=angs))
      }
      print('Minimizing the overlaping area by rotating probes')
      
      #Try top priority targets first and then try swapping most conflicted probe(s) in turn.
      angs_new=angs
      angs_new=wrapMinimizeConflicts(pos=cpos, cegs=cegs, angs=angs, probe_conflicts=conflicts)
      probes=defineprobe(x=cpos[,'x'],y=cpos[,'y'],angs=angs_new, interactive=FALSE)
      
      #Check if could configure:
      conflicts=find_probe_conflicts(probes=probes, pos=cpos, angs=angs_new)
      nconflicts=dim(conflicts)[1]
      
      if (nconflicts >=1) {
      } else{
        print('Woot woot! I managed to find a solution!')
        return(list(pos=cpos,angs=angs_new))
      }
    }
  }
  return(list(pos=cpos,angs=angs_new))
}


#Configuring the standard stars.
configure_stdstars <- function(pos, angs, cegs, stdpos){
  
  print("Adding some standards.")
  
  #Removing the standards that overlap with allocated probes.
  clstdpos=cleanconflictedstds(angs=angs, pos=pos, stdpos=stdpos)
  
  #* Added by Sam- this seems to fix the case where cleanconflictedstds returns a list with NAs in it
  clstdpos=clstdpos[!is.na(clstdpos)[,1],]
  
  if (lengths(clstdpos)[1] < nstdprobes) {
    print(paste('**Seem to only have ', lengths(clstdpos)[1], ' un-clashing standard stars in this field! Exiting and trying again...', sep=''))
    stdflag='StdFail'
    #* Added by Sam so we can catch when we error and when we just can't configure
    stop("HECTOR configuration failed: fewer than two non-conflicting standard stars are available.", call. = FALSE)
  }
  
  draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=clstdpos[,'x'],ys=clstdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  
  #Adding two standards to the probes:
  cs=0 #This will keep track of the configured standards so far.
  std=1 #This will cycle through all available standard stars.
  newpos=pos
  newangs=angs
  newcegs=cegs
  stdflag='none'
  
  while (cs < nstdprobes & std <= length(clstdpos[,1])) {
    nprobes <<- ngalprobes+cs+1
    tempos = rbind(newpos,clstdpos[std,])
    newceg=choose_cegs(pos=tempos)[length(tempos[,'x'])]
    tempcegs=c(newcegs,newceg)
    guessang = first_guess_angs(pos=tempos, cegs=tempcegs, guide=FALSE)[length(newangs)+1]
    tempangs = c(newangs,guessang)
    temprobes = defineprobe(x=tempos[,'x'], y=tempos[,'y'], angs=tempangs, interactive=FALSE)
    tempconflicts = find_probe_conflicts(probes=temprobes, pos=tempos, angs=tempangs)
    nconflicts = dim(tempconflicts)[1]
    
    if (nconflicts > 0) {
      #print(tempos)
      #print(tempangs)
      #print(tempconflicts)
      try(withTimeout(expr=assign('tempangs',MinimizeConflicts(pos=tempos, angs=tempangs, cegs=tempcegs, probe_conflicts=tempconflicts)), timeout=timeout, onTimeout=timeout_behaviour))

      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=clstdpos[,'x'],ys=clstdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      temprobes = defineprobe(x=tempos[,'x'],y=tempos[,'y'],angs=tempangs, interactive=visualise)
      
      tempconflicts = find_probe_conflicts(probes=temprobes, pos=tempos, angs=newangs)
      nconflicts = dim(tempconflicts)[1]
      if (nconflicts > 0){
      }else{
        print(paste0('standard ', as.character(std), ' configured successfully.'))
        newpos=tempos
        newangs=tempangs
        newcegs=tempcegs
        cs = cs+1
        nprobes <<- ngalprobes+cs
      }
    } else{
      print(paste0('standard ', as.character(std), ' configured successfully.'))
      newpos=tempos
      newangs=tempangs
      newcegs=tempcegs
      cs = cs+1
      nprobes <<- ngalprobes+cs
    }
    std = std+1
  }
  if (cs < nstdprobes) {
    print(paste('**Could only configure ',as.character(cs), ' standard stars in this field. Boo! :-('), sep='')
    stdflag='StdFail'
    #* Added by Sam so we can catch when we error and when we just can't configure
    stop("HECTOR configuration failed: fewer than two non-conflicting standard stars could be configured.", call. = FALSE)
  }
  print(list(pos=newpos,angs=newangs, cegs=newcegs, stdflag=stdflag))
  return(list(pos=newpos,angs=newangs,cegs=newcegs, stdflag=stdflag))
}

#Configuring the guide stars.
configure_guidestars <- function(pos, angs, cegs, guidepos){ #In this case, pos must also include the 2 configured standards.
  
  #Removing the standards that overlap with allocated probes.
  clguidepos=cleanconflictedguides(angs=angs, pos=pos, guidepos=guidepos)
  
  #* Added by Sam- this seems to fix the case where cleanconflictedguides returns a list with NAs in it
  clguidepos=clguidepos[!is.na(clguidepos)[,1],]
  guideflag=''

  # Six guides are required for a complete automatic configuration. If the
  # positional cleaning removes some candidates, add all remaining candidates
  # to the queue. They must still be tested by the angular optimiser; being
  # rejected by the initial positional filter is not proof that they cannot
  # be configured.
  if (nrow(guidepos) < ngprobesmax) {
    stop("HECTOR configuration failed: the input contains only ", nrow(guidepos),
         " guide-star candidates; six are required.", call. = FALSE)
  }
  manual_candidate_keys <- character(0)
  if (nrow(clguidepos) < ngprobesmax) {
    print(paste('**Only ', nrow(clguidepos),
                ' guide stars are conflict-free at the initial stage; queueing all remaining candidates for further optimisation.',
                sep=''))

    guide_key <- paste(guidepos[,'x'], guidepos[,'y'], sep='\r')
    clean_key <- paste(clguidepos[,'x'], clguidepos[,'y'], sep='\r')
    remaining <- guidepos[!(guide_key %in% clean_key),,drop=FALSE]
    if (nrow(remaining) > 0) {
      manual_candidate_keys <- paste(remaining[,'x'], remaining[,'y'], sep='\r')
      clguidepos <- rbind(clguidepos, remaining)
    }
  }
  if (nrow(clguidepos) < ngprobesmin) {
    stop("HECTOR configuration failed: fewer than ", ngprobesmin,
         " guide-star candidates are available for manual intervention.", call. = FALSE)
  }

  # The cleaning helper rebuilds a data frame while ordering by quadrant.  Put
  # the original IDs back explicitly so output tables can join guide rows.
  guide_key <- paste(guidepos[,'x'], guidepos[,'y'], sep='\r')
  selected_key <- paste(clguidepos[,'x'], clguidepos[,'y'], sep='\r')
  rownames(clguidepos) <- rownames(guidepos)[match(selected_key, guide_key)]
  
  #Drawing and defining configured hector probes.
  draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=clguidepos[,'x'],yg=clguidepos[,'y'], interactive=visualise)
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  
  #Adding six guide probes. Failed candidates are held aside while all
  #remaining candidates are tested.
  cg=0 #This will keep track of the number of configured guides so far.
  sel_guides=c() #This keeps track of selectec guide star probe numbers.
  gui=1 #This will cycle through all available guide stars.
  tempgangs=c()
  newgangs=c()
  newgpos=clguidepos[0, c('x','y'), drop=FALSE]
  ngprobes <<- 1 #Starting with one probe.
  failed_gpos=clguidepos[0, c('x','y'), drop=FALSE]
  failed_gangs=c()
  failed_guides=c()
  
  while (cg < ngprobesmax & gui <= length(clguidepos[,1])) {
    tempgpos=rbind(newgpos,clguidepos[gui,])
    active_gangs=newgangs
    
    #Redraw:
    draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
    
    probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
    
    gcegs=choose_cegs(pos=tempgpos)
    guessang=first_guess_angs(pos=tempgpos,cegs=gcegs, guide=TRUE)[cg+1]
    tempgangs=c(active_gangs,guessang)
    tempgprobes=defineguideprobe(gxs=tempgpos[,'x'],gys=tempgpos[,'y'],gangs=tempgangs, interactive=visualise)
    
    tempconflicts=find_guide_conflicts(probes=probes, pos=pos, angs=angs, gprobes=tempgprobes, gpos=tempgpos, gangs=tempgangs)
    nconflicts=dim(tempconflicts)[1]

    retained=FALSE
    candidate_key <- paste(clguidepos[gui,'x'], clguidepos[gui,'y'], sep='\r')
    manual_candidate <- candidate_key %in% manual_candidate_keys
    if (nconflicts > 0) {
      active_angs=angs
      newgangs=tempgangs
      solver_result=NULL
      try(withTimeout(expr=assign('solver_result',MinimizeGuideConflicts(pos=pos, angs=angs, cegs=cegs, gpos=tempgpos, gangs=tempgangs, gcegs=gcegs, guide_conflicts=tempconflicts)), timeout=timeout, onTimeout=timeout_behaviour))
      if (is.list(solver_result) &&
          all(c('angs','gangs') %in% names(solver_result))) {
        angs=solver_result$angs
        newgangs=solver_result$gangs
      } else if (!is.null(solver_result)) {
        newgangs=solver_result
      }
      # The joint fallback may have moved target/standard angles as well as
      # guide angles, so rebuild the target polygons before rechecking.
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                         interactive=visualise)
      tempgprobes=defineguideprobe(gxs=tempgpos[,'x'],gys=tempgpos[,'y'],gangs=newgangs, interactive=visualise)
      fov_status=configuration_fov_status(probes, tempgprobes)
      tempconflicts=find_guide_conflicts(probes=probes, pos=pos, angs=angs, gprobes=tempgprobes, gpos=tempgpos, gangs=newgangs)
      nconflicts=dim(tempconflicts)[1]
      retained=(nconflicts == 0 && fov_status$ok)
      if (!retained) {
        # Do not let an unsuccessful joint rescue alter the target field used
        # by subsequent guide candidates.
        angs=active_angs
      }
    } else{
      newgpos=tempgpos
      newgangs=tempgangs
      retained=TRUE
    }

    if (!retained) {
      print(paste('**Could not configure guide star ', as.character(gui),
                  '; trying the next candidate.', sep=''))
      # Keep failed candidates separate from the active configuration. This
      # lets later candidates be tested without an unresolved guide blocking
      # them. Failed candidates can still be appended for manual repair after
      # all candidates have been tried.
      candidate_ang=tempgangs[length(tempgangs)]
      if (length(newgangs) == nrow(tempgpos)) {
        candidate_ang=newgangs[length(newgangs)]
      }
      failed_gpos=rbind(failed_gpos, clguidepos[gui, c('x','y'), drop=FALSE])
      failed_gangs=c(failed_gangs, candidate_ang)
      failed_guides=c(failed_guides, gui)
      newgangs=active_gangs
    } else {
      if (manual_candidate) {
        print(paste('Guide star ', as.character(gui),
                    ' was retained after angular optimisation.', sep=''))
      }
      newgpos=tempgpos
      cg = cg+1
      sel_guides=c(sel_guides,gui)
      ngprobes <<- cg+1
    }
    gui = gui+1
  }

  # If angular optimisation could not place six guides, retain the failed
  # candidates as a clearly flagged manual-intervention configuration. This is
  # deliberately done only after every available candidate has been tested.
  if (cg < ngprobesmax && nrow(failed_gpos) > 0) {
    n_to_fill=min(ngprobesmax-cg, nrow(failed_gpos))
    print(paste('Retaining ', as.character(n_to_fill),
                ' unresolved guide candidate(s) for manual intervention.',
                sep=''))
    fill=seq_len(n_to_fill)
    newgpos=rbind(newgpos, failed_gpos[fill, , drop=FALSE])
    newgangs=c(newgangs, failed_gangs[fill])
    sel_guides=c(sel_guides, failed_guides[fill])
    cg=nrow(newgpos)
    guideflag='GuideFail'
  }
  print(paste('Configured ',as.character(cg), ' guide stars in this field'), sep='')
  if (cg < ngprobesmax) {
    print(paste('**Could only retain ',as.character(cg),
                ' guide stars; manual intervention is required.', sep=''))
    guideflag='GuideFail'
  }
  ngprobes <<- nrow(newgpos)+1

  # Preserve the source IDs after repeated rbind operations.
  selected_key <- paste(newgpos[,'x'], newgpos[,'y'], sep='\r')
  rownames(newgpos) <- rownames(guidepos)[match(selected_key, guide_key)]
  
  return(list(pos=pos, angs=angs, gpos=newgpos, gangs=newgangs,
              guidenum=sel_guides, guideflag=guideflag))
}

configure_HECTOR <- function(pos_master=pos_master){
  
  ngalprobes<<-dngalprobes
  
  if (length(pos_master$pos[,1]) < dngalprobes) {
    ngalprobes<<-length(pos_master$pos[,1])
  }
  
  fieldflag=''
  #Keeping only targets that do not overlap with higher priority targets:
  tpts = SpatialPoints(coords=as.matrix(pos_master$pos[,c('x','y')])) #sp objects for the targets as simple points.
  tokeep=c(1)
  if (length(tpts)>1){
    for (tpt in c(2:length(tpts))) {
      #Following line replaced Sept 2026 by CFoster as used deprecated rgeos.
      #mindist=apply(gDistance(tpts[tpt], tpts[1:tpt-1], byid=TRUE),2,min) #The minimum distance to the nearest higher priority neighbour for each target.
      mindist = min(st_distance(st_as_sf(tpts[tpt]), st_as_sf(tpts[1:(tpt-1)])))
      if (mindist > excl_radius*2.){
        tokeep=c(tokeep,tpt)
      }
      else{
        print("We seem to have a clash before we've even started... Error in Tiling Code??")
        fieldflag='fieldfail'
        #* Added by Sam so we can catch when we error and when we just can't configure
        stop("HECTOR configuration failed: unresolved target conflicts remain after swaps.", call. = FALSE)
      }
    }
  }
  pos_master$pos = pos_master$pos[tokeep,c('x','y')]
  
  #Drawing possible target, standard and guide star positions:
  draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
  swaps=FALSE
  nprobes<<-ngalprobes
  pos=pos_master$pos[1:ngalprobes,]
  
  cegs=choose_cegs(pos=pos)
  angs=first_guess_angs(pos=pos, cegs=cegs)
  angsmin<<-angs #Initialising global variable angsmin.
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
  
  print('Identifying initial conflicts')
  conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
  print(conflicts)
  nconflicts=dim(conflicts)[1]
  if (nconflicts <1) {
    print('Bam! No initial conflicts, this one is done!')
  } else{
    print('Minimizing the overlaping area by rotating probes')
    delta_ang=20 #Start with 20 degrees
    if (nconflicts > 0){
      #Try top priority targets first (for a given time) and then try swapping most conflicted probe(s) in turn.
      #Usually if there is a solution, it is found well before the timer runs out. 
      #This optimal time limit may depend on the density of the field and size of the probes though...?
      
      angs=wrapMinimizeConflicts(pos=pos, angs=angs, cegs=cegs, init_delta_ang=delta_ang, probe_conflicts=conflicts)
      
      draw_all(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'],xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'], interactive=visualise)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=visualise)
      #Check if could configure:
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
    }
  }
  
  
  
  if (nconflicts > 0) {
    print('Could not find a solution, I am going to try swaps. This could take a while, grab a cuppa!')
    swaps=TRUE
    temp=try_swaps(pos_master=pos_master)
    tempos=temp$pos
    angs_new=temp$angs
    print('Looking for remaining conflicts')
    probes=defineprobe(x=tempos[,'x'],y=tempos[,'y'],angs=angs_new, interactive=FALSE)
    conflicts=find_probe_conflicts(probes=probes, pos=tempos, angs=angs_new)
    nconflicts=dim(conflicts)[1]
    if (nconflicts > 0) {
      print('Sorry, I tried my best but really could not configure this field :-(.')
      fieldflag='fieldfail'
      #* Added by Sam so we can catch when we error and when we just can't configure
      stop("HECTOR configuration failed: unresolved target conflicts remain after swaps.", call. = FALSE)
    }else{
      pos=tempos
      angs=angs_new
    }
  } else{
    print('19 galaxy probes are in! Let\'s put some standard stars probes!')
    swaps=FALSE
  }
  
  #Configuring standard probes:
  withstdstars=configure_stdstars(pos=pos, angs=angs, cegs=cegs, stdpos=pos_master$stdpos)
  
  pos=withstdstars$pos #Including the 2 standards in the list of positions.
  angs=withstdstars$angs
  stdflag=withstdstars$stdflag
  cegs=withstdstars$cegs
  print('19 galaxy probes and 2 standard probes are in! Let\'s put some guide probes!')
  
  #Configure six guide stars; retain unresolved candidates for manual repair.
  guides=configure_guidestars(pos=pos, angs=angs, cegs=cegs, guidepos=pos_master$guidepos)
  # The guide-target joint optimiser may have moved a target or standard
  # angle, so carry those updated target angles into the final output.
  pos=guides$pos
  angs=guides$angs
  final_probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs,
                           interactive=FALSE)
  final_gprobes=defineguideprobe(gxs=guides$gpos[,'x'],
                                 gys=guides$gpos[,'y'],
                                 gangs=guides$gangs, interactive=FALSE)
  final_fov_status=configuration_fov_status(final_probes, final_gprobes)
  fovflag='none'
  if (!final_fov_status$ok) {
    print('**At least one configured probe lies outside the HECTOR field of view.')
    print(final_fov_status)
    fovflag='FieldFail'
  }
  if (visualise) {plot_configured_field(pos=pos,angs=angs,gpos=guides$gpos,gangs=guides$gangs, aspdf=FALSE)}
  
  #Saving output to some random global variable for future tweaking in case of a late crash..
  #temp<<-list(pos=pos,angs=angs,gpos=guides$gpos,gangs=guides$gangs, fieldflags=final_config$flags)
  
  #Creating an output file for the robot:
  #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
  azAngs=(pi+atan2(pos[,'y'],pos[,'x']))
  azAngs[azAngs > 3*pi/2]=azAngs[azAngs > 3*pi/2]-2*pi
  azAngs[azAngs < -1*pi/2]=azAngs[azAngs < -1*pi/2]+2*pi
  angs_azAng=angs-azAngs
  angs_azAng[angs_azAng > pi] = angs_azAng[angs_azAng > pi] - 2*pi
  angs_azAng[angs_azAng < 0] = angs_azAng[angs_azAng < 0] + 2*pi
  rads=sqrt(pos[,'x']**2+pos[,'y']**2)
  
  gazAngs=(pi+atan2(guides$gpos[,'y'],guides$gpos[,'x']))
  gazAngs[gazAngs > 3*pi/2]=gazAngs[gazAngs > 3*pi/2]-2*pi
  gazAngs[gazAngs < -1*pi/2]=gazAngs[gazAngs < -1*pi/2]+2*pi
  gangs_gazAng=guides$gangs-gazAngs
  gangs_gazAng[gangs_gazAng > pi] = gangs_gazAng[gangs_gazAng > pi] - 2*pi
  gangs_gazAng[gangs_gazAng < 0] = gangs_gazAng[gangs_gazAng < 0] + 2*pi
  grads=sqrt(guides$gpos[,'x']**2+guides$gpos[,'y']**2)
  
  print('Checking wedge issue.')
  fieldflags=check_wedge_issue(pos=pos, angs=angs, gpos=guides$gpos[1:(ngprobes-1),], gangs=guides$gangs, interactive=visualise)

  manual_intervention <- (fovflag == 'FieldFail' ||
                          guides$guideflag == 'GuideFail' ||
                          fieldflag == 'fieldfail' ||
                          (nzchar(fieldflags) && fieldflags != 'none') ||
                          (nzchar(stdflag) && stdflag != 'none'))
  flag_parts <- c(fovflag, fieldflags, stdflag, guides$guideflag)
  flag_parts <- flag_parts[nzchar(flag_parts) & flag_parts != 'none']
  flags <- if (length(flag_parts)) paste(flag_parts, collapse=' ') else 'none'
  final_config=list(pos=pos, rads=rads, angs=angs, azAngs=azAngs, angs_azAng=angs_azAng, guidenum=guides$guidenum, gpos=guides$gpos, grads=grads, gangs=guides$gangs,gazAngs=gazAngs, gangs_gazAng=gangs_gazAng, swaps=swaps, flags=flags, guide_count=nrow(guides$gpos), guide_required=ngprobesmax, manual_intervention=manual_intervention)
  
  print('Here you go, here\'s your complete configuration!')
  print(final_config)
  
  return (final_config)
}

plot_configured_field <- function(filename='temp.pdf',pos,angs,gpos,gangs, fieldflags='none', aspdf=TRUE){
  if (aspdf) {pdf(file=filename, width=6, height=6)}
  draw_pos(x=pos_master$pos[,'x'],y=pos_master$pos[,'y'], label_probe=FALSE)
  draw_pos(x=pos[,'x'],y=pos[,'y'],label_probe=FALSE,overwrite=FALSE)
  draw_stds(xs=pos_master$stdpos[,'x'],ys=pos_master$stdpos[,'y'],label_probe=FALSE)
  draw_guides(xg=pos_master$guidepos[,'x'],yg=pos_master$guidepos[,'y'])
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=TRUE)
  gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=gangs, interactive=TRUE)
  # Draw labels last so the probe footprints cannot obscure them.  Targets
  # use their historical numeric labels; standards are explicitly marked as
  # S1/S2; guides are marked G1/G2/... to distinguish them from targets.
  n_targets <- min(ngalprobes, nrow(pos))
  if (n_targets > 0) {
    text(pos[seq_len(n_targets), 'x'], pos[seq_len(n_targets), 'y'],
         labels=as.character(seq_len(n_targets)), cex=0.65, font=2)
  }
  if (nrow(pos) > n_targets) {
    std_idx <- seq.int(n_targets + 1, nrow(pos))
    text(pos[std_idx, 'x'], pos[std_idx, 'y'],
         labels=paste0('S', seq_along(std_idx)), cex=0.65, font=2)
  }
  if (nrow(gpos) > 0) {
    text(gpos[,'x'], gpos[,'y'],
         labels=paste0('G', seq_len(nrow(gpos))), cex=0.75, font=2)
  }
  mtext(paste('Flags:', fieldflags), cex=1, side=3)
  # Close only the PDF device opened above. The interactive device must stay
  # open so live configuration plots remain visible to the caller.
  if (aspdf) {dev.off()}
}

check_wedge_issue <- function(pos, angs, gpos, gangs, interactive=FALSE){
  
  ngprobes=length(gangs)
  
  probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=interactive)
  gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=gangs, interactive=interactive)
  
  fieldflags=''
  
  #First check if there are probe heads too close together to let a cable through.
  #If there are pinching probes, check there are no probes within the bend radius 
  #of both pinch ferules. If so, make sure the end of the ferule of wedged probe 
  #is facing out. If not, throw a flag.
  
  #Between target or std probes:
  pinches=c()
  for (probe in c(1:(nprobes))){
    dists=sqrt( (pos[(probe+1):nprobes,'x']- pos[probe,'x'])^2 + (pos[(probe+1):nprobes,'y']- pos[probe,'y'])^2)
    neigh=which(dists < (excl_radius*2+cable_w))+probe
    if (length(neigh) > 0) {
      pinches=rbind(pinches, c(probe,neigh[neigh!=probe],NaN,NaN)) #Identify potential wedges.
    }
  }
  #colnames(pinches)=c('probe1','probe2','guide1','guide2')
  
  #Between guide and target or probes:
  for (probe in c(1:nprobes)){
    dists=sqrt( (gpos[,'x']- pos[probe,'x'])^2 + (gpos[,'y']- pos[probe,'y'])^2)
    neigh=which(dists < (excl_radius*2+cable_w))
    if (length(neigh) > 0) {
      pinches=rbind(pinches, c(probe,NaN,neigh,NaN)) #Identify potential wedges.
    }
  }
  
  #Between guide probes:
  for (gprobe in c(1:(ngprobes-1))){
    dists=sqrt( (gpos[(gprobe+1):ngprobes,'x']- gpos[gprobe,'x'])^2 + (gpos[(gprobe+1):ngprobes,'y']- gpos[gprobe,'y'])^2)
    neigh=which(dists < (excl_radius*2+cable_w))+gprobe
    if (length(neigh) > 0) {
      pinches=rbind(pinches, c(NaN,NaN,gprobe,neigh[neigh!=gprobe])) #Identify potential wedges.
    }
  }
  
  npinches=dim(pinches)[1]
  print(pinches)
  
  #Looking for probes with cable within bend radius (fbr) of both probes in pairs of pinches.
  if (!is.null(npinches)) {
    print(paste0('Found ', as.character(npinches), " potential pinch points."))
    print("Making sure no cable needs to go through those.")
    
    #Calculating ferule tips for probes and guides:
    ftips=calc_ferule_tip_pos(pos[,'x'],pos[,'y'],ang=angs)
    gftips=calc_ferule_tip_pos(x=gpos[,'x'],y=gpos[,'y'],ang=gangs,guide=T)
    
    ftips=rbind(ftips,gftips)
    
    #Cheekily appending the gftips (gpos) positions at the end of the ftips (pos) to simplify coding.
    #First nprobes ftips are for probes, nprobes+1 --> nprobes+ngprobes are guide ferule tip 
    #positions.
    pinches[,3]=pinches[,3]+nprobes
    pinches[,4]=pinches[,4]+nprobes
    
    tempos=rbind(pos,gpos)
    tempangs=c(angs,gangs)
    
    fpinched=c()
    for (pinch in c(1:npinches)){
      pprobes=pinches[pinch,!is.na(pinches[pinch,])]
      dists1=sqrt( (ftips[1:(nprobes+ngprobes),'xt']- ftips[pprobes[1],'xt'])^2 + (ftips[1:(nprobes+ngprobes),'yt']- ftips[pprobes[1],'yt'])^2)
      dists2=sqrt( (ftips[1:(nprobes+ngprobes),'xt']- ftips[pprobes[2],'xt'])^2 + (ftips[1:(nprobes+ngprobes),'yt']- ftips[pprobes[2],'yt'])^2)
      pinched=which((dists1 < fbr) & (dists2 < fbr))
      # Keep only probes that are part of the current pinch pair.  The legacy
      # code used a non-base `%nin%` operator supplied by the old workspace.
      pinched=pinched[!(pinched %in% pprobes)]
      if (length(pinched)>0) {
        #Check if 180 degrees from ang of the pinched probe lands within the broad angle defined between
        #the wedges.
        flipang=tempangs[pinched]+pi
        flipang[flipang > (3*pi/2)]=flipang[flipang > (3*pi/2)]-2.*pi
        flipang[flipang < (-1*pi/2)]=flipang[flipang > (-1*pi/2)]+2.*pi
        pinched=c(
        pinched[((flipang > min(tempangs[pprobes])) & (flipang < max(tempangs[pprobes])) 
                 & (abs(diff(tempangs[pprobes])) < pi))],
        pinched[((flipang < (min(tempangs[pprobes])+2*pi)) & (flipang > max(tempangs[pprobes])) 
                 & (abs(diff(tempangs[pprobes])) > pi))])
        fpinched=c(pinched,fpinched)
      }
    }
    gpinched=fpinched[fpinched > nprobes]-nprobes
    pinched=fpinched[fpinched <= nprobes]
    if (length(pinched>0)){
      fieldflags=paste(fieldflags,'check probe',as.character(pinched))
    }
    if (length(gpinched>0)){
      fieldflags=paste(fieldflags,'check guide',as.character(gpinched))
    }
    if (fieldflags=='') {print('All good :D!')} else {print(fieldflags)}
  }
  
  return(fieldflags)
}

#calculate position of tip of ferule:
calc_ferule_tip_pos <- function(x,y,ang,guide=F) {
  if (guide){
    xt=x-(gtip_l/2.+gprobe_l)*cos(ang)
    yt=y-(gtip_l/2.+gprobe_l)*sin(ang)
  }else{
    xt=x-(tip_l/2.+probe_l)*cos(ang)
    yt=y-(tip_l/2.+probe_l)*sin(ang)
  }
  return(cbind(xt,yt))
}
