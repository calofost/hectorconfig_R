# Dependencies and configuration functions are supplied by the package. This
# file is evaluated inside the private engine environment created by
# hector_checker_app(), so it does not source code or read files relative to
# the user's working directory.

ui <- fluidPage(
  fluidPage(width=8,
            useShinyjs(),
            sidebarLayout(position='left',
                          sidebarPanel(
                            #Input files:
                            fileInput("hexafile", "Choose Hexa file", accept = ".csv"),
                            fileInput("guidefile", "Choose Guide file", accept = ".csv"),
                            hr(),
                            textOutput("selected_probe"),
                            actionButton(inputId='fixProbe', label='Fix Probe.'),
                            actionButton(inputId='ckConflict', label='List conflicted probes.'),
                            textOutput("ck_output"),
                            hr(),
                            downloadButton("downloadHexa", "Save new Hexa file"),
                            downloadButton("downloadGuide", "Save new Guide file"),
                            hr("Only use the following option(s) if you know what you're doing (see instructions):"),
                            h5(""),
                            actionButton(inputId='pReset', label='Reset Target Probes.'),
                            actionButton(inputId='gReset', label='Reset Guide Probes.'),
                            actionButton(inputId='sReset', label='Reset Standard Probes.')
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Configuration",
                                       plotOutput("plotfield", dblclick = "plot_dblclick", click='plot_click')
                              ),
                              tabPanel("Instructions",
                                       h5("1- Upload the output configuration files by clicking the corresponding buttons."), 
                                       h5("2- Double-click on the target/probe-head of the probe you wish to move."), 
                                       h5("3- Single-click towards the direction you wish the tail of the probe to go."),
                                       h5("4- When you're happy with the new position of the probe, click the \"Fix Probe\" button."),
                                       h5("5- Repeat 2-5 for all probes that need moving."), 
                                       h5("6- When you're happy with the field, save the updated hexa and guide files onto your machine
                                          by clicling the relevant buttons."),
                                       h5("Beware,- "),
                                       h5("- the code no longer checks for conflicts, and it does NOT check that you haven't oriented probes 
                                          outside the field-of-view."),
                                       h5("- when starting a new field, be sure to upload both the hexa and guide file 
                                          for that field, otherwise you'll have conflicts ;)."),
                                       h5("Advanced options:-"),
                                       h5("- you can chose to \"Reset Guide Probes\" orientations to the initial guess, but beware of conflicts (not checked)."),
                                       h5("----"),
                                       h6("Questions, issues, comments and compliments (strictly no insults :P)? Contact 
                                          Caroline Foster on c.foster (at) unsw.edu.au or caro.foster (at) gmail.com.")
                                       )
                              )
                            )
                          )
            )
)

server <- function(input, output) {
  
  output$plotfield <- renderPlot(expr={
    magplot(0,0,pch='x',col="red",xlim=c(-fov/2.,fov/2.),ylim=c(-fov/2.,fov/2.), asp=1, xlab='X (mm)', ylab='Y (mm)')
    draw.circle(0,0,radius=fov/2.,border='red')
    draw.circle(0,0,radius=1.,border='black')
    #Draw the cable exit gaps in cyan
    for (cceg in c(1:length(ceg_pos[,1]))){
      dx=wceg/2.* sind(ceg_pos[cceg,'deg'])
      dy=-1.*wceg/2.* cosd(ceg_pos[cceg,'deg'])
      lines(c(ceg_pos[cceg,'x']-dx,ceg_pos[cceg,'x']+dx),c(ceg_pos[cceg,'y']-dy,ceg_pos[cceg,'y']+dy), lwd=3, col='cyan')
    }
    text(as.character(c(1:3)), x=ceg_pos[,'x'],y=ceg_pos[,'y'])
  },width = 650,
  height = 650,
  res = 72)
  observeEvent(list(input$hexafile, input$guidefile), {
    req(input$hexafile)
    req(input$guidefile)
    hfile <<- input$hexafile
    gfile <<- input$guidefile

    hexaconfig <<- read.table(
      file = hfile$datapath, sep = ',', header = TRUE, comment.char = '#',
      colClasses = list("ID" = "character")
    )
    guideconfig <<- read.table(
      file = gfile$datapath, sep = ',', header = TRUE, comment.char = '#',
      colClasses = list("ID" = "character")
    )

    required_hex <- c("fibre_type", "x", "y", "angs")
    required_guide <- c("x", "y", "angs")
    missing_hex <- setdiff(required_hex, names(hexaconfig))
    missing_guide <- setdiff(required_guide, names(guideconfig))
    if (length(missing_hex) || length(missing_guide)) {
      stop(
        paste0(
          "Configured files must contain the geometry columns. Missing Hexa: ",
          paste(missing_hex, collapse = ", "),
          "; missing Guide: ", paste(missing_guide, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Only select the probes, ignoring the sky fibres.
    probes <<- hexaconfig[,'fibre_type'] == 'P'
    pos <<- hexaconfig[probes,c('x','y')]
    angs <<- hexaconfig[probes,'angs']
    gpos <<- guideconfig[,c('x','y')]
    gangs <<- guideconfig[,'angs']
    # The conflict checker and the Fix Probe button operate on these
    # temporary vectors. Initialise them immediately after loading a field
    # so checking a freshly uploaded configuration cannot fail with
    # "object 'tempangs' not found".
    tempangs <<- angs
    tempgangs <<- gangs

    # The exported guide table stores these values as angs/rads/azAngs and
    # angs_azAng.  The interactive app uses the internal g-prefixed names;
    # initialise them on upload so the Guide download button also works
    # immediately, before any probe has been manually edited.
    gazAngs <<- (pi + atan2(gpos[,'y'], gpos[,'x']))
    gazAngs[gazAngs > 3*pi/2] <<- gazAngs[gazAngs > 3*pi/2] - 2*pi
    gazAngs[gazAngs < -1*pi/2] <<- gazAngs[gazAngs < -1*pi/2] + 2*pi
    gangs_gazAng <<- gangs - gazAngs
    gangs_gazAng[gangs_gazAng > pi] <<- gangs_gazAng[gangs_gazAng > pi] - 2*pi
    gangs_gazAng[gangs_gazAng < 0] <<- gangs_gazAng[gangs_gazAng < 0] + 2*pi
    grads <<- sqrt(gpos[,'x']**2 + gpos[,'y']**2)
    guideconfig[,c('gangs','gazAngs','gangs_gazAng','grads')] <<-
      cbind(gangs, gazAngs, gangs_gazAng, grads)

    fieldflags <<- ''
    nhexa <<- length(pos[,1])
    pos_master <<- list(
      pos = pos[c(1:(nhexa-2)),],
      stdpos = pos[c((nhexa-1):nhexa),],
      guidepos = gpos
    )

    output$plotfield <- renderPlot(expr={
      plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
    },width = 650,
    height = 650,
    res = 72)
    
    },
    ignoreNULL = TRUE
  )
  observeEvent(
    input$plot_dblclick,
    {
      cpos=as.numeric(c(input$plot_dblclick$x,input$plot_dblclick$y))
      #Determine which probe is being moved.
      probe1=which(pointDistance(pos,cpos, lonlat=F)==min(pointDistance(pos,cpos, lonlat=F)))
      probe2=which(pointDistance(gpos,cpos, lonlat=F)==min(pointDistance(gpos,cpos, lonlat=F)))
      if(pointDistance(cpos,pos[probe1,], lonlat=F) < pointDistance(cpos,gpos[probe2,], lonlat=F)){
        guide=FALSE
        probe=probe1
        cpos=pos[probe1,]
        if (probe <= ngalprobes){
          output$selected_probe <- renderText({ 
            paste0('Modifying angle for hexa probe ', as.character(probe))
          })
        } else {
          output$selected_probe <- renderText({ 
            paste0('Modifying angle for standard probe ', as.character(probe-ngalprobes))
          })
        }
      }else{
        guide=TRUE
        probe=probe2
        cpos=gpos[probe2,]
        output$selected_probe <- renderText({ 
          paste0('Modifying angle for guide probe ', as.character(probe))
        })
      }
      output$plotfield <- renderPlot(expr={
        plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
        draw.circle(cpos[,'x'],cpos[,'y'],radius=excl_radius, border='red')
        rad=probe_l+cable_l+tip_l/2.
        if (guide){
          points(x=gpos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=5)),y=gpos[probe,'y']+rad*sind(seq(from=0,to=360, by=5)),pch=16,cex=0.2)
          points(x=gpos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=15)),y=gpos[probe,'y']+rad*sind(seq(from=0,to=360, by=15)),pch=16,cex=0.5)
          text(x=gpos[probe,'x']-1*(rad+tip_l)*cosd(seq(from=0,to=359, by=30)),y=gpos[probe,'y']-(rad+tip_l)*sind(seq(from=0,to=359, by=30)), labels=as.character(seq(from=0,to=359, by=30)))
          lines(x=c(gpos[probe,'x'], gpos[probe,'x']-rad), y=c(gpos[probe,'y'], gpos[probe,'y']))
        }else{
          points(x=pos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=5)),y=pos[probe,'y']+rad*sind(seq(from=0,to=360, by=5)),pch=16,cex=0.2)
          points(x=pos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=15)),y=pos[probe,'y']+rad*sind(seq(from=0,to=360, by=15)),pch=16,cex=0.5)
          text(x=pos[probe,'x']-1*(rad+tip_l)*cosd(seq(from=0,to=359, by=30)),y=pos[probe,'y']-(rad+tip_l)*sind(seq(from=0,to=359, by=30)), labels=as.character(seq(from=0,to=359, by=30)))
          lines(x=c(pos[probe,'x'], pos[probe,'x']-rad), y=c(pos[probe,'y'], pos[probe,'y']))
        }
      },width = 650,
      height = 650,
      res = 72)
      #Single click will move the tail towards the click
      observeEvent(
        input$plot_click,
        {
          if (exists('probe')){
            probetail=as.numeric(c(input$plot_click$x,input$plot_click$y))
            x=as.numeric(probetail[1]-cpos[1])
            y=as.numeric(probetail[2]-cpos[2])
            cang=(pi+atan2(y,x))
            if (cang > 3*pi/2) cang=cang-2*pi
            if (cang < -1*pi/2) cang=cang+2*pi
            if (guide){
              tempgangs<<-gangs
              tempgangs[probe]<<-cang
              tempangs<<-angs
            }else{
              tempgangs<<-gangs
              tempangs<<-angs
              tempangs[probe]<<-cang
            }
            output$plotfield <- renderPlot(expr={
              plot_configured_field(pos=pos,angs=tempangs,gpos=gpos,gangs=tempgangs,fieldflags=fieldflags,aspdf=FALSE)
              draw.circle(cpos[,'x'],cpos[,'y'],radius=excl_radius, border='red')
              rad=probe_l+cable_l+tip_l/2.
              if (guide){
                points(x=gpos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=5)),y=gpos[probe,'y']+rad*sind(seq(from=0,to=360, by=5)),pch=16,cex=0.2)
                points(x=gpos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=15)),y=gpos[probe,'y']+rad*sind(seq(from=0,to=360, by=15)),pch=16,cex=0.5)
                text(x=gpos[probe,'x']-1*(rad+tip_l)*cosd(seq(from=0,to=359, by=30)),y=gpos[probe,'y']-(rad+tip_l)*sind(seq(from=0,to=359, by=30)), labels=as.character(seq(from=0,to=359, by=30)))
                lines(x=c(gpos[probe,'x'], gpos[probe,'x']-rad), y=c(gpos[probe,'y'], gpos[probe,'y']))
              }else{
                points(x=pos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=5)),y=pos[probe,'y']+rad*sind(seq(from=0,to=360, by=5)),pch=16,cex=0.2)
                points(x=pos[probe,'x']+1*rad*cosd(seq(from=0,to=360, by=15)),y=pos[probe,'y']+rad*sind(seq(from=0,to=360, by=15)),pch=16,cex=0.5)
                text(x=pos[probe,'x']-1*(rad+tip_l)*cosd(seq(from=0,to=359, by=30)),y=pos[probe,'y']-(rad+tip_l)*sind(seq(from=0,to=359, by=30)), labels=as.character(seq(from=0,to=359, by=30)))
                lines(x=c(pos[probe,'x'], pos[probe,'x']-rad), y=c(pos[probe,'y'], pos[probe,'y']))
              }
            },width = 650,
            height = 650,
            res = 72)
          } else {
            output$selected_probe <- renderText({ 
              paste0('You must double click on a target first.')
            })
          }
        }
      )
    }
  ) 
  
  #Button to fix angle will record new positions if the new angle did not create new conflicts.
  observeEvent(
    input$fixProbe,
    {
      #Checking for created conflicts:
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=tempangs, interactive=FALSE)
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=tempangs)
      nconflicts=dim(conflicts)[1]
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=tempgangs, interactive=FALSE)
      gconflicts=find_guide_conflicts(probes=probes, pos=pos,angs=tempangs,gprobes=gprobes, gpos=gpos, gangs=tempgangs)
      ngconflicts=dim(gconflicts)[1]
      
      output$selected_probe <- renderText({ 
        paste('There are ',nconflicts +ngconflicts,'conflicts in the field', sep=' ')
      })
      output$ck_output <- renderText({ 
        ''
      })
      
      angs<<-tempangs
      gangs<<-tempgangs
      #Creating an output file for the robot:
      #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
      azAngs<<-(pi+atan2(pos[,'y'],pos[,'x']))
      azAngs[azAngs > 3*pi/2]<<-azAngs[azAngs > 3*pi/2]-2*pi
      azAngs[azAngs < -1*pi/2]<<-azAngs[azAngs < -1*pi/2]+2*pi
      angs_azAng<<-angs-azAngs
      angs_azAng[angs_azAng > pi] <<- angs_azAng[angs_azAng > pi] - 2*pi
      angs_azAng[angs_azAng < 0] <<- angs_azAng[angs_azAng < 0] + 2*pi
      rads<<-sqrt(pos[,'x']**2+pos[,'y']**2)

      hexaconfig[hexaconfig[,'fibre_type'] == 'P',c('angs','azAngs','angs_azAng','rads')]<<-cbind(angs,azAngs,angs_azAng,rads)
      #print(hexaconfig)
      
      gazAngs<<-(pi+atan2(gpos[,'y'],gpos[,'x']))
      gazAngs[gazAngs > 3*pi/2]<<-gazAngs[gazAngs > 3*pi/2]-2*pi
      gazAngs[gazAngs < -1*pi/2]<<-gazAngs[gazAngs < -1*pi/2]+2*pi
      gangs_gazAng<<-gangs-gazAngs
      gangs_gazAng[gangs_gazAng > pi] <<- gangs_gazAng[gangs_gazAng > pi] - 2*pi
      gangs_gazAng[gangs_gazAng < 0] <<- gangs_gazAng[gangs_gazAng < 0] + 2*pi
      grads<<-sqrt(gpos[,'x']**2+gpos[,'y']**2)
      guideconfig[,c('gangs','gazAngs','gangs_gazAng','grads')]<<-cbind(gangs,gazAngs,gangs_gazAng,grads)
      #print(guideconfig)
      
      output$plotfield <- renderPlot(expr={
        plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
      },width = 650,
      height = 650,
      res = 72)
    }
  )
  
  observeEvent(
    input$ckConflict,
    { #Checking for created conflicts:
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=tempangs, interactive=FALSE)
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=tempangs)
      nconflicts=dim(conflicts)[1]
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=tempgangs, interactive=FALSE)
      gconflicts=find_guide_conflicts(probes=probes, pos=pos,angs=tempangs,gprobes=gprobes, gpos=gpos, gangs=tempgangs)
      ngconflicts=dim(gconflicts)[1]
      
      if ((nconflicts > 0) | (ngconflicts > 0)){
        conflictedprobes=c(conflicts$Probe_1st,conflicts$Probe_2nd,gconflicts$Probe_1st)
        conflictedstds=conflictedprobes[conflictedprobes>nprobes]-nprobes
        conflictedprobes=conflictedprobes[conflictedprobes<=nprobes]
        conflictedguides=gconflicts$Probe_2nd
        output$ck_output <- renderText({
          paste('CONFLICTED: target probe: ',paste(as.character(unique(conflictedprobes)), collapse=','), '; guide probe:', paste(as.character(unique(conflictedguides)), collapse=','), '; standard probe:', paste(as.character(unique(conflictedstds)), collapse=','), sep=' ')})
      } else {
        output$ck_output <- renderText({ 
          'HURRAY! No conflict detected.'})
      }
    }
  )
  
  observeEvent(
    input$gReset,
    {
      gcegs=choose_cegs(pos=gpos)
      gangs=first_guess_angs(pos=gpos,cegs=gcegs, guide=TRUE)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=FALSE)
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=gangs, interactive=FALSE)
      gconflicts=find_guide_conflicts(probes=probes, pos=pos,angs=tempangs,gprobes=gprobes, gpos=gpos, gangs=tempgangs)
      ngconflicts=dim(gconflicts)[1]
      
      tempgangs<<-gangs
      tempangs<<-angs
      gangs<<-tempgangs
      angs<<-tempangs
      #Creating an output file for the robot:
      #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
      azAngs<<-(pi+atan2(pos[,'y'],pos[,'x']))
      azAngs[azAngs > 3*pi/2]<<-azAngs[azAngs > 3*pi/2]-2*pi
      azAngs[azAngs < -1*pi/2]<<-azAngs[azAngs < -1*pi/2]+2*pi
      angs_azAng<<-angs-azAngs
      angs_azAng[angs_azAng > pi] <<- angs_azAng[angs_azAng > pi] - 2*pi
      angs_azAng[angs_azAng < 0] <<- angs_azAng[angs_azAng < 0] + 2*pi
      rads<<-sqrt(pos[,'x']**2+pos[,'y']**2)
      
      hexaconfig[hexaconfig[,'fibre_type'] == 'P',c('angs','azAngs','angs_azAng','rads')]<<-cbind(angs,azAngs,angs_azAng,rads)
      #print(hexaconfig)
      
      gazAngs<<-(pi+atan2(gpos[,'y'],gpos[,'x']))
      gazAngs[gazAngs > 3*pi/2]<<-gazAngs[gazAngs > 3*pi/2]-2*pi
      gazAngs[gazAngs < -1*pi/2]<<-gazAngs[gazAngs < -1*pi/2]+2*pi
      gangs_gazAng<<-gangs-gazAngs
      gangs_gazAng[gangs_gazAng > pi] <<- gangs_gazAng[gangs_gazAng > pi] - 2*pi
      gangs_gazAng[gangs_gazAng < 0] <<- gangs_gazAng[gangs_gazAng < 0] + 2*pi
      grads<<-sqrt(gpos[,'x']**2+gpos[,'y']**2)
      guideconfig[,c('gangs','gazAngs','gangs_gazAng','grads')]<<-cbind(gangs,gazAngs,gangs_gazAng,grads)
      #print(guideconfig)
      
      output$plotfield <- renderPlot(expr={
        plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
      },width = 650,
      height = 650,
      res = 72)
    }
  )
  
  observeEvent(
    input$pReset,
    {
      cegs=choose_cegs(pos=pos[1:nprobes,])
      tmp=first_guess_angs(pos=pos[1:nprobes,],cegs=cegs, guide=FALSE)
      angs=c(tmp,angs[(nprobes+1):(nprobes+nstdprobes)])
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=FALSE)
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=gangs, interactive=FALSE)
      gconflicts=find_guide_conflicts(probes=probes, pos=pos,angs=tempangs,gprobes=gprobes, gpos=gpos, gangs=tempgangs)
      ngconflicts=dim(gconflicts)[1]
      
      tempgangs<<-gangs
      tempangs<<-angs
      gangs<<-tempgangs
      angs<<-tempangs
      #Creating an output file for the robot:
      #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
      azAngs<<-(pi+atan2(pos[,'y'],pos[,'x']))
      azAngs[azAngs > 3*pi/2]<<-azAngs[azAngs > 3*pi/2]-2*pi
      azAngs[azAngs < -1*pi/2]<<-azAngs[azAngs < -1*pi/2]+2*pi
      angs_azAng<<-angs-azAngs
      angs_azAng[angs_azAng > pi] <<- angs_azAng[angs_azAng > pi] - 2*pi
      angs_azAng[angs_azAng < 0] <<- angs_azAng[angs_azAng < 0] + 2*pi
      rads<<-sqrt(pos[,'x']**2+pos[,'y']**2)
      
      hexaconfig[hexaconfig[,'fibre_type'] == 'P',c('angs','azAngs','angs_azAng','rads')]<<-cbind(angs,azAngs,angs_azAng,rads)
      #print(hexaconfig)
      
      gazAngs<<-(pi+atan2(gpos[,'y'],gpos[,'x']))
      gazAngs[gazAngs > 3*pi/2]<<-gazAngs[gazAngs > 3*pi/2]-2*pi
      gazAngs[gazAngs < -1*pi/2]<<-gazAngs[gazAngs < -1*pi/2]+2*pi
      gangs_gazAng<<-gangs-gazAngs
      gangs_gazAng[gangs_gazAng > pi] <<- gangs_gazAng[gangs_gazAng > pi] - 2*pi
      gangs_gazAng[gangs_gazAng < 0] <<- gangs_gazAng[gangs_gazAng < 0] + 2*pi
      grads<<-sqrt(gpos[,'x']**2+gpos[,'y']**2)
      guideconfig[,c('gangs','gazAngs','gangs_gazAng','grads')]<<-cbind(gangs,gazAngs,gangs_gazAng,grads)
      #print(guideconfig)
      
      output$plotfield <- renderPlot(expr={
        plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
      },width = 650,
      height = 650,
      res = 72)
    }
  )
  
  observeEvent(
    input$sReset,
    {
      cegs=choose_cegs(pos=pos[(nprobes+1):(nprobes+nstdprobes),])
      tmp=first_guess_angs(pos=pos[(nprobes+1):(nprobes+nstdprobes),],cegs=cegs, guide=FALSE)
      angs=c(angs[1:nprobes],tmp)
      probes=defineprobe(x=pos[,'x'],y=pos[,'y'],angs=angs, interactive=FALSE)
      conflicts=find_probe_conflicts(probes=probes, pos=pos, angs=angs)
      nconflicts=dim(conflicts)[1]
      gprobes=defineguideprobe(gxs=gpos[,'x'],gys=gpos[,'y'],gangs=gangs, interactive=FALSE)
      gconflicts=find_guide_conflicts(probes=probes, pos=pos,angs=tempangs,gprobes=gprobes, gpos=gpos, gangs=tempgangs)
      ngconflicts=dim(gconflicts)[1]
      
      tempgangs<<-gangs
      tempangs<<-angs
      gangs<<-tempgangs
      angs<<-tempangs
      #Creating an output file for the robot:
      #Measure the angle corresponding to the shortest distance between the target and the edge of the field.
      azAngs<<-(pi+atan2(pos[,'y'],pos[,'x']))
      azAngs[azAngs > 3*pi/2]<<-azAngs[azAngs > 3*pi/2]-2*pi
      azAngs[azAngs < -1*pi/2]<<-azAngs[azAngs < -1*pi/2]+2*pi
      angs_azAng<<-angs-azAngs
      angs_azAng[angs_azAng > pi] <<- angs_azAng[angs_azAng > pi] - 2*pi
      angs_azAng[angs_azAng < 0] <<- angs_azAng[angs_azAng < 0] + 2*pi
      rads<<-sqrt(pos[,'x']**2+pos[,'y']**2)
      
      hexaconfig[hexaconfig[,'fibre_type'] == 'P',c('angs','azAngs','angs_azAng','rads')]<<-cbind(angs,azAngs,angs_azAng,rads)
      #print(hexaconfig)
      
      gazAngs<<-(pi+atan2(gpos[,'y'],gpos[,'x']))
      gazAngs[gazAngs > 3*pi/2]<<-gazAngs[gazAngs > 3*pi/2]-2*pi
      gazAngs[gazAngs < -1*pi/2]<<-gazAngs[gazAngs < -1*pi/2]+2*pi
      gangs_gazAng<<-gangs-gazAngs
      gangs_gazAng[gangs_gazAng > pi] <<- gangs_gazAng[gangs_gazAng > pi] - 2*pi
      gangs_gazAng[gangs_gazAng < 0] <<- gangs_gazAng[gangs_gazAng < 0] + 2*pi
      grads<<-sqrt(gpos[,'x']**2+gpos[,'y']**2)
      guideconfig[,c('gangs','gazAngs','gangs_gazAng','grads')]<<-cbind(gangs,gazAngs,gangs_gazAng,grads)
      #print(guideconfig)
      
      output$plotfield <- renderPlot(expr={
        plot_configured_field(pos=pos,angs=angs,gpos=gpos,gangs=gangs,fieldflags=fieldflags,aspdf=FALSE)
      },width = 650,
      height = 650,
      res = 72)
    }
  )
  
  output$downloadHexa <- downloadHandler(
    filename = function() {
      paste0(hfile$name)
    },
    content = function(file) {
      # Do not prepend a second probe column when re-saving a file produced
      # by the package or by this app. The old code always added probe from
      # row names, so a second download changed the schema to probe/probe.1
      # and made the file awkward for legacy readers to reload.
      probe_names <- grep("^probe(\\.[0-9]+)?$", names(hexaconfig), value = TRUE)
      if (length(probe_names) == 0) {
        hex_out <- data.frame(probe = rownames(hexaconfig), hexaconfig,
                              row.names = NULL, check.names = FALSE)
      } else {
        keep_names <- c(probe_names[1],
                        setdiff(names(hexaconfig), probe_names))
        hex_out <- hexaconfig[, keep_names, drop = FALSE]
        names(hex_out)[1] <- "probe"
      }
      write.table(hex_out, file, row.names = FALSE, sep = ',',
                  quote = FALSE, append = TRUE)
    }
  )
  output$downloadGuide <- downloadHandler(
    filename = function() {
      paste0(gfile$name)
    },
    content = function(file) {
      guideconfig[,'angs'] = guideconfig[,'gangs']
      guideconfig[,'rads'] = guideconfig[,'grads']
      guideconfig[,'azAngs'] = guideconfig[,'gazAngs']
      guideconfig[,'angs_azAng'] = guideconfig[,'gangs_gazAng']
      # Now get rid of the g columns
      drops <- c("gangs", "grads", "gazAngs", "gangs_gazAng",
                 grep("^probe(\\.[0-9]+)?$", names(guideconfig), value = TRUE))
      final_guide_table <- guideconfig[, !(names(guideconfig) %in% drops),
                                       drop = FALSE]
      write.table(final_guide_table, file, row.names = FALSE, sep = ',',
                  quote = FALSE, append = TRUE)
    }
  )
}

app <- shinyApp(ui, server)
