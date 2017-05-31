## Produces coloured plot

tip.colors <- function (tip.label) {
  dreary <- red <- dull <- lime <- indigo <- ruby <- jade <- straw <- fuschia <- c()
  dreary <- c('Priapulida', 'Ottoia_prolifica', 'Modern_priapulid', 'Cycloneuralia', 'Tubiluchus_Priapulida', 'Cricocosmia')
  
  bronze <- c('Leanchoilia', 'Alalcomenaeus', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Supella_longipalpa')
  straw <- c(
    c('Hadranax', 'Kerygmachela', 'Pambdelurion'),
    c('Siberion', 'Megadictyon', 'Jianshanopodia')
  )
  gold <- c(
    c('Anomalocaris', 'Opabinia', 'Peytoia', 'Laggania', 'Hurdia'),
    c('Eoredlichia', 'Fuxianhuia', 'Euarthropoda', 'Schinderhannes', 'Chengjiangocaris', 'Lyrarapax_unguispinus', 'Pambdelurion_whittingtoni', 'Opabinia_regalis', 'Anomalocaris_canadensis', 'Peytoia_nathorsti', 'Hurdia_victoria','Aegirocassis_benmoulae', 'Lyrarapax_unguispinus', 'Schinderhannes_bartlesi'),
    straw
  )
  racing <- c('Actinarctus_(Heterotardigrada)', 'Macrobiotus_(Eutardigrada)', 'Actinarctus_Heterotardigrada',     'Halobiotus_crispae_Eutardigrada', 'Halobiotus_Eutardigrada', 'Macrobiotus_Eutardigrada', 'Actinarctus_', 'Tardigrada', 'Macrobiotus_', 'Siberian Orsten tardigrade', 'Siberian_Orsten_tardigrade')
  jade <-  'Onychodictyon_ferox'
  galazios <- NULL
  hallucishaniids <- c('Hallucigenia_sparsa', 'Hallucigenia_fortis', 'Hallucigenia_hongmeia', 'Carbotubulus', 'Luolishania', 'Miraluolishania', 'Collins_monster_Emu_Bay', 'Collins_monster_Burgess_Shale', 'Acinocrinus', 'Collinsium', 'Newlishania')
  aquamarine <- c('Euperipatoides_Onychophora', 'Euperipatoides_(Onychophora)', 'Ooperipatellus_Onychophora', 'Plicatoperipatus_Onychophora', 'Tertiapatus_dominicanus')
  onies <- c(aquamarine, 'Antennacanthopodia', 'Cardiodictyon', 'Microdictyon' ,'Paucipodia', 'Xenusion', 'Diania', 'Onychodictyon_gracilis', 'Ilyodes', 'Helenodora', 'Tritonychus_phanerosarkus', 'Tritonychus', 'Orsteny_Long_Legs', 'Orstenotubulus')
  black <- c('Aysheaia', 'Onychodictyon')
  tip.col <- rep('black', length(tip.label))
  tip.col[tip.label %in% straw] <- '#ac9c16'
  tip.col[tip.label %in% gold] <- '#a2751c'
  tip.col[tip.label %in% bronze] <- '#a2751c'
  tip.col[tip.label %in% lime] <- '#99be16'
  tip.col[tip.label %in% racing] <- '#25ac89'
  tip.col[tip.label %in% ruby] <- '#c02349'
  tip.col[tip.label %in% fuschia] <- '#c82298'
  tip.col[tip.label %in% red] <- '#e1001b'
  tip.col[tip.label %in% hallucishaniids] <- '#a62bc5'
  tip.col[tip.label %in% black] <- '#000000'
  tip.col[tip.label %in% dreary] <- '#566666'
  tip.col[tip.label %in% dull] <- '#8e8e8e'
  tip.col[tip.label %in% jade] <- '#25ac89'
  tip.col[tip.label %in% galazios] <- '#009bdd'
  tip.col[tip.label %in% aquamarine] <- '#004d98'
  tip.col[tip.label %in% onies] <- '#aa6cb9'
  tip.col
}

doplot <- function (tr, pdf = FALSE, direction = 'rightwards', font=font, plotw=3, ploth=2.5,
  pts=10, ec=0, bi=FALSE, annot=FALSE, bi.nudge=0, col.factor=1, brightest = 0.9, filename='plot/plot.pdf',
  tip.col, fig=FALSE) {
  #REQUIRE tr, a phylo object.
  if (pdf) {
    pdf(file=paste(filename, "pdf", sep="."), width=plotw, height=ploth, pointsize=pts, colormodel='rgb')
  }
  tip.label <- tr$tip.label
  nTip <- length(tip.label)  

  nTax <- length(tr$tip.label)
  nNode <- tr$Nnode
  if (fig == FALSE && file.exists('lobo.stats.csv')) {
    stats <- loadstats('lobo.stats.csv') 
    tree.stats <- stats[tip.label,]
    inapp<-NULL
    ambig<-NULL
    ambig[1:nTax] <- tree.stats[1:nTax,'amb']
    for (i in (nTax+nNode):(nTax+1)) {
      ambig[i] <- mean(ambig[weightmean(tr, i)])
    }
    ambigs <- ambig[tr$edge[,2]]
    ambigs <- ambigs - min(ambigs)
    #edge.col <- function (x) heat.colors(16^3)[(16^3)*(x)^col.factor*brightest] #, start=0, end=1
    edge.col <- function (x) {X<-x^col.factor*brightest; rgb(X,X,X)} #, start=0, end=1
    ec <- edge.col(ambigs/max(ambigs))
  } else ec <- 'black'
      
  if (fig == FALSE) {
    tip.col = tip.colors(tip.label)
  } else {
    tip.col = 'black';
  }
      
  if (bi) label.offset <- 2*min(tr$edge.length)
  par(cex=0.8)  # Character expansion
  if (direction == 'rightwards') {
    align = 0
  } else {
    align = 0.5
  }
  plot(tr,
    edge.color = ec,
    edge.width = 2,
    font = font,
    cex = 1,
    tip.col = tip.col,
    adj = align,
    label.offset = 0.5,
    use.edge.length = bi,
    direction = direction,
    no.margin = TRUE,
    root.edge=TRUE,
    underscore=TRUE
  )
  if (!is.null(attr(tr, 'pscore'))) legend(1,1,attr(tr, 'pscore'))
   if (direction == 'rightwards' & !bi & fig == FALSE & file.exists('lobo.stats.csv')) {
    ioffset<-nTax-4.5; aoffset<-nTax-0; aheight=5
    max(ambig)
    amarks <- c("100%", "75%", "50%", "25%", "0%")
    for (i in seq(0,1,by=0.005)) {
      rect(xleft = 0, xright=1, ybottom=aoffset-aheight+(aheight*.8*i), ytop = aoffset-aheight+0.025+(aheight*.8*i), col=edge.col(i), border=NA)
    }
    text (0, aoffset,     "Tokens ambiguous or inapplicable:", adj=0)
    text (1.3, seq(aoffset-1, aoffset-aheight, length.out=length(amarks)), amarks, adj=0)
    nudgel <- 1.3
  }
  if (annot) {
    labex <- regexpr("([0-9]+)", tr$node.label); 
    lablen <- attr(labex, 'match.length')
    lab <- ifelse(lablen>0 & lablen <3, substr(tr$node.label, labex, labex+lablen-1), " ")
    nodelabels(lab, adj=c(nudgel + bi.nudge,-0.5), frame='none', cex=0.8)
  }
  
  #ioffset<-nTip-4.5; aoffset<-5
  #for (i in seq(0,1,by=0.01)) {
  #  col = (i^col.factor)*brightest
  #  rect(xleft = 0, xright=1, ybottom=aoffset-3.5+(2.8*i), ytop = aoffset-3.45+(2.8*i), col=rgb(col,col,col), border=NA)
  #}
  #text (0, ioffset,     "Inapplicable data (–)", adj=0)
  #text (1.3, ioffset-1:3, c("90%", "50%", "10%"), adj=0)
  #segments (0, ioffset-0.08-1:3, 1, lwd = c(10,50,90)/width.factor)
  #text (0, aoffset,     "Inapplicable data", adj=0)
  #text (1.3, aoffset-1:3, c("100%", "50%", "0%"), adj=0)
  
  
  if (pdf) {
    dev.off()
  }
}


colplot <- function (tr, taxnames='', direction='rightwards', ec=0, ...) {
  tr1 <- tr
  tip.label <- tr$tip.label
  nTip <- length(tip.label)
  nNode <- tr$Nnode
  # Taxon names and nodes
  roman <- c('Tardigrada', 'Onychophora', 'Priapulida', 'Collins monster', 'Collins', 
  'Collins_monster', 'Collins_monster_Emu_Bay', 'Siberian Orsten tardigrade', 'Siberian_Orsten_tardigrade'
  , 'Modern_priapulid')
  bold <- c('Onychophora', 'Hallucigenia', 'Aysheaia', 'Tardigrada', 'Orstenotubulus',
  'Hallucigenia_sparsa', 'Hallucigenia sparsa', 'Eutardigrada', 'Heterotardigrada',
  'Peripatus_', 'Peripatus_(Onychophora)', 'Actinarctus_', 'Actinarctus_(Heterotardigrada)', 'Macrobiotus_', 'Macrobiotus_(Eutardigrada)',
  'Peripatus ', 'Euperipatoides_Onychophora', 'Peripatus_Onychophora', 'Actinarctus ', 'Actinarctus_Heterotardigrada', 'Macrobiotus ', 'Macrobiotus_Eutardigrada',
  'Peripatus ', 'Peripatus (Onychophora)', 'Actinarctus ', 'Actinarctus (Heterotardigrada)', 'Macrobiotus ', 'Macrobiotus (Eutardigrada)')
  bold <- c('Tritonychus', 'Tritonychus_phanerosarkus')
  extinct <- c('Permochiton', 'Pedanochiton', 'Plumulites', 'Lepidocoleus', 'Turrilepas')
  tip.col <- tip.colors(tip.label)
 
  for (tax in names(taxnames)) {
    taxa <- taxnames[[tax]]
    tr <- drop.tip(tr, which(tr$tip.label%in%taxnames[[tax]]), subtree=TRUE)
    new.clade.name <- paste(tax, " (", length(taxnames[[tax]]), ")", sep="")
    tr$tip.label[length(tr$tip.label)] <- new.clade.name
    roman <- c(roman, new.clade.name)
  }

  nTip <- length(tr$tip.label)
  nNode <- tr$Nnode
  font <- (!(tr$tip.label %in% roman))*2+1+(tr$tip.label %in% bold)
  doplot(tr, direction=direction, font=font, ec=ec, tip.col=tip.col, ...)
}