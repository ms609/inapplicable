data(inapplicable.datasets)
scores <- c(
"Agnarsson2004" =  778 ,
"Aguado2009" =     579 ,
"Aria2015" =       143 ,
"Asher2005" =      345 ,
"Capa2011" =       385 ,
"Conrad2008" =     1761,
"DeAssis2011" =    64  ,
"Dikow2009" =      1611,
"Eklund2004" =     440 ,
"Geisler2001" =    1295,
"Giles2015" =      710 ,
"Griswold1999" =   407 ,
"Liljeblad2008" =  2868,
"Loconte1991" =    539 ,
"Longrich2010" =   131 ,
"OLeary1999" =     508 ,
"OMeara2014" =     273 ,
"Rougier2012" =    1215,
"Rousset2004" =    259 ,
"Sano2011" =       223 ,
"Sansom2010" =     189 ,
"Schulze2007" =    164 ,
"Shultz2007" =     454 ,
"Vinther2008" =    79  ,
"Wetterer2000" =   559 ,
"Wills2012" =      273 ,
"Wilson2003" =     879 ,
"Wortley2006" =    482 ,
"Zanol2014" =      1311,
"Zhu2013" =        638 )

nj.tree <- lapply(inapplicable.phyData, NJTree)


install_github('ms609/inapplicable', rel='cefb5669352aca6425516805f60108063383b6c2')

profvis::profvis(
for (dataset in names(inapplicable.phyData)) {
  cat("\n\n\n ======== NEXT DATASET: ", dataset, "========\n\n\n"
  oTree <- RatchetSearch(nj.tree[[dataset]], inapplicable.phyData[[dataset]], stopAtScore=scores[[dataset]],
  k=1000, maxIt=10000, maxIter=3200, maxHits=12)
}
)
                   