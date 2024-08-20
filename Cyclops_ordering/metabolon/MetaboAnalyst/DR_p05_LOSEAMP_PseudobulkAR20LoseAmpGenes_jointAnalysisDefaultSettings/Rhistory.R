# PID of current job: 4045644
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-SetOrganism(mSet, "hsa")
geneListFile<-"replace_with_your_file_name"
geneList<-readChar(geneListFile, file.info(geneListFile)$size)
mSet<-PerformGeneMapping(mSet, geneList, "hsa", "symbol");
cmpdListFile<-"replace_with_your_file_name"
cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
mSet<-PerformCmpdMapping(mSet, cmpdList, "hsa", "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PrepareIntegData(mSet);
mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "integ", "query");
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA, NA, NA )
mSet<-CreateIntegMatchingTable(mSet);
mSet<-PlotKEGGPath(mSet, "Purine metabolism",566, 490, "png", NULL)
mSet<-RerenderMetPAGraph(mSet, "zoom1721225573242.png",566.0, 490.0, 100.0)
mSet<-SaveTransformedData(mSet)
