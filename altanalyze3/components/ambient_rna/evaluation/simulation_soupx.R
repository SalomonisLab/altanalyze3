# Execute unmodified, pinned official SoupX functions; no Python approximation.
# Matrix is required. Seurat is unnecessary because Matrix Market input and
# externally calculated cellHarmony clusters replace the optional 10x loader.
args = commandArgs(trailingOnly=TRUE)
directory = args[1]; sourceDirectory = args[2]; rho = as.numeric(args[3])
seed = as.integer(args[4])
suppressPackageStartupMessages(library(Matrix))
for (f in c('classFunctions.R','utils.R','estimateSoup.R','setProperties.R',
            'quickMarkers.R','estimateNonExpressingCells.R','autoEstCont.R',
            'adjustCounts.R')) source(file.path(sourceDirectory,f))
meta = read.delim(file.path(directory,'cells.tsv'),check.names=FALSE)
genes = readLines(file.path(directory,'genes.txt'))
toc = as(readMM(file.path(directory,'contaminated.mtx')),'CsparseMatrix')
empty = as(readMM(file.path(directory,'empty.mtx')),'CsparseMatrix')
rownames(toc)=genes; colnames(toc)=meta$cell
rownames(empty)=genes; colnames(empty)=paste0('SIM_EMPTY_',seq_len(ncol(empty)))
start=proc.time()[['elapsed']]
sc=SoupChannel(cbind(toc,empty),toc)
sc=setClusters(sc,setNames(meta$cluster,meta$cell))
write.table(data.frame(gene=rownames(sc$soupProfile),sc$soupProfile),
            file.path(directory,'estimated_profile.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
profileSeconds=proc.time()[['elapsed']]-start
set.seed(seed)
for (method in c('SoupX_known','SoupX_auto')) {
  path=file.path(directory,paste0(method,'.mtx'))
  statusPath=file.path(directory,paste0(method,'_status.tsv'))
  if(file.exists(statusPath)) next
  begin=proc.time()[['elapsed']]
  warnings=character()
  result=tryCatch(withCallingHandlers({
    selected=if(method=='SoupX_known')setContaminationFraction(sc,rho)else
      autoEstCont(sc,doPlot=FALSE,verbose=TRUE)
    selectedRho=mean(selected$metaData$rho)
    corrected=adjustCounts(selected,roundToInt=TRUE,verbose=0)
    computeSeconds=proc.time()[['elapsed']]-begin
    writeMM(corrected,path)
    data.frame(method=method,status='success',rho=selectedRho,
               correction_seconds=computeSeconds,profile_seconds=profileSeconds,error='')
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e)data.frame(method=method,status='failed',rho=NA_real_,
    correction_seconds=proc.time()[['elapsed']]-begin,profile_seconds=profileSeconds,
    error=conditionMessage(e)))
  result$warnings=paste(warnings,collapse=' | ')
  write.table(result,statusPath,sep='\t',row.names=FALSE,quote=TRUE)
}
writeLines(c(capture.output(sessionInfo())),file.path(directory,'R_sessionInfo.txt'))
