#' rrnDBcorrectOTU
#' correct 16S rRNA gene copy numbers by rrnDB database for a OTU table
#' @param otu is a OTU table
#' @param classifer is RDP classifer result file
#' @param rrnDB is rrnDB database by RDP format, which can be downloaded in https://rrndb.umms.med.umich.edu/static/download/


rco = function(otu,classifer,rrnDB){

  hang1 = length(rownames(otu))
  #rrnDB只到class;
  #a为genus,family,order,class
  #b比a低一级
  #从genus到class,即从低到高遍历levels
  a= c(6,5,4,3)
  b= c(7,6,5,4)
  levels= c("genus","family","order","class")
  whole.res = c()

  for (q in 1:4){  #q=1
    spe = classifer[(classifer[,a[q]] != "Unclassified") &  (classifer[,b[q]] == "Unclassified"),]
    rspe = rrnDB[ rrnDB[,2]==levels[q],c(3,9)]

    name.spe = rownames(spe)
    name.otu = rownames(otu)
    name.match = match(name.spe,name.otu)
    match.otu = otu[name.match,]
    #注释信息加到otu中
    otu.spe = cbind(spe[,a[q]],match.otu)

    ##############################
    #用for 循环遍历
    #规定OTU中有而数据库中没有的为-1
    whole = data.frame(matrix(data=NA,
                              nrow=length(rownames(otu.spe)),
                              ncol=length(colnames(otu.spe))),stringsAsFactors=FALSE)

    for (i in 1:length(rownames(otu.spe))){  #i=1
      mat = match(as.character(otu.spe[i,1]),as.character(rspe[,1]))
      if(!is.na(mat)){
        each = cbind(rspe[mat,2],otu.spe[i,-1])
      } else{
        each = cbind(-1,otu.spe[i,-1])
      }
      whole[i,] = each
    }
    #加上物种信息
    spe.res = cbind(taxa=otu.spe[,1],level=rep(levels[q],length(rownames(otu.spe))),whole)
    rownames(spe.res) = rownames(otu.spe)

    #输出到whole.res
    whole.res = rbind(whole.res,spe.res)
    ##############每一个level做完之后都需要将已做过的OTU去掉。加快后面的速度。
    ##############最重要的是避免“有物种信息；unclassifiered；有物种信息”这种情况出现。保证有信息的是最小的level，结果最准确。
    diffotu = setdiff(rownames(otu),rownames(otu.spe))
    otu = otu[diffotu,]
  }
  whole.res = whole.res[complete.cases(whole.res),]
  colnames(whole.res) = c("taxa","level","CopyNumber",colnames(otu))

  #检验一下whole.res是否和OTU表的行数一样
  #结果的行数+属水平数据库空缺的行数+class做完之后还剩下来的OTU
  hang2 = length(rownames(whole.res))+length(rownames(otu))

  if (hang1==hang2){
    print("Well done!")}

  #######################OTU丰度除以拷贝数，得到校正的表格
  #先去掉第三列为-1的值
  simply.res = whole.res[whole.res[,3]!= -1,]
  correct = simply.res[,-c(1:3)] / simply.res[,3]
  correct.table = as.data.frame(cbind(simply.res[,1:2],correct))

  list(whole.res=whole.res,correct.table=correct.table)
}
