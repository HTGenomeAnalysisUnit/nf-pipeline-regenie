#!/usr/bin/env Rscript

# Usage : Rscript --vanilla pheno_chunker.R -f phenotype_file -m model_file -t maximum_proportion_of_missing -s maximum_models_per_chunk



library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="tab separated phenotype file", metavar="xxx.tsv"),
  make_option(c("-m", "--model_file"), type="character", default="models.tsv", 
              help="tab separated model table file [default= %default]", metavar="xxx.tsv"),
  make_option(c("-t", "--tollerance"), type="numeric", default=0.15, 
              help="maximum proportion of missing samples per phenotype", metavar="Number [0-1]"),
  make_option(c("-s", "--maximum_chuck_size"), type="numeric", default=10, 
              help="model table file [default= %default]", metavar="Number [>1]"),
  make_option(c("-a", "--sample_file"), type="character", default=NA, 
              help=".fam/.sample file", metavar="xxx.sample xxx.fam")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


extract.cov=function(x)unlist(strsplit(as.character(x),split="~",fixed=T))[2]|>gsub(pattern = " ",replacement = "")

pheno.chunker=function(pheno_file="data/test_phenotype.tsv" #names of the phenotype file
                      ,model.file="models.tsv"              # model file
                      ,tollerance=0.15                      # maximum number of missing samples per phenotype
                      ,max.chunk.size=10.                   # maximum number of phenotypes per chunck
                      ,fam.file="data/test.fam"
                      ){
   
   require(data.table)

   models=fread(model.file)
   variables=strsplit(models$model,split="~",fixed=T)|>unlist()|>
     strsplit(split="+",fixed=T)|>unlist()|>
     gsub(pattern = " ",replacement = "")|>unique()
   variables=variables[!(variables%in%1)]
      
   pheno=fread(pheno_file,select = c("IID",variables))
   if(all(variables%in%colnames(pheno))){
     print(paste("Loaded variables:",paste(variables,collapse=" ")))
   }else{
     miss=paste(variables[!(variables%in%colnames(pheno))],collapse=",")
     stop(paste("Missing variables in model file:",miss))
     
   }
   row.names(pheno)=pheno$IID
   if(!is.na(fam.file)){
     
     if(length(grep("\\.fam$",fam.file))>0){
       fam=fread(fam.file,header=F)
     }else if(length(grep("\\.sample$",fam.file))>0){
       fam=fread(fam.file,header=F,skip=2)
     }
     names(fam)[1:2]=c("FID","IID")
     fam=fam[match(pheno$IID,fam$IID),]
     pheno=data.table(FID=fam$FID,pheno)
   }else{
     
     pheno=data.table(FID=pheno$IID,pheno)
     
     
   }
   ## Create chunks by model and covariates
   covariates=apply(t(models$model),2,extract.cov) 
   covariates=as.factor(covariates)|>as.numeric()|>as.character()
   models$lab_2=paste(models$trait_type,covariates,models$genetic_model,sep="_")
   models$over.lab=NA
   
   
   
   
   ### Create chunks by missingess
   for(k in unique(models$lab_2)){
     tmp.models=models[models$lab_2==k,]
     pheno.is.na=matrix(nrow=nrow(pheno),ncol=nrow(tmp.models),data = 0)
     rownames(pheno.is.na)=row.names(pheno)
     colnames(pheno.is.na)=tmp.models$model_id
     model.id.list=list()
     for(i in 1:nrow(tmp.models)){
       a=model.frame(as.formula(tmp.models$model[i]),data = pheno)
       pheno.is.na[row.names(a),tmp.models$model_id[i]]=1
       model.id.list[[i]]=row.names(a)
     }
     names(model.id.list)=tmp.models$model_id
     
     overlap.matrix=matrix(ncol=nrow(tmp.models),nrow=nrow(tmp.models))
     rownames(overlap.matrix)=tmp.models$model_id
     colnames(overlap.matrix)=tmp.models$model_id
     sum.vec=colSums(pheno.is.na)
     for(i in 1:ncol(overlap.matrix)){
       for(j in i:ncol(overlap.matrix)){
         val= pheno.is.na[,i]+pheno.is.na[,j]
         val[val==1]=0
         val[val==2]=1
         val=sum(val)/max(c(sum.vec[i],sum.vec[j]))
         overlap.matrix[i,j]=overlap.matrix[j,i]=val
       }
     }   
     if(all(overlap.matrix>(1-tollerance))){
         
         models$over.lab[models$lab_2==k]=1
         
     }else{
        over.perc=0
        over.clust=hclust(as.dist(1-overlap.matrix),method="ward.D2")
        g.n=2
        while(any(over.perc<(1-tollerance)) & g.n<=ncol(overlap.matrix)){
          tmp.gr=cutree(over.clust,k=g.n)   
          over.perc=c()
          for(gr.n in unique(tmp.gr)){
            over.perc=c(over.perc,
                         min(overlap.matrix[names(tmp.gr)[tmp.gr==gr.n],
                                            names(tmp.gr)[tmp.gr==gr.n]]))
               
            
          }
          g.n=g.n+1
          
        }
        n.groups=tmp.gr
        
        models$over.lab[models$lab_2==k]=n.groups[models$model_id[models$lab_2==k]]
       
         
           
         
     }
    }
    models$run_group=paste(models$lab_2,models$over.lab,sep="_") 
    
  ### create maximum trait chunk 
    
    models$run_group=as.numeric(factor(models$run_group,levels=unique(models$run_group)))
    models$size_group=NA
    for(i in models$run_group){
      
      models$size_group[models$run_group==i]=ceiling(1:sum(models$run_group==i)/max.chunk.size)
      
    }
   models$run_group=paste(models$run_group,models$size_group,sep="_")
   models$run_group=as.numeric(factor(models$run_group,levels=unique(models$run_group)))
    
   ##### Create files
   chunk.id.list=list()
   n.miss.chunk=c()
   sample.size=c()
   chunck.model=c()
   chunk.genetic.model=c()
   cov_file=c()
   for(i in unique(models$run_group)){
     tmp.models=models[models$run_group==i,]
     model.id.list=c()
     outcomes=c()
     cova=c()
     chunk.genetic.model=c(chunk.genetic.model,unique(tmp.models$genetic_model))
     
     for(j in 1:nrow(tmp.models)){
       a=model.frame(as.formula(tmp.models$model[j]),data = pheno)
       model.id.list=unique(c(model.id.list,row.names(a)))
       outcomes=c(outcomes,names(a)[1])
       if(any(extract.cov(tmp.models$model[j])!=1)){
         cova=unique(c(cova,names(a)[-1]))
       }
     }
     
     idx.pheno=c("FID","IID",outcomes)
     pheno.matr=pheno[row.names(pheno)%in%model.id.list,..idx.pheno]
     
     n.miss.model=round(apply(pheno.matr[,-c(1,2)],2,FUN = function(x)length(which(is.na(x)))/length(x)),3)
     n.miss.chunk=c(n.miss.chunk,paste(paste(names(n.miss.model),n.miss.model,sep=":"),collapse=","))
     chunck.model=c(chunck.model,paste(tmp.models$model_id,collapse=","))
     dir.create(paste0("chunk_",i))
     
     if(any(extract.cov(tmp.models$model[j])!=1)){
        idx.cova=c("FID","IID",cova)
        cov.matr=pheno[row.names(pheno)%in%model.id.list,..idx.cova]
        fwrite(cov.matr,file=paste("chunk_",i,"/chunk_",i,".cov",sep=""),sep="\t")
        cov_file=c(cov_file,paste("chunk_",i,"/chunk_",i,".cov",sep=""))
     }else{
        cov_file=c(cov_file,"NO_COV_FILE")
     }  
   
     sample.size=c(sample.size,nrow(pheno.matr))
     trait.names=strsplit(tmp.models$model,split="~",fixed=T)
     trait.names=unlist(lapply(trait.names,FUN = function(x)x[1]))|>gsub(pattern = " ",replacement = "")
     names(pheno.matr)[-1]=tmp.models$model_id[match(names(pheno.matr)[-1],trait.names)]
     fwrite(pheno.matr,file=paste("chunk_",i,"/chunk_",i,".pheno",sep=""),sep="\t")
     
   }
   master.table=unique(models[,c("run_group","trait_type")])
   master.table$pheno_file=paste("chunk_",master.table$run_group,"/chunk_",master.table$run_group,".pheno",sep="")
   master.table$cov_file=cov_file
   master.table$genetic_model=chunk.genetic.model
   master.table$run_group=paste0("chunk_",master.table$run_group)
   master.table$sample_size=sample.size
   master.table$model_id=chunck.model
   master.table$perc_miss=n.miss.chunk
   write.table(master.table,file="master_table.tsv",row.names=F,quote=F,sep="\t")
   
}  
      
   
   
pheno.chunker(pheno_file = opt$file,model.file = opt$model_file,
              tollerance = opt$tollerance,
              max.chunk.size = opt$maximum_chuck_size
              , fam.file=opt$sample_file)   
   
   
   












