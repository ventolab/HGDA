# withr::with_libpaths(new = "/home/jovyan/R/lib/", install_github(repo = 'satijalab/seurat', ref = 'release/3.0'))
# source("https://bioconductor.org/biocLite.R")

# HIST1H4C TOP2A CDK1
#ggsave to get svg
library(base)
library(methods)
options(width = 175)
#if (sum(grepl("_64-pc-linux-gnu-library",.libPaths())) == 0){
#	if (base::R.version$minor == "4.2") .libPaths("/nfs/users/nfs_l/lh20/R/x86_64-pc-linux-gnu-library/3.4")
#	if (base::R.version$minor == "5.3") .libPaths("/nfs/users/nfs_l/lh20/R/x86_64-pc-linux-gnu-library/3.5")
#	if (base::R.version$minor == "6.2") .libPaths("/nfs/users/nfs_l/lh20/R/x86_64-pc-linux-gnu-library/3.6")
#}
library(ggplot2)
library(reticulate)
use_python("*")

parselist <-function(string){return(strsplit(gsub(",,",",",gsub("\n", ",",gsub(" ","",string,fixed=T))),",")[[1]])}
printlist <- function(stringlist, append.comma=F, add.quotes=F){
	if (append.comma){
		if (add.quotes == F) for(i in 1:length(stringlist)) cat(paste(stringlist[i],"\n",sep=","))
		else for(i in 1:length(stringlist)) cat(paste("",stringlist[i],",",sep="\""))
	}else{
		for(i in 1:length(stringlist)) cat(paste(stringlist[i],"\n",sep=""))
	}
}

# use readRDS and notify on commandline
noisyreadRDS <- function(basepath, file){
	dapath <- paste(basepath, file, sep="")
	cat(paste("Opening file", dapath,"\n"))
return(readRDS(dapath))}

# use readRDS and notify on commandline
noisysaveRDS <- function(object,basepath, file){
	dapath <- paste(basepath, file, sep="")
	cat(paste("Saving file", dapath,"\n"))
	saveRDS(object,dapath)
}

# use make.names(..., unique=TRUE)
mapFactor <-function(fact, mapping){return(factor(mapping[fact@.Data], levels=unique(mapping)))}
permuteFactor <- function(fact, permutation, extralevels = c() ){
	dad <- setdiff(1:length(levels(fact)), permutation)
	if (length(dad) > 0){
		print("Missing entries in permutation:")
		print(dad)
	}
	if (length(unique(permutation)) != length(permutation)) print("permutation has duplicate entries")
return(factor(as.character(fact), levels=c(levels(fact), extralevels)[permutation]))

}
getSafeLevels <- function(fact.or.characters){
	if (class(fact.or.characters) == "factor") return(levels(fact.or.characters))
	else return(unique(fact.or.characters))
}


filterAndSort <- function(table, flt = c(), ord = c(), decreasing=F){
	if (!is.null(flt) ) {
		filter <- rep(T,nrow(table))
		for(j in 1:length(flt)) {
			if      (flt[[j]][1] == "<")  filter <- filter & (table[[names(flt)[j]]] < as.numeric(flt[[j]][2]))
			else if (flt[[j]][1] == ">")  filter <- filter & (table[[names(flt)[j]]] > as.numeric(flt[[j]][2]))
			else if (flt[[j]][1] == "<=") filter <- filter & (table[[names(flt)[j]]] <=as.numeric(flt[[j]][2]))
			else if (flt[[j]][1] == ">=") filter <- filter & (table[[names(flt)[j]]] >=as.numeric(flt[[j]][2]))
			else filter <- filter & (table[[names(flt)[j]]]  %in% flt[[j]])

		}
		tmptable <- table[filter,,drop=F]
	}
	if (!(ord %in% colnames(tmptable))) stop(paste(ord, "is not a collumn in table"))
return(tmptable[order(tmptable[[ord]],decreasing=decreasing),,drop=F])}


testRGL <- function(data,w=c(),x=c(),y=c(), transform=c(), use.rows=c(), use.cols=c()){
	library(rgl)
	if (class(data) != "list") {
		data <- list(data=data)
		if (is.null(w)) data$w <-matrix(1, dim(data$data)[1], dim(data$data)[2])
		else data$w <- w
	#	if (is.null(x)) data$x <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
	#	else data$x <- x

	#	if (is.null(y)) data$y <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
	#	else data$y <- y
	}else{
		if (!"data" %in% names(data)) stop("input list has 'data' as mandatory field")
		for( i in names(data)){
			if (!i %in% c("data", "x", "y", "w")) print(paste("The nknown field",i," is ignored, valid fields are \"data\", \"x\", \"y\", \"w\" only"))
		}
		
		if (!"w" %in% names(data)) data$w <- matrix(1, dim(data$data)[1], dim(data$data)[2])
		if (!"x" %in% names(data)) data$x <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		if (!"y" %in% names(data)) data$y <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
	}
	if (!is.null(use.rows)){
		if (!is.null(use.cols)){
			for(i in names(data)) data[[i]] <- data[[i]][use.rows,use.cols]
		}else{
			for(i in names(data)) data[[i]] <- data[[i]][use.rows,]
		}
	}else if (!is.null(use.cols)){
		for(i in names(data)) data[[i]] <- data[[i]][,use.cols]
	}	

	for(i in names(transform)){
		if (!i %in% c("data", "x", "y", "w")) stop("allowed transform fiels are \"data\", \"x\", \"y\", \"w\" only")
		if (transform[i] == "pval"){
			data[[i]] <- -0.8 + 0.09 / (0.05 + data[[i]])
			data[[i]][data[[i]] < 0.1] <- 0.1
		}else if (transform[i] == "log10pval"){
			data[[i]] <- 0.05 / (0.05 + exp(data[[i]]/log(10)))
		}else if (transform[i] == "lerfp1"){
			tmp <- as.matrix(log(data[[i]] + 1))
			m <- mean(tmp,na.rm=T)
			v <- var(as.vector(tmp),na.rm=T) ^-0.5
			if (is.na(v)){
				data[[i]] <- 0.5
			}else{
				data[[i]]  <- pnorm((tmp - m)*v)
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "colwise.erf"){
			trformval <- matrix(0, 5, dd[2])
			for(j in 1:dd[2]){
				tmp <- data[[i]][,j]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][,j] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][,j]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "rowwise.erf"){
			trformval <- matrix(0, 5, dd[1])
			for(j in 1:dd[2]){
				tmp <- data[[i]][j,]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][j,] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][j,]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "threshold"){
			flt <- (data[[i]] <= 0.05)
			flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 1;
			data[[i]][!flt] <- 0;
		}else if (transform[i] == "threshold2"){
			flt <- (data[[i]] <= 0.05)
			flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 0;
			data[[i]][!flt] <- 1;
		}

	}

	curdim= dim(data$data)

	#clear3d()
  rgl.init()
  view3d(0,0,scale=c(1/curdim[2],1/curdim[1],1), fov=0)
  winscp <<- list()

     # shade3d(icosahedron3d(), col = "yellow")
      #qmesh3d( matrix(runif(1200), 3,400), 1:400, homogeneous = TRUE, material = rgl.material(color="#FF0000", alpha=0.2), normals = NULL, texcoords = NULL)

	#source("~/Rcode.R"); da <- matrix(runif(16),4,4); colnames(da) <- paste("col",1:4); rownames(da) <- paste("row",1:4); testRGL(da)	

#qmesh3d(vertices, indices, homogeneous = TRUE, material = NULL,              normals = NULL, texcoords = NULL)
	winscp$root <<- currentSubscene3d()
	
	window <- par3d(c("viewport"), dev = rgl.cur(), subscene = winscp$root)

	par3d(viewport = c(window[3]/4,0,window[3]*3/4,3*window[4]/4), subscene = winscp$root);
	
	vdata <- matrix(0,4,4);vdata[4,] <- 1;vdata[3,] <- 1;
	vdata[1,1] <- 0;  vdata[1,2] <- 0; 
	vdata[1,3] <- curdim[2]; vdata[1,4] <- curdim[2];
	vdata[2,1] <- 0;  vdata[2,4] <- 0; 
	vdata[2,2] <- curdim[1]; vdata[2,3] <- curdim[1];


	shade3d(qmesh3d( vdata, matrix(1:4,4,1),material=list(front="culled",back="culled")), normals = NULL, texcoords = NULL)

	pos = (1:curdim[2]) - 0.5
	axes3d('x', labels= colnames(data$data), at=pos,adj=c(0.5,1))

	winscp$Yaxe <- newSubscene3d("replace", "replace", "replace", copyShapes = FALSE, parent = winscp$root,ignoreExtent=T)
	par3d(viewport = c(0,window[4]/4,window[3]*3/4,window[4]*3/4), subscene = winscp$Yaxe);
	winscp$scene <- newSubscene3d("replace", "replace", "replace", copyShapes = FALSE, parent = winscp$root,ignoreExtent=T)
	par3d(viewport = c(window[3]/4,window[4]/4,window[3]*3/4,window[4]*3/4), subscene = winscp$scene);
  useSubscene3d(winscp$Yaxe)

  shade3d(qmesh3d( vdata, matrix(1:4,4,1),material=list(front="culled",back="culled")), normals = NULL, texcoords = NULL)
  pos = (1:curdim[1]) - 0.5

  axes3d('y', labels= rownames(data$data), at=pos,adj=c(1,0.5))

  tmp <- par3d(c("userProjection"), dev = rgl.cur(), subscene = winscp$root)
  par3d(userProjection = translationMatrix(0, 1/3, 0) %*% tmp, dev = rgl.cur(), subscene = winscp$root)
  tmp <- par3d(c("userProjection"), dev = rgl.cur(), subscene = winscp$Yaxe)
  par3d(userProjection = translationMatrix(1/3, 0, 0) %*% tmp, dev = rgl.cur(), subscene = winscp$Yaxe)

	useSubscene3d(winscp$scene)
	
 	shade3d(qmesh3d( vdata, matrix(1:4,4,1),material=list(front="culled",back="culled")), normals = NULL, texcoords = NULL)

	vdata <- matrix(0,4,12*curdim[2]*curdim[1]+4)
	cdata <- rep("#000000",12*curdim[2]*curdim[1]+1)
	vdata[4,] <- 1;vdata[3,] <- 1
	data$data[is.na(data$data)] <- 0
	aur <-range(data$data)
	
	o <- curdim[2] * curdim[1] * 12+1;
	cdata[12*curdim[2]*curdim[1]+1] = "#EEEEEE"
	cdata[12*curdim[2]*curdim[1]+2] = "#EEEEEE"
	cdata[12*curdim[2]*curdim[1]+3] = "#EEEEEE"
	cdata[12*curdim[2]*curdim[1]+4] = "#EEEEEE"
	vdata[1,o  ] <- 0;
	vdata[1,o+1] <- 0;
	vdata[1,o+2] <- curdim[2];
	vdata[1,o+3] <- curdim[2];
	vdata[2,o  ] <- 0;
	vdata[2,o+1] <- curdim[1];
	vdata[2,o+2] <- curdim[1];
	vdata[2,o+3] <- 0;

	for(y in 1:curdim[2]){
		for(x in 1:curdim[1]){
			o <- (x + ((y-1) * curdim[1])) * 12 - 11;
			cdata[o:(o+11)] <- myFireAndIce((data$data[x,y] - aur[1])/(aur[2] - aur[1]))
		
			radata <- c(0.5,0.5) - c(1 - sqrt(data$w[x,y]), 1 - sqrt(data$w[x,y]/2)) /2

			vdata[1,o  ] <- y-radata[1] -0.5;
			vdata[1,o+1] <- y-radata[2] -0.5;
			vdata[1,o+2] <- y -0.5;
			vdata[1,o+3] <- y+radata[2] -0.5;
			vdata[1,o+4] <- y+radata[1] -0.5;
			vdata[1,o+5] <- y+radata[2] -0.5;
			vdata[1,o+6] <- y -0.5;
			vdata[1,o+7] <- y-radata[2] -0.5;
	
			vdata[2,o  ] <- x -0.5;
			vdata[2,o+1] <- x+radata[2] -0.5;
			vdata[2,o+2] <- x+radata[1] -0.5;
			vdata[2,o+3] <- x+radata[2] -0.5;
			vdata[2,o+4] <- x -0.5;
			vdata[2,o+5] <- x-radata[2] -0.5;
			vdata[2,o+6] <- x-radata[1] -0.5;
			vdata[2,o+7] <- x-radata[2] -0.5;

			vdata[,o+8] <- vdata[,o  ]
			vdata[,o+9] <- vdata[,o+3]
			vdata[,o+10] <- vdata[,o+4]
			vdata[,o+11] <- vdata[,o+7]
		}
	}

	shade3d(qmesh3d( vdata, matrix(1:(12*curdim[2]*curdim[1]+4),4,3*curdim[2]*curdim[1]+1), material =list(lit=F, color=cdata)), normals = NULL, texcoords = NULL)
	dev = rgl.cur()

	dapr <- rgl.projection()
	subscene = currentSubscene3d(dev)
        beginL <- function(x, y) {
	    init <- list()
	    print(par3d(c("zoom"), dev = dev, subscene = winscp$root))
            init$root <- par3d(c("userProjection","viewport"), dev = dev, subscene = winscp$root)
            init$Yaxe <- par3d(c("userProjection"), dev = dev, subscene = winscp$Yaxe)
            init$scene <- par3d(c("userProjection"), dev = dev, subscene = winscp$scene)
            init$pos <- c(x/(init$root$userProjection[1,1]*init$root$viewport[3]), 1 - y/(init$root$userProjection[2,2]*init$root$viewport[4]))
            winscp$daval <<- init
        }
        
        updateL <- function(x, y) {
          #for (sub in winscp$listeners) {
            init <- winscp$daval ; #[[as.character(sub)]]
            xlat <- 2*(c(x/(init$root$userProjection[1,1]*init$root$viewport[3]), 1 - y/(init$root$userProjection[2,2]*init$root$viewport[4])) - init$pos)
            mouseMatrix <- translationMatrix(xlat[1], 0, 0)
            par3d(userProjection = mouseMatrix %*% init$root$userProjection, dev = dev, subscene = winscp$root )

            mouseMatrix <- translationMatrix(0, xlat[2], 0)
            par3d(userProjection = mouseMatrix %*% init$Yaxe, dev = dev, subscene = winscp$Yaxe )

            mouseMatrix <- translationMatrix(xlat[1], xlat[2], 0)
            par3d(userProjection = mouseMatrix %*% init$scene, dev = dev, subscene = winscp$scene )
           #}
        }

        beginR <- function(x, y) {
            init <- par3d(c("userProjection","viewport"), dev = dev, subscene = winscp$root)
            init$Yaxe <- par3d(c("userProjection"), dev = dev, subscene = winscp$Yaxe)
            init$scene <- par3d(c("userProjection"), dev = dev, subscene = winscp$scene)

            init$pos <- c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5)
            winscp$daval <<- init
        }

        updateR <- function(x, y) {
            init <- winscp$daval
            xlat <- 2*(c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5) - init$pos)
         
            par3d(userProjection = scaleMatrix(exp(xlat[1]), 1.0, 1.0) %*% init$userProjection, dev = dev, subscene = winscp$root )
            par3d(userProjection = scaleMatrix(1.0, exp(xlat[2]), 1.0) %*% init$Yaxe, dev = dev, subscene = winscp$Yaxe )
            par3d(userProjection = scaleMatrix(exp(xlat[1]), exp(xlat[2]), 1.0) %*% init$scene, dev = dev, subscene = winscp$scene )

        }
	rotaR <- function(x){print(well);}
        rgl.setMouseCallbacks(1, beginL, updateL, dev = dev, subscene = subscene)
        rgl.setMouseCallbacks(2, beginR, updateR, dev = dev, subscene = subscene)
	rgl.setWheelCallback(rotate=rotaR, dev = dev, subscene = subscene)

   #   }
#panzoom2d()

}



cluster <- function(data, prefix){
	library(hclust)
	

}

makeThoseTables <- function(runargs, pair, wpair ){


summary <- list()
summary2 <- list()

geneannot <- read.csv("/lustre/scratch117/cellgen/team218/lh20/geneannotation.tsv",sep="\t")
genecolnames = c("Comparison", "Archtype", "Celltype", "NAME", "DE", "Log2FC", "LogitAuroc", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval","DESeq_basemean", "TPMmean", "FAD_dropout", "Ctrl_dropout", "Alias", "Description", "FullName", "GOslim", "GO")

genetable <- data.frame(row.names=genecolnames)
genetable2 <- data.frame(row.names=genecolnames)

gotable <- data.frame(row.names=c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_dropout", "Ctrl_dropout", "Intersection"))

   for(ss in pair$list.comparisons){
	print(paste("processing", ss))
	
	summary[[ss]] <- matrix(0,12,length(pair$list.celltypes))
	summary2[[ss]] <- matrix(0,12,length(pair$list.celltypes))
	colnames(summary[[ss]]) <- pair$list.celltypes
	rownames(summary[[ss]]) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")
	colnames(summary2[[ss]]) <- pair$list.celltypes
	rownames(summary2[[ss]]) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")

	for(ct in pair$list.celltypes){
	  dn <- gsub(" ", "_",paste(ct,ss))

	  for(i in names(runargs$cel$archtype)){
		if (ct %in% runargs$cel$archtype[[i]]) archt <- i
	  }


	  print(paste("Processing", ct))
	  tmp <- gsub(" ", "_",paste("TOP", ct))
	print(tmp)
	  if (tmp %in% names(pair$signifDE[[ss]])) {
		showlistQ <- pair$signifDE[[ss]][tmp][[1]][,"Corrected_deseq.log10pvalue",drop=F]; 
		showlistQ <- rbind(showlistQ,pair$signifDE[[ss]][ gsub(" ", "_",paste("BOT", ct))][[1]][,"Corrected_deseq.log10pvalue",drop=F])
	  }else {showlistQ <- matrix(0,0,1); colnames(showlistQ) <- "Corrected_deseq.log10pvalue"}
	  if (tmp %in% names(wpair$signifDE[[ss]])) {
		showlistW <- wpair$signifDE[[ss]][tmp][[1]][,"Corrected_wilcox.log10pval",drop=F]
		showlistW <- rbind(showlistW,wpair$signifDE[[ss]][ gsub(" ", "_",paste("BOT", ct))][[1]][,"Corrected_wilcox.log10pval",drop=F])
	  }else {showlistW <- matrix(0,0,1); colnames(showlistW) <- "Corrected_wilcox.log10pval"}
	  tmp <-  gsub(" ", "_",paste("TGO", ct))
	  if (tmp %in% names(pair$signifDE[[ss]])){
		  tmp <- pair$signifDE[[ss]][gsub(" ", "_",paste("TGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){gblock <- cbind(tmp, "+", "DESeq")
		  }else{gblock <- gblock[c(),,drop=F]}
		  colnames(gblock) <- c(colnames(gblock)[1:6], "Tail", "Test")
		  tmp <- pair$signifDE[[ss]][gsub(" ", "_",paste("BGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){
		  	tmp <- cbind(tmp, "-", "DESeq");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test")
		  	gblock <- rbind(gblock, tmp);
		  }
	  }else{
		gblock <- matrix(0,0,8); colnames(gblock) = c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
	  }


	  if ( gsub(" ", "_",paste("TGO", ct)) %in% names(wpair$signifDE[[ss]])){
	  	tmp <- wpair$signifDE[[ss]][ gsub(" ", "_",paste("TGO", ct))][[1]]$list
		 if (!is.null(dim(tmp))) { tmp <- cbind(tmp, "+", "Wilcox");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test"); gblock <- rbind(gblock, tmp);}
		  tmp <- wpair$signifDE[[ss]][ gsub(" ", "_",paste("BGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){tmp <- cbind(tmp, "-", "Wilcox");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test");  gblock <- rbind(gblock, tmp);}	 
	  }

	  print(paste("Processing", ss,ct, "has ", nrow(showlistQ) , nrow(showlistW) ))


	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  ordterm <- rep(0, length(allnames))
	  qmap <- match(allnames , rownames(pair$rawDE$deseq.log2FC))
	  wmap <- match(allnames , rownames(wpair$rawDE$wilcox.log10pval))

	  ordterm <- pair$rawDE$deseq.log2FC[qmap,dn] + wpair$rawDE$wilcox.logitAuroc[wmap,dn] * 3
	  flt <- (allnames %in% rownames(showlistQ))&(allnames %in% rownames(showlistW))
	  ordterm[flt] <- ordterm[flt] + sign(ordterm[flt]) * 10000
	  ordterm <- order(ordterm,decreasing=T)
	  allnames <- allnames[ordterm]
	  qmap <- qmap[ordterm]
	  wmap <- wmap[ordterm]



	  gmap <- match(allnames , geneannot$NAME)
	  print(allnames[is.na(gmap)])

	  if (length(allnames) > 0){
		print(paste("never happened...",length(allnames),length(genecolnames)))
	  block <- data.frame(matrix(0,length(allnames),length(genecolnames)))
	  colnames(block) <-  genecolnames
	
	  block$NAME <- allnames;
	  block$Log2FC <- pair$rawDE$deseq.log2FC[qmap,dn]
	  block$LogitAuroc <- wpair$rawDE$wilcox.logitAuroc[wmap,dn]
	  block$DEseq_Log10pval <- pair$rawDE$deseq.log10pvalue[qmap,dn]
	  block$Wilcox_Log10pval <- wpair$rawDE$wilcox.log10pval[wmap,dn]
	  block$DESeq_basemean <- pair$rawDE$deseq.basemean[qmap,dn]
	  block$TPMmean <- wpair$rawDE$meanTPM[wmap,dn]
	  block$FAD_dropout <- 1 - wpair$rawDE$dropoutPosClass[wmap,dn]
	  block$Ctrl_dropout <- 1 - wpair$rawDE$dropoutNegClass[wmap,dn]
	  block$Alias <- geneannot[gmap, "ALIAS"]
	  block$Description <- geneannot[gmap, "DESCRIPTION"]
  	  block$FullName <- geneannot[gmap, "FULLNAME"]
	  block$GOslim <- geneannot[gmap, "GOSLIM"]
	  block$GO <- geneannot[gmap, "GO"]


	  block$DEseq_adj_Log10pval <- NA;  block$Wilcox_adj_Log10pval <- NA
	  map <- match(block$NAME , rownames(showlistQ))
	  block$DEseq_adj_Log10pval[!is.na(map)] <- showlistQ[map[!is.na(map)],1]
	  map <- match(block$NAME , rownames(showlistW))
	  block$Wilcox_adj_Log10pval[!is.na(map)] <- showlistW[map[!is.na(map)],1]
 	  block$DEseq_adj_Log10pval[block$DEseq_adj_Log10pval > -1.3] <- NA
	  block$Wilcox_adj_Log10pval[block$Wilcox_adj_Log10pval > -1.3] <- NA

	  block$DE <- ""
	  fakeDE <- rep("", length(block$DE))
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q-"
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q+"
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W-",sep="")
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W+",sep="") 
	  
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q-"
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q+"
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W-",sep="")
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W+",sep="")  

	  block$Celltype <- ct
	  block$Comparison <- ss
	  block$Archtype <- archt

	  showlistQ <- showlistQ[showlistQ[,1] <= -1.3 ,,drop=F]; showlistW <- showlistW[showlistW[,1] <= -1.3 ,,drop=F]
	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  flt <- block$NAME %in% allnames
	  print(ct)
	  print(sum(flt))
	  tmp2 <- table(fakeDE)
	   if (sum(flt) > 0){
	  	  tmp <- table(block$DE[flt])
		print(tmp)
	  summary[[ss]]["Q+", ct] = length(grep("Q\\+",block$DE[flt])); summary[[ss]]["Q-", ct] = length(grep("Q\\-",block$DE[flt]))
	  summary[[ss]]["W+", ct] = length(grep("W\\+",block$DE[flt])); summary[[ss]]["W-", ct] = length(grep("W\\-",block$DE[flt]))
	  summary[[ss]]["Q+", ct] = length(grep("Q\\+",block$DE[flt])); summary[[ss]]["Q-", ct] = length(grep("Q\\-",block$DE[flt]))
	  summary[[ss]]["W+", ct] = length(grep("W\\+",block$DE[flt])); summary[[ss]]["W-", ct] = length(grep("W\\-",block$DE[flt]))
 
	  summary[[ss]]["Q+W_", ct] = tmp["Q+"]; summary[[ss]]["Q-W_", ct] = tmp["Q-"]
	  summary[[ss]]["Q_W+", ct] = tmp["W+"]; summary[[ss]]["Q_W-", ct] = tmp["W-"]
	  }
		print(tmp2)	
	  summary2[[ss]]["Q+", ct] = length(grep("Q\\+",fakeDE)); summary2[[ss]]["Q-", ct] = length(grep("Q\\-",fakeDE))
	  summary2[[ss]]["W+", ct] = length(grep("W\\+",fakeDE)); summary2[[ss]]["W-", ct] = length(grep("W\\-",fakeDE))
	  summary2[[ss]]["Q+", ct] = length(grep("Q\\+",fakeDE)); summary2[[ss]]["Q-", ct] = length(grep("Q\\-",fakeDE))
	  summary2[[ss]]["W+", ct] = length(grep("W\\+",fakeDE)); summary2[[ss]]["W-", ct] = length(grep("W\\-",fakeDE))
  
	  summary2[[ss]]["Q+W_", ct] = tmp2["Q+"]; summary2[[ss]]["Q-W_", ct] = tmp2["Q-"]
	  summary2[[ss]]["Q_W+", ct] = tmp2["W+"]; summary2[[ss]]["Q_W-", ct] = tmp2["W-"]
	  for(j in c("Q+W+","Q+W-","Q-W+","Q-W-")) {
		if (sum(flt) > 0)  summary[[ss]][j, ct] = tmp[j]
		summary2[[ss]][j, ct] = tmp2[j]
	
		}
	  genetable <- rbind(genetable, block[block$NAME %in% allnames,,drop=F])
	  genetable2 <- rbind(genetable2, block[!(block$NAME %in% allnames),,drop=F])

	  }
	  if (nrow(gblock) > 0){
	  block2 <- data.frame(matrix(0,nrow(gblock),16))
	  colnames(block2) <-  c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_dropout", "Ctrl_dropout", "Intersection")

	  block2$Domain <- gblock$Domain; block2$ID <- gblock$ID; block2$Term <- gblock$Term; block2$Intersection <- gblock$intersection
	  block2$Test <- gblock$Test; block2$Tail <- gblock$Tail; block2$Comparison <- ss; block2$Celltype <- ct; block2$pvalue <- gblock$pvalue
	  block2$Archtype <- archt;

	  for(j in 1:nrow(gblock)){
		curlist <- parselist(gblock[j,"intersection"])
	  	block2[j,c("MeanLog2FC")] <- mean(pair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"deseq.log2FC"],na.rm=T)
		block2[j,c("MeanLogitAuroc")] <- mean(wpair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"wilcox.logitAuroc"],na.rm=T)
	  	block2[j,c("FAD_dropout")] <- mean(wpair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"dropoutPosClass"],na.rm=T)
		block2[j,c("Ctrl_dropout")] <- mean(wpair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"dropoutNegClass"],na.rm=T)
	  	block2[j,c("DESeq_basemean")] <- mean(pair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"deseq.basemean"],na.rm=T)
		block2[j,c("TPMmean")] <- mean(wpair$signifDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"meanTPM"],na.rm=T)
	 }	
	  gotable<- rbind(gotable,block2) 
	  }

	}  
	}
	# do consensus now
cenetable <- data.frame(row.names=genecolnames)
cenetable2 <- data.frame(row.names=genecolnames)

cotable <- data.frame(row.names=c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_dropout", "Ctrl_dropout", "Intersection"))

	for(ss in names(pair$consensusDE)){
	print(paste("processing", ss))
	summary[[ss]] <- matrix(0,12,length(pair$list.celltypes)); summary2[[ss]] <- matrix(0,12,length(pair$list.celltypes))
	colnames(summary[[ss]]) <- pair$list.celltypes; rownames(summary[[ss]]) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")
	colnames(summary2[[ss]]) <- pair$list.celltypes; rownames(summary2[[ss]]) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")
	for(ct in pair$list.celltypes){

	  for(i in names(runargs$cel$archtype)){if (ct %in% runargs$cel$archtype[[i]]) archt <- i }

	  print(paste("Processing", ct))
	  tmp <- gsub(" ", "_",paste("TOP", ct))
	  if (tmp %in% names(pair$consensusDE[[ss]])) {
		showlistQ <- pair$consensusDE[[ss]][tmp][[1]][,"corr_log10pval",drop=F]; 
		showlistQ <- rbind(showlistQ,pair$consensusDE[[ss]][ gsub(" ", "_",paste("BOT", ct))][[1]][,"corr_log10pval",drop=F])
	  }else {showlistQ <- matrix(0,0,1); colnames(showlistQ) <- "Corrected_deseq.log10pvalue"}
	  if (tmp %in% names(wpair$consensusDE[[ss]])) {
		showlistW <- wpair$consensusDE[[ss]][tmp][[1]][,"corr_log10pval",drop=F]
		showlistW <- rbind(showlistW,wpair$consensusDE[[ss]][ gsub(" ", "_",paste("BOT", ct))][[1]][,"corr_log10pval",drop=F])
	  }else {showlistW <- matrix(0,0,1); colnames(showlistW) <- "Corrected_wilcox.log10pval"}
	  tmp <-  gsub(" ", "_",paste("TGO", ct))
	  if (tmp %in% names(pair$consensusDE[[ss]])){
		  tmp <- pair$consensusDE[[ss]][gsub(" ", "_",paste("TGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){gblock <- cbind(tmp, "+", "DESeq")
		  }else{gblock <- gblock[c(),,drop=F]}
		  colnames(gblock) <- c(colnames(gblock)[1:6], "Tail", "Test")
		  tmp <- pair$consensusDE[[ss]][gsub(" ", "_",paste("BGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){
		  	tmp <- cbind(tmp, "-", "DESeq");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test")
		  	gblock <- rbind(gblock, tmp);
	  	  }
	  }else{
		gblock <- matrix(0,0,8); colnames(gblock) = c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
	  }


	  if ( gsub(" ", "_",paste("TGO", ct)) %in% names(wpair$consensusDE[[ss]])){
	  	tmp <- wpair$consensusDE[[ss]][ gsub(" ", "_",paste("TGO", ct))][[1]]$list
		 if (!is.null(dim(tmp))) { tmp <- cbind(tmp, "+", "Wilcox");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test"); gblock <- rbind(gblock, tmp);}
		  tmp <- wpair$consensusDE[[ss]][ gsub(" ", "_",paste("BGO", ct))][[1]]$list
		  if (!is.null(dim(tmp))){tmp <- cbind(tmp, "-", "Wilcox");  colnames(tmp) <- c(colnames(gblock)[1:6], "Tail", "Test");  gblock <- rbind(gblock, tmp);}	 
	  }

	  print(paste("Processing", ss,ct, "has ", nrow(showlistQ) , nrow(showlistW) ))


	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  ordterm <- rep(0, length(allnames))

	  dn <- gsub(" ", "_",paste(ct,rargs$cel$Consensus[[ss]]))
	  qmap <- match(dn, colnames(pair$rawDE$deseq.log2FC))
	  dn <- dn[!is.na(qmap)]
	  print(dn) 
	  print("here we go")
	  if (length(dn) == 0) print("OMG")
	  qmap <- match(allnames , rownames(pair$rawDE$deseq.log2FC))
	  wmap <- match(allnames , rownames(wpair$rawDE$wilcox.log10pval))

	  ordterm <- rowSums(pair$rawDE$deseq.log2FC[qmap,dn,drop=F]) + rowSums(wpair$rawDE$wilcox.logitAuroc[wmap,dn,drop=F]) * 3
	  flt <- (allnames %in% rownames(showlistQ))&(allnames %in% rownames(showlistW))
	  ordterm[flt] <- ordterm[flt] + sign(ordterm[flt]) * 10000
	  ordterm <- order(ordterm,decreasing=T)
	  allnames <- allnames[ordterm]




	  gmap <- match(allnames , geneannot$NAME)
	  print(allnames[is.na(gmap)])

	  if (length(allnames) > 0){
		print(paste("never happened...",length(allnames),length(genecolnames)))
	  block <- data.frame(matrix(0,length(allnames),length(genecolnames)))
	  colnames(block) <-  genecolnames
	  block$NAME <- allnames;
	

	  qmap <- qmap[ordterm]
	  wmap <- wmap[ordterm]


	  block$Log2FC <- rowSums(pair$rawDE$deseq.log2FC[qmap,dn,drop=F]) / length(dn)
	  block$LogitAuroc <- rowSums(wpair$rawDE$wilcox.logitAuroc[wmap,dn,drop=F])/ length(dn)

	  block$DEseq_Log10pval <- sapply(rowSums(pair$rawDE$deseq.log10pvalue[qmap,dn,drop=F]),function(x){return(pgamma(-x,length(dn),scale=log10(exp(1)),log.p=T,lower.tail=F) / log(10))})
	  block$Wilcox_Log10pval <- sapply(rowSums(wpair$rawDE$wilcox.log10pval[qmap,dn,drop=F]),function(x){return(pgamma(-x,length(dn),scale=log10(exp(1)),log.p=T,lower.tail=F) / log(10))})

#	  block$DEseq_Log10pval <- rowSums(pair$rawDE$deseq.log10pvalue[qmap,dn,drop=F])/ length(dn)
#	  block$Wilcox_Log10pval <- rowSums(wpair$rawDE$wilcox.log10pval[wmap,dn,drop=F])/ length(dn)


	  block$DESeq_basemean <- rowSums(pair$rawDE$deseq.basemean[qmap,dn,drop=F])/ length(dn)
	  block$TPMmean <- rowSums(wpair$rawDE$meanTPM[wmap,dn,drop=F])/ length(dn)
	  block$FAD_dropout <- 1 - ( rowSums(wpair$rawDE$dropoutPosClass[wmap,dn,drop=F])/ length(dn))
	  block$Ctrl_dropout <- 1 - ( rowSums(wpair$rawDE$dropoutNegClass[wmap,dn,drop=F])/ length(dn)) 


	  block$Alias <- geneannot[gmap, "ALIAS"]
	  block$Description <- geneannot[gmap, "DESCRIPTION"]
  	  block$FullName <- geneannot[gmap, "FULLNAME"]
	  block$GOslim <- geneannot[gmap, "GOSLIM"]
	  block$GO <- geneannot[gmap, "GO"]


	  block$DEseq_adj_Log10pval <- NA;  block$Wilcox_adj_Log10pval <- NA
	  map <- match(block$NAME , rownames(showlistQ))
	  block$DEseq_adj_Log10pval[!is.na(map)] <- showlistQ[map[!is.na(map)],1]
	  map <- match(block$NAME , rownames(showlistW))
	  block$Wilcox_adj_Log10pval[!is.na(map)] <- showlistW[map[!is.na(map)],1]
 	  block$DEseq_adj_Log10pval[block$DEseq_adj_Log10pval > -1.3] <- NA
	  block$Wilcox_adj_Log10pval[block$Wilcox_adj_Log10pval > -1.3] <- NA

	  block$DE <- ""
	  fakeDE <- rep("", length(block$DE))
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q-"
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q+"
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W-",sep="")
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W+",sep="") 
	  
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q-"
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q+"
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W-",sep="")
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W+",sep="") 

	  block$Celltype <- ct
	  block$Comparison <- ss
	  block$Archtype <- archt

	  showlistQ <- showlistQ[showlistQ[,1] <= -1.3 ,,drop=F]; showlistW <- showlistW[showlistW[,1] <= -1.3 ,,drop=F]
	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  flt <- block$NAME %in% allnames
	  print(ct)
	  print(sum(flt))
	  tmp2 <- table(fakeDE)
	   if (sum(flt) > 0){
	  	  tmp <- table(block$DE[flt])
	  summary[[ss]]["Q+", ct] = length(grep("Q\\+",block$DE[flt])); summary[[ss]]["Q-", ct] = length(grep("Q\\-",block$DE[flt]))
	  summary[[ss]]["W+", ct] = length(grep("W\\+",block$DE[flt])); summary[[ss]]["W-", ct] = length(grep("W\\-",block$DE[flt]))
	  summary[[ss]]["Q+", ct] = length(grep("Q\\+",block$DE[flt])); summary[[ss]]["Q-", ct] = length(grep("Q\\-",block$DE[flt]))
	  summary[[ss]]["W+", ct] = length(grep("W\\+",block$DE[flt])); summary[[ss]]["W-", ct] = length(grep("W\\-",block$DE[flt]))
 
	  summary[[ss]]["Q+W_", ct] = tmp["Q+"]; summary[[ss]]["Q-W_", ct] = tmp["Q-"]
	  summary[[ss]]["Q_W+", ct] = tmp["W+"]; summary[[ss]]["Q_W-", ct] = tmp["W-"]
	  }
		print(tmp2)	
	  summary2[[ss]]["Q+", ct] = length(grep("Q\\+",fakeDE)); summary2[[ss]]["Q-", ct] = length(grep("Q\\-",fakeDE))
	  summary2[[ss]]["W+", ct] = length(grep("W\\+",fakeDE)); summary2[[ss]]["W-", ct] = length(grep("W\\-",fakeDE))
	  summary2[[ss]]["Q+", ct] = length(grep("Q\\+",fakeDE)); summary2[[ss]]["Q-", ct] = length(grep("Q\\-",fakeDE))
	  summary2[[ss]]["W+", ct] = length(grep("W\\+",fakeDE)); summary2[[ss]]["W-", ct] = length(grep("W\\-",fakeDE))
  
	  summary2[[ss]]["Q+W_", ct] = tmp2["Q+"]; summary2[[ss]]["Q-W_", ct] = tmp2["Q-"]
	  summary2[[ss]]["Q_W+", ct] = tmp2["W+"]; summary2[[ss]]["Q_W-", ct] = tmp2["W-"]
	  for(j in c("Q+W+","Q+W-","Q-W+","Q-W-")) {
		if (sum(flt) > 0)  summary[[ss]][j, ct] = tmp[j]
		summary2[[ss]][j, ct] = tmp2[j]
		}
	  cenetable <- rbind(cenetable, block[block$NAME %in% allnames,,drop=F])
	  cenetable2 <- rbind(cenetable2, block[!(block$NAME %in% allnames),,drop=F])

	  }
	  if (nrow(gblock) > 0){
	  block2 <- data.frame(matrix(0,nrow(gblock),16))
	  colnames(block2) <-  c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_dropout", "Ctrl_dropout", "Intersection")

	  block2$Domain <- gblock$Domain; block2$ID <- gblock$ID; block2$Term <- gblock$Term; block2$Intersection <- gblock$intersection
	  block2$Test <- gblock$Test; block2$Tail <- gblock$Tail; block2$Comparison <- ss; block2$Celltype <- ct; block2$pvalue <- gblock$pvalue
	  block2$Archtype <- archt;

	  for(j in 1:nrow(gblock)){
		curlist <- parselist(gblock[j,"intersection"])
	  	block2[j,c("MeanLog2FC")] <- mean(pair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"deseq.log2FC"],na.rm=T)
		block2[j,c("MeanLogitAuroc")] <- mean(wpair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"wilcox.logitAuroc"],na.rm=T)
	  	block2[j,c("FAD_dropout")] <- mean(wpair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"dropoutPosClass"],na.rm=T)
		block2[j,c("Ctrl_dropout")] <- mean(wpair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"dropoutNegClass"],na.rm=T)
	  	block2[j,c("DESeq_basemean")] <- mean(pair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"deseq.basemean"],na.rm=T)
		block2[j,c("TPMmean")] <- mean(wpair$consensusDE[[ss]][[gsub(" ", "_",paste("ORD", ct))]][curlist,"meanTPM"],na.rm=T)
	 }	
	  cotable<- rbind(cotable,block2) 
	  }
	}

	summary[[ss]][is.na(summary[[ss]])] <- 0
	summary2[[ss]][is.na(summary2[[ss]])] <- 0
}

print("Saving output!")
return(list(cons_gene=cenetable,cons_gene_preFDR=cenetable2,cons_go=cotable, gene=genetable,gene_preFDR=genetable2,go=gotable,summary= summary,summary_preFDR= summary2))}


# reads and compute soup, overwrites raw counts for listed cells
# prefix = c(paste("Ag18",c("WtWt_P","HtWt_N","HtWt_P","WtHt_P","HtHt_N","HtHt_P"),sep=""),paste("Ja19", c("WtWt_K2","WtWt_H9","HtWt_K2","HtWt_H9"),sep=""), paste("Mr19",c("WtWt","HtWt_PSEN", "HmWt_PSEN","HtWt_APP","HmWt_APP"),sep=""))
calcSoupX <- function(sro,bggenes, bggenes2, bggenes3, paths = paste("/warehouse/team218_wh01/lh20/mh/785015", 0:4, "/human93_extended",sep=""), alias = paste("Mr19",c("WtWt","HtWt_PSEN", "HmWt_PSEN","HtWt_APP","HmWt_APP"),sep=""), rho.meta= "SoupX_Rho",flt.cells=c(), meta.cells=c(), meta.cells2=c(),meta.cells3){
	library(SoupX);
	flt.cells2 <- SeuratCellListQuery(sro, meta.cells2, flt.cells)
	flt.cells3 <- SeuratCellListQuery(sro, meta.cells3, flt.cells)
	flt.cells <- SeuratCellListQuery(sro, meta.cells, flt.cells)

	print(sum(flt.cells))
	print(sum(flt.cells2))
	print(sum(flt.cells3))

	if (!rho.meta %in% names(sro@meta.data)) sro@meta.data[[rho.meta]] <- rep(0, nrow(sro@meta.data))
	nrho <- sro@meta.data[[rho.meta]][flt.cells]
	nrho2 <-sro@meta.data[[rho.meta]][flt.cells2]
	nrho3 <- sro@meta.data[[rho.meta]][flt.cells3]
	if (length(paths) != length(alias)) stop("number of path and alias mismatch!")
	for(i in 1:length(paths)){
		orig <- SoupX::load10X(paths[i])
		print("Microglia time?")

		tmp <- adjustCounts(interpolateCellContamination(calculateContaminationFraction(orig, "Channel1", list(IG= bggenes)), "Channel1", interpolationMethod="lowess"))
		colnames(tmp$atoc) <- sub("Channel1__", alias[i], colnames(tmp$atoc))
		map <- match(rownames(sro@meta.data)[flt.cells], colnames(tmp$atoc))
		print(paste("found", sum(!is.na(map)), "cell to overwrite in", alias[i]))
		
		nrho[!is.na(map)] <- tmp$channels$Channel1$rhos[map[!is.na(map)]]
		print(sort(which(flt.cells)[!is.na(map)]))
		sro@assays$RNA@counts[, which(flt.cells)[!is.na(map)]]  <- tmp$atoc[,map[!is.na(map)]]
		
		print("Neuron time?")
		tmp <- adjustCounts(interpolateCellContamination(calculateContaminationFraction(orig, "Channel1", list(IG= bggenes2)), "Channel1", interpolationMethod="lowess"))
		colnames(tmp$atoc) <- sub("Channel1__", alias[i], colnames(tmp$atoc))
		map <- match(rownames(sro@meta.data)[flt.cells2], colnames(tmp$atoc))
		print(paste("found", sum(!is.na(map)), "cell to overwrite in", alias[i]))
		
		nrho2[!is.na(map)] <- tmp$channels$Channel1$rhos[map[!is.na(map)]]
		print(sort(which(flt.cells2)[!is.na(map)]))
		sro@assays$RNA@counts[, which(flt.cells2)[!is.na(map)]]  <- tmp$atoc[,map[!is.na(map)]]
		
		print("Rest time?")

		tmp <- adjustCounts(interpolateCellContamination(calculateContaminationFraction(orig, "Channel1", list(IG= bggenes3)), "Channel1", interpolationMethod="lowess"))
		colnames(tmp$atoc) <- sub("Channel1__", alias[i], colnames(tmp$atoc))
		map <- match(rownames(sro@meta.data)[flt.cells3], colnames(tmp$atoc))
		print(paste("found", sum(!is.na(map)), "cell to overwrite in", alias[i]))
		
		nrho3[!is.na(map)] <- tmp$channels$Channel1$rhos[map[!is.na(map)]]
		print(sort(which(flt.cells3)[!is.na(map)]))
		sro@assays$RNA@counts[, which(flt.cells3)[!is.na(map)]]  <- tmp$atoc[,map[!is.na(map)]]
	}
	sro@meta.data[[rho.meta]][flt.cells] <- nrho 
	sro@meta.data[[rho.meta]][flt.cells2] <- nrho2
	sro@meta.data[[rho.meta]][flt.cells3] <- nrho3

return(sro)}

mergeIntoSparse <- function(spamatrix, datafr, colname, newcolname){
	map <- match(rownames(datafr), rownames(spamatrix))
	print(map)
	ord <- order(map)
	print(map[ord])
	print(datafr[[colnames]][ord])
	tmpsparse <- new("dgCMatrix"
    , i = map[ord]
    , p = c(0L, nrow(datafr))
    , Dim = c(nrow(spamatrix), 1L)
    , Dimnames = list(rownames(spamatrix), NULL )
    , x = datafr[[colnames]][ord]
    , factors = list())
	print("success")
return(tmpsparse)}
argsenv <- function(..., parent=parent.frame()) {
  cl <- match.call(expand.dots=TRUE)

  e <- new.env(parent=parent)
  pf <- parent.frame()
  JJ <- seq_len(length(cl) - 1)
  tagnames <- sprintf(".v%d", JJ)
  for (i in JJ) e[[tagnames[i]]] <- eval(cl[[i+1]], envir=pf)

  attr(e, "tagnames") <- tagnames
  attr(e, "formalnames") <- names(cl)[-1]
  class(e) <- c("environment", "argsenv")
return(e)}

mergeSparseMatrices <- function(...){
	inlist <- argsenv(...)	
	rown <- c()
	coln <- c()
	nbentry <- 0
	for(i in attr(inlist, "tagnames")) {
		rown <- c(rown, rownames(inlist[[i]]))
		coln <- c(coln, colnames(inlist[[i]]))
		nbentry <- nbentry + inlist[[i]]@p[inlist[[i]]@Dim[2]+1]
	}
	rown <- unique(rown)
	out <- new("dgCMatrix", i = rep(0,nbentry), p = rep(0, length(coln)+1)
    , Dim = c(length(rown), length(coln))
    , Dimnames = list(rown, coln )
    , x = rep(0,nbentry))

	for(i in attr(inlist, "tagnames")) {
		for(c in 1:inlist[[i]]@p[inlist[[i]]@Dim[2]]){
			

		} 
	}


return(out)
}

showListsIntersectionSizes <- function(listA, listB, listC=c(), listD=c(), extract=c()){
	if (!is.null(extract)){
		if (length(grep("[Rr]",extract)) != 0){
			listA <- rownames(listA)
			listB <- rownames(listB)
			if (!is.null(listC)) listC <- rownames(listC)
		}else if (length(grep("[Cc]",extract)) != 0){
			listA <- colnames(listA)
			listB <- colnames(listB)
			if (!is.null(listC)) listC <- colnames(listC)

		}
	}

	if (is.null(listC)){
		out <-matrix(0,3,3)
		colnames(out) <- c("!A", "A", "Total")
		rownames(out) <- c("!B", "B", "Total")
		out[1,1] = 0
		out[1,2] = length(setdiff(listA,listB))
		out[1,3] = 0
		out[2,1] = length(setdiff(listB,listA))
		out[2,2] = length(intersect(listA,listB))
		out[2,3] = length(listB)
		out[3,1] = 0
		out[3,2] = length(listA)
		out[3,3] = length(union(listA, listB))
	}else if (is.null(listD)){
		out <-matrix(0,6,3)
		colnames(out) <- c("!A", "A", "Total")
		rownames(out) <- c("!B","B","!C", "C","!B_!C", "B_!C","!B_C", "B_C")
		out[1,1] = 0
		out[1,2] = length(setdiff(listA,listB))
		out[1,3] = 0
		out[2,1] = length(setdiff(listB,listA))
		out[2,2] = length(intersect(listA,listB))
		out[2,3] = length(listB)
		out[3,1] = 0
		out[3,2] = length(listA)
		out[3,3] = length(union(listA, listB))
	}

return(out)}

readDataset <- function(which){
	if (which == "Heart10") return(readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/Heart_10X.rds"));
	if (which == "HeartF") return(readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/Heart_FACS.rds"));
	if (which == "FatF") return(readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/Fat_FACS.rds"));
	stop(paste(which, "is not a known dataset"))
}
tableparse <-function(table, filterpairs=c(), outrownames = c(), expandcol=c(),subsetcol=c(), countuniq =c(), sortcriteria =c(), countuniq.do.filter.zero=T){
	if (is.null(filterpairs)) filter = rep(T, nrow(table))
	else{
		for(i in 1:length(filterpairs)){
			print(length(filterpairs[[i]]))
			if (length(filterpairs[[i]]) == 1){
				newfilter <- table[[names(filterpairs)[i]]] == filterpairs[i]
			}else{
				if (filterpairs[[i]][1] %in% c(">", "<")){
					if (filterpairs[[i]][1] == ">") newfilter <- table[[names(filterpairs)[i]]] > as.numeric(filterpairs[[i]][2]) 
					else newfilter <- table[[names(filterpairs)[i]]] < as.numeric(filterpairs[[i]][2])
				}else{
					tmp <- match(table[[names(filterpairs)[i]]],filterpairs[[i]])
					newfilter <- !is.na(tmp)
				}
			}
			if (i == 1) filter <- newfilter
			else filter <- (filter) & newfilter
		}
	}
	table <-  table[filter,,drop=F]
	if (!is.null(expandcol)){
		mmm <- unique(table[[expandcol[1]]])
		mmm2 <- unique(table[[expandcol[2]]])
		output <- matrix(0,length(mmm), length(mmm2))
		print(dim(output))
		rownames(output) <- mmm
		colnames(output) <- mmm2
		for(i in 1:nrow(table)){
			output[match(table[i,expandcol[1]], mmm),match(table[i,expandcol[2]], mmm2)] <-  table[i,expandcol[3]]
		}
		return(output)
	}
       	if (!is.null(sortcriteria)) {ord <-order(table[[sortcriteria]]);table <- table[ord,]}
	if (!is.null(countuniq)){
		output <- sort(table(table[[countuniq]]))
		if (countuniq.do.filter.zero) output <- output[output!=0] 
		return(output)
	}else{
	   if (!is.null(subsetcol)){
		output <- table[, subsetcol, drop=F]
		}else  output <- table
	}
	if (!is.null(outrownames)) rownames(output) <- table[,outrownames]
	return(output)
}

findRejection <- function(sro, data=c(), truncate.percentile = 0.95, bias= 0 ,meta.partition, pos.cell=c(), pos.meta=c(), use.meta=c(), use.cell=c(), correlation.factor=1.0,bandwidth= 0.25, do.use.adaptive.bias=F, mincell= 100){

	if (length(meta.partition) > 1) {
		mlvl <- meta.partition[2:length(meta.partition)]
		meta.partition <- meta.partition[1]
	}else mlvl <- getSafeLevels(sro@meta.data[[meta.partition]])

	if (is.null(data)){
		data <- cbind(log2(sro@meta.data$nUMI + 1), log2(rowSums(sro@meta.data[,grep("Filtered", colnames(sro@meta.data))])))
	}

	pos.cell <- SeuratCellListQuery(sro, pos.meta,pos.cell)
 	use.cell <- SeuratCellListQuery(sro, use.meta,use.cell)

	if (sum(pos.cell) == nrow(sro@meta.data)) stop("Forgot to define the pos/neg classes")

	weight <- rep(0,nrow(sro@meta.data))
	weight[pos.cell] <- 1

	sampleout <- rep(T, nrow(sro@meta.data))

	densities <- list()
	library(pheatmap)
	library(gridExtra)
	plotlist <- list(A=c(),B=c(),C=c(),D=c())
	imagesize = 512;

	for(i in 1:length(mlvl)){
		curcell <- use.cell & (sro@meta.data[[meta.partition]] == mlvl[i]);
		curdata <- data[curcell,,drop=F]
		print(dim(curdata))
		dacov <- cov(curdata)
		dacov[2,1] <- dacov[2,1] * correlation.factor * correlation.factor;
		dacov[1,2] <- dacov[1,2] * correlation.factor * correlation.factor;


		posout <- InferN0GetDensity(curdata,weight=weight[curcell],cov=dacov, mapsize=512,bandwidth=bandwidth)
		negout <- InferN0GetDensity(curdata,weight=1-weight[curcell],cov=dacov, mapsize=512,bandwidth=bandwidth)

		sselect <- weight[curcell] == 1;
		if ((sum(sselect) < sum(weight[curcell]))||(sum(sselect) > 0)) {

		posout$Query[posout$Query < 0] <- 0
		negout$Query[negout$Query < 0] <- 0
		rejpos <- 1.0 / posout$Query[sselect]
		rejneg <- 1.0 / negout$Query[!sselect]

		if (do.use.adaptive.bias) bias <- 2*(0.5-sum(weight[curcell]) / sum(curcell))
		print(paste("Bias is ", bias))

		if (abs(bias) < 1){
			rejpos <- posout$Query[ sselect] ^ (( bias - 1) / 2) * negout$Query[ sselect] ^ ((1 - bias) / 2) 
			rejneg <- negout$Query[!sselect] ^ ((-bias - 1) / 2) * posout$Query[!sselect] ^ ((1 + bias) / 2) 
		}else if (bias < 0){
			rejpos <- negout$Query[ sselect] / posout$Query[ sselect]
			rejneg <- rep(1, sum(!sselect))
		}else{
			rejpos <- rep(1, sum(sselect))
			rejneg <- posout$Query[!sselect] / negout$Query[!sselect]
		}
		rejpos[is.na(rejpos)|(rejpos<0)] <- 0	
		rejneg[is.na(rejneg)|(rejneg<0)] <- 0
		#print(rejpos[order(rejpos)[length(rejpos)* 0.8]])
		truncp <- rejpos[order(rejpos)[length(rejpos)* truncate.percentile]]
		truncn <- rejneg[order(rejneg)[length(rejneg)* truncate.percentile]]
		
		tmptmp <- rejpos / truncp
		tmptmp[tmptmp > 1] <- 1
		if (sum(tmptmp) < mincell){
			if (length(tmptmp) <= mincell) {tmptmp <- rep(1, length(tmptmp)) ; truncp <- 0.00000001
			}else{
				dif = (mincell - sum(tmptmp)) / (length(tmptmp) - sum(tmptmp))
				tmptmp <- tmptmp + dif * (1 - tmptmp)
			}
		}
		rejpos <- tmptmp

		tmptmp <- rejneg / truncn
		tmptmp[tmptmp > 1] <- 1
		if (sum(tmptmp) < mincell){
			if (length(tmptmp) <= mincell) {
				tmptmp <- rep(1, length(tmptmp))
				truncn <- 0.00000001
			}else{
				dif = (mincell - sum(tmptmp)) / (length(tmptmp) - sum(tmptmp))
				tmptmp <- tmptmp + dif * (1 - tmptmp)
			}
		}
		rejneg <- tmptmp

		



		#rejpos <- 1 / (1 + sqrt(posout$Query/negout$Query))
		#rejneg <- 1 / (1 + sqrt(negout$Query/posout$Query))

		nbpos = sum(weight[curcell])
		nbneg = sum(1-weight[curcell])
		subsample <- rep(0,sum(curcell))
		subsample[sselect] <- rejpos
		subsample[!sselect] <- rejneg
		sampleout[curcell] <- subsample
	

		print(paste("In", mlvl[i],",Would accept", sum(rejpos) / nbpos ,"of pos and",sum(rejneg)  /nbneg ,"neg cells (expects",sum(rejpos), "and",sum(rejneg),")" ))
		densities[[paste(mlvl[i], "Pos_density",sep="_")]] <- posout$Image
		densities[[paste(mlvl[i], "Neg_density",sep="_")]] <- negout$Image

		if (abs(bias) < 1){
			tmp <-  (posout$Image ^ (( bias - 1) / 2) *  negout$Image ^ ((1 - bias) / 2)) / truncp
			tmp2 <-  (negout$Image ^ ((-bias - 1) / 2) *  posout$Image ^ ((1 + bias) / 2)) / truncn
			tmp[tmp > 1] <- 1; tmp <- tmp * posout$Image 
			tmp2[tmp2 > 1] <- 1; tmp2 <- tmp2 * negout$Image
		}else if (bias < 0){
			tmp <- 1.0 / (truncp * negout$Image)
			tmp[tmp > 1] <- 1; tmp <- tmp * posout$Image
			tmp2 <- negout$Image
		}else{
			tmp2 <- 1.0 / (truncn * posout$Image)
			tmp2[tmp2 > 1] <- 1; tmp2 <- tmp2 * negout$Image
			tmp <- posout$Image
		}
		tmp[tmp > 1] <- 1; tmp[is.na(tmp)] <- 0; densities[[paste(mlvl[i], "Pos_sample_density",sep="_")]] <- tmp
		tmp2[tmp2>1] <- 1; tmp2[is.na(tmp2)] <- 0; densities[[paste(mlvl[i], "Neg_sample_density",sep="_")]]  <- tmp2
		
		plotlist[["A"]] <- pheatmap(posout$Image, cluster_rows=F, cluster_cols=F,silent=T)[[4]]
		plotlist[["B"]] <- pheatmap(negout$Image, cluster_rows=F, cluster_cols=F,silent=T)[[4]]
		plotlist[["C"]] <- pheatmap(tmp, cluster_rows=F, cluster_cols=F,silent=T)[[4]]
		plotlist[["D"]] <- pheatmap(tmp2, cluster_rows=F, cluster_cols=F,silent=T)[[4]]
		grid.arrange(arrangeGrob(grobs=plotlist, ncol=2))
	}
	}
return(list(sampleout=sampleout, densities=densities))}



makeFancyFilterDEseqMH <- function(base, basenames = c("PSENNeuNP", "CtrlPSENNeuNP")){
	filter= (base[[basenames[1]]]$deseq.log10pvalue < -1.3)
	for(i in 2: length(basenames)){
		filter = filter | (base[[basenames[2]]]$deseq.log10pvalue < -1.3)
	}
return(filter)}

makeFancyFilterDEseq <- function(base, curlist =c(), basenames = c("PSENNeuNP", "CtrlPSENNeuNP"), curname="", nbtop=5000,selected.colnames=c(),gene.blacklist=c()){
	if (!is.null(selected.colnames)){
		map <- match(selected.colnames, colnames(base[[basenames[1]]]$deseq.log10pvalue))
		map[is.na(map)] <- ncol(base[[basenames[1]]]$deseq.log10pvalue) + 1
		map2 <- match(selected.colnames, colnames(base[[basenames[2]]]$deseq.log10pvalue))
		map2[is.na(map2)] <- ncol(base[[basenames[2]]]$deseq.log10pvalue) + 1
		input <- cbind(base[[basenames[1]]]$deseq.log10pvalue, rep(0, nrow(base[[basenames[1]]]$deseq.log10pvalue)))
		input2 <- cbind(base[[basenames[2]]]$deseq.log10pvalue , rep(0, nrow(base[[basenames[2]]]$deseq.log10pvalue)))
		filter = (input[,map,drop=F] < -1.3)|(input2[,map2,drop=F] < -1.3)
		colnames(filter) <- selected.colnames
	}else filter= ((base[[basenames[1]]]$deseq.log10pvalue < -1.3)|(base[[basenames[2]]]$deseq.log10pvalue < -1.3));
	nbrow <- nrow(base[[basenames[1]]]$deseq.basemean)
	curcolnames <- colnames(filter)
	for(i in 1:length(curcolnames)) {
		tmp <- order(base[[basenames[1]]]$deseq.basemean[,i],decreasing=T);
		tmp2 <- order(base[[basenames[2]]]$deseq.basemean[,i],decreasing=T);
		filter[,i] <- filter[,i] & (!(is.na(match(1:nbrow, tmp[1:nbtop]))&(is.na(match(1:nbrow, tmp2[1:nbtop]))) ) )
	}
	filter[rownames(filter) %in% gene.blacklist,] <- F
return(filter)}





makeFancyFilterWilcox <- function(base, curlist, basenames = c("PSENNeuNP", "CtrlPSENNeuNP"), curname="", dropoutthreshold=0.1,selected.colnames=c(),gene.blacklist=c()){
	if (!is.null(selected.colnames)){
		map <- match(selected.colnames, colnames(base[[basenames[1]]]$deseq.log10pvalue))
		map[is.na(map)] <- ncol(base[[basenames[1]]]$deseq.log10pvalue) + 1
		map2 <- match(selected.colnames, colnames(base[[basenames[2]]]$deseq.log10pvalue))
		map2[is.na(map2)] <- ncol(base[[basenames[2]]]$deseq.log10pvalue) + 1
		input <- cbind(base[[basenames[1]]]$deseq.log10pvalue, rep(0, nrow(base[[basenames[1]]]$deseq.log10pvalue)))
		input2 <- cbind(base[[basenames[2]]]$deseq.log10pvalue , rep(0, nrow(base[[basenames[2]]]$deseq.log10pvalue)))
		filter = (input[,map,drop=F] < -1.3)|(input2[,map2,drop=F] < -1.3)
		colnames(filter) <- selected.colnames
	}else filter= ((base[[basenames[1]]]$deseq.log10pvalue < -1.3)|(base[[basenames[2]]]$deseq.log10pvalue < -1.3));
	curcolnames <- colnames(filter)
	if (curname != "") curlist <- curlist[[curname]]
	filter <- filter[match(rownames(curlist$dropoutPosClass), rownames(filter)),,drop=F]
	for(i in 1:length(curcolnames)) {
		filter[,i] <- filter[,i] & ((curlist$dropoutPosClass[,curcolnames[i]] > dropoutthreshold)|(curlist$dropoutNegClass[,curcolnames[i]] > dropoutthreshold))
	}
	filter[rownames(filter) %in% gene.blacklist,] <- F
return(filter)}

getIntersection <- function(objA, objB){
	if (is.null(dim(objA))){
		map <- match(names(objA), names(objB))
		return(list(first= objA[!is.na(map)], second = objB[map[!is.na(map)]]))
	}else{
		map <- match(rownames(objA), rownames(objB))
		return(list(first= objA[!is.na(map),,drop=F], second = objB[map[!is.na(map)],,drop=F]))
	}
}

multiCbind <- function(listA, listB, prefixA="", prefixB){
	listO <- list()
	for(i in names(listA)){
		colnames(listA[[i]]) <- paste(prefixA, colnames(listA[[i]]),sep="")
		colnames(listB[[i]]) <- paste(prefixB, colnames(listB[[i]]),sep="")
		listO[[i]] <- cbind(listA[[i]], listB[[i]])
	}
return(listO)}


# return a list of n colors, or a named list of n colors if a list of name is provided
myrainbow <- function(n,nbring=c(), huerange = c(0,360) ){
	if (length(nbring) == 0) nbring =0
	if (length(n) != 1){
		nams <- n
		n <- length(n)
	}else{
		nams <- c()
	}
	if (nbring == 0) nbring <- 1 + floor(n / 12);
	fout = c();	
	phase = pi / (nbring);
	term = pi / (nbring);
	for(i in 0:(n-1)){
		fout <- c(fout, hcl((floor(i / nbring) * (huerange[2]- huerange[1]) / floor(n / nbring)) + huerange[1], 100 + 50 *sin(phase + term* (i%%nbring)) , 50 - 25 *cos(phase + term* (i%%nbring))  ))
	}
	if (is.null(nams)) return(fout)
	fout2 = c()
	for(i in 1:n) fout2[[nams[i]]] <- fout[i]
return(fout2)}

myFireAndIce <- function(val){return(colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[1 + (val*40)])}
# use identical hue for subgroups, c(4,5,6) means 15 colors in 3 subgroups with respective sizes
mydoublerainbow <- function(subgroup, huerange = c(0,360), selectedHues=c(), selectedSats=c(), names= c(), permutation=c(), alpha = "",do.optimize.hues=T,useHCL=F){
	nbring = length(subgroup)
	fout = c();
	if (length(permutation) != length(subgroup)) permutation = 1:length(subgroup)
	if (is.null(selectedHues)){
		n = 0
		for(j in 1:nbring) n <- n + 1 + (subgroup[j]/3)

		selectedHues = rep(0,nbring)
		nd <- 0
		for(j in 1:nbring){
			selectedHues[permutation[j]] = (nd  * (huerange[2]- huerange[1]) / n) + huerange[1]
			nd <- nd + 1 + (subgroup[permutation[j]]/3)
		}
	}
		midL <- 20 * (0.25*sin(selectedHues * 0.01745329251994329576923690768489) - cos(selectedHues * 0.05235987755982988730771072305466))
	if (do.optimize.hues){
		selectedHues = selectedHues - 15 * sin(selectedHues *0.05235987755982988730771072305466)
	}

	if (useHCL){
		for(j in 1:nbring){
		phase = pi / subgroup[j];
		term = pi / subgroup[j];
		for(i in 0:(subgroup[j]-1)){
			print(c(selectedHues[j], 250 -3* midL[j] + 50 *sin(phase + term* (i%%subgroup[j])) , midL[j] - 25 * cos(phase + term* (i%%subgroup[j]))))
			fout <- c(fout, hcl(selectedHues[j], 150 - 3 * midL[j] , 50 + midL[j] - 25 * cos(phase + term* (i%%subgroup[j])) ))
		}
		}
	}else{
		sat =1
		for(j in 1:nbring){
		phase = pi / subgroup[j];
		term = pi / subgroup[j];
		if (!is.null(selectedSats)) sat <- selectedSats[j]
		for(i in 0:(subgroup[j]-1)){
			fout <- c(fout, getHCYcolor(selectedHues[j]/360, sat , (0.525 - (0.45 - 1 / (1 + subgroup[j])) * cos( phase + term* (i%%subgroup[j] )))^ 1.5,opt.hue=F ))
		}
		}
	}

	if (alpha != "") fout <- paste(fout, alpha, sep="")

	if (length(names) == 0) return(fout)
	fout2 = c()
	for(i in 1:length(names)) fout2[[names[i]]] <- fout[i]
	return(fout2)
}


mytripplerainbow <- function(subgroup, huerange = c(0,360), selectedHues=c(), names= c(), permutation=c(), alpha = "",do.optimize.hues=T){
	nbring = length(subgroup)
	fout = c();
	if (length(permutation) != length(subgroup)) permutation = 1:length(subgroup)
	if (is.null(selectedHues)){
		n = 0
		for(j in 1:nbring) n <- n + 1 + (subgroup[j]/3)

		selectedHues = rep(0,nbring)
		nd <- 0
		for(j in 1:nbring){
			selectedHues[permutation[j]] = (nd  * (huerange[2]- huerange[1]) / n) + huerange[1]
			nd <- nd + 1 + (subgroup[permutation[j]]/3)
		}
	}
		midL <- 20 * (0.25*sin(selectedHues * 0.01745329251994329576923690768489) - cos(selectedHues * 0.05235987755982988730771072305466))
	if (do.optimize.hues){
		selectedHues = selectedHues + 15 * sin(selectedHues *0.05235987755982988730771072305466)
	}

	for(j in 1:nbring){
		phase = pi / subgroup[j];
		term = pi / subgroup[j];
		for(i in 0:(subgroup[j]-1)){
			print(c(selectedHues[j], 250 -3* midL[j] + 50 *sin(phase + term* (i%%subgroup[j])) , midL[j] - 25 * cos(phase + term* (i%%subgroup[j]))))
			fout <- c(fout, hcl(selectedHues[j], 150 - 3 * midL[j] , 50 + midL[j] - 25 * cos(phase + term* (i%%subgroup[j])) ))
		}
	}
	if (alpha != "") fout <- paste(fout, alpha, sep="")

	if (length(names) == 0) return(fout)
	fout2 = c()
	for(i in 1:length(names)) fout2[[names[i]]] <- fout[i]
	return(fout2)
}



#myPngPrint <- function(name, width=480, height=360){
#    dev.copy(png, paste("/media/sf_SharedFolder/", name, ".png", sep=""), width = width, height = height)
#    dev.off();
#}

equalCheckPalette <- function(sorteddata, cropcol=0){
    colors <- colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(65)

    nbunique = 0
    for(i in 2:length(sorteddata)){
        if (sorteddata[i] != sorteddata[i-1]) nbunique <- nbunique +1
    }


    matout <- 1:length(sorteddata)
    if (nbunique == 0) {
        matout <- matout[33]
        return(matout)
    }
    matout[1] <- colors[1]
    ite=0;
    for(i in 2:length(sorteddata)){
        if (sorteddata[i] != sorteddata[i-1]) ite <- ite +1
        matout[i] <- colors[1 +32 * cropcol+ (64 * (1-cropcol) * ite) / nbunique]
    }
    return(matout)
}


#' Gets any data, spit colors
#'
#' @param list of numbers
#'
#' @export
InferN0Colors <- function(data, filter.na=F, threshold.value=NA, nb.valrange=10, doorder=T, cropcol=0){
	if (class(data) == "factor") {
		for(i in 1:length(levels(data))) cat(paste(i, ":", levels(data)[i],"\n"))
		data <- data@.Data
	}
	valrange <- 1:nb.valrange;
	valcolrs <- c()
	if (!is.na(threshold.value)){
		mat <- 1:length(data);
		mat[data > threshold.value] = "#0088FF";
		mat[data < threshold.value] = "#FF0000";
		mat[is.na(data)] = "#AAAAAA";
		nbnonzero <- sum(!is.na(data))
	}else{
		if (filter.na){

			ndata = data[!is.na(data)]
			nbnonzero = length(ndata)
			mat <- 1:nbnonzero

            mat[order(ndata)] <- equalCheckPalette(sort(ndata), cropcol=cropcol)

			valrange <- 1:nbnonzero;
			tmp <- order(ndata)
			for(i in 1:nbnonzero) valrange[i] <- ndata[tmp[i]]
			valcolrs <- equalCheckPalette(sort(ndata), cropcol=cropcol)
		}else{
			nbnonzero <- sum(!is.na(data))
			mat <- 1:length(data);
            mat[order(data)[1:nbnonzero]] <- equalCheckPalette(sort(data)[1:nbnonzero], cropcol=cropcol)
			mat[is.na(data)] <- "#E0E0E040"
			tmp <- order(data);
			valrange <- 1:nbnonzero;
			tmp <- order(data)
			for(i in 1:nbnonzero) valrange[i] <- data[tmp[i]]
			valcolrs <- equalCheckPalette(sort(data)[1:nbnonzero], cropcol=cropcol )
		}
	}
	
	return(list(col=mat,non.na=nbnonzero,valrange=valrange, valcolrs= valcolrs,hasdist=T))
}

makeColorbar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    #dev.new(width=4, height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

subClusterPCA <- function(sro, meta.partition, nb.var=1000,nb.pc=50, reduction.get= c()  ){ #savxxetxest
	library(Seurat)
	library(M3Drop)
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Seurat object expected!")
	cls <- getSafeLevels(sro@meta.data[[meta.partition]])
	scada <- matrix(0, ncol(sro@raw.data),nb.pc)
	rownames(scada) <- colnames(sro@raw.data)
	tmpgenenames <- 1:nrow(sro@raw.data)
	mainout <- M3DropFeatureSelection(sro@raw.data, mt_method = "fdr",mt_threshold = 0.5)
	if (nrow(mainout) > nb.var) mainout <- rownames(mainout)[order(mainout$q.value)[1:nb.var]]
	else mainout <- rownames(mainout)
	print(paste("got", length(mainout), "main variable genes"))
	rownames(sro@raw.data) <- tmpgenenames
	dakey <- c()
	for(i in 1:length(cls)){
		flt <- sro@meta.data[[meta.partition]] == cls[i]
		print(paste("processing", cls[i], "with",sum(flt),"cells"))

		srom <- CreateSeuratObject(sro@raw.data[,flt,drop=F])
		rownames(srom@raw.data) <- 1:nrow(sro@raw.data)
		mout <- M3DropFeatureSelection(srom@raw.data, mt_method = "fdr",mt_threshold = 0.5)
		if (nrow(mout) > nb.var) mout <- rownames(mout)[order(mout$q.value)[1:nb.var]]
		else mout <- rownames(mout)
		print(paste("got", length(mout), "variable genes"))
		mout <- unique(c(mout,mainout))
		if (length(mout) > nb.var) mout <- mout[1:nb.var]
		srom <- NormalizeData(srom, normalization.method = "LogNormalize",display.progress=F)
		srom <- ScaleData(srom, vars.to.regress=c("nUMI"), genes.use=mout)
		print("Running PCA")
		srom <- tryCatch(RunPCA(srom,pc.genes=mout, pcs.compute = nb.pc,do.print=F), error=function(cond){print("failed! reducing dims");return(tryCatch(RunPCA(srom,pc.genes=mout, pcs.compute = nb.pc/2,do.print=F), error=function(cond){print("failed! reducing dims");return(tryCatch(RunPCA(srom,pc.genes=mout, pcs.compute = nb.pc/10,do.print=F), error=function(cond){print("failed! abort...");return(srom)}))}))})
		if ("pca" %in% names(srom@dr)){
			print("Saving coordinates")
			tmptmp <- dim(scada)[2] - dim(srom@dr$pca@cell.embeddings)[2]
			if (tmptmp != 0) {
				 srom@dr$pca@cell.embeddings <- cbind(srom@dr$pca@cell.embeddings, matrix(0,nrow(srom@dr$pca@cell.embeddings), tmptmp))
				 colnames(srom@dr$pca@cell.embeddings) <- colnames(scada)
			}
			scada[flt,] <- srom@dr$pca@cell.embeddings
			if (is.null(dakey)){
				dakey <- GetDimReduction(srom,reduction.type = "pca", slot = "key") 
				colnames(scada) <- colnames(sro@dr$pca@cell.embeddings)
			}
		}
	}

	if (is.null(reduction.get)) reduction.get <- paste(meta.partition,"subpca",sep='_')
	sro <- SetDimReduction(sro, reduction.type = reduction.get, slot = "cell.embeddings", new.data = scada)
	sro <- SetDimReduction(sro, reduction.type = reduction.get, slot = "key", new.data = dakey)
return(sro)}

subClusterCorrect <- function(sro, meta.partition, meta.batch, reduction.use= c(), reduction.get= "mnnpca", nb.pc=0){ #return matrix since... dlls stuff
	library(scran)
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Seurat object expected!")
	if (is.null(reduction.use)) {
		reduction.use <- paste(meta.partition,"subpca",sep='_')
		if (!reduction.use %in% names(sro@dr)) reduction.use = "pca"		
	}
	if (!reduction.use %in% names(sro@dr)) stop(paste("Reduced dimention", reduction.use, "is missing!"))
	if (nb.pc == 0) nb.pc <- ncol(sro@dr[[reduction.use]]@cell.embeddings)	

	if (is.factor(levels(sro@meta.data[[meta.partition]]))) cls <- levels(sro@meta.data[[meta.partition]])
	else cls <- unique(sro@meta.data[[meta.partition]])
	if (is.factor(levels(sro@meta.data[[meta.batch]])))	blc <- levels(sro@meta.data[[meta.batch]])
	else bls <- unique(sro@meta.data[[meta.batch]])

	fout <- matrix(0,nrow(sro@meta.data),nb.pc);
	rownames(fout) <- rownames(sro@meta.data)
	colnames(fout) <- paste("MNNPC", 1:nb.pc,sep="_")
	mineighbor=20
	for(i in 1:length(cls)){
		mnninput <- list();
		data <- list()
		mflt <- sro@meta.data[[meta.partition]] == cls[i]
		print(paste("processing", cls[i], "with",sum(mflt),"cells"))	
		for(j in 1:length(bls)){
			flt <- mflt & (sro@meta.data[[meta.batch]] == bls[j])
			print(sum(flt))
			newin <- sro@dr[[reduction.use]]@cell.embeddings[flt,,drop=FALSE]
			if (nrow(newin) >= mineighbor) mnninput <- c(mnninput, list(t(newin)))
		}
		if (length(mnninput) > 0) mnnout <- do.call(mnnCorrect, c(mnninput, list(k=mineighbor)))
		newin <- 0
		for(j in 1:length(bls)){
			flt <- mflt & (sro@meta.data[[meta.batch]] == bls[j])
			if (sum(flt) >= mineighbor) fout[flt,] <- mnnout$corrected[[j - newin]]
			else { # not enough cells... no batch correction!
				 newin <- newin+1
				if (sum(flt) > 0) {
					fout[flt,] <- sro@dr[[reduction.use]]@cell.embeddings[flt,,drop=FALSE]
					print(paste("no correction for",cls[i], "in batch", bls[j],"! not enough cells"))
				}
			}
		}
		# re-orthogonalize!
		tmpdata <- fout[mflt,]
		decom <- svd(tmpdata)
		fout[mflt,] <- (decom$u %*% diag(decom$d))
	}
  return(list(sro=sro, cell.embeddings = fout, key = "MNNPC"))	
}
subCluster <- function(sro, meta.partition, reduction.use= "pca", nb.pc=0,resolution=0.8,k=30,reduction_key = "SubTSNE_", todo.index.range=c()){	
	library(Seurat)
	if (class(sro) == "list"){
		if ("sro" %in% names(sro)) {
			reduction.use <- paste(meta.partition,"submnnpca",sep='_') 
			sro$sro <- SetDimReduction(sro$sro, reduction.type = reduction.use, slot = "cell.embeddings", new.data = sro$cell.embeddings)
			sro <- SetDimReduction(sro$sro, reduction.type = reduction.use, slot = "key", new.data = sro$key)
		}
		rt <- paste(meta.partition,"_mnn_subtsne", sep="")
		clt <- paste(meta.partition, "mnn_subcluster", sep="_")

	}else{
		rt <- paste(meta.partition,"_subtsne", sep="")
		clt <- paste(meta.partition, "subcluster", sep="_")
	}
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("not a seurat object!")
	cls <- levels(sro@meta.data[[meta.partition]])
	if (nb.pc == 0) nb.pc <- ncol(sro@dr[[reduction.use]]@cell.embeddings)
	newcl <- data.frame(row.names= rownames(sro@meta.data))
	tsneout <- matrix(0,nrow(sro@dr[[reduction.use]]@cell.embeddings), 2)
	colnames(tsneout) <- paste(reduction_key , 1:2,sep="")
	rownames(tsneout) <- rownames(sro@meta.data)
	if (is.null(todo.index.range)) todo.index.range <- 1:length(cls)
	else{
		if (rt %in% names(sro@dr)) tsneout <- GetDimReduction(sro,reduction.type = rt, slot = "cell.embeddings")
		if (clt %in% colnames(sro@meta.data)) newcl[[1]] <- sro@meta.data[[clt]]
	}
	for( i in todo.index.range){
		flt <- sro@meta.data[[meta.partition]] == cls[i]
		print(paste("processing", cls[i], "with",sum(flt),"cells"))	
		srom <- CreateSeuratObject(sro@raw.data[,flt])
		srom <- SetDimReduction(srom, reduction.type = rt, slot = "cell.embeddings", new.data = sro@dr[[reduction.use]]@cell.embeddings[flt,])
		srom <- SetDimReduction(srom, reduction.type = rt, slot = "key", new.data = GetDimReduction(sro,reduction.type = reduction.use, slot = "key") )
		srom <- FindClusters(srom, reduction.type=rt,dims.use=1:(nb.pc),resolution=resolution,print.output=F, k.param=k,force.recalc=TRUE)
		print("TSNE TIME!")
		srom <- RunTSNE(srom, reduction.use=rt, dims.use=1:(nb.pc),check_duplicates=FALSE,force.recalc=TRUE)
		newcl[rownames(srom@meta.data), 1] <- srom@meta.data$res.0.8
		tsneout[rownames(srom@meta.data) ,] <- srom@dr$tsne@cell.embeddings
	}
	sro <- SetDimReduction(sro, reduction.type = rt, slot = "cell.embeddings", new.data = tsneout)
	sro <- SetDimReduction(sro, reduction.type = rt, slot = "key", new.data = reduction_key)
	sro@meta.data[[clt]] <- newcl[[1]]
return(sro)}

#' Function that return list colors for plot to be used
#'
#' @param list of numbers
#'
#' @export
subTSNE <- function(sro, meta.partition, reduction.use= "pca", nb.pc=0 ,reduction_key = "SubTSNE_", perplexity = 30){
	if (nb.pc == 0) nb.pc <- ncol(sro@dr[[reduction.use]]@cell.embeddings)
	library(Rtsne)
	if (class(sro@meta.data[[meta.partition]]) == "factor") flvl <- levels(sro@meta.data[[meta.partition]])
	else flvl <- unique(sro@meta.data[[meta.partition]])

	if (!reduction.use %in% names(sro@dr)) stop(paste("Missing '", reduction.use, "' in dr slot",sep=""))

	tsneout <- matrix(0,nrow(sro@dr[[reduction.use]]@cell.embeddings), 2)
	colnames(tsneout) <- paste(reduction_key , 1:2,sep="")
	for(i in 1:length(flvl)){
		input <- sro@dr[[reduction.use]]@cell.embeddings[sro@meta.data[[meta.partition]] == flvl[i],,drop =F]
		if (var(diag(var(input))) == 0){
			print(paste("warning, subcluster", flvl[i], "has no dimentions"))
		}else{
			print(paste("Processing",flvl[i]))

			if (nrow(input)  > 1){
				if (ncol(input) >= nrow(input)){ # also prevents to interpret input as distance matrix
					input <- input[, 1:(nrow(input)-1),drop =F] 
				}
				input <- as.matrix(input)
				tmpout <- tryCatch(Rtsne(X = input,dim=2,perplexity = perplexity,check_duplicates=F), error=function(cond){return(
				  tryCatch(Rtsne(X = input,dim=2,perplexity = perplexity/2,check_duplicates=F), error=function(cond){return(
				  tryCatch(Rtsne(X = input,dim=2,perplexity = perplexity/4,check_duplicates=F), error=function(cond){return(
				  tryCatch(Rtsne(X = input,dim=2,perplexity = perplexity/8,check_duplicates=F), error=function(cond){return(
				  tryCatch(Rtsne(X = input,dim=2,perplexity = 1,check_duplicates=F) , error=function(cond){print(paste("Warning, failed to compute TSNE for", flvl[i])); return( list(Y= matrix(0,nrow(input),2))
                                  )}))}))}))}))})
			tsneout[sro@meta.data[[meta.partition]] == flvl[i],] <- tmpout$Y

		}
		}
	}
	rt <- paste(meta.partition,"_subtsne", sep="")
	rownames(tsneout) <- rownames(sro@dr[[reduction.use]]@cell.embeddings)
	sro <- SetDimReduction(sro, reduction.type = rt, slot = "cell.embeddings", new.data = tsneout)
	sro <- SetDimReduction(sro, reduction.type = rt, slot = "key", new.data = reduction_key)
return(sro)}
makeEuclidienTsne <- function(coors){
	library(Rtsne)
	prod <- t(coors) %*% coors
	

}
createEqualSamples <- function(sro, nbparts = 8, ordering, meta.sample, meta.celltype, cell.flt=c(), cell.meta=c()){

	newannot <- rep("NA", nrow(sro@meta.data))
	cell.flt <- SeuratCellListQuery(sro, cell.meta, cell.flt)

	for(i in levels(sro@meta.data[[meta.sample]])){
		for( j in levels(sro@meta.data[[meta.celltype]])){
			flt <- (sro@meta.data[[meta.sample]] == i)&(sro@meta.data[[meta.celltype]] == j)&(cell.flt)
			ord <- order(ordering[flt])
			tag <- paste(i, rep(1:nbparts, (sum(flt)/nbparts)+1)[1:sum(flt)], sep="_subsmp") 
			newannot[flt][ord] <- tag
		}
	}
return(newannot);}

AUCellOverlay <- function(sro, genelistpos=c(), genelistneg=c(), use.cell = c(), use.meta=c() ){

	if (is.null(genelistpos)) genelistpos = c("CADPS2","SATB2","NRGN","SV2B","LINC01378","RORB","ANO3","SLC17A7","LY86.AS1","TESPA1","FOXP2","MEIS2","RMST","IPCEF1","PHYHIP","MAST3","PCDH9.AS2","RXFP1","TYRO3","KLHL2","NGEF","SYT7","LGI4","SLC22A10","PART1","C1orf115","CA11","NPTX1","HS3ST4","HOPX")

	if (is.null(genelistneg)) genelistneg = c("ERBB4","ADARB2","GRIK1","CXCL14","NXPH1","SLC6A1.AS1","GAD1","SOX6","ZNF536","RELN","DLX6.AS1","SLC6A1","GAD2","PLD5","DOCK10","SLC24A3","TAC1","VIP","BTBD11","SST","ALK","KIT","FBXL7","ARL4C","RBMS3","MAF","GRIP2","PROX1","ANK1","IGF1","UBASH3B","ZNF608","CALB2","EGFR","BACH1","TOX3","COL25A1","WLS","MTUS1","GRIN2D","ARX")

	library(AUCell)
	use.cell <- SeuratCellListQuery(sro, use.meta, use.cell)
	
	exprMatrix <- as.matrix(sro@raw.data[c(genelistpos,genelistneg),use.cell])
	cells_rankings <- AUCell_buildRankings(exprMatrix)

	geneSets <- list(geneSet1=genelistpos,geneSet2=genelistneg)
	cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5, nCores=1)

	table <- t(cells_AUC@assays[["AUC"]])
	value <- (table[,1] - table[,2]) / (table[,1] + table[,2])

	return(value)
}

#' Function that return list colors for plot to be used
#'
#' @param list of numbers
#'
#' @export
myPlot <- function(sro,query="",dimred.name="tsne", transform="value", denominator=c(), query2=c(), do.legend=T,do.label=T,pt.size=1,cell.use=c(),size=1, meta.use=c(), cell.set.na=c(), meta.set.na=c(), plot.attribs=c(),color.palette=c("fireice","fire","ice"),na.alpha=0.1,na.color="#888888",do.density=F,do.rgl=c(), bypass.size=c(), bypass.alpha=c(), do.zerocenter.color=F, clamp.value.range =c()){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if (length(query) > 1){
		plotdata <- SeuratCellListQuery(sro, query)
		colortype <- "m"
		use.color <- list(T= "#00FF00", F= "#FF0000")
		use.shape <- rep(16,2)
		query <- paste( query)
	}else{
 		if (query != ""){
 		if (query %in% names(sro@meta.data)) {
		col.lvl <- getSafeLevels(sro@meta.data[[query]])
		if ("z.flt.empty" %in% flags.plot) {mmm <-table(sro@meta.data[[query]]) ; col.lvl <- col.lvl[mmm[match(col.lvl, names(mmm))] > 0]}
		if (("meta.color" %in% names(sro@misc))&&(query %in% names(sro@misc$meta.color))) use.color <- sro@misc$meta.color[[query]]
 		else use.color <- myrainbow(length(col.lvl))
 		if (("meta.shape" %in% names(sro@misc))&&(query %in% names(sro@misc$meta.shape))) use.shape <- sro@misc$meta.shape[[query]]
 		else use.shape <- rep(16,length(col.lvl))
 	 	if (("meta.alpha" %in% names(sro@misc))&&(query %in% names(sro@misc$meta.alpha))) use.alpha <- sro@misc$meta.alpha[[query]]
 		else use.alpha <- rep(1,length(col.lvl))	
		  if (class(sro@meta.data[[query]]) == "factor"){
 		    colortype <- "m"
			if ("z.flt.empty" %in% flags.plot) plotdata <- factor(sro@meta.data[[query]], levels=col.lvl)
	                else plotdata <- sro@meta.data[[query]]
                   }else if (class(sro@meta.data[[query]]) == "character"){
 		    plotdata <- factor(sro@meta.data[[query]], levels=col.lvl)
                     colortype <- "m"
                   }else {
                     plotdata <- sro@meta.data[[query]]
                     colortype <- "i"
 		  }
## 		}else if (que5432k403290ry %in% rownames(sro@data)){
 	##		plotdata <- sro@data[query,]
  	##		plotdata[plotdata == 0] <- NA			
	##		colortype <- "i"
 		}else if (query %in% SeuratRawData(sro, "R")){
			if (transform == "R") plotdata <- SeuratRawData(sro)[query,]
 			else if (transform == "L") plotdata <- log2(SeuratRawData(sro)[query,] +1)
			else plotdata <- SeuratRawData(sro, "N")[query,] # log2(sro@raw.data[query,] +1) 
 			plotdata[plotdata == 0] <- NA
 			colortype <- "i"
                 } else stop(paste("Invalid query Seurat object has no", query ,"in meta.data or data rowname"))
		}else use.color <-myrainbow(length(unique(sro@ident)))
	}
	
	if (!is.null(query2)){
 		if (query2 %in% names(sro@meta.data)) {
		  if (class(sro@meta.data[[query2]]) == "factor"){
                     plotdata2 <- sro@meta.data[[query2]]
                   }else if (class(sro@meta.data[[query2]]) == "character"){
 		    plotdata2 <- as.factor(sro@meta.data[[query2]])
                   }else {
                     plotdata2 <- sro@meta.data[[query2]]
 		  }
 		}else if (query2 %in% SeuratRawData(sro, "R")){
			if (transform == "R") plotdata2 <- SeuratRawData(sro)[query2,]
 			else if (transform == "L") plotdata2 <- log2(SeuratRawData(sro)[query2,] +1)
			else plotdata2 <- SeuratRawData(sro, "N")[query2,] # log2(sro@raw.data[query,] +1) 
 			plotdata2[plotdata2 == 0] <- NA
                 } else stop(paste("Invalid query Seurat object has no", query ,"in meta.data or data rowname"))
	}

	
	
	cell.use <- SeuratCellListQuery(sro, meta.use, cell.use)
	cell.set.na <- SeuratCellListQuery(sro, meta.set.na, cell.set.na,default=F)
	cell.use <- (1:length(cell.use)) %in% sample(which(cell.use))
	
	if ("dr" %in% names(attributes(sro))) drslot = "dr"
	else drslot = "reductions"

	if (colortype != "m"){
		
		plotdata[cell.set.na] <- NA
		plotdata <- plotdata[cell.use]
		if (!is.null(denominator)){
			query <- paste(query, "/",denominator)
			if (denominator %in% names(sro@meta.data)) {
				plotdata <- plotdata / 	sro@meta.data[[denominator]][cell.use]
			}	
		}

		if (transform == "log2p1") plotdata <- log2(plotdata + 1)
		else if (transform == "log2") {plotdata <- log2(plotdata); plotdata[is.infinite(plotdata)] <- NA}
		else if (transform == "log10") {plotdata <- log10(plotdata); plotdata[is.infinite(plotdata)] <- NA}
		else if (transform == "log") {plotdata <- log(plotdata); plotdata[is.infinite(plotdata)] <- NA}

		if (!(is.null(clamp.value.range))){
			plotdata[plotdata < clamp.value.range[1]] <- clamp.value.range[1]
			plotdata[plotdata > clamp.value.range[2]] <- clamp.value.range[2]
		}

	}else{
		plotdata <- plotdata[cell.use]
	}

	if ((!is.null(do.rgl))|(dim(slot(sro,drslot)[[dimred.name]]@cell.embeddings)[2] == 3)){
		library(rgl)
		color = rep("#000000", sum(cell.use))

		if (dim(slot(sro,drslot)[[dimred.name]]@cell.embeddings)[2] == 2) zdim = sro@meta.data$orig.ident[cell.use]@.Data
		else zdim = slot(sro,drslot)[[dimred.name]]@cell.embeddings[cell.use,3];

		if (colortype == "m"){
			color <- use.color[plotdata]
		}else{
			aur <- range(plotdata,na.rm=T)
			color <- (plotdata - aur[1]) / (aur[2] - aur[1])
			flt = !is.na(color)
			color[flt] <- myFireAndIce(color[flt])
			color[!flt] = "#AAAAAA"
		}
		plot3d(x = slot(sro,drslot)[[dimred.name]]@cell.embeddings[cell.use,1], y= slot(sro,drslot)[[dimred.name]]@cell.embeddings[cell.use,2], z= zdim , col=color)
		return;
	}


	if (is.null(plot.attribs)) plot.attribs = list(flags=c())




        gdata <- data.frame(slot(sro,drslot)[[dimred.name]]@cell.embeddings[cell.use,])

	colnames(gdata) <- c("X","Y")



	gdata[["A"]] <- plotdata
	
	if (colortype == "m"){
		p <- ggplot(gdata, aes(x=X,y=Y,color=A,shape=A, alpha=A))
		p <- p + scale_color_manual(name=query,values=use.color,label=col.lvl,drop = FALSE)
		p <- p + scale_shape_manual(name=query,values=use.shape,label=col.lvl,drop = FALSE)
		p <- p + scale_alpha_manual(name=query,values=use.alpha,label=col.lvl,drop = FALSE)
		p <- p + guides(colour = guide_legend(override.aes = list(shape = 15,size = 5,alpha=1)))
	}else{

	

		if (do.density){
			library(InferN0)
			gdata[["A"]][is.na(gdata[["A"]])] <- 0
			output <- InferN0GetDensity(as.matrix(gdata[,1:2]),data=gdata[["A"]], mapsize=256, bandwidth=0.02); 
			return(plotColoredDensity(output$Image,output$Density,isbg=T,plot.attribs=plot.attribs))
		}else{
			p <- ggplot(gdata, aes(x=X,y=Y,color=A, alpha=A))
			if (do.zerocenter.color) {
				daccrange <- range(gdata[["A"]],na.rm=T)
				if (is.infinite(daccrange[1])) daccrange = 1:41
				else{
					if (-daccrange[1] > daccrange[2]) daccrange = 1:(21+ floor(-20 *daccrange[2] /daccrange[1]))
					else daccrange = (21- floor(-20 *daccrange[1] /daccrange[2])):41
				}
			}else daccrange = 1:41

			if ((length(color.palette) != 1)||(color.palette=="fireice")) colpal <- colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]
			else if (color.palette == "fire") colpal <- colorRampPalette(c("#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]
 	       	 else if (color.palette == "ice") colpal <- colorRampPalette(rev(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000")))(41)[daccrange]
			else colpal <- colorRampPalette(c("#FFFFFF", "#000000"))(41)[daccrange]
	
			p <- p + scale_color_gradientn(name=transform,colours=colpal, na.value= na.color)
			p <- p + scale_alpha_continuous(position=NULL,guide="none", na.value=na.alpha, range = c(1, 1))
		}
	}
	if (!"title"  %in% names(plot.attribs)) plot.attribs$title = query
	if (!is.null(bypass.alpha)) {
		p <- p + geom_point(alpha = bypass.alpha[cell.use],size = bypass.size[cell.use])
	}else p <- p + geom_point(size = size)
#	p <- p + guides(colour = guide_legend(override.aes=list(size=5)))
	return(changeStyle(p,plot.attribs)) 	
}

makeSummarySeurat <- function(sro, metaattriblist=c("nUMI", "nGene"), metaattribalias=c(), metafactorlist=c("orig.ident") ){

	if (length(alias) != length( metaattriblist)) alias = metaattriblist

	nbclass <- length(levels(sro@ident))
	nbattrib = length(metaattriblist);
	nbcol = 1 + length(metaattriblist) *4
	
	dacolname = c("nbCell")
	for( j in 1:length(metaattriblist)){
		dacolname <- c(dacolname, paste("Mean", alias[j],sep='_'), paste("Min", alias[j],sep='_'), paste("Median", alias[j],sep='_'), paste("Max", alias[j],sep='_'))
	}

	for(i in 1:length(metafactorlist)){
		nbcol <- nbcol + length(unique(sro@meta.data[[metafactorlist[i]]]))
		dacolname <- c( dacolname, paste(metafactorlist[i],as.character(unique(sro@meta.data[[metafactorlist[i]]])),sep='_'))
		#dacolname <- c( dacolname, as.character(unique(sro@meta.data[[metafactorlist[i]]])))
	}

	data <- matrix(0,nbclass+1, nbcol)
	colnames(data) <- dacolname
	rownames(data) <- c(levels(sro@ident) , "All");

	

	data[(nbclass+1),1] = length(sro@ident@.Data)
	for( j in 1:length(metaattriblist)){
	        data[(nbclass+1),j*4-2] = mean(sro@meta.data[, c(metaattriblist[j])])
                data[(nbclass+1),j*4-1] = min(sro@meta.data[, c(metaattriblist[j])])
                data[(nbclass+1),j*4] = median(sro@meta.data[, c(metaattriblist[j])])
                data[(nbclass+1),j*4+1] = max(sro@meta.data[, c(metaattriblist[j])])
	}
	
	for( i in 1:nbclass){
		data[i,1] = sum(sro@ident@.Data == i)
		for( j in 1:length(metaattriblist)){
			data[i,j*4-2] = mean(sro@meta.data[(sro@ident@.Data == i), c(metaattriblist[j])])
			data[i,j*4-1] = min(sro@meta.data[(sro@ident@.Data == i), c(metaattriblist[j])])
			data[i,j*4] = median(sro@meta.data[(sro@ident@.Data == i), c(metaattriblist[j])])
			data[i,j*4+1] = max(sro@meta.data[(sro@ident@.Data == i), c(metaattriblist[j])])
		}
	}
	curp = 1+4*length(metaattriblist);	
	for( k in 1:length(metafactorlist)){
		dalist = unique(sro@meta.data[[metafactorlist[k]]])
		for(j in 1:length(dalist)){
			if (is.na(dalist[j])) {
			        data[(nbclass+1),curp+j] = sum(is.na(sro@meta.data[[metafactorlist[k]]]))
                                for( i in 1:nbclass){
                                       data[i,curp+j] = sum((sro@ident@.Data == i)&(is.na(sro@meta.data[[metafactorlist[k]]])))
                                }
			}else{
				data[(nbclass+1),curp+j] = sum(sro@meta.data[[metafactorlist[k]]] == dalist[j],na.rm=T) 
				for( i in 1:nbclass){
					data[i,curp+j] = sum((sro@ident@.Data == i)&(sro@meta.data[[metafactorlist[k]]] == dalist[j]) ,na.rm=T)
				}
			}
		}
		curp <- curp + length(dalist)
	}


return(data)}

getHCYcolor <- function(h,c,y,g=2.2,opt.hue=T, alphablend = c()){
    if (length(h) > 1){
        fout <- rep(0, length(h))
	for(i in 1:length(h)) {
		fout[i] <- getHCYcolor(h[i],c[i],y[i],g,opt.hue=opt.hue,alphablend=alphablend)
	}
	return(fout)
    }
    
    if (opt.hue) h <- h - sin(h *18.8495559215) / 24
    
    y <- y ^ g

    res = (h %% 1) * 6
    bas = floor(res)
    res = res - bas
    fout <- switch(as.character(bas), "0"=c(1.0, 1.0 - c * (1-res), 1.0 - c), "1"=c(1.0 - c * res, 1.0, 1.0 - c), "2"=c(1.0 - c, 1.0, 1.0 - c* (1-res) ),"3"=c(1.0-c,1.0 - c*res, 1.0 ), "4"=c(1.0 - c * (1-res), 1.0 - c,1.0),"5"=c(1.0, 1.0 - c, 1.0-c*res)  )
    brscale = fout[1] * 0.3 + fout[2] * 0.59 + fout[3] * 0.11
    if (brscale < y){
	# need to project to be within, fit parabola f[0] = brscale, f[1] = 1 f'[x] =:
	if ((fout[1] < 1)&(fout[1] >0.00001)) {tmp = fout[1]; fout[1] = fout[1] ^ (1/g) ; deriv = 0.3*(1/fout[1]-1)*tmp; }
	else deriv =0
        if ((fout[2] < 1)&(fout[2] >0.00001)) {tmp = fout[2]; fout[2] = fout[2] ^ (1/g) ; deriv = deriv + 0.59*(1/fout[2]-1)*tmp}
	if ((fout[3] < 1)&(fout[3] >0.00001)) {tmp = fout[3]; fout[3] = fout[3] ^ (1/g) ; deriv = deriv + 0.11*(1/fout[3]-1)*tmp}
	deriv = deriv * g
	discr = deriv * deriv - 4.0 * (brscale -y) * (1.0 -brscale - deriv)
	brscale = (-deriv + sqrt(discr)) / (2.0* (1.0 -brscale - deriv))
	#brscale = (y - brscale) / (1- brscale)
	fout[1] = fout[1] + (1 -fout[1]) * brscale
	fout[2] = fout[2] + (1 -fout[2]) * brscale
	fout[3] = fout[3] + (1 -fout[3]) * brscale
	brscale = (fout[1]^g) * 0.3 + (fout[2]^g) * 0.59 + (fout[3]^g) * 0.11
	if (is.null(alphablend)) return(rgb((fout[1]), (fout[2]),(fout[3])))  
	else return(rgb(alphablend[4] * alphablend[1] + (1 - alphablend[4]) * (fout[1]),alphablend[4] * alphablend[2] + (1 - alphablend[4]) *  (fout[2]),alphablend[4] * alphablend[3] + (1 - alphablend[4]) * (fout[3])))
    }else{
	brscale =  y/brscale
        g <- 1 / g
	if (!is.null(alphablend)) return(rgb(alphablend[4] * alphablend[1] + (1 - alphablend[4]) * ((fout[1]* brscale) ^ g), alphablend[4] * alphablend[2] + (1 - alphablend[4]) * ((fout[2]* brscale) ^ g),alphablend[4] * alphablend[3] + (1 - alphablend[4]) * ((fout[3]* brscale) ^ g)))
    	else return(rgb(((fout[1]* brscale) ^ g), ((fout[2]* brscale) ^ g),((fout[3]* brscale) ^ g)))
	}
}

getIntPhaseColor <-function(int,phase){
	grayf <- abs(phase)
	gray <- (1- grayf) * (int * 7 / 10) 
	lum <- int 
	if (phase <0){
		r=0
		g = (7.2*lum - 1) / 5
		if (g < 0) {b = (lum*7.2)^0.4545;g=0}
		else {b=1; g = (0.8*g) ^0.4545;}
	}else{
		b=0
		g = (7 * lum -3) / 4
		if (g < 0) {r = (lum*7/3)^0.4545; g=0}
		else {r=1; g = (0.65*g) ^ 0.4545;}
	}
#	print(0.11*(b^2.2) + 0.59*(g^2.2) + 0.3*(r^2.2))
#	print(int * 7/10)
return(rgb(gray+r*grayf,gray+g*grayf,gray+b*grayf))}

plotColoredDensity <- function(intensity, phase, density=c(), plot.attribs=c(), ispair=F, isbg=F,densitysat=4, pair.rescale=T, trythis=F){
	dd <- dim(intensity)
	gdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 4))
	dacolor <- rep(0,dd[1] * dd[2])
	daalpha <- rep(1,dd[1] * dd[2])
	if (trythis){
		intensity[is.na(intensity)] <- 0
		phase[is.na(phase)] <- 0

		tmp <- intensity + phase
		phase <- (intensity - phase) / tmp
	
	
		intensity <- log2(tmp)
		intensity[intensity < 0]  <- 0
		intensity[is.na(intensity)] <- 0	
		
	
		ir <- range(intensity)
		intensity <- 0.1 + 0.8 * intensity / ir[2]
		phase[is.na(phase)] <- 0


		phase[phase> 1]	<- 1
		phase[phase< -1] <- -1
		print(range(phase))
		print(dim(intensity))	
		print(sum(is.na(phase)))
		print(range(intensity))
		dr <- range(density)
		density <- density*(2 /dr[2])
		density[density > 1] <- 1 
		density <- (1- density) 
		for(j in 1:ncol(intensity)){
			for(i in 1:nrow(intensity)){
				offset <- (i-1+ (j-1) * dd[1])
				if (phase[j,i] <0){
					dacolor[1+(offset)] <- getHCYcolor(0.57, -phase[j,i], intensity[j,i], alphablend=c(0.0,0,0.0, density[j,i]  ))
				}else{
					dacolor[1+(offset)] <- getHCYcolor(0.02, phase[j,i], intensity[j,i], alphablend=c(0.0,0,0.0, density[j,i] ))
				}
			}
		}
	}else{
	if (ispair){
		if (pair.rescale) {
			tmp <- sum(phase,na.rm=T)
			phase <- phase * sum(intensity,na.rm=T)
			intensity <- intensity * tmp
		}
		tmp <- intensity + phase
		phase <- (intensity - phase) / tmp
		phase[is.na(phase)] <- 0
		daval <- range(tmp)
		intensity <- tmp / daval[2]
	}
	if (isbg){
		ir <- range(intensity)
		dr <- range(phase)
		for(j in 1:ncol(intensity)){
			for(i in 1:nrow(intensity)){
				offset <- (i-1+ (j-1) * dd[1])+1
				tmp <- (intensity[j,i]-ir[1]) * 0.6 / (ir[2]-ir[1])
				dacolor[offset] <- getHCYcolor(0 ,1,tmp)
				daalpha[offset] <- densitysat * phase[j,i] / dr[2]
				if (daalpha[offset] > 1) daalpha[offset] =1
			}
		}
	}else{


	for(j in 1:ncol(intensity)){
		for(i in 1:nrow(intensity)){
			offset <- (i-1+ (j-1) * dd[1])
			if (phase[j,i] <0){
				dacolor[1+(offset)] <- getHCYcolor(0.57, -phase[j,i], intensity[j,i] * 0.75)
			}else{
				dacolor[1+(offset)] <- getHCYcolor(0.02, phase[j,i], intensity[j,i] * 0.75)
			}

			}
		}
	}
	}




#	p <- ggplot(data = gdata,mapping=aes(y=Y, x=X,group=I) )	
#	p <- p + geom_polygon(mapping=aes(group = I, y=Y, x=X), fill = rep(dacolor,each=4))
#data = gdata,mapping=aes(y=Y, x=X,group=I) )	
	p <- ggplot(data.frame(expand.grid(1:nrow(intensity),ncol(intensity):1)),aes(Var1,Var2))+ geom_tile(fill = dacolor, alpha= daalpha)
return(changeStyle(p,plot.attribs))}

#' Function that return list colors for plot to be used
#'
#' @param list of numbers
#' @param xform: transformed data rendered default: SeuratNormalized 'T' log10TPM 'R' raw.count 'S' SeuratNormalized&scaled
#'
#' @export
plotSeuratFeat <- function(sro,query, xform= " ", title=NULL,threshold.value=NA,save.pdf=NULL,save.png=NULL, pointscale=1, cropcol= 0, coor=c(), cell.use = c(), colors.use=c(),do.resize=c()){
	#getOption("device")()
	if (!is.null(do.resize)) setPlotSize(width=do.resize[1], height=do.resize[2])

	if (class(sro) == "seurat"){

	if (is.null(title)) title = query

	if (is.null(cell.use)){
		cell.use <- rownames(sro@meta.data);
	} 


	if (length(query) != 1){ # use the query itself!
				dacolors = Infern0Colors(as.vector(query),threshold.value=threshold.value)
        ylab <- "value"
				title <- ""
	}else if (query %in% names(sro@meta.data)){
		if (!is.null(colors.use)){
			dacolors = list(hasdist=F, valcolrs= colors.use)
			dacolors$col = sro@meta.data[[query]] 
		}else dacolors = Infern0Colors(sro@meta.data[cell.use,query],threshold.value=threshold.value)
        ylab <- "value"
	}else{
				data <- switch(substr(xform,1,1),
					T = log10(sro@raw.data[c(query),cell.use]) - log10(sro@meta.data$nUMI) + 6,
					S = sro@scale.data[c(query),cell.use],
					R = sro@raw.data[c(query),cell.use],
					L = log10(sro@raw.data[c(query),cell.use]),
					sro@data[c(query),cell.use]);
				ylab <- switch(substr(xform,1,1),
					T = "Log10 TPM",
					S = "Regressed count",
					R = "UMI count",
					L = "log10 UMI count",
					"Norm.\ncount");
        data[sro@raw.data[c(query),cell.use] == 0] <- NA
        dacolors = Infern0Colors(data,threshold.value=threshold.value, cropcol=cropcol)
	}
    isna <- (col=dacolors$col == "#E0E0E040")
    if (is.null(coor)){
        xcoor <- sro@dr$tsne@cell.embeddings[cell.use,1]
        ycoor <- sro@dr$tsne@cell.embeddings[cell.use,2]
	xtit <- ""
	ytit <- ""
    }else if (length(coor) == 2) {
	xcoor <- sro@meta.data[cell.use,coor[1]]
	ycoor <- sro@meta.data[cell.use,coor[2]]
	xtit <- coor[1]
	ytit <- coor[2]
    }else if (coor == "soup"){
        xcoor <- sro@meta.data[cell.use,"soup_rank"]
	ycoor <- log10(sro@meta.data[cell.use,"nUMI"])
	xtit <- "soup_rank"
	ytit <- "log10_nUMI"
    }else if (coor == "raw.soup") {
        xcoor <- sro@meta.data[cell.use,"soup_fraction"]
	ycoor <- log10(sro@meta.data[cell.use,"nUMI"])
	xtit <- "soup_rank"
	ytit <- "log10_nUMI"
    }else if (coor == "pca"){
        xcoor <- sro@dr$pca@cell.embeddings[cell.use,1]
        ycoor <- sro@dr$pca@cell.embeddings[cell.use,2]
	xtit <- "PC1"
	ytit <- "PC2"
    }else if (coor == "coverage"){
        xcoor <- log10(sro@meta.data[cell.use,"nGene"])
        ycoor <- log10(sro@meta.data[cell.use,"nUMI"])
	xtit <- "log10_nGenes"
	ytit <- "log10_nUMI"
    }
    xlim = range(xcoor)
    ylim = range(ycoor)

    if (dacolors$hasdist){
	layout(matrix(c(1,1,1,2),nrow=1,ncol=4,byrow=T))
	par(mar=c(1,1,1,0) + 0.1)
	p1 <- plot(x=xcoor[isna],y=ycoor[isna],xlim=xlim,ylim=ylim,pch=21
			,xlab=' ', ylab=' ', cex.axis= 1.0, xaxt='n', yaxt='n', cex=1, bg=dacolors$col[isna], col=dacolors$col[isna], main=title)
    isna <- !isna

  par(new=TRUE)
	plot(x=xcoor[isna],y=ycoor[isna],xlim=xlim,ylim=ylim
			,xlab=' ', ylab=' ', cex.axis= 1.0, xaxt='n', yaxt='n', cex=pointscale, col=dacolors$col[isna], main=title)
		par(mar=c(4,2,3,0) + 0.1)

      title(xlab=xtit, line=0.75, cex.lab=1.5)
      title(ylab=ytit, line=0.75, cex.lab=1.5)
	p2 <- plot(x=1:(dacolors$non.na), y=dacolors$valrange, col=dacolors$valcolrs,pch=21,xaxt='n',cex.axis= 1.0, yaxp =c(), bty="l",cex.lab= 1.5, xlab="", main=ylab)

	title(xlab=paste(as.character(dacolors$non.na), "\ncells"), line=3, cex.lab=1.5) # , family="Calibri Light"
    }else{
	library(ggplot2)
	plotdata <- data.frame(cbind(xcoor,ycoor,as.character(dacolors$col)))
	colnames(plotdata) = c(xtit, ytit, query)
	rownames(plotdata) =as.character(1:length(xcoor))
	plotdata[[xtit]] <- as.numeric(as.character(plotdata[[xtit]]))
	plotdata[[ytit]] <- as.numeric(as.character(plotdata[[ytit]]))
	p <- ggplot(plotdata, aes_string(x=xtit,y=ytit,color=query)) +scale_color_manual(values=dacolors$valcolrs) 
        print(p + geom_point()) 
    }
    if (!is.null(save.pdf)){
        dev.print(pdf, paste("/lustre/scratch117/cellgen/team218/lh20/tmppdf/", save.pdf, ".pdf", sep=""))
    }
    if (!is.null(save.png)){
       	dev.copy(png, paste("/lustre/scratch117/cellgen/team218/lh20/tmppdf/", save.png, ".png", sep=""), width = 480, height = 360)
        #else dev.copy(png, paste("/lustre/scratch117/cellgen/team218/lh20/", save.png, ".png", sep=""), width = 480, height = 360)
    #    dev.copy(png, paste("pdfbox/", save.png, ".png", sep=""), width = 480, height = 360)
        dev.off();
    }
	}else if (class(sro) == "SingleCellExperiment"){
		

	}

}
plotSeuratFeatBatchExport <- function(sro, stringtoparse,width = 4, height=3){
	dev.new(width= width,height=height)
	for( i in parselist(stringtoparse) ) plotSeuratFeat(sro, i, save.pdf=i)
	dev.off
}



makeSce_and_mnn  <-function(pathlist, alias=c()){
 library(Seurat)
 library(SingleCellExperiment)
 dalist <- list()
 davar <- list()
 batch <- c()
 aliases <- c()
# cellnames <- c()
 if (is.null(aliases)) aliases <- sub(".*/","", pathlist)

 for( i in 1:length(pathlist)){
   print(paste("processing", alias[i]))
   newite <- SingleCellExperiment(list(counts= Read10X(pathlist[i])))
   colnames(newite) <- paste(aliases[i],colnames(newite),sep="_")
   logcounts(newite) <- counts(newite)
   logcounts(newite)@x <- log2(logcounts(newite)@x + 1)

   dec <- tryCatch(rownames(decomposeVar(newite, trendVar(newite, parametric=T, use.spikes = FALSE))),error=function(cond){return(c())})
   dalist <- c(dalist, newite)    
   davar <- c(davar, dec)
   batch <- c(batch, rep(aliases[i], ncol(counts(newite))))
#    cellnames <- c(cellnames, coddlnames(newite))
 }
print("done") 
#detach(Seurat)
 library(scran)

 davar <- sort(table(as.character(davar)))
 if (length(davar) > 2000) davar <- names(davar)[1:2000]
 else davar <- names(davar)
 mnninput <- list();
 data <- list()
 for( i in 1:length(pathlist)){
   mnninput <- c(mnninput, list(as.matrix(logcounts(dalist[[i]]))[davar,] ))
   data <- c(data, (logcounts(dalist[[i]])))
 }
 print("runningMMN")
 mnnout <- do.call(mnnCorrect, c(mnninput, list(k=20)))
 rownames(omat) <- make.names(rownames(omat),unique=T)
 sce <- SingleCellExperiment(list(logcounts=do.call(cbind,data)))
 tmp <- t(do.call(cbind,mnnout$corrected))
 rownames(tmp) <- colnames(sce)
 reducedDim(sce,"MNN") <- as.matrix(tmp)
 names(batch) <- colnames(sce)
 sce$Batch <- batch
  return(sce)
}



#	set.seed(100)
	# Using irlba to set up the t-SNE, for speed.
#	osce <- runPCA(sce, ntop=Inf, method="irlba")
#	osce <- runTSNE(osce, use_dimred="PCA")
#	ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")set.seed(100)
#	csce <- runTSNE(sce, use_dimred="MNN")
#	ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
	
#	return(sce)
#}

specify_decimal <- function(x,k) trimws(format(round(x,k),nsmall=k))

makeHeatmap <- function(table, rowannot =c(), colannot=c(), plot.attribs=c()){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	library(ggplot2)
	library(stats)
	if (!is.null(colannot)){
		if (is.null(colnames(table))) cdata <- colannot
		else {
			cdata <- colannot[colnames(table),,drop=FALSE]
		}
		if (is.null(colnames(colannot))) colnames(colannot) <- paste("Annotation", 1:ncol(colannot))
	}

	if (!is.null(rowannot)){
		if (is.null(rownames(table))) rdata <- rowannot
		else {
			rdata <- rowannot[rownames(table),,drop=FALSE]
		}
		if (is.null(colnames(rowannot))) colnames(colannot) <- paste("Annotation", 1:ncol(rowannot))
	}

	dtable <- table
	dtable[is.infinite(dtable)] <- NA
	dam <- mean(dtable,na.rm=T)
	if (is.na(dam)) dtable <- matrix(0, nrow(dtable),ncol(dtable))
	else dtable[is.na(dtable)] <- mean(dtable,na.rm=T)
	

	if (nrow(table) > 1){
		row.order <- hclust(dist(dtable), method="complete")$order
		table <- table[row.order,,drop=FALSE]
		if (!is.null(rowannot)) rdata <- rdata[row.order,,drop=FALSE]
	}
	if (ncol(table) > 1){
		col.order <- hclust(dist(t(dtable)), method="complete")$order
		table <- table[,col.order,drop=FALSE]
		if (!is.null(colannot)) cdata <- cdata[col.order,, drop=FALSE]
		
	}
	tt = structure(table, .Dim = c(nrow(table), ncol(table)), .Dimnames = structure(list(Row = rownames(table), Col = colnames(table)), .Names = c("Row", "Col")), class = "table")
	tt <- as.data.frame(tt)
	tt$TiletypeC <- "data"
	tt$TiletypeR <- "data"

	if (!is.null(colannot)) {
		for(i in 1:ncol(cdata)) {
			normin <- cdata[,i]
			isinf  <- is.infinite(normin)
			isnna <- !is.na(cdata[,i])
			normin[isinf] <- NA
			dam <- mean(normin, na.rm=T)
			dav <- sqrt(var(normin,na.rm=T))
			if ((!is.na(dav))&&(dav > 0)) {
				cdata[isnna,i] <- pnorm(cdata[isnna,i], mean =dam, sd = dav)
				colnames(cdata)[i] <- paste(colnames(cdata)[i], "[",specify_decimal(dam,3),"+- ",specify_decimal(dav,3)," stdevs]") 
				
			}else{
				print(paste("failed to normalize", colnames(cdata)[i] ))
				print(dav)
				print(cdata[,i])
				cdata[,i] <- 0
			}
		}
		ctt = as.data.frame(structure(t(cdata), .Dim = c(ncol(cdata), nrow(cdata)), .Dimnames = structure(list(Row = colnames(cdata), Col = rownames(cdata)), .Names = c("Row", "Col")), class = "table"))
		ctt$TiletypeC <- "Annotation"
		ctt$TiletypeR <- "data"

		tt <- rbind(tt,ctt)
	}

	if (!is.null(rowannot)){ 
		for(i in 1:ncol(rdata)) {
			normin <- rdata[,i]
			isinf  <- is.infinite(normin)
			isnna <- !is.na(rdata[,i])
			normin[isinf] <- NA
			dam <- mean(normin, na.rm=T)
			dav <- sqrt(var(normin,na.rm=T))
			if ((!is.na(dav))&&(dav > 0)) {
				rdata[isnna,i] <- pnorm(rdata[isnna,i], mean =dam, sd = dav)
				colnames(rdata)[i] <- paste(colnames(rdata)[i], "[",specify_decimal(dam,3),"+- ",specify_decimal(dav,3)," stdevs]") 
			}else{
				print(paste("failed to normalize", colnames(rdata)[i] ))
				print(dav)
				print(rdata[,i])
				rdata[,i] <- 0
			}
		}
		ctt = as.data.frame(structure(t(rdata), .Dim = c(ncol(rdata), nrow(rdata)), .Dimnames = structure(list(Row = colnames(rdata), Col = rownames(rdata)), .Names = c("Col", "Row")), class = "table"))
		ctt$TiletypeC <- "data"
		ctt$TiletypeR <- "Annotation"
		tt <- rbind(tt,ctt)
	}

	p <- ggplot(as.data.frame(tt),aes(x=Col,y=Row,fill=Freq))
	if (!(("rot.xnames" %in% flags.plot)||("heat.rot.xnames" %in% flags.plot))) p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
	p <- p + geom_tile()
	p <- p + scale_fill_gradientn(colours= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[1:41])

	if (!is.null(colannot)){
		if (!is.null(rowannot))	p <- p + facet_grid(rows = vars(TiletypeC),cols = vars(TiletypeR), scales = "free", space = "free")
		else p <- p + facet_grid(rows = vars(TiletypeC), scales = "free", space = "free")
	}else if (!is.null(rowannot)) p <- p + facet_grid(cols = vars(TiletypeR), scales = "free", space = "free")

return(changeStyle(p, plot.attribs, "heat"))}
# "no.xlabel", "no.jitter", "no.xticks" 
makeSCurvePlot <- function(data, label, metalabel =c(), plot.attribs=c()){
	library(ggplot2)
	llvl <- getSafeLevels(label)
	if (!is.null(metalabel)) mllvl <- getSafeLevels(metalabel)
	
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if ("use.colors" %in% names(plot.attribs)) use.colors <- plot.attribs[["use.colors"]]
	if (class(label) != "factor") label <- factor(label, levels=unique(label))

	library(ggplot2)
	if (length(data) != length(label)) stop("invalid length")


	if ("x.filter.class" %in% names(plot.attribs)){
		filter <- !is.na(match(label, plot.attribs$x.filter.class))
		label <- label[filter]
		data <- data[filter]
	}
	if ("x.remap.class" %in% names(plot.attribs)){
		remap <- match(levels(label), names(plot.attribs$x.remap.class))
		filter <- !is.na(remap)
		levels(label)[filter] <- as.character(plot.attribs$x.remap.class[remap[filter]])
	}

	daorder <- order(label,data)
	tt <- data.frame(matrix(0,length(label),2))
	colnames(tt) <- c("Label", "value")
	tt[,"value"] <- data[daorder];
	tt[,"Label"] <- label[daorder];
	if ("y.log.transform" %in% flags.plot) {
		tt[[2]] <- log2(tt[[2]]+1)
	}	
	if (dim(tt)[2] == 2) {
		p <- ggplot(tt, aes(Label, value, fill=factor(Label)))
		if (("no.xlabel" %in% flags.plot)||("violin.no.xlabel" %in% flags.plot))  p <- p + xlab(NULL)
		if (("no.ylabel" %in% flags.plot)||("violin.no.ylabel" %in% flags.plot))  p <- p + ylab(NULL)

		if ("no.xnamedticks" %in% flags.plot) p <-p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
		else{
			if ("no.xticks" %in% flags.plot) p <-p + theme(axis.ticks.x=element_blank())
			if ("no.xnames" %in% flags.plot) p <-p + theme(axis.text.x=element_blank())
			else if (!(("rot.xnames" %in% flags.plot)||("violin.rot.xnames" %in% flags.plot))) p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
		}
		if ("no.ynamedticks" %in% flags.plot) p <-p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
		else{
			if ("no.yticks" %in% flags.plot) p <-p + theme(axis.ticks.y=element_blank())
			if ("no.ynames" %in% flags.plot) p <-p + theme(axis.text.y=element_blank())
			else if ((("rot.ynames" %in% flags.plot)||("violin.rot.ynames" %in% flags.plot))) p <- p + theme(axis.text.y=element_text(angle = 90, hjust = 1))
		}
		if ("title" %in% names(plot.attribs)) p <- p + ggtitle(plot.attribs[["title"]])
		if (!is.null(use.colors)) p <- p + scale_fill_manual(values=use.colors)
	}else {
		colnames(tt) <- c("Label", "Value", "Type")
		p<-ggplot(tt, aes(factor(label), data,fill=Type));
	}
	
	if ("do.boxplot.instead" %in% flags.plot ) p <- p +  geom_boxplot(color = "0x000000") 
	else p <- p +  geom_violin() 
	if (!"no.jitter" %in% flags.plot) p <- p + geom_jitter(height = 0, width = 0.1)

return(p)}


makeViolinPlot <- function(label,data, metalabel= c(), use.colors=c(), plot.attribs=c(),tmptmp=c()){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if ("use.colors" %in% names(plot.attribs)) use.colors <- plot.attribs[["use.colors"]]
	if (class(label) != "factor") label <- factor(label, levels=unique(label))

	library(ggplot2)
	if (length(data) != length(label)) stop("invalid length")


	if ("x.filter.class" %in% names(plot.attribs)){
		filter <- !is.na(match(label, plot.attribs$x.filter.class))
		label <- label[filter]
		data <- data[filter]
	}
	if ("x.remap.class" %in% names(plot.attribs)){
		remap <- match(levels(label), names(plot.attribs$x.remap.class))
		filter <- !is.na(remap)
		levels(label)[filter] <- as.character(plot.attribs$x.remap.class[remap[filter]])
	}

	daorder <- order(label,data)
	tt <- data.frame(matrix(0,length(label),3))
	colnames(tt) <- c("Label", "value", "rank")
	tt[,"value"] <- data[daorder];
	tt[,"Label"] <- label[daorder];
	for(i in levels(label)){
		tt[(tt[,"Label"] == i),"rank"] <- -0.5 + (((order(tt[(tt[,"Label"] == i),"value"] ))*4) / sum((tt[,"Label"] == i))) 

	}
	tt[tt[,"rank"] > 0.5,"rank"] <- 0.75 - 0.5 *tt[tt[,"rank"] > 0.5,"rank"]
	tt[tt[,"rank"] < -0.5,"rank"]  <- -1.5 - tt[tt[,"rank"] < -0.5,"rank"]*2
	tt[,"rank"] <- tt[,"rank"] * 0.8

	if ("y.log.transform" %in% flags.plot) {
		tt[[2]] <- log2(tt[[2]]+1)
	}	

	if (dim(tt)[2] == 3) {
		p <- ggplot(tt, aes(Label, value, fill=factor(Label), color=factor(Label)))
		if (("no.xlabel" %in% flags.plot)||("violin.no.xlabel" %in% flags.plot))  p <- p + xlab(NULL)
		if (("no.ylabel" %in% flags.plot)||("violin.no.ylabel" %in% flags.plot))  p <- p + ylab(NULL)

		if ("no.xnamedticks" %in% flags.plot) p <-p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
		else{
			if ("no.xticks" %in% flags.plot) p <-p + theme(axis.ticks.x=element_blank())
			if ("no.xnames" %in% flags.plot) p <-p + theme(axis.text.x=element_blank())
			else if (!(("rot.xnames" %in% flags.plot)||("violin.rot.xnames" %in% flags.plot))) p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
		}
		if ("no.ynamedticks" %in% flags.plot) p <-p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
		else{
			if ("no.yticks" %in% flags.plot) p <-p + theme(axis.ticks.y=element_blank())
			if ("no.ynames" %in% flags.plot) p <-p + theme(axis.text.y=element_blank())
			else if ((("rot.ynames" %in% flags.plot)||("violin.rot.ynames" %in% flags.plot))) p <- p + theme(axis.text.y=element_text(angle = 90, hjust = 1))
		}
		if ("title" %in% names(plot.attribs)) p <- p + ggtitle(plot.attribs[["title"]])
		if (!is.null(use.colors)) p <- p + scale_fill_manual(values=use.colors)
	}else {
		colnames(tt) <- c("Label", "Value", "Type")
		p<-ggplot(tt, aes(factor(label), data,fill=Type));
	}
	
	if ("do.boxplot.instead" %in% flags.plot ) p <- p +  geom_boxplot(color = "#000000") 
	else p <- p +  geom_violin(color = "#000000") 
	if (!"no.jitter" %in% flags.plot) p <- p + geom_point(aes(x=as.numeric(Label) + rank, y=value, color=factor(Label)), size=0.25)


return(changeStyle(p,plot.attribs))}

makeScatterPlot <- function(xcoor, ycoor, colorclass, plot.attribs=c(), use.color=c()){
	
	tt <- data.frame(matrix(0,length(label),3))
	colnames(tt) <- c("X","Y", "Col")
	tt$X <- xcoor;
	tt$Y <- ycoor;
	tt$Col <- colorclass
	p <- ggplot(tt, aes(X, Y, fill=factor(Col)))
	p <- p + geom_point()
	p <- p + scale_color_manual(name=query,values=use.color,label=col.lvl,drop = FALSE)

return(changeStyle(p,plot.attribs))}

plotSeuratScatter <- function(sro, xquery, yquery, cquery){
}

changeStyle <- function(p, plot.attribs, classprefix=""){
	library(ggplot2)
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if ("title" %in% names(plot.attribs)) p <- p + ggtitle(plot.attribs[["title"]])
	if ("xlabel" %in% names(plot.attribs)) p <- p + xlab(plot.attribs[["xlabel"]])
	else if (("no.xlabel" %in% flags.plot)||(paste(classprefix,"no.xlabel",sep=".") %in% flags.plot))  p <- p + xlab(NULL)
	if ("ylabel" %in% names(plot.attribs)) p <- p + ylab(plot.attribs[["ylabel"]])
	else if (("no.ylabel" %in% flags.plot)||(paste(classprefix,"no.ylabel",sep=".") %in% flags.plot))  p <- p + ylab(NULL)
	
	if ("ylabel" %in% names(plot.attribs)) p <- p + ylab(plot.attribs[["ylabel"]])
	else if (("no.ylabel" %in% flags.plot)||(paste(classprefix,"no.ylabel",sep=".") %in% flags.plot))  p <- p + ylab(NULL)

	if ("no.legend" %in% flags.plot) p <- p + theme(legend.position="none")
	else{
		if ("legnames.size" %in% names(plot.attribs))  p <- p + theme(legend.text=element_text(size=plot.attribs[["legnames.size"]]))
	}

	# building xticks style
	themearg <- list()
	if (("no.xnamedticks" %in% flags.plot)||(paste(classprefix,"no.xnamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.x=element_blank()) 
		p <- p + theme(axis.text.x=element_blank())
	}else{
		if ("no.xticks" %in% flags.plot) p <- p + theme(axis.ticks.x=element_blank())
		if ("no.xnames" %in% flags.plot) p <- p + theme(axis.text.x=element_blank())
		else if (!(("rot.xnames" %in% flags.plot)||(paste(classprefix,"rot.xnames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("xnames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["xnames.size"]]))
		if ("xnames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["xnames.color"]]))
		if ("xnames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.x= do.call(element_text, themearg))
	themearg <- list()
	if (("no.ynamedticks" %in% flags.plot)||(paste(classprefix,"no.ynamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.y=element_blank()) 
		p <- p + theme(axis.text.y=element_blank())
	}else{
		if ("no.yticks" %in% flags.plot) p <- p + theme(axis.ticks.y=element_blank())
		if ("no.ynames" %in% flags.plot) p <- p + theme(axis.text.y=element_blank())
		else if ((("rot.ynames" %in% flags.plot)||(paste(classprefix,"rot.ynames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("ynames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["ynames.size"]]))
		if ("ynames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["ynames.color"]]))
		if ("ynames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.y= do.call(element_text, themearg))

	if ("scale.xrange" %in% names(plot.attribs)) p <- p + scale_x_continuous(limits = plot.attribs$scale.xrange)
	if ("scale.yrange" %in% names(plot.attribs)) p <- p + scale_y_continuous(limits = plot.attribs$scale.yrange)	

	if ("nb.col.legend" %in% names(plot.attribs)) p <- p + guides(color=guide_legend(ncol = plot.attribs$nb.col.legend),fill=guide_legend(ncol = plot.attribs$nb.col.legend))

	return(p)
}
# first argument is either a ggplot or a list with "plot.attrib" entries,
# returns a list with "plot.attrib" *unless* a ggplot is provided as first argument, in this latter case a ggplot object is returned
plotAttribs <-function(..., ggplot.or.plot.attrib=c(), baseAttribs=c(),no.legend=F,no.jitter=F, no.xticks=F, no.xlabel=F, no.yticks=F, no.ylabel=F, title=c(), xtitle= c(), ytitle=c(),xlabel=c(), ylabel=c(), xnames.bold=F,xnames.color=c(), xnames.size=c(), ynames.bold=F,ynames.color=c(), ynames.size=c(),legnames.size=c(),no.xnamedticks=F,no.ynamedticks=F,nb.col.legend=c()){
	if ((is.null(ggplot.or.plot.attrib))||(inherits(ggplot.or.plot.attrib, "ggplot"))) {
		if (is.null(baseAttribs)) fout <- list(flags=c())
		else fout <- baseAttribs
	}else fout <- ggplot.or.plot.attrib

	fout[["flags"]]  <- c(fout[["flags"]], ...)
	if (no.jitter) fout[["flags"]] <- c(fout[["flags"]], "no.jitter")
	if (no.legend) fout[["flags"]] <- c(fout[["flags"]], "no.legend") 
	if (no.xlabel) fout[["flags"]] <- c(fout[["flags"]], "no.xlabel")
	if (no.xticks) fout[["flags"]] <- c(fout[["flags"]], "no.xticks")
	if (no.xnamedticks) fout[["flags"]] <- c(fout[["flags"]], "no.xnamedticks")
	if (no.ylabel) fout[["flags"]] <- c(fout[["flags"]], "no.ylabel")
	if (no.yticks) fout[["flags"]] <- c(fout[["flags"]], "no.yticks")
	if (no.ynamedticks) fout[["flags"]] <- c(fout[["flags"]], "no.ynamedticks")

	if (xnames.bold) fout[["flags"]] <- c(fout[["flags"]], "xnames.bold")
	if (ynames.bold) fout[["flags"]] <- c(fout[["flags"]], "ynames.bold")


	if (!is.null(xnames.color)) fout[["xnames.color"]] = xnames.color
	if (!is.null(xnames.size)) fout[["xnames.size"]] = xnames.size
	if (!is.null(ynames.color)) fout[["ynames.color"]] = ynames.color
	if (!is.null(ynames.size)) fout[["ynames.size"]] = ynames.size
	if (!is.null(legnames.size)) fout[["legnames.size"]] = legnames.size

	if (!is.null(title)) fout[["title"]] = title
	if (!is.null(xlabel)) fout[["xlabel"]] = xlabel
	if (!is.null(ylabel)) fout[["ylabel"]] = ylabel
	if (!is.null(nb.col.legend)) fout[["nb.col.legend"]] = nb.col.legend


	if (inherits(ggplot.or.plot.attrib, "ggplot")) return(changeStyle(ggplot.or.plot.attrib, fout))
return(fout)}

plotSeuratAUROC <- function(sro, genequery, cellquery, metacl, use.cell=c(), use.meta=c(), pos.cell=c(), pos.meta=c()){
	use.cell <- SeuratCellListQuery(sro, use.meta, use.cell)
	pos.cell <- SeuratCellListQuery(sro, use.meta, use.cell)
	neg.cell <- use.cell & (!pos.cell)
	pos.cell <- use.cell & pos.cell

	
}

plotSeuratBox <-function(sro, query , metacl, count.type="L",use.cell=c(),use.meta=c(), plot.attribs=c(), do.logtransform=F, try.cummul.instead=F,no.jitter=F){
	if (is.null(use.cell)&(is.null(use.meta)))  use.cell = rownames(sro@meta.data)
	else {
		print(sum(SeuratCellListQuery(sro,use.meta,use.cell)))
		use.cell = SeuratCellListQuery(sro,use.meta,use.cell,getnames=T) 
	}
	if (("meta.color" %in% names(sro@misc))&&(metacl %in% names(sro@misc$meta.color))) use.color <- sro@misc$meta.color[[metacl]]
	else use.color <- c()
	if (do.logtransform) plot.attribs$flags <- c(plot.attribs$flags , "y.log.transform")
	if (no.jitter) plot.attribs$flags <- c(plot.attribs$flags , "no.jitter")
	plot.attribs$flags <- c(plot.attribs$flags , "do.boxplot.instead")


	if( length(query) == 1){
		if (query %in% names(sro@meta.data)){
 			return(makeViolinPlot(sro@meta.data[use.cell,metacl], sro@meta.data[use.cell,query],plot.attribs=plot.attribs,use.colors=use.color))
		}else{
			if (query %in% rownames(sro@raw.data)) indind <- match(query, rownames(sro@raw.data))
			else {
			 	indind <- grep(query,rownames(sro@raw.data))
				if (length(indind) == 0) stop("no such gene or meta data")
				if (length(indind) > 1) {
					indind <- indind[1]
					print("more than on gene match query!")
				}
			}
		if (count.type == "L") {
			if (try.cummul.instead) return(plotCummulativeBars(log2(sro@raw.data[indind,use.cell]+1),sro@meta.data[use.cell,metacl],plot.attribs=plot.attribs,use.colors=use.color))
			else return(makeViolinPlot(sro@meta.data[use.cell,metacl], log2(sro@raw.data[indind,use.cell]+1),plot.attribs=plot.attribs,use.colors=use.color))
		}else if (count.type == "N"){
			return(makeViolinPlot(sro@meta.data[use.cell,metacl], sro@data[indind,use.cell],plot.attribs=plot.attribs,use.colors=use.color))
		}else return(makeViolinPlot(sro@meta.data[use.cell,metacl], (log2(sro@raw.data[indind,use.cell]+1) -log2(sro@meta.data[use.cell,"nUMI"]) + 20),plot.attribs=plot.attribs,use.colors=use.color))
		
		}
	}else{
		for(i in 1:length(query)){
			if (query[i] %in% names(sro@meta.data)){
				if (i == 1) data <- cbind(sro@meta.data[use.cell,metacl],sro@meta.data[use.cell,query[i]], rep(i,length(use.cell) ))
				else data <- rbind(data,cbind(sro@meta.data[use.cell,metacl],sro@meta.data[use.cell,query[i]], rep(i,length(use.cell) )))
			}
		}
		tt <- as.data.frame(data)
		colnames(tt) <- c("Label", "Value", "Type")
		print(tt)
		ggplot(tt, aes(x=factor(Label), y=Value,fill=factor(Type))) + geom_violin(position=position_dodge(1)) 
	}
}


plotSeuratViolin <-function(sro, query , metacl, count.type="L",use.cell=c(),use.meta=c(), plot.attribs=c(), do.logtransform=F, try.cummul.instead=F,no.jitter=T,do.box.plot.instead=F){
	if (is.null(use.cell)&(is.null(use.meta)))  use.cell = rownames(sro@meta.data)
	else {
		print(sum(SeuratCellListQuery(sro,use.meta,use.cell)))
		use.cell = SeuratCellListQuery(sro,use.meta,use.cell,getnames=T) 
	}
	if (("meta.color" %in% names(sro@misc))&&(metacl %in% names(sro@misc$meta.color))) use.color <- sro@misc$meta.color[[metacl]]
	else use.color <- c()
	if (do.logtransform) plot.attribs$flags <- c(plot.attribs$flags , "y.log.transform")
	if (no.jitter) plot.attribs$flags <- c(plot.attribs$flags , "no.jitter")


	if( length(query) == 1){
		if (query %in% names(sro@meta.data)){
 			return(makeViolinPlot(sro@meta.data[use.cell,metacl], sro@meta.data[use.cell,query],plot.attribs=plot.attribs,use.colors=use.color))
		}else{
			if (query %in% SeuratRawData(sro, "R")) indind <- match(query, SeuratRawData(sro,"R"))
			else {
			 	indind <- grep(query,SeuratRawData(sro,"R"))
				if (length(indind) == 0) stop("no such gene or meta data")
				if (length(indind) > 1) {
					indind <- indind[1]
					print("more than on gene match query!")
				}
			}
		if (count.type == "L") {
			if (try.cummul.instead) return(plotCummulativeBars(log2(SeuratRawData(sro)[indind,use.cell]+1),sro@meta.data[use.cell,metacl],plot.attribs=plot.attribs,use.colors=use.color))
			else if (do.box.plot.instead) return(makeSCurvePlot(log2(SeuratRawData(sro)[indind,use.cell]+1),sro@meta.data[use.cell,metacl], plot.attribs=plot.attribs,use.colors=use.color))
			else return(makeViolinPlot(sro@meta.data[use.cell,metacl], log2(SeuratRawData(sro)[indind,use.cell]+1),plot.attribs=plot.attribs,use.colors=use.color))
		}else if (count.type == "N"){
			return(makeViolinPlot(sro@meta.data[use.cell,metacl], SeuratRawData(sro,"N")[indind,use.cell],plot.attribs=plot.attribs,use.colors=use.color))
		}else if (do.box.plot.instead) return(makeSCurvePlot((log2(SeuratRawData(sro)[indind,use.cell]+1) -log2(sro@meta.data[use.cell,"nUMI"]) + 20), sro@meta.data[use.cell,metacl],plot.attribs=plot.attribs,use.colors=use.color))
		else return(makeViolinPlot(sro@meta.data[use.cell,metacl], (log2(SeuratRawData(sro)[indind,use.cell]+1) -log2(sro@meta.data[use.cell,"nUMI"]) + 20),plot.attribs=plot.attribs,use.colors=use.color))
		
		}
	}else{
		for(i in 1:length(query)){
			if (query[i] %in% names(sro@meta.data)){
				if (i == 1) data <- cbind(sro@meta.data[use.cell,metacl],sro@meta.data[use.cell,query[i]], rep(i,length(use.cell) ))
				else data <- rbind(data,cbind(sro@meta.data[use.cell,metacl],sro@meta.data[use.cell,query[i]], rep(i,length(use.cell) )))
			}
		}
		tt <- as.data.frame(data)
		colnames(tt) <- c("Label", "Value", "Type")
		print(tt)
		ggplot(tt, aes(x=factor(Label), y=Value,fill=factor(Type))) + geom_violin(position=position_dodge(1)) 
	}
}
computeFreq <- function(sro,xannot,yannot){
	if (length(xannot) == 1) xannot<- sro@meta.data[[xannot]];
	if (length(yannot) == 1) yannot<- sro@meta.data[[yannot]];
	
	if (class(xannot) != "factor") xannot <- factor(xannot, levels = unique(xannot))
	if (class(yannot) != "factor") yannot <- factor(yannot, levels = unique(yannot))

	counts <- matrix(0,length(levels(xannot)), length(levels(yannot)))
	for(i in 1:length(xannot)) counts[ xannot@.Data[i], yannot@.Data[i] ] = counts[xannot@.Data[i], yannot@.Data[i]] +1
	rownames(counts) <- levels(xannot)
	colnames(counts) <- levels(yannot)
return(counts)}

plotCummulativeBars <- function(xannot,hannot, wannot=c(), do.horizscale=T, use.colors=c(), show.label.lowerbound.fraction= 0.05, do.sort.first = F, do.sort.second=F,do.resize=c(), is.vertical=F, plot.attribs=c(),return.counts.instead=F){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if (class(hannot) == "factor") {
		hlabs = levels(hannot)
		if ("x.flt.empty" %in% flags.plot){
			tmp <- table(hannot)
			hlabs <- hlabs[tmp[hlabs] > 0] # remove 0
		}
	}else hlabs = unique(hannot)

	brks <- matrix(0,length(hlabs), 11)
	if ("integer.fix" %in% flags.plot){
		for(i in 1:length(hlabs)){
			curvals <- sort(xannot[hannot == hlabs[i]],na.last=NA)
			print(hlabs[i])
			print(mean(curvals))
			summary <- table(curvals)
			brks[i,2] <- min(curvals);
			tmpint = (length(curvals)-1)
			tmpcoor = 1+(tmpint/16)
			brks[i,3] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(tmpint/8)
			brks[i,4] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(tmpint/4)
			brks[i,5] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(tmpint/2)
			brks[i,6] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(3*tmpint/4)
			brks[i,7] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(7*tmpint/8)
			brks[i,8] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			tmpcoor = 1+(15*tmpint/16)
			brks[i,9] <- curvals[tmpcoor] + (tmpcoor - match(curvals[tmpcoor], curvals))/summary[match(curvals[tmpcoor], names(summary))]
			brks[i,10] <- max(curvals);
		}
	}else{
		for(i in 1:length(hlabs)){
			curvals <- sort(xannot[hannot == hlabs[i]],na.last=NA)
			brks[i,2] <- min(curvals);
			tmpint = (length(curvals)-1)
			brks[i,3] <- curvals[1+(tmpint/8)]
			brks[i,4] <- curvals[1+(tmpint/4)]
			brks[i,5] <- curvals[1+(3*tmpint/8)]
			brks[i,6] <- curvals[1+(tmpint/2)]
			brks[i,7] <- curvals[1+(5*tmpint/8)]
			brks[i,8] <- curvals[1+(3*tmpint/4)]
			brks[i,9] <- curvals[1+(7*tmpint/8)]
			brks[i,10] <- max(curvals);
		}
	}
	brks[1,1] = min(brks[,2])
	brks[1,11] = max(brks[,10])
	if (length(hlabs) > 1){
		for(i in 2:length(hlabs)){
			brks[i,1] <- brks[1,1]
			brks[i,11] <- brks[1,11];
		}
	}
	data <- data.frame(row.names= 1:(length(hlabs) * 10))
	augcol = c("#FFFFFF")
	augnames = c("BG")
	for(i in 1:length(hlabs)){
		for(j in 1:10){
			data[i*10+j-10, "IX"] = (i -1)  / length(hlabs);
			data[i*10+j-10, "AX"] = i  / length(hlabs);
			data[i*10+j-10, "IY"] = brks[i,j];
			data[i*10+j-10, "AY"] = brks[i,j +1];
		}
		data[i*10-9, "COL"] = "BG"
		data[i*10, "COL"] = "BG"
		data[i*10-8, "COL"] = paste(hlabs[i],"1",sep='_')
		data[i*10-1, "COL"] = data[i*10-8, "COL"]
		data[i*10-7, "COL"] = paste(hlabs[i],"2",sep='_')
		data[i*10-2, "COL"] = data[i*10-7, "COL"]
		data[i*10-6, "COL"] = paste(hlabs[i],"3",sep='_')
		data[i*10-3, "COL"] = data[i*10-6, "COL"]
		data[i*10-5, "COL"] = paste(hlabs[i],"4",sep='_')
		data[i*10-4, "COL"] = data[i*10-5, "COL"]

		if (hlabs[i] %in% names(use.colors)){
			augcol = c(augcol, use.color[hlabs[i]], "#000000", "#888888", use.color[hlabs[i]])
		}else{
			augcol = c(augcol, "#CCCCCC", "#888888", "#444444", "#000000")
		}
		augnames = c(augnames, paste(hlabs[i],"1",sep='_'), paste(hlabs[i],"2",sep='_'), paste(hlabs[i],"3",sep='_'), paste(hlabs[i],"4",sep='_'))

	}
	names(augcol) <- augnames

	p <- ggplot(data)
	p <- p + geom_rect(data=data,mapping=aes(xmin=IX, xmax=AX, ymin=IY, ymax=AY, fill=COL), color="black", alpha=0.5)
	p <- p + scale_fill_manual(values=augcol)

	return(p)
}

plotDualBar <- function(xannot,hannot, wannot=c(), do.horizscale=T, use.colors=c(), show.label.lowerbound.fraction= 0.05, do.sort.first = F, do.sort.second=F,do.resize=c(), is.vertical=F, plot.attribs=c(),return.counts.instead=F,make.piechart=F,bypass.counts=c()){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()

	if (!is.null(do.resize)) setPlotSize(do.resize[1],do.resize[2])
	if (length(xannot) != length(hannot)) {
		print(table(xannot))
		print(table(hannot))
		 stop("annotation length mistaches")
	}
	if (class(xannot) == "factor") {
		vlabs = levels(xannot)
		if ("y.flt.empty" %in% flags.plot){
			tmp <- table(xannot)
			print(tmp)
			print(tmp[vlabs])
			vlabs <- vlabs[tmp[vlabs] > 0] # remove 0
		}
	}else vlabs = unique(xannot)
	if (is.vertical == FALSE) vlabs <- rev(vlabs)

	if (class(hannot) == "factor") {
		hlabs = levels(hannot)
		if ("x.flt.empty" %in% flags.plot){
			tmp <- table(hannot)
			hlabs <- hlabs[tmp[hlabs] > 0] # remove 0
		}
	}else hlabs = unique(hannot)



#	if ("x.filter.class" %in% names(plot.attribs)){
#		filter <- !is.na(match(xannot, plot.attribs$x.filter.class))
#		label <- label[filter]
#		data <- data[filter]
#		filter <- !is.na(match(xlabs, plot.attribs$x.filter.class))
#		xlabs <- xlabs[filter]
#	}
	
	if ("x.remap.class" %in% names(plot.attribs)){
	#	print(table(xannot))
	#	remap <- match(as.character(xannot), names(plot.attribs$x.remap.class))
	#	print(table(remap))
	#	filter <- !is.na(remap)
	#	xannot[filter] <- plot.attribs$x.remap.class[remap[filter]]
	#	print(table(xannot))
	#	print(names(plot.attribs$x.remap.class))

		remap <- match(levels(xannot), names(plot.attribs$x.remap.class))
		filter <- !is.na(remap)
		levels(xannot)[filter] <- as.character(plot.attribs$x.remap.class[remap[filter]])
	}

#	if ("x.filter.empty" %in% flags.plot){
#		tmp <- table(xannot)
#		tmp <- tmp[tmp !=0]
#		remap <- match(xlabs, names(tmp))
#		xlabs <- xlabs[!is.na(remap)]
#		print(tmp)
#		print(remap)	
#	}

	if (do.sort.first) {
		#if (is.na(as.numeric(levels(vlabs)[1]))) vlabs <- sort(vlabs)
		vlabs <- sort(vlabs)
	}
	if (do.sort.second) {
		#if (is.na(as.numeric(levels(nlabs)[1]))) nlabs <- sort(nlabs)
		#else nlabs <- sort(nlabs)
		hlabs <- sort(hlabs)
	}

	if ("x.rev.order" %in% flags.plot) hlabs <- rev(hlabs)
	if ("y.rev.order" %in% flags.plot) vlabs <- rev(vlabs)
	
	if (is.null(use.colors)) use.colors <- myrainbow(length(hlabs))
	counts <- matrix(0,length(hlabs), length(vlabs))


	if (length(wannot) == length(xannot)){
		for(i in 1:length(xannot)) {
			pa <- match(hannot[i], hlabs);
			pb <- match(xannot[i], vlabs);
			counts[pa, pb] = counts[ pa, pb] + wannot[i]
		}
	}else{
		for(i in 1:length(xannot)) {
			pa <- match(hannot[i], hlabs);
			pb <- match(xannot[i], vlabs);
			counts[pa, pb] = counts[ pa, pb] + 1
		}
	}
	if (!is.null(bypass.counts)) {
		counts <- bypass.counts
		hlabs <- rownames(counts)
		vlabs <- colnames(counts)
	}
	if (return.counts.instead){
		rownames(counts) <- hlabs
		colnames(counts) <- vlabs
		return(counts)
	}
	print("alive")
	tt = structure(counts, .Dim = c(length(hlabs), length(vlabs)), .Dimnames = structure(list(Cluster = hlabs, verti = vlabs), .Names = c("Cluster", "verti")), class = "table")
	tt <- as.data.frame(tt)
	tt[["Frequency"]] <- tt[["Freq"]]
	if (do.horizscale) {
		denums <- colSums(counts)
		names(denums) <- vlabs
		for(i in 1:nrow(tt)) tt[i, "Frequency"] <- tt[i, "Frequency"] / denums[tt[i, "verti"]]
	}
	tt[["Label"]] <- as.character(tt[["Freq"]])
	for(i in 1:nrow(tt)) {
		if (is.na(tt[i,"Frequency"])||(tt[i,"Frequency"] <= show.label.lowerbound.fraction)) tt[[i,"Label"]]  <- ""
	}

	library(ggplot2)
	p <- ggplot(as.data.frame(tt),aes(x=factor(verti),y=Frequency,fill=Cluster))
	if ("no.xnamedticks" %in% flags.plot) p <-p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	else{
		if ("no.xticks" %in% flags.plot) p <-p + theme(axis.ticks.x=element_blank())
		if ("no.xnames" %in% flags.plot) p <-p + theme(axis.text.x=element_blank())
		else if (!(("rot.xnames" %in% flags.plot)||("bars.rot.xnames" %in% flags.plot))) p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
	}

	if (is.vertical){
	#if ("x.rev.order" %in% flags.plot)  p <- p + geom_bar(stat="identity",position=position_stack(reverse = TRUE))
	#else
	p <- p + geom_bar(stat="identity",position=position_stack(reverse = FALSE))
	if ("make.pie" %in% flags.plot) p <- p + coord_polar(theta="y")
	else{
	p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = FALSE))
	}
	
	p <- p +  scale_fill_manual(values=use.colors)
	}else{
	#if ("x.rev.order" %in% flags.plot) p <- p + geom_bar(stat="identity",position=position_stack(reverse = FALSE))
	#else
	p <- p + geom_bar(stat="identity",position=position_stack(reverse = TRUE))
	if ("make.pie" %in% flags.plot ) {
		p <- p + coord_polar(theta="y")
	 	p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = TRUE))
		p <- p +  scale_fill_manual(values=use.colors)
	}else{
  p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = TRUE))+
  scale_fill_manual(values=use.colors)+
  coord_flip()
	}
	}
	p <- p + guides(fill=guide_legend(ncol=1))
return(changeStyle(p,plot.attribs,"bar"))}

plotSeuratBars <- function(sro, vannot, hannot, wannot = c(), do.horizscale=T, meta.use=c(), cell.use=c(), use.colors=c(), show.label.lowerbound.fraction= 0.05, do.sort.first = F, do.sort.second=F, do.resize=c(), is.vertical=T,plot.attribs=c(), return.counts.instead=F){
	cell.use <- SeuratCellListQuery(sro,meta.use, cell.use)
	if (("meta.color" %in% names(sro@misc))&&(hannot %in% names(sro@misc$meta.color))) use.colors <- sro@misc$meta.color[[hannot]]
	if (!is.null(wannot)) return(plotDualBar(sro@meta.data[cell.use, c(vannot)], sro@meta.data[cell.use,c(hannot)], wannot= as.numeric(sro@meta.data[cell.use(wannot)]),  do.horizscale=do.horizscale, use.colors=use.colors, show.label.lowerbound.fraction=show.label.lowerbound.fraction,do.sort.first =do.sort.first,do.sort.second =do.sort.second,do.resize=do.resize,is.vertical=is.vertical,plot.attribs=plot.attribs,return.counts.instead=return.counts.instead))
	else return(plotDualBar(sro@meta.data[cell.use, c(vannot)], sro@meta.data[cell.use,c(hannot)], do.horizscale=do.horizscale, use.colors=use.colors, show.label.lowerbound.fraction=show.label.lowerbound.fraction, do.sort.first =do.sort.first, do.sort.second=do.sort.second,do.resize=do.resize,is.vertical=is.vertical,plot.attribs=plot.attribs,return.counts.instead=return.counts.instead))
}
setPlotSize <-function(width, height){
	while (!is.null(dev.list()))  dev.off()
	dev.new(width=width, height=height)
}
setColorSlot <-function(sro, meta,colors){
	clnames <- levels(sro@meta.data[[meta]])
	if (!(meta %in% colnames(sro@meta.data))) stop(paste(meta,"is not a meta.data entry in seurat object"))
	if (length(clnames) != length(colors)) stop(paste("Expected", length(clnames) , "colors but got", length(colors), "instead!"))
	fout = c()
	for(i in 1:length(clnames)) fout[[clnames[i]]] <- colors[i]
	sro@misc$meta.color[[meta]] <- fout
return(sro)}




makegraph <- function(sro, which, cells.use=c(),plot.attribs=c()){
return(switch(which,
"T"=TSNEPlot(sro,group.by="MTG_cluster", colors.use=mydoublerainbow(c(1,5,4,3,4,2,2))),
"C"= plotSeuratBars(sro, "orig.ident", "MTG_cluster",cells.use=cells.use, use.colors=mydoublerainbow(c(1,5,4,3,4,2,2)),do.resize=c(12,3), is.vertical=T),
"c"= plotSeuratBars(sro, "orig.ident", "MTG_cluster",cells.use=cells.use, is.vertical=T,plot.attribs=plot.attribs),
"d"= plotSeuratBars(sro, "orig.ident", "FED_cluster",cells.use=cells.use, is.vertical=T,plot.attribs=plot.attribs),

"S"=plotSeuratFeat(sro, "MTG_simil", do.resize=c(4,3)),
"Q"= plotSeuratViolin(sro,"MTG_simil", "orig.ident",cells.use=cells.use,plot.attribs=plot.attribs),
"P"= plotSeuratViolin(sro,"MTG_simil", "MTG_cluster",cells.use=cells.use,plot.attribs=plot.attribs),
"W"=plotSeuratViolin(sro,"MTG_simil", "orig.ident",cells.use=cells.use,plot.attribs=plot.attribs),
"i"= plotSeuratBars(sro, "orig.ident", "res.0.8",cells.use=cells.use, is.vertical=T,plot.attribs=plot.attribs),
"q"=plotSeuratFeat(sro, "orig.ident", coor=c("MTG_simil", "nUMI"),do.resize=c(4,3),colors.use=myrainbow(length(unique(sro@meta.data$orig.ident))) ),
"r"= plotSeuratViolin(sro,"MTG_simil", "res.0.8",cells.use=cells.use,plot.attribs=plot.attribs),
stop("unknown graph")))}

makeFencyLayout <- function(sro, which){
	library(grid)
	p1 <- makegraphs(sro, which[1])
	par(new=TRUE)
	p2 <- makegraphs(sro, which[2])

	vplayout <- function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(2,1)))
	print(p1, vp=vplayout(1,1))
	print(p2, vp=vplayout(2,1))
}

# s for split, w for weight, xyz  
makeMultiSeuratGraph  <- function(sro, which, sannot, sremap = c(),  nrow = 1, swaprowcol=F , plot.attribs=c(), plot.attribs.edge=c() ){
	if (length(sannot) == 1){
		if (!sannot %in% colnames( sro@meta.data)) stop("no such meta data")
		sannot <- sro@meta.data[[sannot]]
	}
	if (class(sannot) != "factor") sannot <- factor(sannot, levels= unique(sannot))
	if (is.null(plot.attribs.edge)) plot.attribs.edge <- plot.attribs 
	
	#vplayout <- function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)
	#grid.newpage()
	#pushViewport(viewport(layout = grid.layout(length(annots),1)))

	if (!is.null(sremap)){ # remap annotations found in list, remapping to NA effectively filters
		remap <- match(levels(sannot), names(sremap))
		filter <- !is.na(remap)
		levels(sannot)[filter] <- as.character(sremap[remap[filter]])
	}
	if (length(levels(sannot)) > 10) print("Warning, trying to make a large ammount of plots")
	graphlist <- list()
	if (swaprowcol) {
		isedge <- rep(F,length(levels(sannot)))
		isedge[length(levels(sannot))] <- T
	}else{
		isedge <- rep(F,length(levels(sannot)))
		isedge[length(levels(sannot))] <- T
	}
	for(i in 1:length(levels(sannot))){
		if (isedge[i]) graphlist[[i]] <- makegraph(sro, which, rownames(sro@meta.data)[sannot@.Data == i],plot.attribs=c(plot.attribs.edge, list(title=levels(sannot)[i] )) )
		else graphlist[[i]] <- makegraph(sro, which, rownames(sro@meta.data)[sannot@.Data == i],plot.attribs=c(plot.attribs, list(title=levels(sannot)[i] )) )
	}
	if (swaprowcol) grid_arrange_shared_legend(graphlist, position="right",nrow = length(levels(sannot))/nrow ,ncol = nrow)
	else grid_arrange_shared_legend(graphlist, position="right",ncol = length(levels(sannot))/nrow ,nrow = nrow)
}
posnegremap <- function(sro, sannot){
	if (length(sannot) == 1){
		if (!sannot %in% colnames( sro@meta.data)) stop("no such meta data")
		sannot <- sro@meta.data[[sannot]]
	}
	if (class(sannot) == "factor") slabs <- levels(sannot)
	else slabs <- unique(sannot)

	ispos <- grepl("\\+", slabs)
	mapped <- rep("", length(ispos))
	mapped[ispos] <- "NeuN+"
	mapped[!ispos] <- "NeuN-"
return(setNames(as.list(mapped), slabs))}

posnegremap2 <- function(sro, sannot){
	if (length(sannot) == 1){
		if (!sannot %in% colnames( sro@meta.data)) stop("no such meta data")
		sannot <- sro@meta.data[[sannot]]
	}
	if (class(sannot) == "factor") slabs <- levels(sannot)
	else slabs <- unique(sannot)

	mapped <- gsub("\\+","",slabs)
	mapped <- gsub("\\-","",mapped)
	print(mapped)
return(setNames(as.list(mapped), slabs))}

makeMultiGraph <- function(sro, which){
return(switch(which,
	"M"=makeMultiSeuratGraph(sro, "P","orig.ident",swaprowcol=T,plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"))),
	"N"=makeMultiSeuratGraph(sro, "r","orig.ident",swaprowcol=T,plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"))),
	"O"=makeMultiSeuratGraph(sro, "r","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"))),
	"F"=makeMultiSeuratGraph(sro, "c","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,
plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident")),
plot.attribs.edge=list(flags=c("no.jitter","no.xlabel", "no.ylabel"),x.remap.class=posnegremap2(sro,"orig.ident"))
),
	"d"=makeMultiSeuratGraph(sro, "d","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,
plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident")),
plot.attribs.edge=list(flags=c("no.jitter","no.xlabel", "no.ylabel"),x.remap.class=posnegremap2(sro,"orig.ident"))
),


	"m"=makeMultiSeuratGraph(sro, "W","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,
plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident")),
plot.attribs.edge=list(flags=c("no.jitter","no.xlabel", "no.ylabel"),x.remap.class=posnegremap2(sro,"orig.ident"))
),
	"p"=makeMultiSeuratGraph(sro, "P","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,
plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident")),
plot.attribs.edge=list(flags=c("no.jitter","no.xlabel", "no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident"))
),
		"n"=makeMultiSeuratGraph(sro, "i","orig.ident",sremap=posnegremap(sro,"orig.ident"),swaprowcol=T,
plot.attribs=list(flags=c("no.jitter", "no.xlabel","no.ylabel","no.xnamedticks"),x.remap.class=posnegremap2(sro,"orig.ident")),
plot.attribs.edge=list(flags=c("no.jitter","no.xlabel", "no.ylabel"),x.remap.class=posnegremap2(sro,"orig.ident"))
),

stop("unknown graph")))}


getGGlegend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
return(tmp$grobs[[leg]])}

grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right"), do.share.legend=T,do.newpage=T) {
library(ggplot2)
library(gridExtra)
library(grid)
  if (do.share.legend){
  position <- match.arg(position)
  cur_legend <- getGGlegend(plots[[1]])
  lheight <- sum(cur_legend$height)
  lwidth <- sum(cur_legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            cur_legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           cur_legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth))
		)
  }else{
  gl <- lapply(plots, function(x) x )
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- arrangeGrob(gl)
  }
  if (do.newpage) grid.newpage()
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}


EdgeRcalc <- function(sro, cell.listPos, cell.listNeg){
	library(edgeR)
	dge <- DGEList(sro@raw.data[, c(cell.listPos, cell.listNeg)])
	group_edgeR <- factor(c(rep("P", length(cell.listPos)), rep("N",length(cell.listNeg))))
	dge <- estimateDisp(dge, design = model.matrix(~ group_edgeR), trend.method = "none")
	fit <- glmFit(dge, design)
	res <- glmLRT(fit)
	pVals <- res$table[,4]
	names(pVals) <- rownames(res$table)
	pVals <- p.adjust(pVals, method = "fdr")
	DE_Quality_AUC(pVals)	
}

plainDEvolcano <- function(mat, bulk_pos, bulk_neg){
	cold <- data.frame(rep("P",length(bulk_pos)))
	gdata <- data.frame(rownames(mat))
	colnames(gdata) <- c("Label")
	gdata$X <- 0
	gdata$Y <- 0
	for(i in 1:nrow(mat)){
		res <- t.test(mat[i,bulk_pos],mat[i,bulk_neg] )
		gdata$Y[i] = -log(res$p.value) / log(10)
		gdata$X[i] = (log(res$estimate[1]) - log(res$estimate[2])) / log(2)
	}
return(gdata)}
plainDESeq <- function(mat, bulk_pos, bulk_neg){
	library(DESeq2)
	flt <- bulk_pos | bulk_neg
	
	cold <- data.frame(rep("P",length(bulk_pos)))
	cold$X1[bulk_neg] = "N"
	print(table(cold$X1))
	colnames(cold) <- "Sample"
	cold <- cold[flt,]
	mat <- mat[,flt]
	deobj <- DESeqDataSetFromMatrix(mat,cold, design=~ Sample)	
	deobj <- tryCatch({DESeq(deobj)},  error = function(e) {
		print("Attempt 1 failed!")
		return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "iterate") ;DESeq(deobj)},error=function(e){
			print("Attempt 2 failed!")
			return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "ratio"); DESeq(deobj);},error= function(e) {
				print("all failed")
				return(c())}))}))})
return(deobj)}

DEcalc <- function(sro, cell.listPos, cell.listNeg, nbreplicates, use.meta.as.replicate= c(), cell.weight =c(), gene.use=c(), min.nbcell.threshold=10, bypass.posclass=c(), bypass.negclass=c(),do.fdr.correction=F, do.permutation.test="none"){
	library("DESeq2")
	library("Matrix")
#	if (length(cell.listPos) == ncol(sro@data)) cell.listPos <- colnames(sro@data)[cell.listPos]
#	if (length(cell.listNeg) == ncol(sro@data)) cell.listNeg <- colnames(sro@data)[cell.listNeg]


	if (is.null(gene.use)) gene.use <- rep(T,SeuratRawData(sro, "D")[1])
	if (is.null(use.meta.as.replicate)){
		if ((length(cell.listPos) < nbreplicates)||(length(cell.listNeg) < nbreplicates)){
			print("Warning, not enought cells")
			return(c())
		}
		dataN <- matrix(0,sum(gene.use), nbreplicates)
		dataP <- matrix(0,sum(gene.use), nbreplicates)
		nbp <- length(cell.listPos)
		nbn <- length(cell.listNeg)
		for(i in 1:nbreplicates){
			flt <- (1:nbp %% nbreplicates) == (i-1)
			if (sum(flt) == 1) dataP[,i] <- SeuratRawData(sro)[gene.use,cell.listPos[flt],drop=F ]
			else dataP[,i] <- (Matrix::rowSums(SeuratRawData(sro)[gene.use,cell.listPos[flt],drop=F]))  
			flt <- (1:nbn %% nbreplicates) == (i-1)
			if (sum(flt) == 1) dataN[,i] <- SeuratRawData(sro)[gene.use,cell.listNeg[flt],drop=F]
			else dataN[,i] <- (Matrix::rowSums(SeuratRawData(sro)[gene.use,cell.listNeg[flt],drop=F]))  
		}
		if (!is.null(bypass.negclass)) {dataN <- bypass.negclass; rownames(dataN) <- 1:sum(gene.use)}

		cold <- data.frame(c(rep("P",ncol(dataP)),rep("N",ncol(dataN))))
		data <- as.matrix(cbind(dataP,dataN))
	}else{ # uses cells associated with meta to build replicates
		mlvl <- getSafeLevels(sro@meta.data[[use.meta.as.replicate]])
		print(table(mlvl))
		print(paste("minimu is", min.nbcell.threshold))	
		print(paste("nbgene", sum(gene.use)))
		print("hello")
		if (do.permutation.test == "cell"){
			dataP <- data.frame(row.names = 1:sum(gene.use))
			dataN <- data.frame(row.names = 1:sum(gene.use))
			for(i in 1:length(mlvl)){
				flt <- (sro@meta.data[[use.meta.as.replicate]] == mlvl[i])
				nbp <- sum(flt & cell.listPos)
				nbn <- sum(flt & cell.listNeg)
				dalist <- sample(c(which(flt & cell.listPos), which(flt & cell.listNeg)))
				if (sum(nbp) >= min.nbcell.threshold) dataP[[paste(mlvl[i],"+",sep="")]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,dalist[1:nbp],drop=F]) 
				if (sum(nbn) >= min.nbcell.threshold) dataN[[paste(mlvl[i],"-",sep="")]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,dalist[(nbp+1):length(dalist)],drop=F])
			}
		}else{
		if (!is.null(bypass.posclass)){
			if (nrow(bypass.posclass) == sum(gene.use)) dataP <- bypass.posclass
			else dataP <- bypass.posclass[gene.use,]
			rownames(dataP) <- paste("Pos",1:sum(gene.use),sep="")
		}else{
			dataP <- data.frame(row.names = 1:sum(gene.use))
			for(i in 1:length(mlvl)){
				flt <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[i]) & cell.listPos
				flt2 <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[i]) & cell.listNeg 
				if ((sum(flt) >= min.nbcell.threshold)&(sum(flt2) >= min.nbcell.threshold)) dataP[[paste(mlvl[i],"+",sep="")]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt,drop=F]) 
			}
		}
		if (!is.null(bypass.negclass)) {
			if (nrow(bypass.negclass) == sum(gene.use)) dataN <- bypass.negclass
			else dataN <- bypass.negclass[gene.use,]
			rownames(dataN) <- paste("Neg",1:sum(gene.use),sep="")
		}else{
			dataN <- data.frame(row.names = 1:sum(gene.use))
			for(i in 1:length(mlvl)){
				flt <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[i]) & cell.listNeg 
				flt2 <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[i]) & cell.listPos 
				if ((sum(flt) >= min.nbcell.threshold)&(sum(flt2) >= min.nbcell.threshold)) dataN[[paste(mlvl[i],"-",sep="")]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt,drop=F])
			}
		}
		}
		

		cold <- data.frame(c(rep("P",ncol(dataP)),rep("N",ncol(dataN))))	
		data <- as.matrix(cbind(dataP,dataN))
		rownames(data) <- SeuratRawData(sro,"R")[gene.use]
	}
	print(paste(ncol(dataP),"pos and",ncol(dataN), "neg samples"))
	
	if ((ncol(dataP) == 0)||(ncol(dataN) == 0)) return(c())	
	if (do.permutation.test == "sample") cold  <- cold[sample(1:nrow(cold)),,drop=F]

	colnames(cold) <- "Sample"
	cold$Sample <- as.factor(cold$Sample)
	deobj <- DESeqDataSetFromMatrix(data,cold, design=~ Sample)
	print("Running DEseq")

	deobj <- tryCatch({DESeq(deobj)},  error = function(e) {
		print("Attempt 1 failed!")
		return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "iterate") ;DESeq(deobj)},error=function(e){
			print("Attempt 2 failed!")
			return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "ratio"); DESeq(deobj);},error= function(e) {
				print("all failed")
				return(c())}))}))})
	#	print("Attempt 1 failed!")
	#	return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "iterate") ;return(DESeq(deobj))},error=function(e){
	#		print("Attempt 2 failed!")
	#		return(tryCatch({deobj <- estimateSizeFactors(deobj,type= "ratio"); return(DESeq(deobj));},error= function(e) {
	#			print("all failed")
	#			return(c())}))}))})
	print("DESEQ done")
	if (is.null(deobj)) return(c())
	if (do.fdr.correction) out <-results(deobj, contrast=c("Sample", "P", "N"))
	else out <- results(deobj, contrast=c("Sample", "P", "N"),pAdjustMethod = "none")
	print("result done")

	out$pvalue[is.na(out$pvalue)] <- 1
	rownames(out) <- make.names(SeuratRawData(sro,"R")[gene.use],unique=T)
return(list(deobj = deobj, res=out))}

SROcalc <- function(sro, cell.listPos, cell.listNeg, test= "NB",do.fdr.correction=F){
	#out <- FindMarkers(sro, cells.1 = rownames(sro@meta.data)[cell.listPos], cells.2 = rownames(sro@meta.data)[cell.listNeg], test.use = test
	if ((length(cell.listPos) == 0)|(length(cell.listNeg) == 0)) return(c())
	sub <- SeuratRawData(sro)[,c(cell.listPos,cell.listNeg)]
	if (test == "MAST"){
		out <- Seurat:::MASTDETest(data.use = sub, cells.1=cell.listPos , cells.2=cell.listNeg,latent.vars=NULL,verbose = T )
	}else if (test == "NB"){
		flt <- (Matrix::rowSums(SeuratRawData(sro)[,cell.listPos] != 0) > 2)&(Matrix::rowSums(SeuratRawData(sro)[,cell.listNeg] != 0) > 2)
		print(sum(flt))
		out <- Seurat:::GLMDETest(data.use =sub[flt,], cells.1=cell.listPos , cells.2=cell.listNeg, test.use="negbinom")
	}else{
		out <- Seurat:::LRDETest(data.use =sub, cells.1=cell.listPos , cells.2=cell.listNeg)
	}
return(out)}

SeuratRawData <- function(input, getAttrib=c()){	
	if (is.null(getAttrib)){
		if (class(input) == "Seurat") {
			if (class(input@assays[[names(input@assays)[1]]]) == "Assay") return(input@assays[[names(input@assays)[1]]]@counts)
			else return(input@assays[[names(input@assays)[1]]])
		}else if (class(input) == "seurat") return(input@raw.data)
		else stop("Unknown class for input")
	}else if (grepl("[Rr]", getAttrib)){
		if (class(input) == "Seurat") {
			print(rownames(input@assays[[names(input@assays)[1]]]@counts)[1:2]) # not usre why... but that's silly print is needed
			if (class(input@assays[[names(input@assays)[1]]]) == "Assay") return(rownames(input@assays[[names(input@assays)[1]]]@counts))
			else return(rownames(input@assays[[names(input@assays)[1]]]))		
		}else if (class(input) == "seurat") return(rownames(input@raw.data))
		else stop("Unknown class for input")
	}else if (grepl("[Cc]", getAttrib)){
		if (class(input) == "Seurat") {
			if (class(input@assays[[names(input@assays)[1]]]) == "Assay") return(colnames(input@assays[[names(input@assays)[1]]]@counts))
			else return(colnames(input@assays[[names(input@assays)[1]]]))			
		}else if (class(input) == "seurat") return(colnames(input@raw.data))
		else stop("Unknown class for input")
	}else if (grepl("[Dd]", getAttrib)){
		if (class(input) == "Seurat") {
			if (class(input@assays[[names(input@assays)[1]]]) == "Assay") return(dim(input@assays[[names(input@assays)[1]]]@counts))
			else return(dim(input@assays[[names(input@assays)[1]]]))
		}else if (class(input) == "seurat") return(dim(input@raw.data))
		else stop("Unknown class for input")
	}else if (grepl("[Nn]", getAttrib)){
		if (class(input) == "Seurat") {
			if (class(input@assays[[names(input@assays)[1]]]) == "Assay") return(input@assays[[names(input@assays)[1]]]@data)
			else return(input@assays[[names(input@assays)[1]]]@data)
		}else if (class(input) == "seurat") return(input@data)
		else stop("Unknown class for input")
	}
}

SeuratGeneListQuery <- function(sro, metaquery = c(), listquery= c(),negate=F,default=T){
	srow <- SeuratRawData(sro,"R")
	if (!is.null(listquery)){
		if (length(listquery) != length(srow)) {
			listquery <- !is.null(match(listquery, srow))
		}else{
			listquery[is.na(listquery)] <- F
			listquery <- listquery
		}
	}else{
		if (is.null(metaquery)) return(rep(default,length(srow)))
		listquery <- rep(T,length(srow))
		print(length(srow))
		print(sum(listquery))
	}

	if (is.null(metaquery)) return(rep(T,length(srow)))
	if (length(metaquery) == 1) {
		if (!metaquery %in% colnames(sro@misc$meta.gene)) stop(paste("Missing", metaquery, "in seurat meta.data"))
		return(sro@misc$meta.gene[[metaquery]])
	}else{
		if (!metaquery[1] %in% colnames(sro@misc$meta.gene)) {
			if (substr(metaquery[1],1,1) %in% c(">", "<")){
				if (substr(metaquery[1],2,2) == "=") {
					if (substr(metaquery[1],1,1) == "<") dalist <- sro@misc$meta.gene[[substr(metaquery[1],3,10000)]] <= as.numeric(metaquery[2])
					else dalist <- sro@misc$meta.gene[[substr(metaquery[1],3,10000)]] >= as.numeric(metaquery[2])
				}else{
					if (substr(metaquery[1],1,1) == "<") dalist <- sro@misc$meta.gene[[substr(metaquery[1],2,10000)]] < as.numeric(metaquery[2])
					else dalist <- sro@misc$meta.gene[[substr(metaquery[1],2,10000)]] > as.numeric(metaquery[2])
				}
				dalist[is.na(dalist)] <- F
			}else{
				if (substr(metaquery[1],1,1) == "!"){
					negate = T
					metaquery[1] <- substr(metaquery[1],2,10000)
				}
				if (!metaquery[1] %in% colnames(sro@misc$meta.gene)) stop(paste("Missing", metaquery[1], "in seurat meta.data"))
				dalist <- !is.na(match(sro@misc$meta.gene[[metaquery[1]]] , metaquery[2:length(metaquery)]))
				
			}
		}else dalist <- !is.na(match(sro@misc$meta.gene[[metaquery[1]]] , metaquery[2:length(metaquery)]))
	}
	if (negate) dalist <- !dalist
	print("OUTOUT")
	print(metaquery)
	print(sum(dalist))
	print(sum(listquery))
	print(sum(dalist & listquery))
	return(dalist & listquery)
}





SeuratCellListQuery <- function(sro, metaquery=c(), listquery=c(),default=T,getnames=F,negate=F){
	if (!is.null(listquery)){
		if (length(listquery) != nrow(sro@meta.data)) {
			listquery <- !is.null(match(listquery, rownames(sro@meta.data)))
		}else{
			listquery[is.na(listquery)] <- F
			listquery <- listquery
		}
		if (is.null(metaquery)) return(listquery)
	}else{
		if (is.null(metaquery)) return(rep(default,nrow(sro@meta.data)))
		listquery <- rep(T,nrow(sro@meta.data))
	}
	if (class(metaquery) == "list"){
		dalist <- rep(T,nrow(sro@meta.data))
		for(i in names(metaquery)){
			if (!i %in% colnames(sro@meta.data)) stop(paste("Missing", i, "in seurat meta.data"))
			listquery = listquery & (sro@meta.data[[i]] %in% c(metaquery[[i]]))
		}
	}else{
	if (length(metaquery) == 1) {
		if (!metaquery %in% colnames(sro@meta.data)) stop(paste("Missing", metaquery, "in seurat meta.data"))
		dalist <- sro@meta.data[[metaquery]]
	}else{
		if (!metaquery[1] %in% colnames(sro@meta.data)) {
			if (substr(metaquery[1],1,1) %in% c(">", "<")){
				if (substr(metaquery[1],2,2) == "=") {
					if (substr(metaquery[1],1,1) == "<") dalist <- sro@meta.data[[substr(metaquery[1],3,10000)]] <= as.numeric(metaquery[2])
					else dalist <- sro@meta.data[[substr(metaquery[1],3,10000)]] >= as.numeric(metaquery[2])
				}else{
					if (substr(metaquery[1],1,1) == "<") dalist <- sro@meta.data[[substr(metaquery[1],2,10000)]] < as.numeric(metaquery[2])
					else dalist <- sro@meta.data[[substr(metaquery[1],2,10000)]] > as.numeric(metaquery[2])
				}
				dalist[is.na(dalist)] <- F
			}else{
				if (substr(metaquery[1],1,1) == "!"){
					negate = T
					metaquery[1] <- substr(metaquery[1],2,10000)
				}
				if (!metaquery[1] %in% colnames(sro@meta.data)) stop(paste("Missing", metaquery[1], "in seurat meta.data"))
				dalist <- !is.na(match(sro@meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
			}
		}else dalist <- !is.na(match(sro@meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
	}
        }
        if (negate) dalist <- !dalist
        if (getnames) return(rownames(sro@meta.data)[dalist & listquery]) else return(dalist & listquery)
}
SeuratCellMetaPermutation <- function(sro, nb.permutations,meta.partition, part.factor){
        print("Generating permutations")
        fout <-  matrix(0,nrow(sro@meta.data),nb.permutations)
        for(i in 1:length(levels(part.factor))){
                tmptab <- table(sro@meta.data[[meta.partition]]@.Data[part.factor@.Data == i])
                tab <- c();
                for(j in 1:length(tmptab)) tab <- c(tab, rep(names(tmptab)[j], tmptab[j]))
                for(j in 1:nb.permutations) fout[part.factor@.Data == i,j] <- sample(tab)
        }
return(fout)}

DEcalcMarkers <- function(sro, meta.partmarkers, meta.partition=c(),meta.use=c(), cell.use=c(), use.meta.as.replicate= c(),nb.replicates = 4, bypass.rownames=c(), gene.use=c(), min.nbcell.threshold=10){
        cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	pmarklvl <- levels(sro@meta.data[[meta.partmarkers]])
        out <- list()
	print(pmarklvl)
        for(i in 1:length(levels(sro@meta.data[[meta.partmarkers]]))){
		print(paste(levels(sro@meta.data[[meta.partmarkers]])[i], " got ", sum(sro@meta.data[[meta.partmarkers]]@.Data[cell.use] == i), " cells to comparte to the rest"))
                if (sum(sro@meta.data[[meta.partmarkers]]@.Data[cell.use] == i) != 0){
			print(paste(pmarklvl[i], ": ", sum(sro@meta.data[[meta.partmarkers]]@.Data[cell.use] == i)))
                        out[[gsub(" ", "_",pmarklvl[i])]] <- DEcalcIntersection(sro,meta.partition=meta.partition,cell.positive= (sro@meta.data[[meta.partmarkers]]@.Data == i), cell.use=cell.use,use.meta.as.replicate=use.meta.as.replicate,bypass.rownames=bypass.rownames, gene.use=gene.use, min.nbcell.threshold= min.nbcell.threshold)

                }
        }
return(out)}

SROcalcIntersection <- function(sro, meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), use.meta.as.replicate= c(),nb.replicates = 0, cell.weight =c(), bypass.rownames=c(), do.gprofile=T, gene.use=c(), gene.metause=c(), enqueue.go=F, namask.pvalue.threshold=0, bypass.negclass=c(),do.fdr.correction=F, do.permutation.test=c("none","cell","sample"),nb.partition.permute=0, sampling.meta=c(),min.nbcell.threshold=10){

        library(Seurat)
        cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
        cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)

        if (is.null(meta.partition)) lvl <- "All"
        else if (length(meta.partition) == 1) {
                if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
                else lvl <- unique(sro@meta.data[[meta.partition]])
        }else{
                lvl <- meta.partition[2:length(meta.partition)]
                cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
                meta.partition <- meta.partition[1]
        }

        tmpout <- list()
        for(i in 1:length(lvl)){
                if (is.null(meta.partition)) flt <- cell.use
                else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
                print(paste(lvl[i], "has", sum(cell.positive & flt), "to", sum(!cell.positive & flt), "to compare"))
                tmpout[[lvl[i]]] <- SROcalc(sro, which(cell.positive & flt) , which(flt &(!cell.positive)), do.fdr.correction=do.fdr.correction)
        }
return(tmpout)}


DEcalcIntersection <- function(sro, meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), use.meta.as.replicate= c(),nb.replicates = 0, cell.weight =c(), bypass.rownames=c(), do.gprofile=T, gene.use=c(), gene.metause=c(), enqueue.go=F, namask.pvalue.threshold=0, bypass.negclass=c(),do.fdr.correction=F, do.permutation.test=c("none","cell","sample"),nb.partition.permute=0, sampling.meta=c(),min.nbcell.threshold=10){
	if (nb.replicates == 0){
		if (is.null(use.meta.as.replicate)) nb.replicates = 4
		else nb.replicates = 1
	}

	if (length(do.permutation.test) == 3) do.permutation.test = "none"
	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	if (is.null(meta.partition)) lvl <- "All"
	else if (length(meta.partition) == 1) {
		if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
		else lvl <- unique(sro@meta.data[[meta.partition]])
	}else{
		lvl <- meta.partition[2:length(meta.partition)]
		cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
		meta.partition <- meta.partition[1]
	}

	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)
	gene.use <- SeuratGeneListQuery(sro, gene.metause, gene.use)
	print(paste("using",sum(gene.use), "genes"))
	
	if (is.null(gene.use)) stop("give gene list")
	if (is.null(cell.weight)) cell.weight <- rep(1,nrow(sro@meta.data))
	if (is.null(bypass.rownames)) bypass.rownames <- make.names(SeuratRawData(sro,"R")[gene.use],unique=T)
	pval <- data.frame(row.names=bypass.rownames)
	log2FC <- data.frame(row.names=bypass.rownames)
	meanExp <- data.frame(row.names=bypass.rownames)
	dropoutPosClass <- data.frame(row.names=bypass.rownames)
	dropoutNegClass <- data.frame(row.names=bypass.rownames)
	print(dim(pval))

	if (!is.null(sampling.meta)) cell.use <- cell.use & (runif(length(cell.use)) < sro@meta.data[[sampling.meta]]);


	if (nb.partition.permute != 0){
		permute.pval <- data.frame(row.names=bypass.rownames)
		print("Counting cell in celltypes:")	
		if (!is.null(use.meta.as.replicate)) dafactor <- as.factor(paste(sro@meta.data[[use.meta.as.replicate]],cell.use))
		else dafactor <- as.factor(cell.use)
		permflt <- SeuratCellMetaPermutation(sro,nb.partition.permute,meta.partition,part.fact= dafactor);
		permutation <- data.frame(row.names=bypass.rownames)
	}
	for(i in 1:length(lvl)){
		if (is.null(meta.partition)) flt <- cell.use
		else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
		print(paste("DE for ", sum(flt & cell.positive), "against", sum(flt & !cell.positive) ,"cells for class", lvl[i]));
             		out <- DEcalc(sro, flt & cell.positive, flt & !cell.positive, nb.replicates, use.meta.as.replicate=use.meta.as.replicate,gene.use=gene.use,bypass.negclass=bypass.negclass,do.fdr.correction=do.fdr.correction,do.permutation.test=do.permutation.test,min.nbcell.threshold=min.nbcell.threshold)
	#	out = tryCatch({DEcalc(sro, flt & cell.positive, flt & !cell.positive, nb.replicates, use.meta.as.replicate=use.meta.as.replicate,gene.use=gene.use)
        #	},  error = function(e) {print(paste("error detected for",lvl[i]));return(NULL)})
		print(paste("got", is.null(out)))	
		if (!is.null(out)){	
		pval[[ lvl[i] ]] <- out$res$pvalue
		log2FC[[ lvl[i] ]] <- out$res$log2FoldChange
		meanExp[[ lvl[i] ]] <- out$res$baseMean
		tmptmp <- Matrix::rowSums(SeuratRawData(sro)[gene.use,(flt & cell.positive)] != 0) / sum(flt & cell.positive)
		names(tmptmp) <-  make.names(rownames(tmptmp),unique=T)
		dropoutPosClass[[ lvl[i] ]] <- tmptmp
		tmptmp <- Matrix::rowSums(SeuratRawData(sro)[gene.use,(flt & (!cell.positive))] != 0) / sum(flt & (!cell.positive))
		names(tmptmp) <-  make.names(rownames(tmptmp),unique=T)
		dropoutNegClass[[ lvl[i] ]] <- tmptmp

		if (nb.partition.permute != 0){
			cmore <- rep(0.5,length(out$res$pvalue))
			tmptmp <- match(lvl[i], levels(sro@meta.data[[meta.partition]]))
			for(j in 1:nb.partition.permute){
	         		print(paste("Permutation no", j, "with", sum(cell.use & (permflt[,j] == tmptmp) & cell.positive), "against", sum(cell.use & (permflt[,j] == tmptmp) & !cell.positive) ,"cells for class", lvl[i]) )
				pout <- DEcalc(sro, (permflt[,j] == tmptmp) & cell.positive & cell.use, cell.use &(permflt[,j] == tmptmp) & !cell.positive, nb.replicates, use.meta.as.replicate=use.meta.as.replicate,gene.use=gene.use,bypass.negclass=bypass.negclass,do.fdr.correction=do.fdr.correction,do.permutation.test=do.permutation.test, min.nbcell.threshold=min.nbcell.threshold)
				reor <- match(rownames(out$res), rownames(pout$res) )
				
				cmpf <- pout$res$pvalue[reor] < out$res$pvalue
				cmore[cmpf] <- cmore[cmpf] + 1
				print(paste("outclassed by", sum(cmpf), "reord", sum(!is.na(reor))))
			}
			permute.pval[[lvl[i]]]	<- cmore / (nb.partition.permute + 1)

		}

		}
	}

#	tosort <- log2FC
#	tosort[sapply(tosort,is.infinite)] <- NA
#	if (namask.pvalue.threshold != 0) tosort[pval >= namask.pvalue.threshold] <- NA

#	daval <- rowSums(tosort,na.rm=T)
#	if ((length(daval)- sum(is.na(daval))) != 0){
#		reord <- order(daval,decreasing=T,na.last=NA)
#		log2FC <- log2FC[reord,,drop=F]
#		meanExp <- meanExp[reord,,drop=F]
#		pval <- pval[reord,,drop=F]
#		dropoutPosClass <- dropoutPosClass[reord,,drop=F]	
#		dropoutNegClass <- dropoutNegClass[reord,,drop=F]
#		if (nb.partition.permute != 0) permute.pval <- permute.pval[reord,,drop=F]
#	}
	if (nb.partition.permute != 0) extra <- list(permute.pval = permute.pval)
	else extra <- list()
	print("would exit...")
return(c(list(deseq.log10pvalue = log10(pval), deseq.log2FC= log2FC, deseq.basemean=meanExp,dropoutPosClass=dropoutPosClass, dropoutNegClass = dropoutNegClass, runargs=list(meta.partition=meta.partition, meta.positive=meta.positive,cell.positive=cell.positive,meta.use=meta.use,cell.use=cell.use) ), extra ))}

WilcHyperMarkers <- function(sro,  meta.partmarkers,meta.use=c(), cell.use=c(), use.meta.as.replicate= c(),nb.size.partition = 4, size.ordering= "nCount_RNA", bypass.rownames=c(), gene.use=c() , nb.best.mark=10){
        library(InferN0)
	pmarklvl <- getSafeLevels(sro@meta.data[[meta.partmarkers]])
        cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	scp <- InferN0Init(t(SeuratRawData(sro)))
        tmp <- as.factor(sro@meta.data[[meta.partmarkers]])
	if (length(levels(tmp)) > 1000) return(c())
	out <- Infern0FindMarkers(scp,nb.size.partition, tmp@.Data, order(sro@meta.data[[size.ordering]]))
#        for(i in 1:length(pmarklvl)){
 #               if (sum(sro@meta.data[[meta.partmarkers]]@.Data[cell.use] == i) != 0){
  #                      out[[pmarklvl[i]]] <- WilcHyperIntersection(sro,scp,meta.partition=meta.partition,cell.positive= (sro@meta.data[[meta.partmarkers]]@.Data == i), cell.use=cell.use,bypass.rownames=bypass.rownames, gene.use=gene.use)
                        
   #             }
  #      }


	nbmax <- rep(0, ncol(out$Zscore))
	modm <- out$Zscore
	dabestmark <- rep("tofill", ncol(out$Zscore) * nb.best.mark)
	nbdone=0;
	curord = order(modm,decreasing=T)
	k=1
	while(nbdone < ncol(out$Zscore)){
		i = floor((curord[k]-1) / nrow(out$Zscore))+1
		newone <- rownames(out$Zscore)[((curord[k] -1) %% nrow(out$Zscore))+1]
		if (!(newone %in% dabestmark)){
			nbmax[i] <- nbmax[i] + 1
			dabestmark[(i-1) * nb.best.mark + nbmax[i]] <- newone
			j <- ((curord[k] -1) %% nrow(out$Zscore))+1
			print(paste(i,"and", rownames(out$Zscore)[j], "with",modm[curord[k]]))
			if (nbmax[i] <  nb.best.mark) {
				k <- k +1;
				modm[ j, ] <- -1000
			}else{
				modm[, i] <- -1000
				curord = order(modm,decreasing=T)
				k <- 1
				nbdone <- nbdone + 1
			}
		}else{
			k <- k +1;
		}
	}
return(list(output=out, bestmarkers=dabestmark))}







WilcHyperIntersection <- function(sro, scp, sampling.prob=c(), meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), meta.partition.to.combine=c(),do.partition.Rside=F, cell.weight=c(),bypass.rownames=c(),namask.pvalue.threshold=0,do.plotZscores=0,do.quartile.average=T, do.permutation.test=F, nb.threads=4,do.downsample=F,do.reorder.bysignif=F,nbpart=4){
	library(InferN0)
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Expects a seurat object!")
	if (class(scp) != "InferN0Scope") stop("Expects an InferN0 scope!")
	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	if (is.null(meta.partition)) lvl <- "All"
	else if (length(meta.partition) == 1) {
		if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
		else lvl <- unique(sro@meta.data[[meta.partition]])
	}else{
		lvl <- meta.partition[2:length(meta.partition)]
		cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
		meta.partition <- meta.partition[1]
	}

	if (!is.null(meta.partition.to.combine)) cplvl <- unique(sro@meta.data[[meta.partition.to.combine]])

	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)
	#gene.use <- SeuratGeneListQuery(sro, gene.use, gene.metause)
	if (is.null(cell.weight)) cell.weight <- rep(1,nrow(sro@meta.data))
	if (is.null(bypass.rownames)) {bypass.rownames <- make.names(SeuratRawData(sro,"R")[gene.use],unique=T)}
	log10pval <- data.frame(row.names=bypass.rownames)
	Zscore <- data.frame(row.names=bypass.rownames)
	log2FC <- data.frame(row.names=bypass.rownames)
	meanTPM <- data.frame(row.names=bypass.rownames)
	LogitAuroc <- data.frame(row.names=bypass.rownames)
	dropoutPosClass <- data.frame(row.names=bypass.rownames)
	dropoutNegClass <- data.frame(row.names=bypass.rownames)
	CoverageEnrichment <- data.frame(row.names=bypass.rownames)
	Weight <- data.frame(row.names=bypass.rownames)
	gene.use <- match(scp$cell.names, make.names(SeuratRawData(sro,"R"),unique=T) )
	size <- Matrix::colSums(SeuratRawData(sro)[gene.use,,drop=F])
	print(paste("could not find", sum(is.na(gene.use))))
	print(scp$cell.names[is.na(gene.use)])
	if (sum(is.na(gene.use)) > 0) gene.use[is.na(gene.use)] <- F
	
	if (!is.null(sampling.prob)) sample.succ <- (runif(nrow(sro@meta.data)) < sampling.prob)
	else sample.succ <- rep(T, nrow(sro@meta.data))

	for(i in 1:length(lvl)){
		if (is.null(meta.partition)) flt <- cell.use
		else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
		if (is.null(sampling.prob)) print(paste("DE for ", sum(flt & cell.positive), "against", sum(flt & !cell.positive) ,"cells for class", lvl[i]))
		else {
			print(paste("DE for ", sum(sample.succ & flt & cell.positive), "against", sum(sample.succ & flt & !cell.positive) ,"cells for class", lvl[i], "(sampled from",sum(flt & cell.positive),"and",sum(flt & !cell.positive),")") );
			flt <- flt & sample.succ;
		}

		if ((sum(flt & cell.positive) != 0)&&(sum(flt & !cell.positive) != 0)){	

		if (is.null(meta.partition.to.combine)||(do.partition.Rside == F)){
			if (!is.null(meta.partition.to.combine)) {

				if (do.permutation.test){
					nb <- sum(flt)
					nbp <- sum(flt & cell.positive)
					cell.positive[flt] <- sample(c(rep(T,nbp), rep(F,nb-nbp)))
				}
				out <- InferN0ZeroCorrectedWilcox2(scp, flt & cell.positive, flt & !cell.positive, list.partition = sro@meta.data[[meta.partition.to.combine]]@.Data, total_count=size,do.plotZscores=do.plotZscores,do.quartile.average=do.quartile.average, nb.threads= nb.threads,do.downsample=do.downsample,nbpart=nbpart)
			
			}else{
				if (do.permutation.test){
					nb <- sum(flt)
					nbp <- sum(flt & cell.positive)
					print(nbp)
					print(nb-nbp)
					cell.positive[flt] <- sample(c(rep(T,nbp), rep(F,nb-nbp)))
				}
	            		out <- InferN0ZeroCorrectedWilcox2(scp, flt & cell.positive, flt & !cell.positive, total_count=size,do.plotZscores=do.plotZscores,do.quartile.average=do.quartile.average, nb.threads=nb.threads,do.downsample=do.downsample,nbpart=nbpart)
			}
		}else{
			nbin =0
			for( j in 1:length(cplvl)){
				print(paste("within ", cplvl[j], ", which is ", sum( (sro@meta.data[[meta.partition.to.combine]] == cplvl[j]) &  flt & cell.positive), "against", sum( (sro@meta.data[[meta.partition.to.combine]] == cplvl[j]) & flt & !cell.positive) ));
				nbcells = sum((sro@meta.data[[meta.partition.to.combine]] == cplvl[j]) &  flt)
				if (nbcells != 0){ 
			out2 <- InferN0ZeroCorrectedWilcox2(scp, (sro@meta.data[[meta.partition.to.combine]] == cplvl[j]) & flt & cell.positive , (sro@meta.data[[meta.partition.to.combine]] == cplvl[j]) & flt & !cell.positive, total_count=size,do.plotZscores=do.plotZscores, nb.threads= nb.threads,nbpart=nbpart)
			if (!is.null(out2)){
				out2$Zscore <- out2$Zscore * sqrt(nbcells)
				out2$log2FC <- out2$log2FC * nbcells
				out2$MeanTPM <- out2$MeanTPM * nbcells
				out2$LogitAuroc <- out2$LogitAuroc *nbcells
				out2$CoverageEnrichment <- out2$CoverageEnrichment * nbcells
				if (nbin == 0) {out <- out2}
				else{
				cpout$Zscore <- out$Zscore + out2$Zscore
					out$log2FC <- out$log2FC + out2$log2FC
					out$MeanTPM <- out$MeanTPM + out2$MeanTPM
					out$LogitAuroc <- out$LogitAuroc + out2$LogitAuroc
					out$CoverageEnrichment <- out$CoverageEnrichment + out2$CoverageEnrichment

				}
				nbin <- nbin + nbcells 
				}
			}
			}
			if (nbin == 0) out <- NULL
			else{
				out$Zscore <- out$Zscore /sqrt(nbin)
				out$log2FC <- out$log2FC / nbin
				out$MeanTPM <- out$MeanTPM / nbin
				out$LogitAuroc <- out$LogitAuroc / nbin
				out$CoverageEnrichment <- out$CoverageEnrichment / nbin 
			}

		}

		if (!is.null(out)){
			if (is.null(dim(out$Zscore))){
			Zscore[[ lvl[i] ]] <- out$Zscore;
			log10pval[[ lvl[i] ]] <- ( pnorm(-out$Zscore,log.p=T) + log(2)) / log(10); 					
			flippval <- out$Zscore < 0 ; 				
			flippval[is.na(flippval)] <- F ; 				
			log10pval[[ lvl[i] ]][flippval] <- (pnorm(out$Zscore[flippval], log.p=T) + log(2)) / log(10);	
			log2FC[[ lvl[i] ]] <- out$log2FC;
			LogitAuroc[[ lvl[i] ]] <- out$LogitAuroc;					
			meanTPM[[ lvl[i] ]] <- out$MeanTPM;
			CoverageEnrichment[[ lvl[i] ]] <- out$CoverageEnrichment;
			Weight[[ lvl[i] ]] <- out$Weight;	
			print(paste(sum(log10pval[[lvl[i]]] < -1.3),"out of",length(out$Zscore), "are significant"));
			print("Zscore range:")
			print(range(out$Zscore))
			dropoutPosClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & cell.positive,drop=F] != 0) / sum(flt & cell.positive)
			dropoutNegClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & (!cell.positive),drop=F] != 0) / sum(flt & (!cell.positive))
		}else{
			print("out$Zscore is null!!!")
				for(j in 1:ncol(out$Zscore)){
					Zscore[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- out$Zscore
					log10pval[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- ( pnorm(-out$Zscore,log.p=T) + log(2)) / log(10);	
			flippval <- out$Zscore < 0
			flippval[is.na(flippval)] <- F
			log10pval[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]][flippval] <- (pnorm(out$Zscore[flippval], log.p=T) + log(2)) / log(10);	
			log2FC[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- out$log2FC;
			LogitAuroc[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- out$LogitAuroc;	
			meanTPM[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- out$MeanTPM;
			CoverageEnrichment[[ paste(lvl[i],colnames(out$Zscore)[j], sep="_") ]] <- out$CoverageEnrichment;
				
			}
			
			dropoutPosClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & cell.positive,drop=F] != 0) / sum(flt & cell.positive)
			dropoutNegClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & (!cell.positive),drop=F] != 0) / sum(flt & (!cell.positive))
			}	
		}else{print("out is null");}

		}
	}

	if (do.reorder.bysignif){
		tosort <- log2FC
		tosort[sapply(tosort,is.infinite)] <- NA
		if (namask.pvalue.threshold != 0) tosort[log10pval >= log(namask.pvalue.threshold) / log(10)] <- NA

		daval <- rowSums(tosort,na.rm=T)
		if ((length(daval)- sum(is.na(daval))) != 0){
			reord <- order(daval,decreasing=T,na.last=NA)
			log2FC <- log2FC[reord,,drop=F]
			meanTPM <- meanTPM[reord,,drop=F]
			Zscore <- Zscore[reord,,drop=F]
			log10pval <- log10pval[reord,,drop=F]
			LogitAuroc <- LogitAuroc[reord,,drop=F]
			Weight <- Weight[reord,,drop=F]
	
			CoverageEnrichment <- CoverageEnrichment[reord,,drop=F]
			dropoutPosClass <- dropoutPosClass[reord,,drop=F]
			dropoutNegClass <- dropoutNegClass[reord,,drop=F]
		}
	}
	Weight <- Weight[1,]
return(list(wilcox.log10pval = log10pval, wilcox.zscore= Zscore, wilcox.log2FC= log2FC,wilcox.logitAuroc= LogitAuroc,  meanTPM=meanTPM, CoverageEnrichment= CoverageEnrichment, dropoutPosClass=dropoutPosClass,dropoutNegClass=dropoutNegClass, Weight=Weight))}

WilcHyperCrossPermutation <- function(sro, scp, sampling.prob=c(), nb.permute=399, meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), meta.partition.to.combine=c(),bypass.rownames=c(),do.plotZscores=0,nb.threads=4,do.downsample=F,nbpart=4){
	if (nb.permute == 0) return(WilcHyperIntersection(sro, scp, sampling.prob=sampling.prob, meta.partition=meta.partition, meta.positive=meta.positive,  cell.positive =cell.positive,meta.use=meta.use, cell.use=cell.use, meta.partition.to.combine=meta.partition.to.combine,bypass.rownames=bypass.rownames,do.plotZscores=do.plotZscores, nb.threads=nb.threads,do.downsample=do.downsample,nbpart=nbpart ))


	library(InferN0)
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Expects a seurat object!")
	if (class(scp) != "InferN0Scope") stop("Expects an InferN0 scope!")
	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	if (is.null(meta.partition)) lvl <- "All"
	else if (length(meta.partition) == 1) {
		if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
		else lvl <- unique(sro@meta.data[[meta.partition]])
	}else{
		lvl <- meta.partition[2:length(meta.partition)]
		cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
		meta.partition <- meta.partition[1]
	}

	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)
	permute.pval <- data.frame(row.names=bypass.rownames)

	log10pval <- data.frame(row.names=bypass.rownames)
	Zscore <- data.frame(row.names=bypass.rownames)
	log2FC <- data.frame(row.names=bypass.rownames)
	meanTPM <- data.frame(row.names=bypass.rownames)
	LogitAuroc <- data.frame(row.names=bypass.rownames)
	dropoutPosClass <- data.frame(row.names=bypass.rownames)
	dropoutNegClass <- data.frame(row.names=bypass.rownames)
	CoverageEnrichment <- data.frame(row.names=bypass.rownames)
	
	Weight <- data.frame(row.names=bypass.rownames)
	gene.use <- match(scp$cell.names, make.names(SeuratRawData(sro,"R"),unique=T) )
	size <- Matrix::colSums(SeuratRawData(sro)[gene.use,,drop=F])
	print(paste("could not find", sum(is.na(gene.use))))
	print(scp$cell.names[is.na(gene.use)])
	print("Counting cell in celltypes:")	

	if (!is.null(sampling.prob)) cell.use <- cell.use & (runif(length(cell.use)) < sampling.prob)

	tmptab <- table(sro@meta.data[[meta.partition]]@.Data[cell.positive & cell.use])
	ptab <- c();
	for(i in 1:length(tmptab)) ptab <- c(ptab, rep(names(tmptab)[i], tmptab[i]))
	tmptab <- table(sro@meta.data[[meta.partition]]@.Data[(!cell.positive) & cell.use])
	ntab <- c();
	for(i in 1:length(tmptab)) ntab <- c(ntab, rep(names(tmptab)[i], tmptab[i]))


	print("Generating permutations")	
	permflt <- matrix(0,length(cell.use),nb.permute)
	for(j in 1:nb.permute){
		permflt[cell.positive & cell.use, j] <- sample(ptab)
		permflt[(!cell.positive) & cell.use, j] <- sample(ntab)	
	}

	for(i in 1:length(lvl)){
		if (is.null(meta.partition)) flt <- cell.use
		else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
		print(paste("DE for ", sum(flt & cell.positive), "against", sum(flt & !cell.positive) ,"cells for class", lvl[i]));
		
		if ((sum(flt & cell.positive) != 0)&&(sum(flt & !cell.positive) != 0)){	
		
		

		out <- InferN0ZeroCorrectedWilcox2(scp, flt & cell.positive, flt & !cell.positive, total_count=size,do.plotZscores=do.plotZscores,do.quartile.average=T,nb.threads=nb.threads,do.downsample=do.downsample)
		cmore <- rep(0.5,length(out$Zscore))
		if (!is.null(out)){
		
		Zsng <- rep(1, length(out$Zscore))
		dflt <- (out$Zscore < 0)
		dflt[is.na(dflt)] <- F
		Zsng[dflt] <- rep(-1,sum(dflt))
		

		for( j in 1:nb.permute){
			tmptmp <- match(lvl[i], levels(sro@meta.data[[meta.partition]]))
			print(paste("Permuted DE no",j, "for ", sum((permflt[,j] == tmptmp) & cell.positive), "against", sum((permflt[,j] == tmptmp) & (!cell.positive)) ,"cells for class", lvl[i]));
			pout <- InferN0ZeroCorrectedWilcox2(scp, (permflt[,j] == tmptmp) & cell.positive, (permflt[,j] == tmptmp) & (!cell.positive), total_count=size,do.plotZscores=do.plotZscores,do.quartile.average=T,nb.threads=nb.threads,do.downsample=do.downsample)
			#hist(pout$Zscore, breaks=100)
			cmpf <- (pout$Zscore - out$Zscore) * Zsng > 0
			cmpf[is.na(cmpf)] <- T
			cmore[cmpf] <- cmore[cmpf] + 1
		}

			permute.pval[[lvl[i] ]] <- cmore / (nb.permute+1)
			Zscore[[ lvl[i] ]] <- out$Zscore;
			log10pval[[ lvl[i] ]] <- ( pnorm(-out$Zscore,log.p=T) + log(2)) / log(10); 					
			flippval <- out$Zscore < 0 ; 				
			flippval[is.na(flippval)] <- F ; 				
			log10pval[[ lvl[i] ]][flippval] <- (pnorm(out$Zscore[flippval], log.p=T) + log(2)) / log(10);	
			log2FC[[ lvl[i] ]] <- out$log2FC;
			LogitAuroc[[ lvl[i] ]] <- out$LogitAuroc;					
			meanTPM[[ lvl[i] ]] <- out$MeanTPM;
			CoverageEnrichment[[ lvl[i] ]] <- out$CoverageEnrichment;
			Weight[[ lvl[i] ]] <- out$Weight;	

			dropoutPosClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & cell.positive,drop=F] != 0) / sum(flt & cell.positive)
			dropoutNegClass[[ lvl[i] ]] <- Matrix::rowSums(SeuratRawData(sro)[gene.use,flt & (!cell.positive),drop=F] != 0) / sum(flt & (!cell.positive))
		}
	}}
	
	return(list(wilcox.log10pval = log10pval, wilcox.zscore= Zscore, wilcox.log2FC= log2FC,wilcox.logitAuroc= LogitAuroc,  meanTPM=meanTPM, CoverageEnrichment= CoverageEnrichment, dropoutPosClass=dropoutPosClass,dropoutNegClass=dropoutNegClass, Weight=Weight,permute.pval=permute.pval))}


WilcHyperOverwrite <- function(sro, scp, sampling.prob=c(), meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), meta.partition.to.combine=c(),do.partition.Rside=F, cell.weight=c(),namask.pvalue.threshold=0,nb.threads=4,do.downsample=F,do.reorder.bysignif=F,do.fdr.correct=T){
	library(InferN0)
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Expects a seurat object!")
	if (class(scp) != "InferN0Scope") stop("Expects an InferN0 scope!")
	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	if (is.null(meta.partition)) lvl <- "All"
	else if (length(meta.partition) == 1) {
		if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
		else lvl <- unique(sro@meta.data[[meta.partition]])
	}else{
		lvl <- meta.partition[2:length(meta.partition)]
		cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
		meta.partition <- meta.partition[1]
	}
	if (!is.null(meta.partition.to.combine)) cplvl <- unique(sro@meta.data[[meta.partition.to.combine]])
	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)
	if (is.null(cell.weight)) cell.weight <- rep(1,nrow(sro@meta.data))
	gene.use <- match(scp$cell.names, make.names(SeuratRawData(sro,"R"),unique=T) )
	size <- Matrix::colSums(SeuratRawData(sro)[gene.use,,drop=F])
	if (sum(is.na(gene.use)) > 0) gene.use[is.na(gene.use)] <- F
	if (!is.null(sampling.prob)) sample.succ <- (runif(nrow(sro@meta.data)) < sampling.prob)
	else sample.succ <- rep(T, nrow(sro@meta.data))
	for(i in 1:length(lvl)){
		if (is.null(meta.partition)) flt <- cell.use
		else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
		if (is.null(sampling.prob)) print(paste("DE for ", sum(flt & cell.positive), "against", sum(flt & !cell.positive) ,"cells for class", lvl[i]))
		else {
			print(paste("DE for ", sum(sample.succ & flt & cell.positive), "against", sum(sample.succ & flt & !cell.positive) ,"cells for class", lvl[i], "(sampled from",sum(flt & cell.positive),"and",sum(flt & !cell.positive),")") );
			flt <- flt & sample.succ;
		}
		if ((sum(flt & cell.positive) != 0)&&(sum(flt & !cell.positive) != 0)){	
	       		out <- InferN0ZeroCorrectedWilcox2(scp, flt & cell.positive, flt & !cell.positive, total_count=size,do.plotZscores=F,nbpart=1,do.quartile.average=F,nb.threads=nb.threads,do.downsample=do.downsample, do.output.overwritematrix=T,do.fdr.correct=do.fdr.correct)
		}
#		if (!is.null(out)){
#			# ooverwrite non-zero with DE
#			dapart <- t(sro@assays$RNA@data[,flt,drop=F])
#			for(i in 1:length(out$LogitAuroc)){
#				y = match(names(out$LogitAuroc)[i], colnames(dapart))
#				if (!is.na(y)){
#					l = dapart@p[y+1] - dapart@p[y]
#					if (l > 0) dapart@x[(dapart@p[y]+1):dapart@p[y+1]] <- rep(2.0 - 1.0 / (1.0 +  exp(out$LogitAuroc[i])) , l)
#				}
#			}
#			sro@assays$RNA@data[,flt] <- t(dapart)
#		}else{
#
#		}
	}
	sro@assays$RNA@data <- t(scp$data)
return(list(sro=sro,out=out))}




# Augmented probablistic index
# 
# @description If the model includes nuisance factors that are not of primary interest (for example 
# library size, batch, cell cycle ets), then the estimated probablistic index will be marginalized 
# over these factors to increase efficiency and statistical power. At the current version, it is 
# implemeted only for a two group comparison study. See Vermeulen et al 2015 for further.
# 
# @param pim.fit PIM output for a particular gene
# @param y a vector of the response variable 
# @param a grouping factor 
# @param x a matrix of covariates. A variable is presented in a column
# @param link a character vector with a single value that determines the used link function. 
# Possible values are "logit", "probit" and "identity". The default is "logit".
# @param y.pobs pseudo observations from a pim result
# @param alpha siginificance level
# 
# @return a list object
# 
# @references Vermeulen, K., Thas, O., & Vansteelandt, S. (2015). Increasing the power of the Mann‐Whitney test in randomized experiments through flexible covariate adjustment. \emph{Statistics in medicine}, 34(6), 1012-1030.
# @examples
# #example
#  
# @keywords internal
augmentedPI <- function(pim.fit, y, a, x,link="logit", y.pobs=NULL, alpha=0.05){
  
  a <- as.factor(a)
  if(length(unique(a))!=2){
    #message("Augmented marginal PI is not applicable for a grouping factor with level morethan 2.")
    return(c(Augmented.MPI=NA, 
             Std.Error=NA, 
             Wald.statistic=NA, 
             p.value=NA))
  }
  else if(is.null(x)){
   # message("Augmented marginal PI is not applicable onl if there is additional covariate.")
    return(c(Augmented.MPI=NA, 
             Std.Error=NA, 
             Wald.statistic=NA, 
             p.value=NA))
  }
  else{
    a <- ifelse(a==levels(a)[1], 0, 1)
    n <- length(a)
    p <- sum(a)/n
    
    
    x.dat    <- x
    if(is.null(colnames(x.dat))) {colnames(x.dat) <- paste0("U", 1:ncol(x))}
    
    x <- as.matrix(do.call(cbind, lapply(x, function(xx){
      if(is.factor(xx)) {
        as.numeric(xx)}
      ##as.numeric(ifelse(xx==levels(xx)[1], 0, 1))}
      else if(is.character(xx)){
        xx=as.factor(xx)
        as.numeric(xx)
        #as.numeric(ifelse(xx==levels(xx)[1], 0, 1))
        #message("One of the character variables in the U matrix is converted to a factor and a dummy variable is created")
      }
      else if(is.numeric(xx)){as.numeric(xx)}
      else{stop("Unable to identify the type of one of the variables in the U matrix")}
    })))
    #x <- as.matrix(x)
    
    if(is.null(colnames(x))) {colnames(x) <- colnames(x.dat)}
    data.sub <- as.data.frame(cbind(y,a,x)) 
    #data.sub <- as.data.frame(do.call(data.frame, data.sub))
    
    # Point estimate:
    # pim.fit <- pim::pim(formula=as.formula(paste("y~a+",paste(colnames(x.dat), collapse="+"))),
    #                     link=link, data=data.sub)
    coef <- pim::coef(pim.fit) 
    MW   <- 1-wilcox.test(y~a,exact=FALSE)$statistic/(sum(a)*sum(1-a))
    
    if(link=="probit"){
      augMW<-MW
      for(i in 1:n){
        augMW<-augMW+sum((1/(n*(n-1))-(1-a[i])*a[-i]/
                            (sum(a)*(n-sum(a))))*pnorm(coef[1]+t(-x[i,]+t(x[-i,]))
                                                       %*%coef[2:(dim(x)[2]+1)]))}
    }
    else if(link=="logit"){
      augMW<-0
      for(i in 1:n){
        augMW<-augMW+sum(expit(coef[1]+t(-x[i,]+t(x[-i,])) %*%coef[2:(dim(x)[2]+1)]))/(n*(n-1))
      }
    }
    else{stop("not a valid link function")}
    
    # Standard error:
    pseudo.y <- pseudo(y) 
    a1.hat<-sapply(1:length(a), A=a, pseudo.Y=pseudo.y, FUN=vec.a1hat)
    a2.hat<-sapply(1:length(a), A=a, pseudo.Y=pseudo.y, FUN=vec.a2hat)
    phi0<-vec.phi0(a,a1.hat,a2.hat,augMW)
    
    if(link=="probit"){
      pred<-t(sapply(coef=coef,x=x,1:dim(x)[1],vec.pred.probit))
      alphahat<-sapply(A=a,pred=pred,1:length(a),vec.alphahat)
      phiest<-phi0+(a-p)*(alphahat-mean(alphahat)+mean((1-a)*a1.hat/(1-p)^2-a*a2.hat/p^2))
    }
    else if(link=="logit"){
      pred<-t(sapply(1:dim(x)[1], coef=coef,x=x,vec.pred.logit))
      alphahat<-sapply(1:length(a),A=a,pred=pred,vec.alphahat)
      phiest<-phi0+(a-p)*alphahat
    }
    else{print("not a valid link function")}
    
    
    se.augMW<-sqrt(mean(phiest^2)/n) 
    W <-(augMW-0.5)/se.augMW # Wald test statistic:
    p.value<-2*pnorm(abs(W),lower.tail=FALSE) # p-value Wald test:
    # CI<-augMW+c(-1,1)*se.augMW*qnorm(1-alpha/2)
    
    
    return(c(Augmented.MPI=augMW, 
             Std.Error=se.augMW, 
             Wald.statistic=W, 
             p.value=p.value))
  }  
  
}

PimCalc <- function(expres.mat,X.mat, BPPARAM=BiocParallel::SerialParam(), do.approx=T){ 

#	for(i in 1:length(lvl)){
	
	#	y <- SeuratRawData(sro)[,cell.use]
	#	y@x <- y@x+1)
	
	#	mod <- pim.fit( )
	#	design <- pim::new.pim.formula(as.formula(paste("y", model.formula)), data=d)
	#	tag.pim.env <- pim::new.pim.env(d, ...)
	#	Y.pobs <- pim::response(design)
	#	model.tag <- pim.fit(x=pim::model.matrix(design), y=Y.pobs, penv = tag.pim.env, ...)
	
	if (!do.approx) expres.mat@x <- log2(expres.mat@x+1)

	colnames(X.mat) <- c("Positive","Rdepth")
	model.formula <- "~Positive+Rdepth" 
	
	pim.models.list <- BiocParallel::bplapply(seq_len(nrow(expres.mat)), function(tag){
      #pim.models.list <- lapply(seq_len(nrow(expres.mat)), function(tag){ 
      y <- expres.mat[tag,]
      d <- as.data.frame(cbind.data.frame(y,X.mat)) 
#	print(sum(y != 0))
	if (sum(y !=0) <= 5){
	}else{
      # Fit models

      tag.pim.env <- pim::new.pim.env(d) # Creat PIM enviroment 

      

	if (do.approx){
		d <- as.data.frame(do.call(data.frame, d))
		Cns <- rep(1, nrow(d))
		fml <- as.formula(paste("survival::Surv(y, Cns)", model.formula))
		model.tag <- try(coxph(fml, data = d, ties="efron"), silent=TRUE)


		if(model.tag != "error"){
			b=-1*coef(model.tag)
			names(b)      <- colnames(pim::vcov(model.tag))
			augmented.MPI <- c(Augmented.MPI=NA,  Std.Error=NA,  Wald.statistic=NA, p.value=NA)
			tag.pim.res   <- list(b=b, v=pim::vcov(model.tag), tag.name = rownames(expres.mat)[tag], augmented.MPI=augmented.MPI)
			#tag.pim.res
		#	print(model.tag)
		#	print(tag.pim.res)
			return(c(coef(summary(model.tag))[,5],model.tag$coefficients, sum(y[X.mat[,1] == 0] !=0), sum(y[X.mat[,1] != 0] !=0), sum(y[X.mat[,1] == 0]), sum(y[X.mat[,1] != 0])))
		}

	}else{
	      design      <- pim::new.pim.formula(as.formula(paste("y",model.formula,sep="")), data = d) # design
	      Y.pobs      <- pim::response(design) # Calculate pseudo observations  
		model.tag   <- tryCatch(suppressWarnings(pim.fit(x=pim::model.matrix(design), y=Y.pobs,penv = tag.pim.env)), error = function(e){"error"})
		if(model.tag != "error"){
        	        b=coef(model.tag)
		        names(b)      <- colnames(pim::vcov(model.tag))
		        augmented.MPI <- augmentedPI(pim.fit=model.tag, y=y, a=X.mat[,1,drop=F], x=X.mat[,2,drop=F], link=link, y.pobs=NULL)
		        tag.pim.res   <- list(b=b, v=pim::vcov(model.tag), tag.name = rownames(expres.mat)[tag],augmented.MPI=augmented.MPI)
		       tag.pim.res
		}
        }



	if (model.tag == "error"){
        print("got error") 
	   return("FAILED")
#     i   if(!is.null(U.mat)){
#         r <- sum(sapply(1:ncol(X.mat[,1,drop=F]), function(x){
#            if(is.factor(X.mat[,x])){length(levels(X.mat[,x]))-1}
#            else{1}
#          }))+ sum(sapply(1:ncol(X.mat[,2,drop=F]), function(x){
#            if(is.factor(U.mat[,x])){length(levels(U.mat[,x]))-1}
#            else{1}
#          })) 
 ##         
 #         names <- c(as.character(sapply(1:ncol(X.mat[,1,drop=F]), function(x){
 #           if(is.factor(X.mat[,x])){paste0(colnames(X.mat)[x],levels(X.mat[,x])[-1])}
 #           else{colnames(X.mat)[x]}
 #         })),
  #        as.character(sapply(1:ncol(U.mat), function(x){
  #          if(is.factor(U.mat[,x])){paste0(colnames(U.mat)[x],levels(U.mat[,x])[-1])}
  #          else{colnames(U.mat)[x]}
  ##        })))
  #      }else{
  #        r <- sum(sapply(1:ncol(X.mat), function(x){
  #          if(is.factor(X.mat[,x])){length(levels(X.mat[,x]))-1}
  #          else{1}
   #       }))
   #       
   #       names <- c(as.character(sapply(1:ncol(X.mat), function(x){
   #         if(is.factor(X.mat[,x])){paste0(colnames(X.mat)[x],levels(X.mat[,x])[-1])}
   #         else{colnames(X.mat)[x]}
   #       })))
  #      }
        
         
        # main.test <- gdata::unmatrix(main.test,byrow=FALSE)
        # main.test
        
  #      b=rep(NA, r) ; names(b) <- names
  #      v=matrix(rep(NA, r*r), r, r) ; rownames(v) <- colnames(v) <- names
  #      augmented.MPI=c(Augmented.MPI=NA,  Std.Error=NA,  Wald.statistic=NA, p.value=NA)
  #      tag.pim.res <- list(b=b, v=v, tag.name = rownames(expres.mat)[tag], augmented.MPI=augmented.MPI)
  #      tag.pim.res
      }
      }
    },
    BPPARAM = BPPARAM)
    #parallel::stopCluster(cl)

    out <- matrix(0,nrow(expres.mat),8)
    rownames(out) <- rownames(expres.mat)
    colnames(out) <- c("Positive","Rdepth", "Positive_coef", "Rcoef_coef", "geneCov_pos", "geneCov_neg", "nbUMI_pos", "nbUMI_neg")

    print(paste("got", length(pim.models.list),"and",nrow(expres.mat)))
    for(i in 1:length(pim.models.list)){
	if (is.null(pim.models.list[[i]])) out[i,] = NA
	else if (pim.models.list[[i]] == "FAILED") out[i,] = NA
	else out[i,] = pim.models.list[[i]] 

    }
    out <- out[!is.na(out[,1]),,drop=F]
return(out[order(out[,1]),,drop=F]);}

PimIntersection <- function(sro,meta.partition=c(), meta.positive=c(),  cell.positive =c(),meta.use=c(), cell.use=c(), BPPARAM=BiocParallel::SerialParam(), do.approx=T,nb.partition.permute=0){
	library(pim)
	library(survival)

	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	if (is.null(meta.partition)) lvl <- "All"
	else if (length(meta.partition) == 1) {
		if (class(sro@meta.data[[meta.partition]]) == "factor") lvl <- levels(sro@meta.data[[meta.partition]])
		else lvl <- unique(sro@meta.data[[meta.partition]])
	}else{
		lvl <- meta.partition[2:length(meta.partition)]
		cell.use <- cell.use & SeuratCellListQuery(sro, meta.partition)
		meta.partition <- meta.partition[1]
	}
	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)

	pval <- data.frame(row.names=SeuratRawData(sro, "R"))
        val <- data.frame(row.names=SeuratRawData(sro, "R"))

	for(i in 1:length(lvl)){
		print(paste("Processing", lvl[i]))
		if (is.null(meta.partition)) flt <- cell.use
		else flt <- (sro@meta.data[[meta.partition]] == lvl[i]) & cell.use
		expres.mat <- SeuratRawData(sro)[,flt]
		X.mat <- cbind(cell.positive[flt], log2(sro@meta.data$nCount_RNA[flt]))
		out <- PimCalc(expres.mat,X.mat,do.approx=do.approx, BPPARAM =BiocParallel::bpparam())
	
		newcol <- rep(NA, nrow(pval))
		map <- match(SeuratRawData(sro, "R"), rownames(out))
		newcol[!is.na(map)] <- out$Positive[map[!is.na(map)]]
		pval[[lvl[i]]] <- newcol
		newcol[!is.na(map)] <- out$Positive_coef[map[!is.na(map)]]
		val <- newcol
	}
	pval <- log10(pval)
	pval[is.na(pval)] <- 0
return(list(coxph.log10pval = pval, coxph.coef = val))}



sandwich.estimator<-function(U, U.diff, 
                             g1, g2, 
                             shared.factor=1, switched.factor=1,
                             self.factor=1)
{
  #note: this estimator is greatly optimized, based on T Lumley's code!
  
  if(any(fullyequal <- g1==g2 ))
  {
    #the contribution of these terms to variance and covariance is zero anyway, because
    #I_{ii} = 0.5 always (ie: constant), so we drop them out!
    U<-U[!fullyequal,,drop=FALSE]
    g1<-g1[!fullyequal]
    g2<-g2[!fullyequal]
  }
  
  usum1.tmp <- rowsum(U,g1,reorder=FALSE)
  usum2.tmp <- rowsum(U,g2,reorder=FALSE)
  
  ng <- length(union(g1,g2))
  
  usum1  <- matrix(nrow = ng, ncol = ncol(usum1.tmp),0)
  usum2  <- matrix(nrow = ng, ncol = ncol(usum2.tmp),0)
  usum1[unique(g1),] <- usum1.tmp
  usum2[unique(g2),] <- usum2.tmp
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1,usum2) + crossprod(usum2,usum1)
  uDiag<-crossprod(U) #Is counted twice as shared.factor, but needs to be counted as self.factor
  
  utu<-shared.factor*utushared  + 
    (switched.factor)*utuswitched +
    (self.factor-2*shared.factor)*(uDiag)
  
  #if the inverses occur (ij, ji), they are counted doubly as switched!!
  #However: they should be counted as - self
  
  mx<-ng+1
  
  uids<-g1*mx+g2
  invuids<-g2*mx+g1
  invrowperrow<-match(uids, invuids)
  rowsWithInv<-(!is.na(invrowperrow)) 
  if(any(rowsWithInv)) {
    rowsWithInv<-which(rowsWithInv)
    tmp<-sapply(rowsWithInv, function(ri){
      res<-U[ri,] %*% t(U[invrowperrow[ri],]) 
      return(res)
    })

    if(is.null(dim(tmp))){
      if(length(rowsWithInv)==1) tmp<-matrix(tmp, ncol=1)
      else tmp<-matrix(tmp, nrow=1)
    }
    
    uijji<-matrix(rowSums(tmp), ncol=ncol(usum1.tmp))
    
    utu<-utu-(2 *  switched.factor + self.factor)*uijji
  }
  
  #skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if(is.null(U.diff)) return (utu)
  
  U.diff.inv<-solve(U.diff)
  return(U.diff.inv%*%utu%*%U.diff.inv)
}



LigerIntersection <- function(sro, partition, genelist, use.meta.as.replicate=c(),  cell.positive=c(),meta.positive=c(), cell.use=c(), meta.use=c(), meta.gene=c(), use.gene=c(),bypass.posclass=c(), bypass.negclass=c() ,min.nbcell.threshold=10){
        partlvl <- getSafeLevels(sro@meta.data[[partition]])
	if (!is.null(use.meta.as.replicate)) mlvl <- getSafeLevels(sro@meta.data[[use.meta.as.replicate]])

	library(liger)
	cell.positive <- SeuratCellListQuery(sro, meta.positive, cell.positive)
	cell.use <- SeuratCellListQuery(sro, meta.use , cell.use)
	use.gene <- SeuratGeneListQuery(sro, meta.gene,use.gene)
	input <- rep(0, sum(use.gene))

	fout <- list()
	print(paste("Starting running liger on up to",sum(use.gene),"genes"))
	for(i in 1:length(partlvl)){
		flt <- sro@meta.data[[partition]]@.Data == i
		#inputpos <- sro@raw.data[use.gene,flt & cell.positive ,drop=F]
		#inputneg <- sro@raw.data[use.gene,flt & (!cell.positive) ,drop=F]
		#inputpos@x <- log2(inputpos@x + 1)
		#inputneg@x <- log2(inputneg@x + 1)
		#sp <- Matrix::rowSums(inputpos)
		#sn <- Matrix::rowSums(inputneg)
		#inputpos@x <- inputpos@x * inputpos@x
		#inputneg@x <- inputneg@x * inputneg@x
		#ssq <- Matrix::rowSums(inputpos) + Matrix::rowSums(inputneg)
		#n <- sum(flt)
		#cor = E((nx - sum x)*(ny- sum y)) / sqrt((ssq + (sp - sn)^2)  * (n - (ncol(inputpos) - ncol(inputneg))^2 ))
		#input = 
		if (!is.null(use.meta.as.replicate)){
			if (!is.null(bypass.posclass)){dataP <- bypass.posclass; rownames(dataP) <- paste("Pos",1:sum(use.gene),sep="")
			}else{
				dataP <- data.frame(row.names = 1:sum(use.gene))
				for(j in 1:length(mlvl)){
					lflt <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[j]) & flt & cell.positive 
					if (sum(lflt) >= min.nbcell.threshold) dataP[[paste(mlvl[j],"+",sep="")]] <- Matrix::rowSums(sro@raw.data[use.gene,lflt,drop=F]) 
				}
			}
			if (!is.null(bypass.negclass)) {dataN <- bypass.negclass; rownames(dataN) <- paste("Neg",1:sum(use.gene),sep="")
			}else{
				dataN <- data.frame(row.names = 1:sum(use.gene))
				for(j in 1:length(mlvl)){
					lflt <- ( sro@meta.data[[use.meta.as.replicate]] == mlvl[j]) & flt & (!cell.positive)
					if (sum(lflt) >= min.nbcell.threshold) dataN[[paste(mlvl[j],"-",sep="")]] <- Matrix::rowSums(sro@raw.data[use.gene,lflt,drop=F])
				}
			}
			objective <- c(rep(1,ncol(dataP)),rep(-1,ncol(dataN)))
			subinput <- t(log2(as.matrix(cbind(dataP,dataN))+1))
		}else{
			subinput <- sro@raw.data[flt, ]
			objective <- rep(-1, sum(flt)) 
			objective[cell.positive[flt]] <- 1
			subinput@x <- log2(subinput@x + 1)
			subinput <- log2(sro@raw.data[use.gene,flt,drop=F]+1)
		}
		print(paste("Processing", partlvl[i], "with", sum(flt),"samples/cells"))
		print(table(objective))
		input <- cor(subinput, objective)
		names(input) <- rownames(sro@raw.data)[use.gene]
		input <- input[!is.na(input)]
		if (class(genelist) == "list") {
			print("running bulk liger")
			fout[[partlvl[i]]] <- bulk.gsea(input, genelist)
		}else{
			print("running liger")
			 fout[[partlvl[i]]] <- gsea(input, genelist)
	}	}
return(fout)}

plotAgreement <- function(lista, listb, do.plot.instead=c(), plot.x.axis = "log2FC", plot.y.axis = "log2FC", pvalue.color.threshold= c(0.05, 0.05) , plot.attribs=c(), plot.grid.genelist=c()){
	#genelist <- rownames(lista$

	slota_FC <- grep( plot.x.axis, names(lista))
	slotb_FC <- grep( plot.y.axis, names(listb))

	genelist <- unique(c(rownames(lista[[slota_FC]]), rownames(listb[[slotb_FC]])))

	concat <- data.frame(row.names=genelist)
	for(i in colnames(lista[[slota_FC]])){
		concat[[paste("A_FC", i,sep="_")]] <- lista[[slota_FC]][[i]]
	}
	for(i in colnames(listb[[slota_FC]])){
		concat[[paste("B_FC", i,sep="_")]] <- listb[[slotb_FC]][[i]]
	}	

	slota_FC <- grep("pval", names(lista))
	slotb_FC <- grep("pval", names(listb))


	if (length(pvalue.color.threshold) == 1) pvalue.color.threshold <- c(pvalue.color.threshold,pvalue.color.threshold)

	for(i in colnames(lista[[slota_FC]])){
		selector <- (lista[[slota_FC]][[i]]/pvalue.color.threshold[1] < listb[[slota_FC]][[i]]/pvalue.color.threshold[2] )
		isna <- is.na(selector)
		selector[isna] <- F
		concat[[paste("minpvalue", i,sep="_")]][selector] <- lista[[slota_FC]][[i]][selector] /pvalue.color.threshold[1]
		concat[[paste("minpvalue", i,sep="_")]][!selector] <- listb[[slotb_FC]][[i]][!selector] /pvalue.color.threshold[2]
		concat[[paste("minpvalue", i,sep="_")]][isna] <- NA
		
		concat[[paste("a_pval", i,sep="_")]] <- lista[[slota_FC]][[i]]
		concat[[paste("b_pval", i,sep="_")]] <- listb[[slota_FC]][[i]]
	}


	if (!is.null(do.plot.instead)) {
		if (!("title" %in% names(plot.attribs))) plot.attribs$title <- do.plot.instead
		if (!paste("A_FC", do.plot.instead,sep="_") %in% names(concat)) stop(paste(do.plot.instead,"is not an colname in the input data"))

		alpha <- sapply(concat[[paste("minpvalue", do.plot.instead,sep="_")]], function(x){if (is.na(x)) return(0); if (x < 1) return(1); if (x > 2) return(0); return(2 - x)})
		
		nalist <- (is.na(concat[[paste("A_FC", do.plot.instead,sep="_")]]) | is.na( concat[[paste("B_FC", do.plot.instead,sep="_")]])|(alpha == 0) )

		if (!is.null(plot.grid.genelist)) {
			slota_FC <- grep( plot.x.axis, names(lista))
			slotb_FC <- grep( plot.y.axis, names(listb))
			cdata <- matrix(0,length(plot.grid.genelist),ncol(lista[[slota_FC]]) + ncol(listb[[slotb_FC]]))
			rownames(wdata) <- plot.grid.genelist
                        for(i in 1:ncol(lista[[slota_FC]])) cdata[,i] <- lista[[slota_FC]][plot.grid.genelist,i]
                        for(i in 1:ncol(listb[[slotb_FC]])) cdata[,i+ncol(lista[[slota_FC]])] <- listb[[slotb_FC]][plot.grid.genelist,i]
                        slota_FC <- grep( "pval", names(lista))
                        slotb_FC <- grep( "pval", names(listb))
			wdata <- matrix(0,length(plot.grid.genelist),ncol(lista[[slota_FC]]) + ncol(listb[[slotb_FC]]))
			if (length(grep("log",names(lista)[grep( "pval", names(lista))])) == 0) {
				for(i in 1:ncol(lista[[slota_FC]])) wdata[,i] <- 0.05 / (0.05 + lista[[slota_FC]][plot.grid.genelist,i])
                        }else{
				for(i in 1:ncol(lista[[slota_FC]])) wdata[,i] <- 0.05 / (0.05 + exp(lista[[slota_FC]][plot.grid.genelist,i] * log(10)))
			}
			if (length(grep("log",names(listb)[grep( "pval", names(listb))])) == 0) {
				for(i in 1:ncol(listb[[slotb_FC]])) wdata[,i+ncol(lista[[slota_FC]])] <- 0.05 / (0.05 + listb[[slotb_FC]][plot.grid.genelist,i])
                        }else{
				for(i in 1:ncol(listb[[slotb_FC]])) wdata[,i+ncol(lista[[slota_FC]])] <- 0.05 / (0.05 + exp(listb[[slotb_FC]][plot.grid.genelist,i] * log(10)))
			}
                        p <- plotDataGrid(cdata,wdata)
		}else{
				
			p = plotLabels(concat[[paste("A_FC", do.plot.instead,sep="_")]][!nalist], concat[[paste("B_FC", do.plot.instead,sep="_")]][!nalist], names=rownames(concat)[!nalist], color = concat[[paste("minpvalue", do.plot.instead,sep="_")]][!nalist],alpha=alpha[!nalist],plot.attribs=plot.attribs)
			p <- p + scale_color_gradientn(values= c(0.0,0.5,0.51,1),colors=c("#000000FF","#000000FF","#88888844","#88888844"))
		}
		return(p)

	}
return(concat)}






subSampledSeurat<- function(sro, fraction = 0.2, remove.sparse.zeros=F){
	datacur=1;
	datavector <- sro@raw.data@x
 	pb <- txtProgressBar(max = ncol(sro@raw.data), title="Subsampling SeuratObject", style = 3)
	for(i in 1:ncol(sro@raw.data)){
		setTxtProgressBar(pb, i)
		surrow <- datavector[(sro@raw.data@p[i]+1):sro@raw.data@p[i+1]]
		names(surrow) <- (sro@raw.data@p[i]+1):sro@raw.data@p[i+1]
		surrow <- sort(surrow);
		curgroup <- table(surrow)
		nbumi = sum(surrow)
		whtorem <- sort(sample(nbumi, floor(nbumi * fraction), replace=F))
		wh <- 1
		prior <- 0
		priorname <- 1
		for(j in 1:length(curgroup)){
			while(wh <= length(whtorem)){
				if (whtorem[wh] > prior + curgroup[j] * as.integer(names(curgroup[j]))) break;
				daname <- as.integer(names(surrow)[priorname + floor((whtorem[wh] - prior)/as.integer(names(curgroup[j])))])
				datavector[daname]  <- datavector[daname] - 1
				wh <- wh + 1
			}
			priorname <- priorname + curgroup[j]
			prior <- prior + curgroup[j] * as.integer(names(curgroup[j]))
		}
		#print(paste(nbumi, "->", sum(datavector[(sro@raw.data@p[i]+1):sro@raw.data@p[i+1]]) , sum(datavector[(sro@raw.data@p[i]+1):sro@raw.data@p[i+1]]) /nbumi))
	}
	print("done!")
	sro@raw.data@x <- datavector
	if (remove.sparse.zeros){
		print("Removing Zeroes")
		sro@raw.data <- Matrix(sro@raw.data, sparse=T)
	}
return(sro)}



parseGprofiler <-function(gprofilerout, gprofile.nbterms=0, bg.size,query, nbterm=0, do.bitmap=T, do.heatmap=T, plot.attribs=c(), gene.annot = c(), filter.lone.annot=T,filter.no.annot=T, heatmapquery=c(), allowed.classes = c("tf","CC","BP","MF","keg","REA")){
	if (is.null(heatmapquery)) heatmapquery <- query
	if (do.heatmap) do.bitmap=T
	if (class(gprofilerout) == "list"){
	}else{
		if (nbterm == 0) {
			if ((is.null(gprofilerout))||(nrow(gprofilerout) == 0)) return(list())
			nbterm <- nrow(gprofilerout)
		}
		fout <- data.frame(row.names=paste("Annot", 1:nbterm,sep=""))
		ord <- order(gprofilerout$p.value)
		
		if (nrow(gprofilerout) < nbterm){
			return(list())
		}else{
			fout[["Term"]] <- gprofilerout$term.name[ord]
			fout[["Enrich"]] <- gprofilerout$overlap.size * bg.size /  (gprofilerout$term.size * gprofilerout$query.size)
			fout[["pvalue"]] <- gprofilerout$p.value
		}
		fout <-  list(list=fout)
		if (do.bitmap){
			fout[["bitmap"]] <- matrix(F, nrow(gprofilerout), length(heatmapquery));
			rownames(fout[["bitmap"]]) <- make.names(gprofilerout$term.name[ord], unique=T)
			rownames(fout[["bitmap"]])[gprofilerout$domain[ord] == "tf"] <- gprofilerout$term.id[ord][gprofilerout$domain[ord] == "tf"]
			colnames(fout[["bitmap"]]) <- heatmapquery
			for(i in 1:nrow(gprofilerout)){
				parsed <-parselist( gprofilerout$intersection[ord[i]])
				for( j in parsed) {
				  tmpind <- match(j,heatmapquery) 
				  if (!is.na(tmpind)) fout[["bitmap"]][i, tmpind ] <- T
				}
			}
			if (do.heatmap) {
				curatt <- plotAttribs(plot.attribs, ylabel="Annotations", xlabel="Genes");
				map <- !is.na( match(gprofilerout$domain, allowed.classes))
				if ((sum(map) > gprofile.nbterms)&(gprofile.nbterms !=0)){
					input <- fout[["bitmap"]][map,,drop=F][order(gprofilerout$p.value[map])[1:gprofile.nbterms],, drop=F]
				}else input <- fout[["bitmap"]][map,,drop=F]
				if (filter.lone.annot) input <- input[rowSums(input) > 1,,drop=F] 
				if (filter.no.annot) input <- input[,colSums(input) > 0,drop=F]

				#input <- input * range(gene.annot)[2]
				rowannot <- data.frame(row.names=  rownames(gprofilerout));
				rowannot$log10pvalue <- -log10(gprofilerout$p.value)
				rowannot$category <- 0
				rowannot$category[gprofilerout$domain == "tf"] <- 1
				rowannot$category[gprofilerout$domain == "CC"] <- 2
				rowannot$category[gprofilerout$domain == "BP"] <- 3
				rowannot$category[gprofilerout$domain == "MF"] <- 4
				rowannot$category[gprofilerout$domain == "keg"] <- 5
				rowannot$category[gprofilerout$domain == "rea"] <- 6
				rowannot <- rowannot[ord,]
				rownames(rowannot) <- rownames(fout[["bitmap"]])
				fout[["heatmap"]] <- makeHeatmap(input, colannot=gene.annot, rowannot= rowannot, plot.attribs= curatt)
			}
		}
		return(fout);
	}
}



parseGprofiler2 <-function(gprofilerout, gprofile.nbterms=0, bg.size,query, nbterm=0, do.bitmap=T, plot.attribs=c(), gene.annot = c(), plot.maxlist=25, allowed.classes = c("tf","CC","BP","MF","keg","REA"), make.ggplots=T){
	library(stats)

	if (class(gprofilerout) == "list"){
	}else{
		if (nbterm == 0) {
			if ((is.null(gprofilerout))||(nrow(gprofilerout) == 0)) return(list())
			nbterm <- nrow(gprofilerout)
		}
		fout <- data.frame(row.names=paste("Annot", 1:nbterm,sep=""))
		ord <- order(gprofilerout$p.value)
		
		if (nrow(gprofilerout) < nbterm){
			return(list())
		}else{
			fout[["pvalue"]] <- gprofilerout$p.value[ord]
			fout[["Domain"]] <- gprofilerout$domain[ord]
			fout[["ID"]] <- gprofilerout$term.id[ord]
			fout[["Term"]] <- gprofilerout$term.name[ord]
			fout[["Enrich"]] <- gprofilerout$overlap.size[ord] * bg.size /  (gprofilerout$term.size[ord] * gprofilerout$query.size[ord])
			fout[["intersection"]] <- gprofilerout$intersection[ord]
		}
		fout <-  list(list=fout)
		
		bitmap <- matrix(F, nrow(gprofilerout), length(query));
		rownames(bitmap) <- make.names(gprofilerout$term.name, unique=T)
		rownames(bitmap)[gprofilerout$domain == "tf"] <- gprofilerout$term.id[gprofilerout$domain == "tf"]
		colnames(bitmap) <- query
		dom <- as.factor(gprofilerout$domain)
		pgeneannot <- gene.annot[match(query, rownames(gene.annot)),,drop=F]
		annot.annot <- matrix(F, nrow(gprofilerout) , ncol(gene.annot))

		for(i in 1:nrow(gprofilerout)){
			parsed <-parselist( gprofilerout$intersection[i])
			for( j in parsed) {
				tmpind <- match(j,query) 
				if (!is.na(tmpind)) bitmap[i, tmpind ] <- T
			}
			annot.annot[i,] <- colSums(pgeneannot[bitmap[i,],,drop=F]) / sum(bitmap[i,])
		}

		annot.annot <- cbind(annot.annot, -log10(gprofilerout$p.value))
		rownames(annot.annot) <- rownames(bitmap)
		colnames(annot.annot) <- c(colnames(gene.annot), "Enrich-log10pvalue")


		if (length(query) > plot.maxlist) heatmapquery <- query[1:plot.maxlist]
		else heatmapquery <- query
		if (length(heatmapquery) > 1) gene.order <- hclust(dist(t(bitmap[,1:length(heatmapquery),drop=F])), method="complete")$order
		else gene.order <- c(1)

		if (nrow(bitmap) > 1) {
			if (nrow(bitmap) > plot.maxlist) {
				annot.order <- ord[hclust(dist(bitmap[ord[1:plot.maxlist],,drop=F]), method="complete")$order]
			}else annot.order <- hclust(dist(bitmap), method="complete")$order
		}else annot.order <- c(1)

		rowflt <- !is.na(match(dom@.Data, match(c("CC", "BP", "MF"), levels(dom)) ))
		
		if (sum(rowflt) > 1){
			annotord <- ord[rowflt[ord]]
			if (length(annotord) > plot.maxlist) {
				go.order <- annotord[hclust(dist(bitmap[annotord[1:plot.maxlist],,drop=F]), method="complete")$order]
			}else go.order <- annotord[hclust(dist(bitmap[annotord,,drop=F]), method="complete")$order]
		}else go.order <- which(rowflt)

		rowflt <- !is.na(match(dom@.Data, match(c("keg", "rea"), levels(dom)) ))
		if (sum(rowflt) > 1){
			annotord <- ord[rowflt[ord]]
			if (length(annotord) > plot.maxlist) pw.order <- annotord[hclust(dist(bitmap[annotord[1:plot.maxlist],,drop=F]), method="complete")$order]
			else pw.order <- annotord[hclust(dist(bitmap[annotord,,drop=F]), method="complete")$order]
		}else pw.order <- which(rowflt)

		rowflt <- !is.na(match(dom@.Data, setdiff(1:length(levels(dom)),match(dom@.Data, c("CC","keg", "rea", "BP", "MF" ))), levels(dom)) )
		if (sum(rowflt) > 1){
			annotord <- ord[rowflt[ord]]
			if (length(annotord) > plot.maxlist) misc.order <- annotord[hclust(dist(bitmap[annotord[1:plot.maxlist],,drop=F]), method="complete")$order]
			else misc.order <- annotord[hclust(dist(bitmap[annotord,,drop=F]), method="complete")$order]
		}else misc.order <- which(rowflt)

		pgeneannot <- gene.annot[match(heatmapquery, rownames(gene.annot)),,drop=F]

		if (make.ggplots) {
			if (length(gene.order) > 0) fout[["plot.genes"]] <- plotDataGrid(pgeneannot[gene.order,,drop=F],plot.attribs=plot.attribs,transform=list(data = "colwise.erf"))
			if (length(annot.order) > 0) fout[["plot.annot"]] <- plotDataGrid(annot.annot[annot.order,,drop=F],plot.attribs=plot.attribs,transform=list(data = "colwise.erf"))
			if (length(go.order) > 0) fout[["plot.go"]] <- plotDataGrid(annot.annot[go.order,,drop=F],plot.attribs=plot.attribs,transform=list(data = "colwise.erf"))
			if (length(pw.order) > 0) fout[["plot.pathway"]] <- plotDataGrid(annot.annot[pw.order,,drop=F],plot.attribs=plot.attribs,transform=list(data = "colwise.erf"))
			if (length(misc.order) > 0) fout[["plot.misc.annot"]] <- plotDataGrid(annot.annot[misc.order,,drop=F],plot.attribs=plot.attribs,transform=list(data = "colwise.erf"))
		}

		#curatt <- plotAttribs(plot.attribs, ylabel="Annotations", xlabel="Genes");
		#map <- !is.na( match(gprofilerout$domain, allowed.classes))
		#if ((sum(map) > gprofile.nbterms)&(gprofile.nbterms !=0)){
		#	input <- bitmap[map,,drop=F][order(gprofilerout$p.value[map])[1:gprofile.nbterms],, drop=F]
		#}else input <- fout[["bitmap"]][map,,drop=F]

		#input <- input * range(gene.annot)[2]
		#rowannot <- data.frame(row.names=  rownames(gprofilerout));
		#rowannot$log10pvalue <- -log10(gprofilerout$p.value)
		#rowannot <- rowannot[ord,]
		#rownames(rowannot) <- rownames(fout[["bitmap"]])
		#fout[["heatmap"]] <- makeHeatmap(input, colannot=gene.annot, rowannot= rowannot, plot.attribs= curatt)
		
		return(fout);
	}
}

prepareExtremeFilter <- function(datalist, soupde, thresh.pval.tosoup = -1.3, thresh.expr.rank=0, rankgrep = "basemean", pickmax=T, thresh.dropout=0){
	filt <- (soupde$deseq.log10pvalue < thresh.pval.tosoup)
	if (thresh.expr.rank != 0){
		slot <- grep(rankgrep, names(datalist))
		if (is.na(slot)) stop(paste("Could not find entry in list containing", rankgrep))
		for(i in 1:ncol(datalist[[1]])) {tmp <- order(datalist[[slot]],decreasing=pickmax); filt[,i] <- filt[,i] & (!(is.na(match(1:nrow(base$FadNeuNN$deseq.basemean), tmp[1:thresh.expr.rank]))))  }
	}
	if (thresh.dropout != 0){
		for(i in 1:ncol(datalist[[1]])) filt[,i] <- filt[,i] & (!(is.na(match(1:nrow(base$FadNeuNN$deseq.basemean), tmp[1:thresh.expr.rank]))))
	}

return(filt)}

# just changed nbterm max...
getExtremes <- function(datalist, maxlist.toreport= 1000, plot.maxlist=25, value.name= c(), meta.names=c(), pval.grep.string= "pval", pval.threshold=0.05, filter.name.and.val=c(), filter.table=c(), do.fdr.correction=F,gprofile.args=list(hier_filtering="strong",ordered_query=T),gprofile.nbterms=100,plot.attribs=c(),selected.colnames=c(), extrafilter=c(),make.ggplots=F,funcfilt=c(), use.funky.threshold=c()){

        slot_pval <- grep( pval.grep.string, names(datalist))
        islog10pvalue <- (length(grep("log", names(datalist)[slot_pval])) == 1)
        pval.threshold = log10(pval.threshold)
        if (is.null(value.name)) slot_val <- 1
        else {
		slot_val <- grep( value.name, names(datalist))
		if (length(slot_val) == 0) stop(paste("Did not find attribute", value.name))
	}
	if (!is.null(filter.table)) {
		mflt <- matrix(F, nrow(datalist[[slot_pval]]), ncol(filter.table))

		rownames(mflt) <- rownames(datalist[[slot_pval]])
		colnames(mflt) <- colnames(filter.table)
#		print("yes")
#		print(colnames(filter.table))
		map <- match(rownames(datalist[[slot_pval]]), rownames(filter.table))
#		print(dim(mflt))
#		print(dim(filter.table))
#		print(sum(!is.na(map)))
#		print(colnames(mflt))
#		print("nsdewsize")

#		print(class(filter.table[map[!is.na(map)],]))
		mflt <- as.data.frame(mflt)
#		print(class(mflt))
		mflt[(!is.na(map)),] <- filter.table[map[!is.na(map)],]
#		print(class(mflt))
#		print("newsize")
#		print(colnames(mflt))
#		print(dim(mflt))
	}else {
		mflt <- c()
	}
        if (!is.null(filter.name.and.val)){
                slot_filter <- grep(filter.name.and.val[1], names(datalist))
		if (is.null(mflt)) mflt <- datalist[[slot_filter]] > as.numeric(filter.name.and.val[2])
		mflt <- mflt & (datalist[[slot_filter]] > as.numeric(filter.name.and.val[2]))
        }
	fout <- list()
	slot_meta <- rep(1,length(meta.names))
	if (gprofile.nbterms > 0) library(gProfileR)
		print(dim(mflt))

	for(i in 1:length(meta.names)) slot_meta[i] <- grep( meta.names[i], names(datalist))
	if (is.null(selected.colnames)) selected.colnames=colnames(datalist[[slot_val]])
	if (length(selected.colnames) != 0){
        for(i in 1:length(selected.colnames)){
		if (selected.colnames[i] %in% colnames(datalist[[slot_pval]])){
		if (length(slot_meta) == 0) topass <- c()
		else{
			for(j in 1:length(slot_meta)){
				if (j == 1) topass <- datalist[[slot_meta[j]]][,selected.colnames[i],drop=F]
				else topass <- cbind(topass,datalist[[slot_meta[j]]][,selected.colnames[i],drop=F])
			}
			colnames(topass) <- meta.names
		}
		if (!is.null(extrafilter)) curf <- (datalist[[extrafilter]][,selected.colnames[i],drop=T] < 0.05) #use dropout threshold?
                else curf <- rep(T, nrow(datalist[[slot_val]]))
                if (!is.null(mflt)) {
			#print(colnames(mflt))
			#print(selected.colnames)
			#print(i)
			#print(dim(mflt))
			#print(length(curf))
			curf <- curf & mflt[,selected.colnames[i]]

		}
		
		if (!is.null(topass)) toorder <- cbind(datalist[[slot_val]][curf,selected.colnames[i],drop=F], datalist[[slot_pval]][curf,selected.colnames[i],drop=F], datalist[[slot_pval]][curf,selected.colnames[i],drop=F], topass[curf,,drop=F])
		else toorder <- cbind(datalist[[slot_val]][curf,selected.colnames[i],drop=F], datalist[[slot_pval]][curf,selected.colnames[i],drop=F], , datalist[[slot_pval]][curf,selected.colnames[i],drop=F] )

		print(paste("Processing", selected.colnames[i], "with" , nrow(toorder), "unfiltered genes"))

		colnames(toorder)[1] <- names(datalist)[slot_val]
		colnames(toorder)[2] <- names(datalist)[slot_pval]
		colnames(toorder)[3] <- paste("Corrected",names(datalist)[slot_pval],sep='_')
                if (!islog10pvalue) toorder[,2] <- log10(toorder[,2])
                toorder <- toorder[ (!is.na(toorder[,1]))&(!is.na(toorder[,2])),,drop=F]
		if (nrow(toorder) > 0) {
    			if (!"custom_bg" %in% names(gprofile.args)) thisgarg <- c(list(custom_bg=intersect(rownames(toorder), funcfilt)), gprofile.args)
			else thisgarg <- gprofile.args
		
		if (do.fdr.correction){
                        birange <- c(1,nrow(toorder))
			if (birange[2] > 1) {
                        pord <- order(toorder[,2])
			curthreshold = pval.threshold - log10(nrow(toorder))
                        while( birange[1] != birange[2]){
                                test = floor((birange[1] + birange[2])/2)
		#		print(paste(test, "is ", toorder[pord[test],2] - log10(test),  "vs", curthreshold ))
				if (toorder[pord[test],2] - log10(test) < curthreshold) birange[1] = test + 1
				else birange[2] = test
                        }
		#	print(paste(i,"", nrow(toorder), birange[1], "gives", (log10(nrow(toorder)) - log10(birange[1])) ))
		#	print(paste("fdr correction ", (log10(nrow(toorder)) - log10(birange[1]))))
			}
			daffflllt<-toorder[,2] <pval.threshold
			toorder[,3] <- toorder[,3] + (log10(nrow(toorder)) - log10(birange[1]))
			toorder[toorder[,3] >= pval.threshold,3] <- 0
			fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]] <- toorder[order(toorder[,1]),, drop=F]
			toorder <- toorder[ daffflllt,,drop=F]
			

		#	print(dim(toorder))
                }else{
			 fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]] <- toorder[order(toorder[,1]),, drop=F]
			 toorder <- toorder[ toorder[,2] <pval.threshold,,drop=F]
		}

		


		
		print(paste("Function enrichment computed on", sum(toorder[,1]> 0), "upregulated and" , sum(toorder[,1]< 0), "downregulated genes", length(thisgarg$custom_bg) - nrow(toorder) , "bg genes"))

		locflt <- sum(toorder[,1]> 0)
                if (locflt != 0) {
			ord <- order(toorder[,1],decreasing=T)[1:locflt]
			fout[[gsub(" ", "_",paste("TOP", selected.colnames[i],sep="_"))]] <- toorder[ord,]
			if (!is.null(use.funky.threshold)){
				if (use.funky.threshold == "deseq"){
					locflt <- (toorder[,1]> 0)&&(toorder[,2] < -1.3 - (0.3 / ( 0.0625+ toorder[,1]^2)))
				}else{
					locflt <- (toorder[,1]> 0)&&(toorder[,2] < -1.3 - (0.3 / ( 0.0625+ 16 * toorder[,1]^2)))
				}
			}else{locflt <- toorder[,1]> 0}
			ord <- setdiff(order(toorder[,1],decreasing=T), which(!locflt))
			if (length(ord) >= maxlist.toreport) ord <- ord[1:maxlist.toreport]

			if (gprofile.nbterms > 0 ){
			query <- rownames(toorder)[ord]
			if (length(query) > 0){	
			if (!is.null(funcfilt)) query <- query[query %in% funcfilt]

			if (length(query) > 1000) query <- query[1:1000] #do not kill Gprofiler... please

			outputH <- do.call(gprofiler,c(list(query=query) ,thisgarg))
			#outputH <- parseGprofiler2(outputH, bg.size = length(fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]]),query=query,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,])
			outputH <- parseGprofiler2(outputH, bg.size = length(thisgarg$custom_bg),query=query,plot.maxlist=plot.maxlist,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,],make.ggplots=make.ggplots)
			fout[[gsub(" ", "_",paste("TGO", selected.colnames[i],sep="_"))]] <- outputH
			}}
		}
		locflt <- sum(toorder[,1]< 0)
		if (locflt != 0) {
			ord <- order(toorder[,1],decreasing=F)[1:locflt]
			fout[[gsub(" ", "_",paste("BOT", selected.colnames[i],sep="_"))]] <- toorder[ord,]
			if (!is.null(use.funky.threshold)){
				if (use.funky.threshold == "deseq"){
					locflt <- (toorder[,1]< 0)&&(toorder[,2] < -1.3 - (0.3 / ( 0.0625+ toorder[,1]^2)))
				}else{
					locflt <- (toorder[,1]< 0)&&(toorder[,2] < -1.3 - (0.3 / ( 0.0625+ 16 * toorder[,1]^2)))
				}
			}else{locflt <- toorder[,1]< 0}
			ord <- setdiff(order(toorder[,1],decreasing=F), which(!locflt))
			if (length(ord) >= maxlist.toreport) ord <- ord[1:maxlist.toreport]

			if (sum(toorder[,1]< 0) < maxlist.toreport) ord <-order(toorder[,1])[1:(sum(toorder[,1]< 0))]
			else ord <- order(toorder[,1])[1:maxlist.toreport];
	 		if (gprofile.nbterms > 0 ){
			query <- rownames(toorder)[ord]
			if (!is.null(funcfilt)) query <- query[query %in% funcfilt]
			if (length(query) > 0){
			if (length(query) > 1000) query <- query[1:1000] #do not kill Gprofiler... please

			outputH <- do.call(gprofiler,c(list(query=query) ,thisgarg))
			#outputH <- parseGprofiler2(outputH, bg.size = length(fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]]),query=query,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("Low in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,])
			outputH <- parseGprofiler2(outputH, bg.size = length(thisgarg$custom_bg),query=query,plot.maxlist=plot.maxlist,plot.attribs = plotAttribs(plot.attribs,title=paste("Low in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,],make.ggplots=make.ggplots)
			fout[[gsub(" ", "_",paste("BGO", selected.colnames[i],sep="_"))]] <- outputH
			}}
		}
		print(paste("Found", sum((toorder[,1]> 0)&(toorder[,3] <pval.threshold)) , "upregulated and" , sum((toorder[,1]< 0)&(toorder[,3] <pval.threshold)), "downregulated genes after FDR correction"))
		
		}
        }}}
return(fout)}

getMultiExtremes <- function(datalist, consensus, maxlist.toreport= 1000, plot.maxlist=25, value.name= c(), meta.names=c(), pval.grep.string= "pval", pval.threshold=0.05, filter.table=c(), do.fdr.correction=F,gprofile.args=list(hier_filtering="strong",ordered_query=T),gprofile.nbterms=100,plot.attribs=c(),make.ggplots=F,funcfilt=c(), use.funky.threshold=c()){
	
        slot_pval <- grep( pval.grep.string, names(datalist[[1]]))
        islog10pvalue <- (length(grep("log", names(datalist[[1]])[slot_pval])) == 1)
        pval.threshold = log10(pval.threshold)
        if (is.null(value.name)) slot_val <- 1
        else {
		slot_val <- grep( value.name, names(datalist[[1]]))
		if (length(slot_val) == 0) stop(paste("Did not find attribute", value.name))
	}
	fout <- list()
	slot_meta <- rep(1,length(meta.names))
	if (gprofile.nbterms > 0) library(gProfileR)
	for(i in 1:length(meta.names)) slot_meta[i] <- grep( meta.names[i], names(datalist[[1]]))

	
	selected.colnames=c()
	for(i in consensus){
		selected.colnames <- unique(c(selected.colnames, colnames(datalist[[i]][[1]])))
	}

	all_meta <- c(names(datalist[[1]])[c(slot_pval,slot_val)] ,meta.names)
	print(all_meta)
	if (length(selected.colnames) != 0){
        for(i in 1:length(selected.colnames)){
		# averaging time
		conscount <- rep(0,nrow(filter.table[[1]]))
		topass <- data.frame(row.names = rownames(filter.table[[1]]))
		for(k in consensus){
			thatmap <- match(rownames(filter.table[[1]]), rownames(datalist[[k]][[all_meta[1]]]))
	#		if (sum(is.na(thatmap))
			if ((selected.colnames[i] %in% colnames(filter.table[[k]]))&(selected.colnames[i] %in% colnames(datalist[[k]][[all_meta[1]]]))) {
				if (ncol(topass) == 0){
					topass$nb_unfiltered <- rep(0,nrow(topass))
					topass$corr_log10pval <- rep(0,nrow(topass))
					topass[[all_meta[1]]] <- rep(0,nrow(topass))
					flt <- filter.table[[k]][,selected.colnames[i]]
					print(paste("nas2", sum(is.na(flt))))
					flt[is.na(flt)] <- F
					topass$nb_unfiltered[flt] <- 1
					for(j in 2:length(all_meta)){
						topass[[all_meta[j]]] <- rep(0,nrow(topass))
						topass[[all_meta[j]]][flt] <- datalist[[k]][[all_meta[j]]][thatmap,selected.colnames[i]][flt]
					}
					ispos = (datalist[[k]][[all_meta[2]]][thatmap,selected.colnames[i]] > 0)
					isneg = !ispos
					ispos[is.na(ispos)] <- F
					isneg[is.na(isneg)] <- F
					topass$corr_log10pval[flt & ispos] <- datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][flt & ispos]
					topass[[all_meta[1]]][flt & isneg] <- datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][flt & isneg]
	
				}else{
					flt <- filter.table[[k]][,selected.colnames[i]]
					print(paste("nas", sum(is.na(flt))))
					flt[is.na(flt)] <- F
					topass$nb_unfiltered[flt] <- topass$nb_unfiltered[flt] + 1
					for(j in 2:length(all_meta)){
						topass[[all_meta[j]]][flt] <- topass[[all_meta[j]]][flt] + datalist[[k]][[all_meta[j]]][thatmap,selected.colnames[i]][flt]
					}
					ispos = (datalist[[k]][[all_meta[2]]][thatmap,selected.colnames[i]] > 0)
					isneg = !ispos
					ispos[is.na(ispos)] <- F
					isneg[is.na(isneg)] <- F
					topass$corr_log10pval[flt & ispos] <- topass$corr_log10pval[flt & ispos] + datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][flt & ispos]
					topass[[all_meta[1]]][flt & isneg] <- topass[[all_meta[1]]][flt & isneg] + datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][flt & isneg]
				}
			}
		}
		if (ncol(topass) > 0){
		print("Averaging")
		for(j in 4:length(all_meta)){
			topass[[all_meta[j]]] <- topass[[all_meta[j]]] / topass$nb_unfiltered
		}
		print("Gamma-correct")
		for(j in 1:nrow(topass)){
			if (is.na(topass[[all_meta[2]]][j])){
				topass[[all_meta[1]]][j] <- 0
			}else if ((topass[[all_meta[2]]][j]) >0){
				if (topass$nb_unfiltered[j] >1){
					#print(paste("Cons for ",topass$nb_unfiltered[j]))
					#print(topass[[all_meta[1]]][j])
					topass[[all_meta[1]]][j] <- pgamma(-topass$corr_log10pval[j],topass$nb_unfiltered[j],scale=log10(exp(1)),log.p=T,lower.tail=F) / log(10)
					#print(topass[[all_meta[1]]][j])
				}else{topass[[all_meta[1]]][j] <- topass$corr_log10pval[j]}
			}else{
				if (topass$nb_unfiltered[j] >1){
					#print(paste("Cons for ",topass$nb_unfiltered[j]))
					#print(topass[[all_meta[1]]][j])
					topass[[all_meta[1]]][j] <- pgamma(-topass[[all_meta[1]]][j],topass$nb_unfiltered[j],scale=log10(exp(1)),log.p=T,lower.tail=F) / log(10)
					#print(topass[[all_meta[1]]][j])
				}
			}
		}
		toorder <- topass
                toorder <- toorder[ (!is.na(toorder[,4]))&(!is.na(toorder[,3]))&(!is.na(toorder[,1]))&(toorder[,1]>0),,drop=F]
		print(paste("Processing", selected.colnames[i], "with" , nrow(toorder), "unfiltered genes"))
		if (nrow(toorder) > 0) {
    			if (!"custom_bg" %in% names(gprofile.args)) thisgarg <- c(list(custom_bg=intersect(rownames(toorder),funcfilt)), gprofile.args)
			else thisgarg <- gprofile.args
		print("making the correction")

		if (do.fdr.correction){
                        birange <- c(1,nrow(toorder))
			if (birange[2] > 1) {
                        pord <- order(toorder[,3])
			curthreshold = pval.threshold - log10(nrow(toorder))
                        while( birange[1] != birange[2]){
                                test = floor((birange[1] + birange[2])/2)
		#		print(paste(test, "is ", toorder[pord[test],2] - log10(test),  "vs", curthreshold ))
				if (toorder[pord[test],3] - log10(test) < curthreshold) birange[1] = test + 1
				else birange[2] = test
                        }
		#	print(paste(i,"", nrow(toorder), birange[1], "gives", (log10(nrow(toorder)) - log10(birange[1])) ))
		#	print(paste("fdr correction ", (log10(nrow(toorder)) - log10(birange[1]))))
			}
			daffflllt<-toorder[,3] <pval.threshold
			toorder[,2] <- toorder[,3] + (log10(nrow(toorder)) - log10(birange[1]))
			toorder[toorder[,2] >= pval.threshold,2] <- 0
			fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]] <- toorder[order(toorder[,4]),, drop=F]
			toorder <- toorder[ daffflllt,,drop=F]
		#	print(dim(toorder))
                }else{
			 fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]] <- toorder[order(toorder[,4]),, drop=F]
			 toorder <- toorder[ toorder[,3] <pval.threshold,,drop=F]
		}
	#	if (!is.null(funcfilt)) toorder <- toorder[rownames(toorder) %in% funcfilt,,drop=F]

		print("making the extreemes")
		locflt <- sum(toorder[,4]> 0)
                if (locflt != 0) {
			ord <- order(toorder[,4],decreasing=T)[1:locflt]
			fout[[gsub(" ", "_",paste("TOP", selected.colnames[i],sep="_"))]] <- toorder[ord,]
			if (!is.null(use.funky.threshold)){
				if (use.funky.threshold == "deseq"){
					locflt <- (toorder[,4]> 0)&&(toorder[,3] < -1.3 - (0.3 / ( 0.0625+ toorder[,4]^2)))
				}else{
					locflt <- (toorder[,4]> 0)&&(toorder[,3] < -1.3 - (0.3 / ( 0.0625+ 16 * toorder[,4]^2)))
				}
			}else{locflt <- toorder[,4]> 0}
			ord <- setdiff(order(toorder[,4],decreasing=T), which(!locflt))
			if (length(ord) >= maxlist.toreport) ord <- ord[1:maxlist.toreport]

			if (do.fdr.correction) fout[[gsub(" ", "_",paste("TOP", selected.colnames[i],sep="_"))]] <- toorder[ord,]
			else fout[[gsub(" ", "_",paste("TOP", selected.colnames[i],sep="_"))]] <- toorder[ord, setdiff(1:ncol(toorder),2)]
			if (gprofile.nbterms > 0 ){
			query <- rownames(toorder)[ord]
			if (!is.null(funcfilt)) query <- query[query %in% funcfilt]
			if (length(query) > 0){
			if (length(query) > 1000) query <- query[1:1000] #do not kill Gprofiler... please
		#	print("doGpro")
			outputH <- do.call(gprofiler,c(list(query=query) ,thisgarg))
			#outputH <- parseGprofiler2(outputH, bg.size = length(fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]]),query=query,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,])
		#	print("doGproparse")

			outputH <- parseGprofiler2(outputH, bg.size = length(thisgarg$custom_bg),query=query,plot.maxlist=plot.maxlist,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,],make.ggplots=make.ggplots)
			fout[[gsub(" ", "_",paste("TGO", selected.colnames[i],sep="_"))]] <- outputH
			}}
		}
		locflt <- sum(toorder[,4]< 0)
		if (locflt != 0) {
			ord <- order(toorder[,4],decreasing=F)[1:locflt]
			fout[[gsub(" ", "_",paste("BOT", selected.colnames[i],sep="_"))]] <- toorder[ord,]
			if (!is.null(use.funky.threshold)){
				if (use.funky.threshold == "deseq"){
					locflt <- (toorder[,4]< 0)&&(toorder[,3] < -1.3 - (0.3 / ( 0.0625+ toorder[,4]^2)))
				}else{
					locflt <- (toorder[,4]< 0)&&(toorder[,3] < -1.3 - (0.3 / ( 0.0625+ 16 * toorder[,4]^2)))
				}
			}else{locflt <- toorder[,4]< 0}
			ord <- setdiff(order(toorder[,4],decreasing=F), which(!locflt))
			if (length(ord) >= maxlist.toreport) ord <- ord[1:maxlist.toreport]
			
			if (sum(toorder[,4]< 0) < maxlist.toreport) ord <-order(toorder[,4])[1:(sum(toorder[,4]< 0))]
			else ord <- order(toorder[,4])[1:maxlist.toreport];
	 		if (gprofile.nbterms > 0 ){
			query <- rownames(toorder)[ord]
			if (!is.null(funcfilt)) query <- query[query %in% funcfilt]
			if (length(query) > 0){
			if (length(query) > 1000) query <- query[1:1000] #do not kill Gprofiler... please
		#	print("doGpro")

			outputH <- do.call(gprofiler,c(list(query=query) ,thisgarg))
		#	print("doGproparse")

			#outputH <- parseGprofiler2(outputH, bg.size = length(fout[[gsub(" ", "_",paste("ORD", selected.colnames[i],sep="_"))]]),query=query,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,])
			outputH <- parseGprofiler2(outputH, bg.size = length(thisgarg$custom_bg),query=query,plot.maxlist=plot.maxlist,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", selected.colnames[i])),gprofile.nbterms=gprofile.nbterms,gene.annot=toorder[ord,],make.ggplots=make.ggplots)
			fout[[gsub(" ", "_",paste("BGO", selected.colnames[i],sep="_"))]] <- outputH
			}}
		}
		print(paste("Enrichment computed in", sum(toorder[,4]> 0), "upregulated and" , sum(toorder[,4]< 0), "downregulated",length(thisgarg$custom_bg) - nrow(toorder) , "bg genes"))
		print(paste("Found", sum((toorder[,4]> 0)&(toorder[,2]<pval.threshold )), "upregulated and" , sum((toorder[,4]< 0)&(toorder[,2]<pval.threshold )), "downregulated gene after FDR correction"))

        }}}}
return(fout)}

# log in base "e" for pvalue
getIDR <- function(logPvalmatrix, em.steps = 10,em.substep=100, epsilon = 0.001){ 
	print(dim(	logPvalmatrix))
	mu <- rep(-3,ncol(logPvalmatrix))
	sigmas <- rep(1,ncol(logPvalmatrix))	
	p  <- 0.001
	sqprec <- diag(ncol(logPvalmatrix))

	z <- matrix(0, nrow(logPvalmatrix),ncol(logPvalmatrix))
	trz <- matrix(0, nrow(logPvalmatrix),ncol(logPvalmatrix))
	LLtrace <- rep(0, em.steps)	
	for( i in 1:em.steps){
		# resample
		for(j in 1:ncol(logPvalmatrix)){
			#range <- c(min(-3,mu[j] - 3 * sigmas[j] ), max(3,mu[j] + 3 * sigmas[j] ))
			# w <- qnorm()
			w <- seq(min(-6,mu[j] - 6 * sigmas[j] ),max(6,mu[j] + 6 * sigmas[j] ), length=1000)
			wcdf <- log(p * pnorm(w,mu[j],sigmas[j]) + (1.0 -p ) * pnorm(w,0,1))
			for(k in 1:999){
				whi <- which(logPvalmatrix[,j] >= wcdf[k] & logPvalmatrix[,j] < wcdf[k+1])
				z[whi,j] <- w[k] + (logPvalmatrix[whi,j] - wcdf[k]) * (w[k+1]-w[k])/ (wcdf[k+1] - wcdf[k])
			}
			z[(logPvalmatrix[,j] < wcdf[1])  ,j] <- w[1]
			z[(logPvalmatrix[,j] > wcdf[1000])  ,j] <- w[1000]
		}
		print(plot(z[,1], z[,2])) 
	
		# e-step
			
		for(j in 1:em.substep){	
			ll0 <- -0.5 * rowSums(z ^ 2)
			trz <- ((z - rep(mu, each = nrow(logPvalmatrix))) %*% sqprec)
			ll1 <-  determinant(sqprec, logarithm=T)$modulus -0.5 * rowSums(trz ^ 2)
			ez <- 1 / (1 + (1-p) * exp(ll0 - ll1)/p)
		
			tmp <- sum(ll0 + log( (1 - ez) * (1-p) + ez * p * exp(ll1 - ll0)))
		
			if (j > 1){
				if (abs(tmp - LLtrace[i]) < epsilon * (1 + abs(tmp))) break;
			}
	
			LLtrace[i] <- tmp
			print(paste("SubE step ", LLtrace[i]))
			# m-step
			p <- sum(ez)
			for(k in 1:ncol(trz)) trz[,k] <- z[,k] * ez
			mu <- colSums(trz) / p
			Sigma <- matrix(0,2,2)
			for(k in 1:nrow(trz)) Sigma <- Sigma + ((z[k,] - mu) %*% (ez[k] * t(z[k,] - mu)))
			Sigma <- Sigma / p
			dec <- svd(Sigma)
			sqprec <- dec$u %*% ((dec$d ^ -0.5) * dec$v)
			p <- p / nrow(logPvalmatrix)
			trz <- ((z - rep(mu, each = nrow(logPvalmatrix))) %*% sqprec)
			ll1 <-  determinant(sqprec, logarithm=T)$modulus -0.5 * rowSums(trz ^ 2)
			ez <- 1 / (1 + (1-p) * exp(ll0 - ll1)/p)
			print(paste("SubM step ", sum(ll0 + log( (1 - ez) * (1-p) + ez * p * exp(ll1 - ll0)))))
		}
		LLtrace[i] <- tmp
		print(paste("Substep end LL ", LLtrace[i]))
		sigmas = diag(sqprec) ^ -1
	}
	print(LLtrace)
	return( - log(1 + (1-p) * exp(ll0 - ll1)/p) )
}

getMultiExtremesIDR <- function(datalist, consensus, maxlist.toreport= 1000, plot.maxlist=25, value.name= c(), meta.names=c(), pval.grep.string= "pval", pval.threshold=0.05, gene.list=c(), gprofile.args=list(hier_filtering="strong",ordered_query=T),gprofile.nbterms=100,plot.attribs=c(),make.ggplots=F,funcfilt=c()){
        slot_pval <- grep( pval.grep.string, names(datalist[[1]]))
        islog10pvalue <- (length(grep("log", names(datalist[[1]])[slot_pval])) == 1)
        pval.threshold = log10(pval.threshold)
        if (is.null(value.name)) slot_val <- 1
        else {
		slot_val <- grep( value.name, names(datalist[[1]]))
		if (length(slot_val) == 0) stop(paste("Did not find attribute", value.name))
	}
	fout <- list()
	slot_meta <- rep(1,length(meta.names))
	if (gprofile.nbterms > 0) library(gProfileR)
	for(i in 1:length(meta.names)) slot_meta[i] <- grep( meta.names[i], names(datalist[[1]]))

	
	selected.colnames=c()
	for(i in consensus){
		selected.colnames <- unique(c(selected.colnames, colnames(datalist[[i]][[1]])))
	}

	all_meta <- c(names(datalist[[1]])[c(slot_pval,slot_val)] ,meta.names)
	if (length(selected.colnames) != 0){
	        for(i in 1:length(selected.colnames)){
			conscount <- rep(0,length(gene.list))
			topasspos <- data.frame(row.names = gene.list)
			topassneg <- data.frame(row.names = gene.list)
			topass <- data.frame(row.names = gene.list)

			# define common filter...
			for(k in consensus){
				thatmap <- match(gene.list, rownames(datalist[[k]][[all_meta[1]]]))
				if (ncol(topass) == 0){
					topass$corr_log10pval <- rep(0,nrow(topass))
					for(j in 2:length(all_meta)){
						topass[[all_meta[j]]] <- datalist[[k]][[all_meta[j]]][thatmap,selected.colnames[i]]
					}
				}
				ispos = (datalist[[k]][[all_meta[2]]][thatmap,selected.colnames[i]] > 0)
				isneg = !ispos
				isna = is.na(ispos)
				ispos[isna] <- F
				isneg[isna] <- F
				topasspos[[k]][ispos] <- datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][ispos] * log(10)
				topassneg[[k]][isneg] <- datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][isneg] * log(10)
				topasspos[[k]][isna] <- log(0.05 + 0.9 * runif(sum(isna))) 
				topassneg[[k]][isna] <- log(0.05 + 0.9 * runif(sum(isna)))
				topasspos[[k]][isneg] <- log(1.0 - exp(log(10) * datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][isneg]))
				topassneg[[k]][ispos] <- log(1.0 - exp(log(10) * datalist[[k]][[all_meta[1]]][thatmap,selected.colnames[i]][ispos]))
			}	
		
			pos <- getIDR(topasspos)
			return()
			neg <- getIDR(topassneg)



		}

		






	}

}

getExtremesPathfinder <- function(inputlist, celltypenames = c(), filter.table=c(),pval.grep.string=c(),value.name=c(),meta.names=c(),visualize_pathways=F, gmtfile=c(),gene_sets = "BioCarta"){
	library(pathfindR)
	slotpval <- grep(pval.grep.string, names(inputlist))
 	slotfc <- match(value.name, names(inputlist))
	out <- list()
	if (!is.null(gmtfile)) {
		gene_sets = "Custom"
		custom_gene = rownames(inputlist[[slotpval]])
	}
	if (is.null(celltypenames)) {
		if (is.null(filter.table)) celltypenames = colnames(inputlist[[slotfc]])
		else celltypenames =  colnames(filter.table)
	}
	for(i in 1:length(celltypenames)){
		print(paste("Processing", celltypenames[i]))
		input <- cbind(rownames(inputlist[[slotpval]]), inputlist[[slotfc]][,match(celltypenames[i],colnames(inputlist[[slotfc]])),drop=F], 10^inputlist[[slotpval]][,match(celltypenames[i],colnames(inputlist[[slotpval]])),drop=F])
		if (!is.null(filter.table)){
			j <- match(celltypenames[i], colnames(filter.table))
			map <- match(rownames(input), rownames(filter.table))
			print(paste("found filter col", j, "and", sum(is.na(map))))
			input <- input[filter.table[rownames(input) ,j] ,,drop=F]
			if (gene_sets == "Custom"){
				custom_gene = input[,1]
			}
		}
		print(paste("input to pathfinder has",dim(input)[1], "rows."))	

		out[[celltypenames[i]]] <- pathfindR::run_pathfindR(input,output_dir=paste("pathfindR_Results", celltypenames[i], sep="_" ),iterations=4, visualize_pathways=visualize_pathways,custom_gene=custom_gene,custom_pathways =gmtfile,gene_sets=gene_sets)
	}
return(out)
}

getExtremesFast <- function(inputlist, gmtfile, nbperm= 10000,celltypenames = c(), filter.table=c(),pval.grep.string=c(),value.name=c(),meta.names=c(),visualize_pathways=F,gene_sets = "BioCarta",extrafilter=c()){
	library(fgsea)
	slotpval <- grep(pval.grep.string, names(inputlist))
 	slotfc <- match(value.name, names(inputlist))
	out <- list()
	gset <- gmtPathways(gmt.file=gmtfile)
	if (is.null(celltypenames)) {
		if (is.null(filter.table)) celltypenames = colnames(inputlist[[slotfc]])
		else celltypenames =  colnames(filter.table)
	}
	for(i in 1:length(celltypenames)){
		print(paste("Processing", celltypenames[i]))
		if (!is.null(extrafilter)) curf <- (inputlist[[extrafilter]][,celltypenames[i],drop=F] < 0.05)
		else curf <- rep(T, nrow(inputlist[[slotpval]]))
 

		#input <- cbind(rownames(inputlist[[slotpval]]), inputlist[[slotfc]][,match(celltypenames[i],colnames(inputlist[[slotfc]])),drop=F], 10^inputlist[[slotpval]][,match(celltypenames[i],colnames(inputlist[[slotpval]])),drop=F])
		input <- inputlist[[slotfc]][,match(celltypenames[i],colnames(inputlist[[slotfc]]))]
		names(input) <- rownames(inputlist[[slotpval]])
		input <- input[curf]
		if (!is.null(filter.table)){
			j <- match(celltypenames[i], colnames(filter.table))
			map <- match(names(input), rownames(filter.table))
			print(paste("found filter col", j, "and", sum(is.na(map))))
			input <- input[filter.table[names(input) ,j]]
		}
		input <- input[order(input)]

		tmpout <- fgsea(gset,input,nperm= nbperm)

		out[[celltypenames[i]]] <- tmpout[tmpout$padj < 0.05,,drop=F]
	}
return(out)
}

getExtremePathfinder <- function(inputlist, celltypenames = c(), filter.table=c(),visualize_pathways=F){
	library(pathfindR)
	slotpval <- grep("pval", names(inputlist))
 	slotfc <- grep("FC", names(inputlist))
	out <- list()
	if (is.null(celltypenames)) celltypenames = colnames(inputlist[[slotfc]][i])
	for(i in 1:ddlength(celltypenames)){
		print(paste("Processing", celltypenames[i]))
		input <- cbind(rownames(inputlist[[slotpval]]), inputlist[[slotfc]][,match(celltypenames[i],colnames(inputlist[[slotfc]])),drop=F], 10^inputlist[[slotpval]][,match(celltypenames[i],colnames(inputlist[[slotpval]])),drop=F])
		if (!is.null(filter.table)){
			j <- match(celltypenames[i], colnames(filter.table))
			map <- match(rownames(input), rownames(filter.table))
			print(paste("found filter col", j, "and", sum(is.na(map))))
			input <- input[filter.table[rownames(input) ,j] ,,drop=F]
		}
		out[[celltypenames[i]]] <- pathfindR::run_pathfindR(input, output_dir=paste("pathfindR_Results", celltypenames[i], sep="_" ),visualize_pathways=visualize_pathways)
	}
return(out)
}


getExtremeWhich <- function(value=c(), pval=c(), metavalues=c(),  nb.top.bot = c(50,50), usewholeranking=T, pval.threshold= 0.05, is.ordering.in.collumns=T, bypass.rownames=c(), save.as.csv= c(), gprofile.nbterms=25, gprofile.args=list(hier_filtering="strong",ordered_query=T),plot.attribs=c(), label.value="log2FC", label.meta="meanExp", label.pvalue= "Pval", outstr=c(),gene.use=c(),filter.no.annot.inheatmap=F,maxquery=1000, bgfilter =T, do.flip.pvalue=F){
	if (!is.null(outstr)) { # assumes it is a structure outputed my DEcalcIntersection or WilcHyperIntersection
		slot_FC <- grep("log2FC", names(outstr))
		slot_pval <- grep("pval", names(outstr))
		if ("deseq.log2FC" %in% names(outstr)){
			value <- outstr$deseq.log2FC
			value[is.na(value)] <- 0
			pval <- outstr$deseq.pvalue
			metavalues <- log2(outstr$deseq.basemean+1)
			label.pvalue= "log10Pval"
			label.meta <- "log2basemean"
		}else{
			value <- outstr$wilcox.logitAuroc
			value[is.na(value)] <- 0

			pval <- -abs(outstr$wilcox.pvalue -0.5) + 0.5
			metavalues <- log2(outstr$meanTPM+1)
			label.pvalue= "log10Pval"
			label.meta <- "log2basemean"
		}
	}



	if (gprofile.nbterms > 0) library(gProfileR)
	if (length(nb.top.bot) == 1){nbpos = nb.top.bot;nbneg = nb.top.bot
	}else{	nbpos = nb.top.bot[1];	nbneg = nb.top.bot[2]}
	if (is.null(dim(value))){
		if (is.null(gene.use)) gene.use <- rep(T, length(value))
		if (is.null(bypass.rownames)) bypass.rownames <- make.names(names(value), unique=T)
		if (!is.null(gprofile.args)) {
			if (!"custom_bg" %in% names(gprofile.args)) gprofile.args <- c(list(custom_bg=bypass.rownames), gprofile.args)
		}

		fout <- data.frame(row.names = c(paste("Top", 1:nbpos,sep=""), paste("Bot", nbneg:1,sep="")))
		flt <- (pval < pval.threshold) & gene.use
		flt[is.na(flt)] <- F
		ord <- order(value[flt])
		tmp <- length(ord);
		if (tmp < nbpos) wh <- c(ord[tmp:1], rep(tmp+1, nbpos - tmp))
		else wh <- ord[tmp:(tmp-nbpos+1)]
		if (tmp < nbneg) wh <- c(wh, ord[tmp:1], rep(tmp+1, nbneg - tmp))
		else wh <- c(wh, ord[nbneg:1])
		fout[[ "name" ]] <- c(bypass.rownames[flt], NA)[wh]
		fout[[ "value" ]] <- c(value[flt],NA)[wh]
		fout[[ "log10pval" ]] <- log10(c(pval[flt],NA)[wh])

		if (!is.null(metavalues)){
			if (is.null(dim(metavalues))) fout[["metavalue"]] <- c(metavalues[flt], NA)[wh]
			# else TODO
		}
		if (gprofile.nbterms > 0) {
			cdata <- cbind(pval,value,metavalues)
			rownames(cdata) <- bypass.rownames	
			outputH <- do.call(gprofiler,c(list(query=fout$name[1:nbpos]),gprofile.args))
			outputH <- parseGprofiler(outputH, bg.size = length(value),query=fout$name[1:nbpos],do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", colnames(value)[i])),gene.annot=cdata)
			outputL <- do.call(gprofiler,c(list(query=fout$name[(nbpos+nbneg):(nbpos+1)]),gprofile.args))
			outputL <- parseGprofiler(outputL, bg.size = length(value),query=fout$name[(nbpos+nbneg):(nbpos+1)],do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("Low in", colnames(value)[i])), gene.annot=cdata)
		}else return(fout)
		return(list(names=fout,low.annot = outputL, high.annot = outputH))
	}else{

		if (is.null(gene.use)) gene.use <- rep(T, nrow(value))
		fout <- list()
	if (is.ordering.in.collumns){
		if (is.null(bypass.rownames)) bypass.rownames <- make.names(rownames(value),unique=T)
		fout <- list()
		thisgarg <- gprofile.args
		for(i in 1:ncol(value)){
			print(paste("Processing", colnames(value)[i]))

			if (!"custom_bg" %in% names(thisgarg)) {
				if(!bgfilter) thisgarg <- c(list(custom_bg=bypass.rownames[gene.use]), thisgarg)
				else{
					flt <- (pval[,i] < pval.threshold) & gene.use

					if (is.null(dim(metavalues))) {
						tmpthr <- min(metavalues[flt])
						thisgarg <- c(list((custom_bg=bypass.rownames[(metavalues >= tmpthr) &gene.use])), thisgarg)
					}else{
						tmpthr <- min(metavalues[flt,i])
						thisgarg <- c(list((custom_bg=bypass.rownames[(metavalues[,i] >= tmpthr) &gene.use])), thisgarg)
						print(paste("bg has", sum((metavalues[,i] >= tmpthr) &gene.use), "gene from the",sum(gene.use), "conscidered"));	
					}
				}
			}
			slotname <- gsub("-","_",gsub(" ", "_", colnames(value)[i])) 
			flt <- (pval[,i] < pval.threshold) & gene.use
			flt[is.na(flt)] <- F
			ord <- order(value[,i][flt])
			tmp <- length(ord);
			curnbneg <- nbneg
			curnbpos <- nbpos

			if (usewholeranking){
				curlist <- data.frame(row.names = paste("Rank", 1:tmp,sep=""))
				wh <- rev(ord)
				curnbneg <- length(wh)
				curnbpos <- length(wh)
			}else{
				curlist <- data.frame(row.names = c(paste("Top", 1:nbpos,sep=""), paste("Bot", nbneg:1,sep="")))
				if (tmp < nbpos) wh <- c(ord[tmp:1], rep(tmp+1, nbpos - tmp))
				else wh <- ord[tmp:(tmp-nbpos+1)]
				if (tmp < nbneg) wh <- c(wh, ord[tmp:1], rep(tmp+1, nbneg - tmp))
				else wh <- c(wh, ord[nbneg:1])
			}
			curlist[["name"]] <- c(bypass.rownames[flt], NA)[wh]
			curlist[["value"]] <- c(value[flt,i],NA)[wh]
			curlist[["log10pval"]] <- log10( c(pval[flt,i],NA)[wh])
		
			if (!is.null(metavalues)){
				if (is.null(dim(metavalues))) curlist[["metavalue"]] <- c(metavalues[flt], NA)[wh]
				else curlist[["metavalue"]] <- c(metavalues[flt,i], NA)[wh]
			}
	
			if (gprofile.nbterms > 0) {
				cdata <- cbind(log10(pval[,i]),value[,i],metavalues[,i])
				rownames(cdata) <- bypass.rownames
				colnames(cdata) <- c(label.pvalue, label.value, label.meta)
				query = c(bypass.rownames[flt], NA)[wh[1:curnbpos]]
				print(paste("Gprofiler query has", length(query), "genes!"))
				if (length(query) > maxquery) query <- query[1:maxquery]
				outputH <- do.call(gprofiler,c(list(query=query),gprofile.args))

				if (curnbpos > 100) tmphistquery <- query
				else tmphistquery <- c(bypass.rownames[flt], NA)[wh[1:100] ]
				outputH <- parseGprofiler(outputH, bg.size = length(value),query=query,heatmapquery=tmphistquery,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("High in", colnames(value)[i] )),gene.annot=cdata,filter.no.annot=filter.no.annot.inheatmap)
				query = c(bypass.rownames[flt], NA)[wh[length(wh):(length(wh)-curnbneg)] ]	
				if (curnbneg > 100) tmphistquery <- query
				else tmphistquery <- c(bypass.rownames[flt], NA)[wh[length(wh):(length(wh)-curnbneg)] ]

				query = c(bypass.rownames[flt], NA)[wh[length(wh):(length(wh)-curnbneg)] ]
				
				print(paste("Gprofiler query has", length(query), "genes!"))
				if (length(query) > maxquery) query <- query[1:maxquery]

				outputL <- do.call(gprofiler,c(list(query=query),gprofile.args))
				outputL <- parseGprofiler(outputL, bg.size = length(value),query=query,heatmapquery=tmphistquery,do.heatmap=T,plot.attribs = plotAttribs(plot.attribs,title=paste("Low in", colnames(value)[i] )),gene.annot=cdata,filter.no.annot=filter.no.annot.inheatmap)
				fout[[ slotname ]] <- list(list=curlist, low.annot= outputL, high.annot=outputH)
			}else fout[[ slotname ]] <- curlist
		}
	}else{
		if (is.null(bypass.rownames)) bypass.rownames <- make.names(colnames(value),unique=T)
		for(i in 1:nrow(value)){
			curlist <- data.frame(row.names = c(paste("Top", 1:nbpos,sep=""), paste("Bot", nbneg:1,sep="")))
			flt <- (pval[i,] < pval.threshold)
			flt[is.na(flt)] <- F
			ord <- order(value[i,][flt])
			tmp <- length(ord);
			if (tmp < nbpos) wh <- c(ord[tmp:1], rep(tmp+1, nbpos - tmp))
			else wh <- ord[tmp:(tmp-nbpos+1)]
			if (tmp < nbneg) wh <- c(wh, ord[tmp:1], rep(tmp+1, nbneg - tmp))
			else wh <- c(wh, ord[nbneg:1])
			curlist[[ "name" ]] <- c(bypass.rownames[flt], NA)[wh]
			curlist[[ "value" ]] <- c(value[i,flt],NA)[wh]
			curlist[[ "pval" ]] <- c(pval[i,flt],NA)[wh]

			if (!is.null(metavalues)){
				if (is.null(dim(metavalues))) curlist[["metavalue"]] <- c(metavalues[flt], NA)[wh]
				# else TODO
			}
			if (gprofile.nbterms > 0) {
				cdata <- t(rbind(pval[i,],value[i,],metavalues[i,]))
				rownames(cdata) <- bypass.rownames
				outputL <- do.call(gprofiler,c(list(query=curlist$name[1:nbpos]),gprofile.args))
				outputL <- parseGprofiler(outputL, bg.size = length(value),query=curlist$name[1:nbpos],do.heatmap=T,  plot.attribs = plotAttribs(plot.attribs,title=paste("Low in", colnames(value)[i] )))
				outputH <- do.call(gprofiler,c(list(query=curlist$name[(nbpos+nbneg):(nbpos+1)]),gprofile.args))
				outputH <- parseGprofiler(outputH, bg.size = length(value),query=curlist$name[(nbpos+nbneg):(nbpos+1)],do.heatmap=T, plot.attribs = plotAttribs(plot.attribs,title=paste("High in", colnames(value)[i] )))

				fout[[ colnames(value)[i] ]] <- list(list=curlist, low.annot= outputL, high.annot=outputH)
			}else fout[[ colnames(value)[i] ]] <- curlist
		}
	}
	}
return(fout)} 
saveAsCSV <-function(object){
	if (class(object) == "list"){
		

	}

}


# posneg is a vector of T/F
plotAUROCcomparison <- function(data, posneg, cluster, do.plot=T,plot.zscore.denum=1,do.return.fig=F){
	library(ggplot2)
	d <- dim(data)
	clsize <- table(sort(cluster))
	psize <- clsize
	ssize <- rep(1,length(clsize)+1)
	for(j in 1:length(clsize)){
		psize[j] <- sum(cluster[posneg] == names(clsize)[j])
		ssize[j+1] = ssize[j] + clsize[j]
	}
	fout <- list();
	fout$auroc <- matrix(NA,d[1],length(clsize))
	fout$Zscore <- matrix(NA,d[1],length(clsize))
	for(i in 1:d[1]){
		currow <- data[i,]
		ord <- order(cluster,currow)
		for(j in 1:length(clsize)){
			auroc=0; last= c(0,0);
			tie =0;
			if ((psize[j] == clsize[j])||(psize[j] == 0)){
				fout$Zscore[i,j] = 0
				fout$auroc[i,j] = 0.5 # no curve no data for one or two class
			}else { # at least 2 data point because of above
				if (posneg[ord[ssize[j]]]) cur <- c(1,0)
				else cur <- c(0,1)
				lastval <- currow[ord[ssize[j]]]
				for(k in ord[(ssize[j]+1):(ssize[j+1]-1)]){
					if (lastval != currow[k]){
						auroc = auroc + (cur[1] - last[1]) * (last[2] + cur[2]);
						tmp = cur[1] + cur[2] - last[1] - last[2]
						tie <- tie + tmp * (tmp * tmp - 1.0)
						lastval = currow[k];
						last <- cur
					}
					if (posneg[k]) cur[1] <- cur[1] +1
					else cur[2] <- cur[2] +1
				}
				auroc = auroc + (cur[1] - last[1]) * (last[2] + cur[2]);
				tmp = cur[1] + cur[2] - last[1] - last[2]
				tie <- tie + tmp * (tmp * tmp - 1.0)
				var = cur[1] * cur[2] * ((cur[1] +cur[2] + 1) - tie / ((cur[1] +cur[2]) * (cur[1] +cur[2]-1)))/3
				if (var == 0){ # all equal... meaningless test
					fout$Zscore[i,j] = 0
					fout$auroc[i,j] = 0.5 # no curve no data for one or two class
				}else{
					fout$auroc[i,j] <- auroc / (2.0 * cur[1] *cur[2])
					fout$Zscore[i,j] <- (auroc - cur[1] * cur[2]) / sqrt(cur[1] * cur[2] * ((cur[1] +cur[2] + 1) - tie / ((cur[1] +cur[2]) * (cur[1] +cur[2]-1)))/3)
				}
			}
		}
	}
	rownames(fout$auroc) <- rownames(data)
	colnames(fout$auroc) <- names(clsize)
	rownames(fout$Zscore) <- rownames(data)
	colnames(fout$Zscore) <- names(clsize)

	if (do.plot){
	#	library(pheatmap)
	#	if (do.estimate.Zscore){
	#		library(pracma)
	#		erftr <- matrix(sapply(fout, function(x){return(0.5 + 0.5 * erf(x * sqrt(0.5)))}),d[1],length(clsize))
	#		rownames(erftr) <- rownames(data)
	#		colnames(erftr) <- names(clsize)
	#		pheatmap( erftr, color= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41), breaks=seq(0,1,by=1/41) )	
	#	}else{

	#		pheatmap( fout$auroc, color= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41), breaks=seq(0,1,by=1/41) )
	#	}
		aurange <- range(as.vector(fout$auroc));
		if (1 - aurange[2] < aurange[1]) aurange <- c(1 - aurange[2], aurange[2])
		else aurange <- c(aurange[1], 1 - aurange[1])

		curd <- dim(fout$auroc)
		ggdata <- as.data.frame(cbind( rep(1:curd[1], curd[2]), rep(1:curd[2],each=curd[1]),as.vector(fout$auroc) , sapply(as.vector(fout$Zscore), function(x){return(1.0 / ( 1.0 + ((plot.zscore.denum / x)^2)))})))
		colnames(ggdata) <- c("Y", "X", "Auroc", "Significance")
		ggdata$X <- as.factor(ggdata$X)
		ggdata$Y <- as.factor(ggdata$Y)
		p <- ggplot(ggdata, aes(X,Y,color=Auroc, size=Significance)) + theme(axis.text.x=element_text(angle=90,vjust=0.5)) + scale_x_discrete(breaks= 1:curd[2], labels=names(clsize) ) + xlab(NULL) + scale_y_discrete(breaks= 1:curd[1], labels=rownames(data) ) + ylab(NULL) + theme(axis.text=element_text(color="#000000",face="bold",size=10)) + geom_point() + scale_colour_gradientn(colours=c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"),limits=aurange) + scale_size(limits=c(0,1))
		 
	if (do.return.fig) return(p)
	else print(p)
	}
return(fout)
}
showMultiAUROC <- function(data, posneg, partition, ignorezero=F){
	library(ggplot2)
	if (class(partition) == "factor") plvl <- levels(partition)
	else plvl <- unique(partition);
	for(i in 1:length(plvl)){
		steps <- table(data[partition == plvl[i]]);
		steps <- steps[order(names(steps))]
		psteps <- table(data[(partition == plvl[i])&(posneg)]);
		gdata <- matrix(0,length(steps)+1,3);
		nsize <- sum(!posneg[partition == plvl[i]]);
		psize <- sum(posneg[partition == plvl[i]]);
		cur <- c(0,0)
		gdata[1,3] = i
		for(j in 1:length(steps)){
			hehe <- match(names(steps)[j], names(psteps))
			if (!is.na(hehe)){
				cur[2] <- cur[2] + psteps[hehe];	
				cur[1] <- cur[1] + steps[j] - psteps[hehe];
			}else{
				cur[1] <- cur[1] + steps[j];
			}
			gdata[(j+1),1] <- cur[1] / nsize;
			gdata[(j+1),2] <- cur[2] / psize;
			gdata[(j+1),3] <- i
		}
		if (i == 1) adata <- gdata
		else adata <- rbind(adata,gdata)
	}
	colnames(adata) <- c("X", "Y", "Z")
	p <- ggplot(data.frame(adata),mapping = aes(x = X, y = Y, group=Z)) + geom_line(size=.75)
return(p);
}

myggplotSeurat <- function(sro,pt.shape=NULL,pt.size=1){# , shape = factor(x = pt.shape)
		id <- cbind(sro@dr$tsne@cell.embeddings, as.numeric(sro@meta.data$res.0.8))
		colnames(id) <- c("x","y","i")
		pl <- ggplot(data.frame(id),mapping = aes(x = x, y = y) ) + geom_point(mapping = aes(colour = factor(x = i)), size = pt.size)
return(pl)}


#    if (do.label) {
#        centers <- data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(x = x), 
#            y = median(x = y))
#        p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
#            y = y), size = 0, alpha = 0) + geom_text(data = centers, 
#            mapping = aes(label = ident), size = label.size)
#    }

#' Function that return argument, really pointless
#'
#' @param some number
#'
#' @export
#setGeneric("show", function(object,arg){return("well")})

InferN0Heatmap <- function(data, permute,show_rownames=F,show_colnames=F,cluster_cols=T, ...){
	library(pheatmap)
	pheatmap(data[permute,]
		,color= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)
		,show_rownames=show_rownames,show_colnames=show_colnames,
		,cluster_rows=F,cluster_cols=cluster_cols,
		,...
		)
}

PlotCosine <- function(data, genename){
	pos = match(genename,rownames(data))
	cosval <- matrix(0,nrow(data),1)
	rownames(cosval) <- rownames(data)
	colnames(cosval) <- c("Cosine")

	thattranspose <- t(data)
	thatrow <- data[pos,]
	for(i in 1:nrow(data)){
		cosval[i] <- (thatrow %*% thattranspose[,i])[1,1]
		if ((i %% 100) == 0) print(i)
	}
	cosval = cosval / cosval[pos]
	cosval[pos] <- 0
	return(cosval)
}

QueryOverlay <- function(sro, genename, genename2, doheatmap=F){
#	library(gridExtra)

	if (doheatmap){
		rowdata <- sro@raw.data[genename,]
		rowdata2 <- sro@raw.data[genename2,]
		dim1 = log2(max(sro@raw.data[genename,])+1)+1
		dim2 = log2(max(sro@raw.data[genename2,])+1)+1
		counts <- matrix(0,dim1,dim2)

		for(i in 1:length(rowdata)){
			coor <- c(log2(rowdata[i]+1)+1,log2(rowdata2[i]+1)+1)
			if ((coor[1] != 1)&&(coor[2] !=1)) {
				counts[coor[1],coor[2]] = counts[coor[1],coor[2]] + 1
			}
		}
		print(counts)
		pheatmap(counts, color=colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(17)
		,cluster_rows=F,cluster_cols=F,labels_row=0:dim1,labels_col=0:dim2)
	}else{
		getOption("device")()
				rowdata <- sro@data[genename,]
			rowdata2 <- sro@data[genename2,]
		par(mfrow = c(1,2))
		#rowdata <- sro@dr$tsne@cell.embeddings[,1]
		rowdata[rowdata ==0] <- NA

		mmm <- Infern0Colors(rowdata)
		p1 <- plot(sro@dr$tsne@cell.embeddings,col=mmm$col, main=genename)

		#rowdata <- sro@dr$tsne@cell.embeddings[,1]
		rowdata2[rowdata2 ==0] <- NA
		mmm2 <- Infern0Colors(rowdata2)
		p2 <- plot(sro@dr$tsne@cell.embeddings,col=mmm2$col, main=genename2)
	}
	#plot(sro@raw.data[genename,], sro@raw.data[genename2,], main="plot 1",cluster_rows=F,cluster_cols=F)
	#do.call(grid.arrange,list(p1,p2,ph[[4]]),layout_matrix=rbind(c(1,1,1),c(1,2,3)))
}


# predecated?
makeSeuratObjectFromCellRanger <- function(path, mito_grep = "", ribo_grep = "", min_genes = 200, min_cell = 3) {
    # ex: mito_grep = "^MRP[SL][0123456789]\|^MT-" ribo_grep = "^RP[SL][0123456789]"
    library("Seurat")
    data <- Read10X(path);
    srt <- CreateSeuratObject(raw.data = data,min.cell = min_cell, min.genes = min_genes)
    if (nchar(mito_grep) != 0){
        mito.genes <- grep(pattern = mito_grep, x = rownames(x= srt@data), value = TRUE)
        percent.mito <- Matrix::colSums(srt@data[mito.genes,]) / Matrix::colSums(srt@data)
        srt <- AddMetaData(object = srt, metadata = percent.mito,col.name = "percent.mito")
    }
    if (nchar(ribo_grep) != 0){
        names <- grep(pattern = ribo_grep, x = rownames(x= srt@data), value = TRUE)
        percent <- Matrix::colSums(srt@data[names,]) / Matrix::colSums(srt@data)
        srt <- AddMetaData(object = srt, metadata = percent,col.name = "percent.ribo")
    }
    return(srt)
}



#' Construct SingleCellExperiement from Cellranger data path
#'
#' @param directory path to a cellranger output folder, which expected that contains 3 files named barcodes.tsv genes.tsv matrix.mtx respectively. (ex: "./data/")
#'
#' @export
makeSingleCellExperimentFromCellRanger <- function(path) {
    library("SingleCellExperiment")
    data <- readRangerMatrix(path)
    sce <- SingleCellExperiment(assays = list(counts = data))
    rowData(sce)$feature_symbol <- rownames(data)
    colData(sce)$feature_symbol <- colnames(data)
    return(sce)
}



#' Construct SingleCellExperiement from Cellranger data path
#'
#' @param directory path to a cellranger output folder, which expected that contains 3 files named barcodes.tsv genes.tsv matrix.mtx respectively. (ex: "./data/")
#'
#' @export
makeSingleCellExperimentFromSeurat <- function(sro, celllist =c(), genelist=c(), fetch.dr=T, fetch.metadata=T, silly.matrix.convert=F) {
    library("SingleCellExperiment");
  	if (length(celllist) == 0) celllist <- rep(T,SeuratRawData(sro,"D")[2])
  	if (length(genelist) == 0) genelist <- rep(T,SeuratRawData(sro,"D")[1])	
	else if (length(genelist) != SeuratRawData(sro,"D")[1]) {
		genelist <- !is.na(match(SeuratRawData(sro,"R"), genelist))
	}
    		sce <- SingleCellExperiment(assays = list(counts = SeuratRawData(sro)[genelist,celllist,drop=F], logcounts= SeuratRawData(sro)[genelist,celllist,drop=F]));
		if (class(logcounts(sce)) == "dgCMatrix") {
			logcounts(sce)@x <- log2(logcounts(sce)@x + 1);
			if (silly.matrix.convert) logcounts(sce) <- as.matrix(logcounts(sce))
		}else logcounts(sce) <- log2(logcounts(sce) + 1);
		
		rowData(sce)$feature_symbol <- SeuratRawData(sro,"R")[genelist]
		colData(sce)$feature_symbol <- SeuratRawData(sro,"C")[celllist]
		if (fetch.metadata){
			for(i in colnames(sro@meta.data)){
				colData(sce)[[i]] <- sro@meta.data[celllist,c(i)]
			}
		}
		if (fetch.dr){
			if ("tsne" %in% names(sro@dr)){
				reducedDim(sce, "tSNE") <- sro@dr$tsne@cell.embeddings[celllist,,drop=F] 
			}
			if ("pca" %in% names(sro@dr)){
				reducedDim(sce, "PCA") <- sro@dr$pca@cell.embeddings[celllist,,drop=F]
			}
		}
return(sce)}
#' Construct Seurat from SingleCellExperiement object
#'
#' @param sce single cell experiment object
#'
#' @export
makeSeuratFromSingleCellExperiment <- function(sce, celllist =c(), assay.raw.data.name="counts", makeSparse=T) {
    library("SingleCellExperiment");
    library("Seurat");
    library("Matrix");
  	if (length(celllist) == 0){
		rawcounts <- Matrix(assays(sce)[[assay.raw.data.name]], sparse=makeSparse)
		if (is.null(rownames(rawcounts))){
			rownames(rawcounts) <- rowData(sce)$feature_symbol
		}else{
			rownames(rawcounts) <- rownames(sce)
		}
		if (is.null(colnames(rawcounts))){
			colnames(rawcounts) <- colData(sce)$feature_symbol
		}else{
			colnames(rawcounts) <- colnames(sce)
		}
	}
	sro <- CreateSeuratObject( rawcounts )
	toadd <-data.frame(row.names = colnames(sce))
	for(i in names(colData(sce))){
		toadd[[i]] <- colData(sce)[[i]]
	}
	if (ncol(toadd) > 0) sro <- AddMetaData(sro, toadd)

	for( drn in reducedDimNames(sce)){
		tmp <- reducedDim(sce,drn)
		colnames(tmp) <- paste(drn,1:ncol(tmp),sep="_")
		rownames(tmp) <- rownames(sro@meta.data)
		sro <- SetDimReduction(sro, reduction.type = drn, slot = "cell.embeddings", new.data = tmp)
	        sro <- SetDimReduction(sro, reduction.type = drn, slot = "key", new.data = drn)	
	}
    return(sro)
}

#sce <- prepareSCMapReference(srop); saveRDS(sce, "MDS.sce."
# run all markers, and save in misc slot
DoFindAllMarkers <- function(sro, do.not.use.seurat=T){
	if (do.not.use.seurat){
		library(InferN0)
		scptr <- InferN0Init(t(sro@raw.data))
		output <- InferN0FindMarkers(scptr,nbpart=4,ordering=order(sro@meta.data$nUMI)-1,partition= (sro@ident@.Data - 1))
		sro@misc$mymarkers <- output
	}else{
		tmp <- FindAllMarkers(sro);
		sro@misc$allmarkers <- tmp
	}
return(sro)}

makeSCFindIndex <- function(sce_or_sro, clustername){
	require(scfind)
	require(SingleCellExperiment)
	if (class(sce_or_sro) == "seurat") sce_or_sro <- makeSingleCellExperimentFromSeurat(sce_or_sro)
	colData(sce_or_sro)$cell_type1 <- colData(sce_or_sro)[[clustername]]
	return(buildCellTypeIndex(sce_or_sro, "Cell_type"))
}
scfindroutine <- function(t){
	require(ggplot2)
	#Find the marker genes for the T cells in the thymus
	d <- cellTypeMarkers(t, "Thymus.T cell", top.k=5)
	a <- t@index$findCellTypes(as.character(d$genes), "Thymus")
	ns <- paste0("Thymus.", rownames(t@metadata$Thymus.umap))
	e <- rep(F, length(ns))
	inds <- which(ns==names(a)[2])
	e[inds[a[[2]]]] <- T
	data <- data.frame(type=ns, x=t@metadata$Thymus.umap[,1], y=t@metadata$Thymus.umap[,2], exp=e)
	ggplot(data, aes(x=x, color=e, y=y, shape=type)) + geom_point() + theme_minimal() + xlab("UMAP1") + ylab("UMAP2") + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
}



#sce <- prepareSCMapReference(srop); saveRDS(sce, "MDS.sce."

# run all markers, and save in misc slot
MyFindAllMarkers <- function(sro){
return(sro)}
getMappableList <- function(sro, soup.threshold=1){
	if (!("allmarkers" %in% names(sro@misc))) sro <- DoGetAllMarkers(sro)
	damap <- match(sro@misc$allmarkers[[7]], rownames(sro@data))
	soupcol <- grep("soupfraction", colnames(sro@misc$genemetadata))
	print(soupcol)
	plot(rowSums(sro@misc$genemetadata[damap,soupcol]) / length(soupcol), log10(sro@misc$allmarkers[[5]]))
}
countPairOccurence <- function(datafr, colA,colB){

	rows <- unique(datafr[,c(colA)])
	cols <- unique(datafr[,c(colB)])
	fout <- matrix(0,length(rows),length(cols))
	rownames(fout) <- rows
	colnames(fout) <- cols
	for(i in 1:length(rows)){
		subr <- datafr[ (datafr[,c(colA)] == rows[i])  , c(colB)]
		for(j in 1:length(cols)) fout[i,j] <- sum(subr == cols[j])
	}
return(fout)}

scmapindexClusterVariant <- function(object, cluster_col) {
    if (!cluster_col %in% colnames(colData(object))) {
        stop("Please define an existing cluster column of the `colData` slot of the input object using the `cluster_col` parameter!")
        return(object)
    }
    tmp <- object[rowData(object)$scmap_features, ]
    gene <- cell_class <- exprs <- NULL
    exprs_mat <- as.matrix(logcounts(tmp))
#    rownames(exprs_mat) <- as.datdda.frame(rowData(tmp))$feature_symbol
#    colnames(exprs_mat) <- as.data.frame(colData(tmp))[[cluster_col]]
    
    # calculate median feature expression in every cell class of object
 #   exprs_mat <- reshape2::melt(exprs_mat)
#
  #  colnames(exprs_mat) <- c("gene", "cell_class", "exprs")

   # exprs_mat <- exprs_mat %>% group_by(gene, cell_class) %>% summarise(med_exprs = median(exprs))
   # exprs_mat <- reshape2::dcast(exprs_mat, gene ~ cell_class, value.var = "med_exprs")
   # rownames(exprs_mat) <- exprs_mat$gene
   # index <- exprs_mat[, 2:ncol(exprs_mat), drop = FALSE]
    label <- as.data.frame(colData(tmp))[[cluster_col]]
    if (class(label) == "factor") lvls <- levels(label)
    else lvls <- unique(lvls)
    nbcells <- rep(0,length(lvls))
    index <- matrix(0, sum(rowData(object)$scmap_features)  ,length(lvls))

    colnames(index) <- lvls
    rownames(index) <- sort(rowData(object)$feature_symbol[rowData(object)$scmap_features])
    geneorder <- match(rowData(object)$feature_symbol[rowData(object)$scmap_features], rownames(index))
    for(i in 1:length(lvls)){
	nbcells[i] <- sum(label == lvls[i])
	if (nbcells[i] != 0) {
	subset <- exprs_mat[, label == lvls[i]]
       	qr <- c(1 + floor(nbcells[i]/20), 1 + floor(nbcells[i]/2), 1 + floor((19*nbcells[i]-1)/20)) # middle 80% range
#	for(j in 1:sum(rowData(object)$scmap_features)) index[geneorder[j],i] <- sort(subset[j,])[tmptmp] # computes median
	for(j in 1:sum(rowData(object)$scmap_features)) index[geneorder[j],i] <- mean(sort(subset[j,])) # computed mean (original scmap... but does not crash ^^')
#	for(j in 1:sum(rowData(object)$scmap_features)) index[geneorder[j],i] <- mean(sort(subset[j,])[qr[1]:qr[3]]) # computes average of middle 50% range
#	for(j in 1:sum(rowData(object)$scmap_features)) index[geneorder[j],i] <- mean(sort(subset[j,])[qr]) # computes average of middle 50% range

        }
    }   

    metadata(object)$scmap_cluster_index <- as.data.frame(index[, (nbcells != 0)])
    return(object)
}

prepareSCMapReference <- function(sce_or_sro, clname="", nbfeature=500, selected.features= c(),cell.use=c(), cell.meta=c(), output.index.only=F, use_scmap_variant=F, do.indexcells.k=0){
	library(scmap)
	if ((class(sce_or_sro) == "seurat")||(class(sce_or_sro) == "Seurat")){
		cell.use <- SeuratCellListQuery(sce_or_sro, cell.meta,cell.use)	
		if (sum(cell.use) == 0) stop("Got an empty cell list")
		if (clname == "") clname = "res.0.8"
		if (length(selected.features) == 0){
			sce <- makeSingleCellExperimentFromSeurat(sce_or_sro, celllist=cell.use,fetch.dr=F);
			sce <- selectFeatures(sce, nbfeature)
		}else{

			sce <- makeSingleCellExperimentFromSeurat(sce_or_sro, celllist=cell.use, genelist=!is.na(match(SeuratRawData(sce_or_sro, "R"),selected.features)),fetch.dr=F);
			sce <- setFeatures(sce, rownames(sce))
		}
	}else if (class(sce_or_sro) == "SingleCellExperiment"){
		sce <- selectFeatures(sce_or_sro, nbfeature)
	}else stop("expected a seurat object or single cell object")
	print("indextipe")
	if (use_scmap_variant) sce <- scmapindexClusterVariant(sce, clname)
	else sce <- indexCluster(sce, clname)
	print("indexdo")
	if (do.indexcells.k != 0) {
		 print("making Cell to Cell index")
		 sce <- indexCell(sce, k=do.indexcells.k)
	}

	if (output.index.only) {
		print("ga")
		if (do.indexcells.k != 0) sce <- list(index = metadata(sce)$scmap_cluster_index, scindex = metadata(sce)$scmap_cell_index)
		else sce <- list(index = metadata(sce)$scmap_cluster_index)
		if (class(sce_or_sro) == "seurat"){
			if (class(sce_or_sro@meta.data[[clname]]) == "factor"){ #permuting 
				map <- match(levels(sce_or_sro@meta.data[[clname]]), colnames(sce$index))
				map <- map[!is.na(map)]
				map <- c(map, which(is.na(match(colnames(sce$index),levels(sce_or_sro@meta.data[[clname]])))))
				sce$index <- sce$index[,map]
			}
			if ("meta.color" %in% names(sce_or_sro@misc)){
				if (clname %in% names(sce_or_sro@misc$meta.color)){
					sce <- c(sce, list(colors=  sce_or_sro@misc$meta.color[[clname]]))
					if (is.null(names(sce_or_sro@misc$meta.color[[clname]]))) names(sce$colors) <- levels(sce_or_sro@meta.data[[clname]])[1:length(sce$colors)]
				}
			}
			if (do.indexcells.k != 0) sce$scindex$metadata <- sce_or_sro@meta.data[colnames(sce$scindex$subclusters),clname]
		}
		return(sce)
	}
return(sce)}
# sro<- readRDS("neuroep2.sro.rds"); sce<- readRDS("MDS.sce.rds"); source("~/work/readFromNew.R");  runSCMapFromSeurat(sro,sce)

# sce is Single cell Experiment *or* the index only (latter works for cell.NN == 0 only)
runSCMapFromSeurat <-  function(sro, sce, clustername="scmap",celllist =c(), cluster.threshold= 0, cell.NN=0, clname="", rename.cluster="",rename.do.reorder=T, color.use=c(), make.names.unique=F, return.table.instead=F) {
	library(scmap)
	output = list();
	if (class(sce) == "list"){
		genelist = rownames(sce$index)
	}else{
		genelist = rownames(sce)
	}
	print(genelist)

	if (class(sro) != "SingleCellExperiment") tmpsce <- makeSingleCellExperimentFromSeurat(sro, celllist,genelist=genelist,silly.matrix.convert=T,fetch.dr=F);
	if (make.names.unique) colData(tmpsce)$feature_symbol <- make.names(colData(tmpsce)$feature_symbol, unique=T)
	if (class(sce) == "list"){
		if ((clustername == "scmap")&("alias" %in% names(sce))) clustername <- sce$alias
		if (!"meta.color" %in% names(sro@misc) ) sro@misc$meta.color <- list()
		classnames <- colnames(sce$index)
		if (make.names.unique) rownames(sce$index) <- make.names(rownames(sce$index),unique=T)
		#subind <- sce$index[!is.na(match(rownames(sce$index), rownames(sro@raw.data))),]
		print("Mapping cell to clusters")
		output$cluster <- scmapCluster(tmpsce, list(sce$index),threshold=cluster.threshold)
	}else if (class(sce) == "data.frame") {
		cell.NN = 0
		classnames <- colnames(sce)
		if (make.names.unique) rownames(sce) <- make.names(rownames(sce),unique=T)
		output$cluster <- scmapCluster(tmpsce, list(sce),threshold=cluster.threshold)
	}else{
		if (make.names.unique) rownames(metadata(sce)$scmap_cluster_index) <- make.names(rownames(metadata(sce)$scmap_cluster_index),unique=T)
		classnames <-  colnames(metadata(sce)$scmap_cluster_index)
		output$cluster <- scmapCluster(tmpsce, list(metadata(sce)$scmap_cluster_index),threshold=cluster.threshold)
	}
	print("done")
	annot <- data.frame(output$cluster$scmap_cluster_labs, output$cluster$scmap_cluster_siml)
	if (sum(annot[[1]] == "unassigned") > 0) classnames <- c(classnames,"unassigned")
	annot[[1]] <- factor(annot[[1]], levels=classnames)
	rownames(annot) <- rownames(sro@meta.data)
	colnames(annot) <- c(paste(clustername,c("cluster", "simil"), sep='_'))
	if (return.table.instead) return(annot)	

	library(Seurat)
	output$sro <- AddMetaData(sro,annot)
	if (!"meta.color" %in% names(output$sro@misc) ) print("wtf")
	if (cell.NN > 0){
		print(class(sce))
		if (class(sce) == "list"){
			print("Running cell-cell clustering")
			tmptmp <- scmapCell(tmpsce, list(sce$scindex),w=cell.NN)
			outout <- rep(0, ncol(tmptmp[[1]]$cells))
			for(i in 1:length(outout)) {
				outout[i] <- sum(sce$scindex$metadata[tmptmp[[1]]$cells[,i]] * tmptmp[[1]]$similarities[,i])
			}
			outout <- outout / colSums(tmptmp[[1]]$similarities)
			return(list(a=tmptmp, b=outout))
		}else{
			tmptmp <- scmapCell(tmpsce, list(metadata(sce)$scmap_cell_index),w=cell.NN)
			output$cell <- matrix(NA,dim(tmptmp[[1]]$cell)[1],dim(tmptmp[[1]]$cell)[2])
			for(i in 1:dim(tmptmp[[1]]$cell)[2]) output$cell[,i] <- colData(sce)[[clname]][tmptmp[[1]]$cell[,i]]
			annot2 <- data.frame(output$cell[1,])
			rownames(annot2) <- rownames(sro@meta.data)
			makeSingleCellExperimentFromSeuratcolnames(annot2) <- c(paste(clustername,"bestcell", sep='_'))
			output$sro <- AddMetaData(output$sro,annot2)
			for(i in 1:dim(tmptmp[[1]]$cell)[2]) annot2[i,1] = names(sort(table(output$cell[,i]),decreasing=T))[1]
			colnames(annot2) <- c(paste(clustername,"majoritycell", sep='_'))
			output$sro <- AddMetaData(output$sro,annot2)
		}
	}
	if (rename.cluster != "") {
		output$sro <- DoRenameClusters(output$sro,rename.cluster, colnames(annot)[1],colnames(annot)[2],rename.do.reorder)
		#if (class(sce) == "list"){
		#			#	output$sro@misc$meta.color[[rename.cluster]] <- mydoublerainbow(length(levels(output$sro@meta.data[[rename.cluster]])))
		#	names(output$sro@misc$meta.color[[rename.cluster]]) <- levels(output$sro@meta.data[[rename.cluster]])
		#}
	}
	if (class(sce) == "list"){
		fout = c()
		for(i in 1:length(sce$colors)) fout[[classnames[i]]] <- sce$colors[classnames[i]]
		fout[["unassigned"]] <- "#FFFFFF88"
		output$sro@misc$meta.color[[colnames(annot)[1]]] <- fout
	}else if (!is.null(color.use)){
		fout = c()
		for(i in 1:length(color.use)) fout[[classnames[i]]] <- color.use[i]
		fout[["unassigned"]] <- "#FFFFFF88"
		output$sro@misc$meta.color[[colnames(annot)[1]]] <- fout
	}
return(output)}

DoTransferMeta<-function(sro.target, sro.source, target.metaname, source.metaname){
	map <- match(rownames(sro.target@meta.data), rownames(sro.source@meta.data))
	if (class(sro.source@meta.data[[source.metaname]]) == "factor")	dalvls <- levels(sro.source@meta.data[[source.metaname]])
	else dalvls <- unique(sro.source@meta.data[[source.metaname]])
	data <- matrix("Missing", dim(sro.target@meta.data)[1],1)
	if (sum(is.na(map)) != 0) dalvls <- c(dalvls, "Missing")
	data[!is.na(map) ,1] <- as.character(sro.source@meta.data[[source.metaname]])[map[!is.na(map)]]
	sro.target@meta.data[[target.metaname]] <- factor(data, levels = dalvls)
	if ("meta.color" %in% names(sro.source@misc)){
		if (sum(is.na(map)) != 0) sro.target@misc$meta.color[[target.metaname]] <- c(sro.source@misc$meta.color[[source.metaname]], "#FFFFFF88")
		else sro.target@misc$meta.color[[target.metaname]] <- sro.source@misc$meta.color[[source.metaname]]
	}
return(sro.target)}

DoRenameClusters <- function(sro, meta.to.rename, meta.mapnames, meta.mapweight=c(), puritythreshold=0.5,do.reorder=T,do.merge=F){
	if (class(sro@meta.data[[meta.to.rename]]) != "factor") sro@meta.data[[meta.to.rename]] <- as.factor(sro@meta.data[[meta.to.rename]])
	vlabs = levels(sro@meta.data[[meta.to.rename]])
	if (class(sro@meta.data[[meta.mapnames]])  == "factor") hlabs = levels(sro@meta.data[[meta.mapnames]])
	else hlabs = unique(sro@meta.data[[meta.mapnames]])
	nbcell = dim(sro@meta.data)[1]
	counts <- matrix(0,length(hlabs), length(vlabs))
	if (!is.null(meta.mapweight)){
		for(i in 1:nbcell) counts[ match(sro@meta.data[i, meta.mapnames], hlabs), match(sro@meta.data[i, meta.to.rename], vlabs)] = counts[ match(sro@meta.data[i, meta.mapnames], hlabs), match(sro@meta.data[i, meta.to.rename], vlabs)] + sro@meta.data[i, meta.mapweight]
	}else{
		for(i in 1:nbcell) counts[ match(sro@meta.data[i, meta.mapnames], hlabs), match(sro@meta.data[i, meta.to.rename], vlabs)] = counts[ match(sro@meta.data[i, meta.mapnames], hlabs), match(sro@meta.data[i, meta.to.rename], vlabs)] + 1
	}
	bestname <- paste("Unknown",vlabs,sep='_')
	bestmap <- rep(65535,length(vlabs))
	purety <- rep(0,length(vlabs))
	print(vlabs)
	print(counts)
	for(j in 1:length(vlabs)) {
		print(purety)
		daord <- order(counts[,j],decreasing=T)
		bestmap[j] <- daord[1]
		counts[,j] <- counts[,j] / sum(counts[,j],na.rm=T)
		purety[j] <- counts[daord[1],j]
		print(daord)
		print(puritythreshold)
		print(purety[j])
		if (puritythreshold > purety[j]) bestmap[j] <- 65535
	}
	print(bestmap)
	hehe <- table(bestmap)
	print(table)
	reorder <- rep(0,length(vlabs))
	cur <- 1
	for(i in 1:length(hehe)){

		if (hehe[i] == 1) {
			tmpname = match(as.numeric(names(hehe)[i]) ,bestmap)
			if (names(hehe)[i] == 65535) bestname[tmpname] <- "Mixture"
			else bestname[tmpname] <- hlabs[as.numeric(names(hehe)[i])]
			reorder[tmpname] <- cur
			cur <- cur + 1
		}else{
			which <- bestmap == as.numeric(names(hehe)[i])
			daord <- order(purety[which])
			if (names(hehe)[i] == 65535) tmpname <- "Mixture"
			else tmpname <- hlabs[as.numeric(names(hehe)[i])]
			
			for(j in 1:length(daord)){
				bestname[which][daord[j]] = paste(tmpname," Type", rawToChar(as.raw(64+j)),sep="")
				reorder[which][daord[j]] = cur
				cur <- cur + 1
			}
		}
	}
	levels(sro@meta.data[[meta.to.rename]]) <- bestname 
	if (do.reorder) {
		daord <- match(1:length(vlabs),reorder)
		bestname <- bestname[daord]
		sro@meta.data[[meta.to.rename]] <- factor(sro@meta.data[[meta.to.rename]], levels = bestname)
		counts <- counts[,daord]
	}
	rownames(counts) <- hlabs;
	colnames(counts) <- bestname;
	sro@misc$rename.fingerprint[[meta.mapnames]] <- counts
return(sro)}
#' Run SC3
#'
#' @param directory path to a cellranger output folder, which expected that contains 3 files named barcodes.tsv genes.tsv matrix.mtx respectively. (ex: "./data/")
#'
#' @export
runSC3 <- function(sce){
	if (sce@assays[["counts"]] != "matrix") sce@assays[["counts"]] <- as.matrix(sce@assays[["counts"]])
	if (sce@assays[["logcounts"]] != "matrix") sce@assays[["logcounts"]] <- as.matrix(sce@assays[["logcounts"]])
	sce <- SC3::sc3_prepare(sce);
}

#' Run SC3
#'
#' @param directory path to a cellranger output folder, which expected that contains 3 files named barcodes.tsv genes.tsv matrix.mtx respectively. (ex: "./data/")
#'
#' @export
DoMnnCorrectPCA<- function(sro, meta, vargenelist, nbneighbor = 32, nb.pc = 50){
	len <- length(srolist)
	if (class(sro@meta.data[[meta]]) == "factor") blvl = levels(sro@meta.data[[meta]])
	else blvl = unique(sro@meta.data[[meta]])
	input <- list();
	for( i in 1:length(blvl)) input <- c(input , sro@data[,sro@meta.data[[meta]] == blvl[i] ] )
	library(scran)	
	out <- mnnCorrect(input, k=nbneighbor)

	data <- out$corrected[1];
	for( i in 2:length(blvl)) cbind(data,out$corrected[i])

	tmpdata <- sro@data
	sro@data <- data

	sro <- ScaleData(sro, vars.to.regress=c("nUMI"), genes.use=vargenelist)
	if (class(vargenelist) != "character") vargenelist <- sro@var.genes
	else sro@var.genes <-  vargenelist
	
	if (is.null(cell.weight)) sro <- RunPCA(sro, pc.genes=vargenelist, pcs.compute = nb.pc, reduction.name = "mnnpca")
	sro@data <- tmpdata
	outdata <- sro@data
	sro@scale.data <- NULL
return(list(sro= sro, data = outdata))}

CreateSparse  <- function(path_s, prefix = "", gene_metalist = c(""), ...){
	if (class(path_s) != "character"){
		stop("Expects a path or a list of paths for populating seurat objects")
	}
	library("Seurat")
	if (length(path_s) != length(prefix)){
		prefix <- paste(as.character(1:length(path_s)), "_", sep="")
	}
	metacol <- c()
	for(i in 1:length(path_s)){
		tdata <- Read10X(path_s[i]);
		tdata <- tdata[,(Matrix::colSums(tdata) != 0)]
		colnames(tdata) <- paste(prefix[i], colnames(tdata), sep="_")
		if (i == 1) data <- tdata
		else data = cbind(data,tdata)
	}
	return(data)
}

# source("~/work/readFromNew.R");sro <- CreateSeurat("/warehouse/team218_wh01/lh20/references/Fan2018Spatial/GSE103723_", is.floatingpoint=T)
# c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")
CreateSeurat  <- function(path_s, prefix = "", gene_metalist = c(), do.gunzip=F, is.floatingpoint=F, ...){
	if (class(path_s) != "character"){
		stop("Expects a path or a list of paths for populating seurat objects")
	}
	library("Seurat")
	library("InferN0")
	if (length(path_s) != length(prefix)){
		prefix <- paste(as.character(1:length(path_s)), "_", sep="")
	}
	metacol <- c()

	for(i in 1:length(path_s)){
#		tdata <- InferN0ReadRangerMatrixCPP(path_s[i], do.soup.regression=F ,threshold.genes.ignored=c("thisgenedoesnotexistsurely"),do.plot.soup=F, is.floatingpoint=is.floatingpoint)$data
	
		tdata <- Read10X(path_s[i]);
		if (length(prefix[i]) != 0) colnames(tdata) <- paste(prefix[i], colnames(tdata), sep="_")
		if (i == 1) data <- tdata
		else data = cbind(data,tdata)
		metacol <- c(metacol , rep(prefix[i], ncol(tdata)))
	}
	if (length(gene_metalist) != 0){
		tometa = match(gene_metalist, rownames(data));
		metadata <- t(data[tometa,]);
		summy <- as.matrix(rowSums(as.matrix(metadata)));
		colnames(summy) <- c("sumMetaUMI")
		tmp <- CreateSeuratObject(data[setdiff(1:(dim(data)[1]),tometa),], ...)
		tmp <- AddMetaData(tmp,  data.frame(cbind(summy, as.matrix(metadata))))
	}else{
		tmp <- CreateSeuratObject(data, ...)
	}
	#tmp <- AddMetaData(tmp, metacol, col.name="orig.ident")
	return(tmp)
}

ListSeuratMeta <- function(sro){
	print(names(sro@dr))
	print(colnames(sro@meta.data))
}

CreateSeuratSimpler <- function(paths, prefix){
	library(Seurat)
	for(i in 1:length(prefix)) {
		print(paste("Processing", prefix[i]))
		out <- Read10X(paths[i])
		colnames(out) <- paste(prefix[i],colnames(out))
		if (i == 1) {
			concat <- out
			sample <- rep(prefix[i], ncol(out))
		}else{
			concat <- cbind(concat, out)
			sample <- c(sample, rep(prefix[i], ncol(out)))
		}
	}
	print(dim(concat))
	print(length(sample))
	sro <- CreateSeuratObject(concat)
	sro@meta.data$orig.ident =as.factor(sample)
return(sro)}

# CreateSeuratCustomfilter( paste("/lustre/scratch117/cellgen/team218/lh20/outoutSIGAXX75605", as.character(37:42), "_I-0/outs/raw_gene_bc_matrices/humanintrons/" ,sep=""), c("WW+","HW-","HW+","WH+","HH-","HH+"), paste("/lustre/scratch117/cellgen/team218/lh20/out75605", as.character(37:42),"_E-0/outs/raw_gene_bc_matrices/GRCh38/",sep=""), threshold.UMI=200, threshold.GENE=75,threshold.soup= c(0.08,0.03,0.02,0.02,0.03,0.03), threshold.soup.slope=c(0.05,0.04,0.05,0.05,0.04,0.05)  )
# c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")
# genemeta <- readRDS("/lustre/scratch117/cellgen/team218/lh20/human_genetypes.rds")
CreateSeuratCustomfilter  <- function(path_s, prefix = "", exononlypath= c(), genemeta = c(),  gene_metalist = c(), custom.celllist.import=c(), threshold.UMI=50, threshold.GENE=25, threshold.soup=1.0,threshold.soup.slope=0.0,do.plot.soup=F,do.print.soup.pdf=F,do.soup.regression=T,do.rename=F,add.meta.data=c(), nb.soup.samples=32, do.gzip=T){
	library(InferN0)
	library(Matrix)
	if (do.print.soup.pdf) do.plot.soup<-T
	if (class(path_s) != "character"){
		stop("Expects a path or a list of paths for populating seurat objects")
	}
	library("Seurat")
	if (length(path_s) != length(prefix)){
		prefix <- paste(as.character(1:length(path_s)), "_", sep="")
	}
	if (length(threshold.UMI) != length(path_s)) threshold.UMI = rep(threshold.UMI,  length(path_s))
	if (length(threshold.GENE) != length(path_s)) threshold.GENE = rep(threshold.GENE,  length(path_s))
	if (length(threshold.soup) != length(path_s)) threshold.soup = rep(threshold.soup,  length(path_s))
	if (length(threshold.soup.slope) != length(path_s)) threshold.soup.slope = rep(threshold.soup.slope,  length(path_s))
	tedata = c()
	data=c()
	soupy=c()
	
	daorigid <- c()
	if (is.null(gene_metalist)) gene_metalist <- c("thisgenedoesnotexistsurely")
	metaannot <- NULL
	if (!is.null(custom.celllist.import)) customcelllist  <- c()
	if (do.gzip) library(R.utils)

	for(i in 1:length(path_s)){
		print(paste("processing", prefix[i]))
		if (do.plot.soup) souptitle <- prefix[i]
		else souptitle  <- c()
		if (do.gzip) {
			for( j in c( paste(path_s[i], c("barcodes.tsv.gz", "features.tsv.gz","matrix.mtx.gz"),sep="") ,  paste(exononlypath[i], c("barcodes.tsv.gz", "features.tsv.gz","matrix.mtx.gz"),sep=""))) {
				if (file.exists(j)){
					print(paste("uncompressing", j))
					gunzip(j)
				}
			}
		}
		if (threshold.UMI[i] == 0){
			if (do.plot.soup) souptitle <- prefix[i]
			else souptitle  <- c()
			idata <- InferN0ReadRangerMatrixCPP(path_s[i], do.soup.regression=do.soup.regression,threshold.genes.ignored=gene_metalist,do.soup.regression=T,plot.soup.title=souptitle,nb.soup.samples=nb.soup.samples)
		}else{

			idata <- InferN0ReadRangerMatrixCPP(path_s[i], threshold.low.UMI=threshold.UMI[i], threshold.low.GENE=threshold.GENE[i],threshold.soup.fraction=threshold.soup[i],threshold.soup.logUMIslope=threshold.soup.slope[i] ,threshold.genes.ignored=gene_metalist,do.soup.regression=do.soup.regression,plot.soup.title=souptitle, nb.soup.samples=nb.soup.samples)
		}
		if (!is.null(custom.celllist.import)){
			datlist <- as.character(read.csv(custom.celllist.import[i],sep='\t',header=F)[,1])
			customcelllist <- c(customcelllist, !is.na(match(colnames(idata$data), datlist)))
		}
		#print(paste(prefix[i],colnames(idata$data),sep="_"))
		#print(ncol(idata$data))
		#print(length(colnames(idata$data)))
		idata$data@Dim[1] <- length(idata$data@Dimnames[[1]])
		idata$soup.samples@Dim[1] <- length(idata$soup.samples@Dimnames[[1]])
		colnames(idata$data) <- paste(prefix[i],colnames(idata$data),sep="_")
		daorigid <- c(daorigid, rep(prefix[i], ncol(idata$data)))
		if (do.print.soup.pdf) dev.print(pdf, paste("tmppdf/", prefix[i] , "_soup.pdf",sep=""))
		if ("data" %in% names(idata)){
			if (length(exononlypath) != 0){
				edata <- InferN0ReadRangerMatrixCPP(exononlypath[i], threshold.low.UMI=threshold.UMI[i]*0.25, threshold.low.GENE=threshold.GENE[i]*0.25,do.soup.regression=F,do.output.data.only=T,nb.soup.samples=1)
				
				edata$data@Dim[1] <- length(edata$data@Dimnames[[1]])
				colnames(edata$data) <- paste(prefix[i],colnames(edata$data),sep="_")
				damap <- match(colnames(idata$data), colnames(edata$data))
				damap[is.na(damap)] <- dim(edata$data)[2]
				edata$data = cbind(edata$data, Matrix(0,dim(edata$data)[1],1, sparse=T))
				edata$data <- edata$data[,damap]
				if (is.null(tedata)) tedata <- edata$data
				else tedata <- cbind(tedata, edata$data)
			}
			if (do.soup.regression){
				metaantmp <-  cbind(idata$soupfraction_cell,(match(1:length(idata$soupfraction_cell),order(as.numeric(as.character(idata$soupfraction_cell)))) -0.5)/length(idata$soupfraction_cell), rep(prefix[i], length(idata$soupfraction_cell)))
				

				idata$hardfiltered.UMIs@Dim[1] <- length(idata$hardfiltered.UMIs@Dimnames[[1]])
				idata$filtered.data@Dim[1] <- length(idata$filtered.data@Dimnames[[1]])
				tmptmp <-as.matrix(idata$soupfraction_gene)
				
				print(dim(tmptmp))
				print(dim(idata$filtered.data))
				tmptmp <- cbind(tmptmp, (-0.5 + order(tmptmp[,1]))/dim(tmptmp)[1], Matrix::rowSums(idata$hardfiltered.UMIs[,1:4]),  Matrix::rowSums(idata$hardfiltered.UMIs[,5:ncol(idata$hardfiltered.UMIs)]), Matrix::rowSums(idata$filtered.data) )
				colnames(tmptmp)<- c(paste(prefix[i],c("soupfactor","souprank", "errorUMIs","hardfiltedUMI", "softfilteredUMI" ),sep="_"))
			}else{
				metaantmp = cbind(rep(prefix[i], dim(idata$data)[2]))
				tmptmp = NULL
			}
			print("aliveE")
			rownames(metaantmp) <- colnames(idata$data) 
			if (is.null(data)) {
				accnames <- idata$gene_accession
				gdata <- tmptmp
				data <- idata$data
				if ("soup.samples" %in% names(idata)) {colnames(idata$soup.samples) <- paste(prefix[i],colnames(idata$soup.samples), sep="_"); soupy <- idata$soup.samples}
			}else{
				gdata <- cbind(gdata,  tmptmp)
				data = cbind(data,idata$data)
				if ("soup.samples" %in% names(idata)) {colnames(idata$soup.samples) <- paste(prefix[i],colnames(idata$soup.samples), sep="_"); soupy <- cbind(soupy,idata$soup.samples)}

			}
			print("aliveF")
			if (is.null(metaannot)) metaannot <- metaantmp
			else metaannot <-rbind(metaannot ,metaantmp)
		}

		if (do.gzip) {
			
			for( j in c( paste(path_s[i], c("barcodes.tsv", "features.tsv","matrix.mtx"),sep="") ,  paste(exononlypath[i], c("barcodes.tsv", "features.tsv","matrix.mtx"),sep=""))) {
				if (substr(j,2,3) != "wa"){
				if (file.exists(j)){
					print(paste("compressing", j))
					gzip(j)
				}
			}}
		}

	}
	if (do.rename){
		tmp <- gsub("\tGene Expression","",rownames(data))
		tmpens <- gsub("\t.*","", tmp)
		rownames(data) <- gsub(".*\t","", tmp)
	}
	if (length(gene_metalist) != 1){
		tometa = match(gene_metalist, rownames(data));
		tometa <- tometa[!is.na(tometa)]
		metadata <- data[tometa,,drop=F];
		melvl <- levels(genemeta$biotype)
		memap <- match(rownames(metadata), genemeta$name)
		if (sum(is.na(memap)) > 0){
			print(rownames(metadata)[is.na(memap)])
		}
		memap[!is.na(memap)] <- genemeta$biotype@.Data[memap[!is.na(memap)]]
		summy <- data.frame(row.names = colnames(data))
		for(i in 1:length(melvl)){
			entry <- colSums(metadata[memap == i,,drop=F])
			if (sum(entry) > 0) summy[[paste("Filtered",melvl[i], sep="_")]] <- entry
		}
		melvl <- rownames(data)
		datab <- table(melvl);datab <- names(datab)[datab > 1]; datab <- !is.na(match(rownames(data), datab))
		rownames(data)[datab] <- gsub("\\.","-",make.names(melvl[datab],unique=T))
		tmp <- CreateSeuratObject(data)
		
		data <- data[setdiff(1:(dim(data)[1]), tometa),,drop= F]
		tmp@meta.data$rawnGene <- tmp@meta.data$nGene
		tmp@meta.data$rawnUMI <- tmp@meta.data$nUMI
		tmp@meta.data$nGene <- Matrix::colSums(data != 0)
		tmp@meta.data$nUMI <- Matrix::colSums(data)
		tmp <- AddMetaData(tmp, summy)
	}else{
		melvl <- rownames(data)
		datab <- table(melvl);datab <- names(datab)[datab > 1]; datab <- !is.na(match(rownames(data), datab))
		rownames(data)[datab] <- gsub("\\.","-",make.names(melvl[datab],unique=T))
		tmp <- CreateSeuratObject(data)
	#	print(plot(match(colnames(tmp@assays$RNA@counts),  colnames(data)), 1:ncol(tmp@assays$RNA@counts)))

	}
	tmp@meta.data$orig.ident <- daorigid
	if (!is.null(add.meta.data)){
		projmeta <- add.meta.data[as.character(tmp@meta.data$orig.ident),,drop =F]
		rownames(projmeta) <- rownames(tmp@meta.data)
		for(i in colnames(projmeta)) tmp@meta.data[[i]]	<- as.factor(projmeta[,c(i)])
	}
	
	if (!is.null(genemeta)){
		memap <- match(accnames, genemeta$ID)
		tmp@misc$meta.gene <- genemeta[memap ,, drop=F]
		tmp@misc$meta.gene$nUMI <- Matrix::rowSums(SeuratRawData(tmp))
		tmp@misc$meta.gene$nCell <- Matrix::rowSums(SeuratRawData(tmp) != 0)
		rownames(tmp@misc$meta.gene) <- SeuratRawData(tmp, "R")
		return(tmp)
	}
	if (!is.null(custom.celllist.import)) tmp@meta.data$CellrangerList <- customcelllist

	if (length(exononlypath) != 0) {
		metaannot <- cbind(as.data.frame( Matrix::colSums(tedata) / Matrix::colSums(data))  , metaannot)
		colnames(metaannot) <- c("exon_fraction", rep(c("soup_fraction", "soup_rank"),ifelse(do.soup.regression,1,0)), "orig.ident")
		tmp@misc$exondata <- tedata 
		
	}else colnames(metaannot) <- c(rep(c("soup_fraction", "soup_rank"),ifelse(do.soup.regression,1,0)), "orig.ident")
	if (!is.null(gdata)) tmp@misc$genemetadata <- gdata
	return(tmp)
	print("adding Meta data:")
	rownames(metaannot) <- colnames(data)
	tmp <- AddMetaData(tmp, data.frame(metaannot))
	if (do.soup.regression){
		tmp@meta.data[["soup_fraction"]] <- as.numeric(as.character(tmp@meta.data[["soup_fraction"]]))
		tmp@meta.data[["soup_rank"]] <- as.numeric(as.character(tmp@meta.data[["soup_rank"]]))
	}
	if ("exon_fraction" %in% names(tmp@meta.data)) tmp@meta.data[["exon_fraction"]] <- as.numeric(as.character(tmp@meta.data[["exon_fraction"]]))	

	if (!is.null(soupy)) tmp@misc$soup.samples = soupy
	if ("soup_fraction" %in% colnames(tmp@meta.data)){
		soupcoor <- cbind(tmp@meta.data$soup_fraction, log2(tmp@meta.data$nCount_RNA))
		colnames(soupcoor) <- c("S_1","S_2")
		rownames(soupcoor) <- rownames(tmp@meta.data)
		tmp@reductions$soup <- CreateDimReducObject(embeddings=soupcoor, key = "S_",assay="RNA")
		
#		tmp <- SetDimReduction(tmp, reduction.type = "soup", slot = "cell.embeddings", new.data = soupcoor)
#		tmp <- SetDimReduction(tmp, reduction.type = "soup", slot = "key", new.data = "S")

		soupcoor <- cbind(tmp@meta.data$soup_fraction, tmp@meta.data$exon_fraction , log2(tmp@meta.data$nCount_RNA))
		colnames(soupcoor) <- c("S_1","S_2","S_3")
		rownames(soupcoor) <- rownames(tmp@meta.data)

		tmp@reductions$soupexon <- CreateDimReducObject(embeddings=soupcoor, key = "S_",assay="RNA")
#		tmp <- SetDimReduction(tmp, reduction.type = "soupexon", slot = "cell.embeddings", new.data = soupcoor)
#		tmp <- SetDimReduction(tmp, reduction.type = "soupexon", slot = "key", new.data = "S")
	}
	
	return(tmp)
}

inferCiteSeqIds <- function(sro, pathlist, unmap.subtring = "unmap", threshold= 0.05, metadata = list(), dbl.names=c(), partition=c()){
	library(Seurat)
	sro@misc$CiteSeqData <- list()
	sro@meta.data$CiteSeq_EnrichmentPvalue <- 0
	sro@meta.data$CiteSeq_Alt_EnrichmentPvalue <- 0
	sro@meta.data$CiteSeq_Count <- 0
	sro@meta.data$CiteSeq_Alt_Count <- 0
	sro@meta.data$CiteSeq_Total_Count <- 0
	sro@meta.data$i7.ident <- sro@meta.data$orig.ident
	daord <- setdiff(levels(sro@meta.data$orig.ident), names(pathlist))
	sro@meta.data$orig.ident <- as.character(sro@meta.data$orig.ident)
	tetra <- c()
	
	if (length(metadata) != 0){
		for(i in 1:length(metadata)) sro@meta.data[[names(metadata)[i]]] <- ""
	}
	for(i in 1:length(pathlist)){
		flt <- sro@meta.data$orig.ident == names(pathlist)[i]
		print(paste("processing", names(pathlist)[i], "with", sum(flt), "cells"))
		if (nchar(pathlist[i]) == 0){
			daord <- c(daord, names(pathlist)[i])
			if (length(metadata) != 0){
				for(j in 1:length(metadata)) {
					print(names(metadata)[j])
					sro@meta.data[[names(metadata)[j]]][flt] <-metadata[[j]][[names(pathlist)[i]]]
				}
			}
			print("alive")
		}else{
		hehe <- as.matrix(Seurat::Read10X(as.character(pathlist[i])))
		print(class(hehe))
		print(dim(hehe))
		if (nrow(hehe) == 4){
			print("fixin!")
			print(rownames(hehe))
		#	hehe <- rbind(hehe[1:3,], hehe[4,,drop=F] * 0 ,hehe[4,,drop=F])
			hehe <- rbind(hehe[1:3,], hehe[4,,drop=F] * 0)
			rownames(hehe)[4] <- paste(names(pathlist)[i], "unused", sep = "_")
			print(dim(hehe))
			print(rownames(hehe))
		}else{
			hehe <- hehe[1:4,]
		}
		daunknownname <- paste(names(pathlist)[i], "Unknown",sep="_")	
		sro@meta.data$orig.ident[flt] <- daunknownname
		cellnames <- sub("^.*_", "" , sub("-[1234567890]*", "", rownames(sro@meta.data)[flt]))
		print(paste("barcode with Ns:", sum(grepl("N",colnames(hehe)))))
		map <- match(cellnames, colnames(hehe))
		datags <- t(hehe[,map[!is.na(map)]])
		tmp <- grepl(unmap.subtring, colnames(datags)) ; datags <- datags[, !tmp]

		if (is.null(partition)) {
			daoutput <- funTagEval(datags)
		}else{
			tmpl <- levels(partition)
			recap <- partition[flt]
			daoutput <- list(best = rep(0, nrow(datags)), pvalue= rep(0, nrow(datags)),second.best= rep(0, nrow(datags)),second.pvalue= rep(0, nrow(datags)), normalized.logpvalue= rep(0, nrow(datags)), second.normalized.logpvalue= rep(0, nrow(datags)))
			for(k in 1:length(tmpl)){
				print(table(recap[!is.na(map)]))
				sub <- which(recap[!is.na(map)] == tmpl[k])
				print(length(sub))
				daoutput2 <- funTagEval(datags[sub,])
				for(m in names(daoutput2)){
					daoutput[[m]][sub] <- daoutput2[[m]]
				}
			}
		}
		print("alive!")
		print(table(daoutput$best))
		if (sum(tmp) != 0) daunknownname <- colnames(datags)[which(tmp)[1]]


		print(ncol(datags))
		curthresh <- threshold / ncol(datags)
		sflt <- which(flt)[!is.na(map)]
		sro@meta.data$CiteSeq_EnrichmentPvalue[sflt] <-	daoutput$pvalue
		sro@meta.data$CiteSeq_Alt_EnrichmentPvalue[sflt] <- daoutput$second.pvalue
		sro@meta.data$CiteSeq_Total_Count[sflt] <- rowSums(datags)
		for(j in 1:length(daoutput$best)){
			sro@meta.data$CiteSeq_Count[sflt[j]] <- datags[j,daoutput$best[j]]
			sro@meta.data$CiteSeq_Alt_Count[sflt[j]] <- datags[j,daoutput$second.best[j]]
		}
		
		if (ncol(datags) == 4){
			print("TETRA!!!!!")
			if (is.null(tetra)) tetra <- matrix(0,nrow(sro@meta.data),3)
			tetra[sflt, 1] <- rowSums(datags[,2:3]) /  rowSums(datags)
			tetra[sflt, 2] <- rowSums(datags[,3:4]) / rowSums(datags)
			tetra[sflt, 3] <- rowSums(datags[,c(2,4)]) /  rowSums(datags)
		}
		daoutput$best[daoutput$pvalue >= curthresh] <- length(colnames(datags))+1
		daoutput$best[is.na(daoutput$best)] <- length(colnames(datags))+1
		altbest <- daoutput$best
		tmpbest <- (!is.na(daoutput$pvalue)) &  (!is.na(daoutput$second.pvalue)) & (daoutput$second.pvalue < curthresh)
		print(sum(tmpbest))
		altbest[tmpbest] <- daoutput$second.best[tmpbest]
		if (length(metadata) != 0){
			for(j in 1:length(metadata)) {
				tmpbest <- altbest
				tmpbest[ (metadata[[j]][[names(pathlist)[i]]][altbest] != metadata[[j]][[names(pathlist)[i]]][daoutput$best])] <- length(colnames(datags))+1
				print(table(altbest,daoutput$best))
				sro@meta.data[[names(metadata)[j]]][flt] <-metadata[[j]][[names(pathlist)[i]]][length(colnames(datags))+1]
				sro@meta.data[[names(metadata)[j]]][sflt] <- metadata[[j]][[names(pathlist)[i]]][tmpbest]
			}
		}
		daoutput$best[daoutput$second.pvalue < curthresh] <- length(colnames(datags))+2
		if (is.null(dbl.names)) daunknownname2 <- paste(names(pathlist)[i], "_Dblt", sep="")
		else daunknownname2 <- dbl.names[[names(pathlist)[i]]]
		print(sum(is.na(sro@meta.data$orig.ident)))
		sro@meta.data$orig.ident[flt] <- rep(daunknownname, sum(flt))
		print(sum(is.na(sro@meta.data$orig.ident)))
		
		daoutput$best[is.na(daoutput$best)] <- length(colnames(datags))+1
		print(daoutput$best)
		sro@meta.data$orig.ident[sflt] <- c(colnames(datags), daunknownname, daunknownname2)[daoutput$best] 
		print(sum(is.na(sro@meta.data$orig.ident)))			

		daord <- c(daord, colnames(datags), daunknownname,daunknownname2)
		sro@misc$CiteSeqData[[names(pathlist)[i]]] <- hehe
		}
	}
	print(table(sro@meta.data$i7.ident, sro@meta.data$orig.ident))
	print(daord)
	sro@meta.data$orig.ident <- factor(sro@meta.data$orig.ident, levels = daord)
	
	print(table(sro@meta.data$orig.ident))
	if (is.null(tetra)) return(sro)
		colnames(tetra) <- c("CS_1","CS_2", "CS_3")
		rownames(tetra) <- rownames(sro@meta.data)
		sro@reductions$citeseq <- CreateDimReducObject(embeddings=tetra, key = "CS_",assay=sro@active.assay)
	return(sro)
}

exportSeuratToTables <- function(sro, pathprefix){
	library(R.utils)	
	library(InferN0)
	scp <- InferN0Init(sro@assays$RNA@counts)
	InferN0ExportDataAsMtx(scp,paste(pathprefix, "_", sep=""),do.gzip=T)
	reduction <- c()
	for(i in names(sro@reductions)){
		if (is.null(reduction)) reduction <- sro@reductions[[i]]@cell.embeddings
		else reduction <- cbind(reduction, sro@reductions[[i]]@cell.embeddings)
	}
	filename =  paste(pathprefix, "_coor.tsv",sep="")
	write.table(reduction, filename, sep="\t", col.names=NA, quote=F)
	gzip(filename)
	filename = paste(pathprefix, "_meta.tsv",sep="")
	write.table(sro@meta.data, filename, sep="\t", col.names=NA,quote=F)
	gzip(filename)


}

exportDEToTables <- function(tab, pathprefix){
	library(R.utils)
	filename =  paste(pathprefix, "_DE_genes.tsv",sep="")
	write.table(tab$gene, filename, sep="\t", col.names=NA,quote=F); gzip(filename)
	filename =  paste(pathprefix, "_DE_genes_preFDR.tsv",sep="")
	write.table(tab$gene_preFDR, filename, sep="\t", col.names=NA,quote=F); gzip(filename)
	filename =  paste(pathprefix, "_DE_terms_and_pathways.tsv",sep="")
	write.table(tab$go, filename, sep="\t", col.names=NA,quote=F); gzip(filename)
	filename =  paste(pathprefix, "_consensusDE_genes.tsv",sep="")
	write.table(tab$consensus.gene, filename, sep="\t", col.names=NA,quote=F); gzip(filename)
	filename =  paste(pathprefix, "_consensusDE_genes_preFDR.tsv",sep="")
	write.table(tab$consensus.gene_preFDR, filename, sep="\t", col.names=NA,quote=F); gzip(filename)
	filename =  paste(pathprefix, "_consensusDE_terms_and_pathways.tsv",sep="")
	write.table(tab$consensus.go, filename, sep="\t", col.names=NA,quote=F); gzip(filename)

	filename =  paste(pathprefix, "_tables.rds",sep="")
	saveRDS(list(gene=tab$gene, gene_preFDR=tab$gene_preFDR, go=tab$go, consensus.gene =tab$consensus.gene, consensus.gene_preFDR= tab$consensus.gene_preFDR,consensus.go=tab$consensus.go, summary=tab$summary, summary_preFDR=tab$summary_preFDR), filename);
	
	filename =  paste(pathprefix, "_raw_tables.rds",sep="")
	saveRDS(list(DEseq=tab$DEseq, Wilcox=tab$Wilcox, log2FC_FAD_to_Empty= tab$log2FC_FAD_to_Empty, log2FC_CTRL_to_Empty= tab$log2FC_CTRL_to_Empty), filename);
}


addExonMetaDataSeurat <- function(sro, paths_all, prefix, low.UMI.threshold=50, cellnamelenght=0){
	library("InferN0")
  if (class(paths_all) != "character"){
		stop("Expects 'paths_exon' to be path or a list of paths")
	}
	if (class(prefix) != "character"){
		stop("Expects 'prefix' to a list of prefixes to cellnames appended to default cell names, in each sample ('_' not included)")
	}
	if (length(paths_all) != length(prefix)){
		stop("Expects equal size for path lists, (and coherent ordering of elements)")
	}
	d <- dim(sro@raw.data)
	if (!("raw.exon.data" %in% names(sro@misc))){
		sro@misc$raw.exon.data <- Matrix(0,nrow= d[1], ncol=d[2],sparse=T)
		rownames(sro@misc$raw.exon.data) <- rownames(sro@raw.data)
		colnames(sro@misc$raw.exon.data) <- colnames(sro@raw.data)
	}
	if (!("gene.meta.data" %in% names(sro@misc))){
		sro@misc$gene.meta.data <- Matrix(0,nrow= d[1], ncol=0,sparse=T)
		rownames(sro@misc$gene.meta.data) <- rownames(sro@raw.data)
	}
	for(i in 1:length(paths_all)){
		edata <- InferN0ReadRangerMatrixCPP(paths_all[i], threshold.low.UMI=low.UMI.threshold,threshold.genes.ignored=c("thisgenedoesnotexistsurely"))
		damap <- match(rownames(edata$filtered_data),rownames(sro@misc$gene.meta.data))
		curcol <- paste("Filtered_UMIs",prefix,sep="_")
		if (!(curcol %in% colnames(sro@misc$gene.meta.data))){
			tmpcolname <- colnames(sro@misc$gene.meta.data)
			sro@misc$gene.meta.data <- cbind(sro@misc$gene.meta.data, Matrix(0,nrow= d[1], ncol=1,sparse=T))
			colnames(sro@misc$gene.meta.data) <- c(tmpcolname, curcol)
		}
		sro@misc$gene.meta.data[damap[!is.na(damap)],c(curcol)] <- Matrix::rowSums(edata$filtered_data)[(1:length(damap))[!is.na(damap)]]
	}
return(sro)
}

addIntronMetaDataSeurat <- function(sro, paths_exon, prefix, low.UMI.threshold=50, cellnamelenght=0){
	library("InferN0")
  if (class(paths_exon) != "character"){
		stop("Expects 'paths_exon' to be path or a list of paths")
	}
	if (class(prefix) != "character"){
		stop("Expects 'prefix' to a list of prefixes to cellnames appended to default cell names, in each sample ('_' not included)")
	}
	if (length(paths_exon) != length(prefix)){
		stop("Expects equal size for path lists, (and coherent ordering of elements)")
	}
	d <- dim(sro@raw.data)

	if (!("raw.exon.data" %in% names(sro@misc))){
		sro@misc$raw.exon.data <- Matrix(0,nrow= d[1], ncol=d[2],sparse=T)
		rownames(sro@misc$raw.exon.data) <- rownames(sro@raw.data)
		colnames(sro@misc$raw.exon.data) <- colnames(sro@raw.data)
	}

	for(i in 1:length(paths_exon)){
		#edata <- Read10X(paths_exon[i]);
#		if ("package:InferN0" %in% search()) 
		edata <- InferN0ReadRangerMatrixCPP(paths_exon[i], threshold.low.UMI=low.UMI.threshold,threshold.genes.ignored=c("thisgenedoesnotexistsurely"))$data
	#	else edata <- InferN0ReadRangerMatrix(paths_exon[i], low.UMI.threshold=low.UMI.threshold)

		plen <- nchar(prefix)
		subset <- (substr(colnames(sro@raw.data),0,plen) == prefix)
		if (cellnamelenght == 0) damap <- match( substr(colnames(sro@raw.data)[subset], plen+2, 9999999), colnames(edata))	
		else damap <- match( substr(colnames(sro@raw.data)[subset], plen+2, plen+1+cellnamelenght), substr(colnames(edata), 1,16))

		
		rowmap <- match(rownames(sro@misc$raw.exon.data), rownames(edata))
		if (length(rownames(edata)) != sum(!is.na(rowmap))) print(paste("Warning only ", sum(!is.na(rowmap)), " genes out of ", length(rownames(edata)), " were found within ", paths_exon[i],sep=''))

		
		

		subset <- which(!is.na(damap))
		if (length(subset) != length(damap)) print(paste("Warning only ", length(subset), " cells out of ", length(damap), " were found within ", paths_exon[i],sep=''))

		damap <- damap[!is.na(damap)]
		edata <- edata[,damap]

		#edata[,subset] <- 0;

		narowmap <- which(!is.na(rowmap))
		rowmap <- rowmap[!is.na(rowmap)]
	  edata <- edata[rowmap,]
		rownames(edata) <- rownames(sro@misc$raw.exon.data)[narowmap];
		colnames(edata) <- colnames(sro@misc$raw.exon.data)[subset];

		#print(max(rowmap))
		#print(length(rowmap))
		#print(dim(sro@misc$raw.exon.data[,subset]))
		#print(dim(edata))
		#print("time to skrink...")

	  #edata[,subset] <- 0
		print("time to merge...")

		sro@misc$raw.exon.data <- InferN0SparseMatrixMerge(sro@misc$raw.exon.data, edata)

		print("time to zero...")		


#		colnames(tdata) <- paste(prefix[i], colnames(tdata), sep="_")
#		if (i == 1) data <- tdata
#		else data = cbind(data,tdata)
#		metacol <- c(metacol , rep(prefix[i], ncol(tdata)))
	}
	print("updating metadata")
	meta <- data.frame(1.0 - (Matrix::colSums(sro@misc$raw.exon.data)/ Matrix::colSums(sro@raw.data)));
	colnames(meta) <- c("Intronic_Fraction")
	library(Seurat)
	sro <- AddMetaData(sro,meta)
return(sro)}

SingleCellExpFromSeurat <- function(sro, convert.to.dense=F){
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("not a seurat object")

	library(SingleCellExperiment)
	if (convert.to.dense){
		assay <- list(counts=as.matrix(sro@raw.data))
		if (!is.null(sro@data)) assay$normcounts <- as.matrix(sro@data)
		if (!is.null(sro@scale.data))	assay$logcounts <- as.matrix(sro@scale.data)
	}else{
		assay <- list(counts=sro@raw.data)
		if (!is.null(sro@data)) assay$normcounts <- sro@data
		if (!is.null(sro@scale.data))	assay$logcounts <- sro@scale.data
	}

	sce <- SingleCellExperiment(assays= assay, rowData = data.frame(gene_names = rownames(sro@raw.data)), colData = data.frame(cell_names = colnames(sro@raw.data)) )
	rowData(sce)$feature_symbol <- rownames(sro@data)
	colData(sce)$feature_symbol <- colnames(sro@data)
	return(sce)
}


DoFindMarkerGenes <- function(input, mt_threshold=0.01, filtersoup=1.0, do.NB.fit=F, use.genes =c(), meta.genes=c(), use.cells =c(), meta.cells=c(), min.count.per.gene=10,suppress.plot=F){
	  library(M3Drop)
		#if (is.null(input@calc.params$CreateSeuratObject$normalization.method)){
		#	stop("seurat object is not normalized")
		#}
		#output <- 3DropDifferentialExpression(as.matrix(sro@raw.data),mt_method = "fdr",mt_threshold = mt_threshold)

		if ((class(input) == "seurat")|(class(input) == "Seurat")){
			 use.genes <- SeuratGeneListQuery(input, meta.genes,use.genes)
		 	 use.cells <- SeuratCellListQuery(input, meta.cells,use.cells)
			 rawdata <-SeuratRawData(input)[use.genes,use.cells, drop=F]
		}else{
			 if (is.null(use.cells)) rawdata = input
			 else rawdata = input[,use.cells]
		}
		tmp <- Matrix::rowSums(rawdata) >= min.count.per.gene 
		rawdata <- rawdata[tmp,]
		use.genes <- use.genes & tmp

		if (filtersoup == 1.0) {
			if (do.NB.fit) m3input <- rawdata / input@meta.data@nUMI
			else m3input <- rawdata
		}else{	# assume metadata is available!
			soupcols <- grep("souprank", colnames(input@misc$genemetadata))
			makeit <- (rowSums(input@misc$genemetadata[,soupcols]) < (filtersoup * length(soupcols)))
			filter <- match(rownames(rawdata), rownames(input@misc$genemetadata))
			filter[!is.na(filter)] <- makeit[filter[!is.na(filter)]]
			filter[is.na(filter)] <- F
			print(paste(sum(filter)," gene passed soup threshold"))
			if (do.NB.fit) m3input <- (rawdata / input@meta.data@nUMI[use.cells])[filter,]
			else m3input <- rawdata[filter,]
		}
		if (do.NB.fit){
			output <- NBumiFeatureSelectionHighVar(NBumiFitModel(m3input))			
		}else{
			print(dim(m3input))
			if (exists("M3DropFeatureSelection")) output <- M3DropFeatureSelection(m3input, mt_method = "fdr",mt_threshold = mt_threshold,suppress.plot=suppress.plot)
			else {
				output <- M3DropDifferentialExpression(as.matrix(m3input), mt_method = "fdr",mt_threshold = mt_threshold,suppress.plot=suppress.plot)
			}
			return(output)
		}
		if ("genemetadata" %in% names(input@misc)){
			soupcols <- grep("souprank", colnames(input@misc$genemetadata))
			map <- match(rownames(output),rownames(input@misc$genemetadata))
			tomap <- c((rowSums(input@misc$genemetadata[,soupcols]) / length(soupcols)), 1)
			map[is.na(map)] <- dim(input@misc$genemetadata)[1]
			output$soupfraction_average <- tomap[map]	
		}
	return(output)
}

MyFindMarkerGenes <- function(input, xlim = NA){
	if (class(input) == "seurat"){
		gene_info <- bg__calc_variables(sro@raw.data);

		      xes <- log(gene_info$s)/log(10);
		      put_in_order <- order(xes);
		      dens.col <- densCols(xes, gene_info$p, colramp=colorRampPalette(c("grey75","black")))

		if (!suppress.plot) {
		      	par(fg="black")
			if (!(sum(is.na(xlim)))) {
			      	plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
			} else {
			      	plot(xes,gene_info$p, main="", ylab="", xlab="", col = dens.col,pch=16, ylim=c(0,1))
			}
			title(ylab="Dropout Rate", line=2)
			title(xlab="log10(expression)", line=2)
		}
		invisible(list(gene_info = gene_info, xes=xes, order=put_in_order));
	}
}
bg__calc_variables <- function(expr_mat) {
		if  (class(expr_mat) == "dgCMatrix"){
			library(Matrix)
			p <- 1 - Matrix::rowSums(expr_mat > 0)/ncol(expr_mat);
			s <- Matrix::rowSums(expr_mat) / ncol(expr_mat)
			return(list(p = p, s_stderr = sqrt( (rowMeans(expr_mat^2) - s^2)/ncol(expr_mat)) , p_stderr = sqrt(p*(1-p)/ncol(expr_mat)) )); #sparsep_stderr = sqrt(p*(1-p)/ncol(expr_mat)))
		}
    if (class(expr_mat) != "matrix" & class(expr_mat) != "dgCMatrix") {
	warning("Warning: not a recognized matrix class, coercing to 'matrix'.")
	expr_mat <- as.matrix(expr_mat)
    }
    if (sum(is.na(expr_mat)) > 0) {
	stop("Error: Expression matrix contain NA values.");
    }
     sum_zero <- prod(dim(expr_mat)) - sum(expr_mat > 0)
    lowest <- min(expr_mat)
    if (lowest < 0) {stop("Error: Expression matrix cannot contains negative values! Has the matrix been log-transformed?")}
    
    # Deal with strangely normalized data
    if (lowest > 0) {
        warning("Warning: No zero values (dropouts) detected will use minimum expression value instead.")
        # If no zeros in expression matrix convert minimum value into zero
        #expr_mat <- round(expr_mat, digits=2) # Round to accomodate errors
        min_val <- lowest+0.05 # Instead of rounding to accomodate errors
        expr_mat[expr_mat == min_val] <- 0;
    }
    if (sum_zero < 0.1*prod(dim(expr_mat))) {
        # Less than 10% zeros
        warning("Warning: Expression matrix contains few zero values (dropouts) this may lead to poor performance.")
    }
    
    p <- 1 - rowSums(expr_mat > 0)/ncol(expr_mat)
    if (sum(p == 1) > 0) {
	warning(paste("Warning: Removing", sum(p==1),"undetected genes."))
	expr_mat <- expr_mat[p < 1,]
        p <- 1-rowSums(expr_mat > 0)/ncol(expr_mat)
    }

    p_stderr <- sqrt(p*(1-p)/ncol(expr_mat))
    s <- rowMeans(expr_mat)
    s_stderr <- sqrt( (rowMeans(expr_mat^2) - s^2)/ncol(expr_mat) ) #sparse matrix friendly
#    s_stderr <- rowSds(expr_mat)/sqrt(ncol(expr_mat))
    names(s_stderr) <- rownames(expr_mat)
    return(list(s = s, p = p, s_stderr = s_stderr, p_stderr = p_stderr))
}
BenchDESeq <- function(datamatrix, posneg){
	library("DESeq2")
	cold <- data.frame(posneg)
	rownames(cold) <- colnames(datamatrix)
	colnames(cold) <- c("condition")
	cold$condition <- as.factor(cold$condition)
	destr <- DESeqDataSetFromMatrix(datamatrix,cold, ~ condition)
	destr <- estimateSizeFactors(destr,type= "iterate")
	DEseq(destr);
}

sumDuplicateRows <- function(mat){
	entries <- table(rownames(mat)) 
	entries <-entries[entries>1]
	if (length(entries) == 0) return(mat);
	entries <- names(entries)
	todrop <- c()
	toflt <- rep(F, nrow(mat))
	for(i in 1:length(entries)){
		flt <- (rownames(mat) == entries[i])
		print(sum(is.na(flt)))
		print(sum(flt))
		curval <- Matrix::colSums(mat[flt,])
		flt <- which(flt)
		toflip <- match(entries[i],rownames(mat))
		mat[toflip,] <- curval
		flt[toflip] <- F
		toflt <- toflt | flt
	}
	return(mat[!toflt,]);
}

# vargene: list of gene or number of verable genes to find
NormalizeSeuratNew <- function(sro, metabatch, vargenes, nbdims=30,nbpc=100, use.cells= c(), meta.cells=c(),do.UMAP=T,suppress.plot=F, vargene.filter=c(),resolution=0.8, vars.to.regress= NULL, mincells=200){

	use.cells <- SeuratCellListQuery(sro, meta.cells,use.cells)
	tmpmeta <- sro@meta.data
	tmpmisc <- sro@misc
	library(Seurat)
	if (length(vargenes) == 1){
		if (is.null(vargene.filter)) vargene.filter = length(SeuratRawData(sro, "R"))
	} 

	if (!is.null(metabatch)){
		sron <- CreateSeuratObject(SeuratRawData(sro))
		sron@meta.data$Batch= sro@meta.data[, metabatch]
		ssron <- SplitObject(object = sron, split.by = "Batch")
		tflt <- rep(T,length(ssron))
		for(i in 1:length(ssron)){
			tflt[i] = (nrow(ssron[[i]]@meta.data) >= mincells) 
		}
		ssron <- ssron[tflt]
		for (i in 1:length(ssron)) {
			print(i)
			ssron[[i]] <- NormalizeData(object = ssron[[i]], verbose = FALSE)
			print("ok")
			if (length(vargenes) > 1) ssron[[i]]@assays$RNA@var.features = vargenes
			else{
				dalist <- rownames(DoFindMarkerGenes(ssron[[i]]@assays$RNA@counts,suppress.plot=suppress.plot,use.genes=vargene.filter ))
				if (length(dalist) > vargenes) dalist <- dalist[1:vargenes]
				ssron[[i]]@assays$RNA@var.features <- dalist
			}
		}
		
		print("Done")
		dadims <- 1:nbdims
		sro <- IntegrateData(FindIntegrationAnchors(ssron[names(ssron)], dims = dadims), dims =dadims, verbose =F)
		DefaultAssay(sro) <- "integrated"
		if (length(vargenes) == 1){
			vargenes <- c();
			for (i in 1:length(ssron)) vargenes <- c(vargenes , ssron[[i]]@assays$RNA@var.features)
			vargenes <- unique(vargenes)
			print(vargenes)
		}
	}else if (length(vargenes) == 1){
		dalist <- rownames(DoFindMarkerGenes(SeuratRawData(sro),suppress.plot=suppress.plot,use.genes=vargene.filter))
		if (length(dalist) > vargenes) dalist <- dalist[1:vargenes]
		vargenes <- dalist
	}
	print("Adding Metadata")
	map <- match(rownames(sro@meta.data),rownames(tmpmeta))
	for(i in setdiff(colnames(tmpmeta), colnames(sro@meta.data))) sro@meta.data[[i]] <- tmpmeta[[i]][map]
	sro@misc <- tmpmisc
	if ("exondata" %in% names(sro@misc)) sro@misc$exondata <- sro@misc$exondata[,map]
	sro@misc$permutation <- map
	
	if (sum(!use.cells) != 0){

		mro <- SubsetData(sro, cells = use.cells)
		print(paste("Scaling Data", sum(use.cells)) )
		mro <- tryCatch(ScaleData(mro,features = vargenes,vars.to.regress=vars.to.regress), error=function(e){print("Scaling failed");return(sro)})

		print("Computing PCA")
		mro <- tryCatch(RunPCA(mro, npcs = nbpc,features=vargenes), error=function(e){print("PCA failed");return(sro)})
	
		print("Computing KNN")
		mro <- tryCatch(FindNeighbors(mro, reduction = "pca", dims = 1:nbpc), error=function(e){print("KNN failed");return(sro)})
		print("Computing Clustering")
		mro <- tryCatch(FindClusters(mro, ,resolution=resolution), error=function(e){print("CLustering failed");return(sro)})

		clustername = paste(mro@active.assay, "_snn_res.0.8", sep="")
		print(paste("Import", clustername))
		sro@meta.data[[clustername]] = "filtered"
		sro@meta.data[[clustername]][use.cells] =  mro@meta.data[[clustername]] 
		if (do.UMAP) mro <- tryCatch(RunUMAP(mro, reduction = "pca", dims = 1:nbpc), error=function(e){print("UMAP failed");return(sro)})
		print("Computing TSNE")
		mro <- tryCatch(RunTSNE(mro, reduction = "pca", dims = 1:nbpc), error=function(e){print("TSNE failed");return(sro)})
		for(redf in names(mro@reductions)){
			print(paste( redf ,"transfer"))
			damat <- matrix(0, nrow(sro@meta.data), ncol(mro@reductions[[redf]]@cell.embeddings))
			damat[use.cells,] = mro@reductions[[redf]]@cell.embeddings
			sro@reductions[[redf]] <- CreateDimReducObject(damat, key = redf, assay = mro@reductions[[redf]]@assay.used)
		}

	}else{
	


		print("Scaling Data")
		sro <- tryCatch(ScaleData(sro,features = vargenes,vars.to.regress=vars.to.regress), error=function(e){print("Scaling failed");return(sro)})

		print("Computing PCA")
		sro <- tryCatch(RunPCA(sro, npcs = nbpc,features=vargenes), error=function(e){print("PCA failed");return(sro)})
	
		print("Computing KNN")
		sro <- tryCatch(FindNeighbors(sro, reduction = "pca", dims = 1:nbpc), error=function(e){print("KNN failed");return(sro)})
		print("Computing Clustering")
		sro <- tryCatch(FindClusters(sro, ,resolution=resolution), error=function(e){print("CLustering failed");return(sro)})

		#if (do.UMAP) sro <- tryCatch(RunUMAP(sro, reduction = "pca", dims = 1:nbpc), error=function(e){print("UMAP failed");return(sro)})
		print("Computing TSNE")
		sro <- tryCatch(RunTSNE(sro, reduction = "pca", dims = 1:nbpc), error=function(e){print("TSNE failed");return(sro)})
	}
	# reinsert filtered data!


#	print("Computing TSNE")
#	sro <- tryCatch(RunTSNE(sro, reduction = "pca", dims = 1:nbpc,dim.embed=3,reduction.name="tsne3d"), error=function(e){print("TSNE3D failed");return(sro)})

return(sro)
}

# vargene: list of gene or number of verable genes to find
NormalizeSeuratPCs <- function(sro, metabatch, vargenes, nbdims=30,nbpc=100, use.cells= c(), meta.cells=c(),suppress.plot=F, vargene.filter=c()){
	use.cells <- SeuratCellListQuery(sro, meta.cells,use.cells)
	tmpmeta <- sro@meta.data
	tmpmisc <- sro@misc
	library(Seurat)
	if (is.null(vargene.filter)) vargene.filter = 1:length(SeuratRawData(sro, "R"))
	else vargene.filter = !(SeuratRawData(sro, "R") %in% vargene.filter) 
	print(sum(vargene.filter))

	sron <- CreateSeuratObject(SeuratRawData(sro)[,use.cells])
	sron@meta.data$Batch= sro@meta.data[use.cells, metabatch]
	ssron <- SplitObject(object = sron, split.by = "Batch")
	print(length(ssron))
	for (i in 1:length(ssron)) {
		ssron[[i]] <- NormalizeData(object = ssron[[i]], verbose = FALSE)
		if (length(vargenes) > 1) ssron[[i]]@assays$RNA@var.features = vargenes
		else{
			dalist <- rownames(DoFindMarkerGenes(ssron[[i]]@assays$RNA@counts,suppress.plot=suppress.plot,use.genes=vargene.filter ))
			print(paste("var g", length(dalist)))
			if (length(dalist) > vargenes) dalist <- dalist[1:vargenes]
			ssron[[i]]@assays$RNA@var.features <- dalist
		}
	}
	dadims <- 1:nbdims
	sros <- IntegrateData(FindIntegrationAnchors(ssron[names(ssron)], dims = dadims), dims =dadims, verbose =F)
	DefaultAssay(sros) <- "integrated"
	if (length(vargenes) == 1){
		vargenes <- c();
		for (i in 1:length(ssron)) vargenes <- c(vargenes , ssron[[i]]@assays$RNA@var.features)
		vargenes <- unique(vargenes)
	}
	sros <- tryCatch(ScaleData(sros,features = vargenes,vars.to.regress=vars.to.regress), error=function(e){print("Scaling failed");return(sros)})
	sros <- tryCatch(RunPCA(sros, npcs = nbpc,features=vargenes), error=function(e){print("PCA failed");return(sros)})
return(sros)}


RunSlingshoti <- function(sro, clusters){


}


# now does DR 
NormalizeSeurat <- function(sro,vargenelist,nb.pc = 20, cell.weight=c(),purge.scale=T,do.show.elbow=T){
        library("Seurat")
	if (!is.null(cell.weight)) {
		if (length(cell.weight) != ncol(sro@raw.data)) stop("Invalid length for cell.weight")
		sro@meta.data$cell.weight <- cell.weight
	}
	sro <- NormalizeData(sro, normalization.method = "LogNormalize",display.progress=F)
	sro <- ScaleData(sro, vars.to.regress=c("nUMI"), genes.use=vargenelist)
	if (class(vargenelist) != "character") vargenelist <- sro@var.genes
	else sro@var.genes <-  vargenelist
	print("Done Scaling data")
#	saveRDS(sro,paste("/lustre/scratch117/cellgen/team218/lh20/frozen", dim(sro@raw.data)[2],".sro.rds",sep=""))
	
	if (is.null(cell.weight)) {
		sro <- RunPCA(sro, pc.genes=vargenelist, pcs.compute = nb.pc)
		pcname = "pca"
	}else{
		var <- matrix(0, length(vargenelist), length(vargenelist));
		mean <- matrix(0, length(vargenelist));
		submatrix <- sro@scale.data[vargenelist,]
		pb <- txtProgressBar()
		for(i in 1:ncol(submatrix)){
			setTxtProgressBar(pb, i / ncol(submatrix))
			tmp <- submatrix[,i] *cell.weight[i];
			mean <- mean + tmp
			var <- var + (tmp  %*% t(submatrix[,i])) 
		}
		print("Done Suffitient Stats")
		tmp <- sum(cell.weight)
		mean <- mean / tmp;
		var <- (var / tmp) - mean %*% t(mean)
		tmp <- svd(var)
		tmp$u <- tmp$u[,1:nb.pc]
		tmp$d <- tmp$d[1:nb.pc]
		var <- matrix(0, ncol(submatrix), nb.pc);
		pb <- txtProgressBar()
		for(i in 1:ncol(submatrix)){
			setTxtProgressBar(pb, i / ncol(submatrix))
			var[i,] <- t(submatrix[,i]) %*% tmp$u 
		}
		rownames(var) <- colnames(sro@scale.data)
		colnames(var) <- paste("PC", 1:nb.pc,sep="")
		rownames(tmp$u) <- vargenelist
		colnames(tmp$u) <- paste("PC", 1:nb.pc,sep="")
		pcname = "cwpca"
		sro@dr[[pcname]] <- new("dim.reduction",cell.embeddings=var,gene.loadings=tmp$u , sdev=tmp$d,key = "PC")
		print("Dona Saving PCA object")
	}

	if (purge.scale) sro@scale.data <- NULL# remove scaledata for saving memory...
	if ( (do.show.elbow)&&(pcname == "pca")) PCElbowPlot(sro,num.pc=nb.pc)
return(sro)}
DoHarmonize <- function(sro, orig.ident.to.batch=c(), metaslot= "orig.ident", group.by.secondary = NULL, theta = 1, theta2 = 1, sigma = 0.1, alpha = 0.1, nclust = 100, 
    tau = 0, block.size = 0.05, max.iter.harmony = 10, max.iter.cluster = 200, epsilon.cluster = 1e-05, epsilon.harmony = 1e-04, burn.in.time = 10, plot_convergence = F, include_projection=F){
	library(harmony)
	library(dplyr)
	dameta <-  as.factor(sro@meta.data[[metaslot]])
	if (is.null(orig.ident.to.batch)){
		cellsample <- dameta
	}else{
		permute <- match(levels(dameta), names(orig.ident.to.batch))
		if (sum(is.na(permute)) != 0) {
			print(setdiff(levels(dameta), names(orig.ident.to.batch)))
			stop(paste("missing id to batch in mapping for the above names"))
		}
		permute <- orig.ident.to.batch[levels(dameta)]
		cellsample <- permute[sro@meta.data[[metaslot]]@.Data]
	}
    	if (!is.null(group.by.secondary)) batches_secondary <- sro@meta.data[[group.by.secondary]]
    	else batches_secondary <- NULL

	harmonyEmbed <- HarmonyMatrix(sro@dr$pca@cell.embeddings, 
        cellsample, batches_secondary, theta, 
        theta2, sigma, alpha, nclust, tau, block.size, max.iter.harmony, 
        max.iter.cluster, epsilon.cluster, epsilon.harmony, burn.in.time, 
        plot_convergence)
    	rownames(harmonyEmbed) <- row.names(sro@meta.data)
    	colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
    	
	sro <- SetDimReduction(sro, reduction.type = "harmony", slot = "cell.embeddings", new.data = harmonyEmbed)
        sro <- SetDimReduction(sro, reduction.type = "harmony", slot = "key", new.data = "Harmony")
 	if (include_projection){
        	sro <- ProjectDim(sro,reduction.type = "harmony", replace.dim = T, do.print = F)
	}
return(sro)}
	
	#rotmatrix = diag(0,length(vr))
	#for(i in 2:length(vr)){
	#	c = sqrt(1.0 - vr[i] * vr[i])
	#	trow <- c * rotmatrix[,1] - s *rotmatrix[,i]
	#	rotmatrix[,i] <-  c * rotmatrix[,i] + s *rotmatrix[,i]
 	#	rotmatrix[,1] <- trow
	#}

RegressOut <- function(sro, posneg){
	posnegvec <- (2 * posneg - 1)
	vr <- posnegvec %*% sro@dr$pca@cell.embeddings
	vr <- vr / norm(vr, "2")
	out <- sro@dr$pca@cell.embeddings
	for(i in 1:dim(out)[1]){
		k <- vr %*% out[i,];
		out[i,] <- out[i,] - t(vr) * k[1,1]
	}
	sro@dr$pca@cell.embeddings <- out
return(sro)}


# assumes NormalizeSeurat was run first
DoItAllSeurat <- function(sro, reduction.type= "pca", resolution=0.8, k=30, nb.pc=0,  do.diffusion =F){
	library("Seurat")
	if (sum(dim(sro@data) == dim(sro@raw.data)) != 2) stop("Run 'NormalizeSeurat' first for Normalization and Dimentionnality reduction")
	if ((reduction.type %in% names(sro@dr))&&(dim(sro@dr[[reduction.type]]@cell.embeddings)[2] >= nb.pc)) { #PCA already there! good
	}else stop(paste("Missing \"",reduction.type ,"\"Dimention reduction, Run 'NormalizeSeurat' or RunPCA first",sep=""))
	if (nb.pc == 0) nb.pc <- dim(sro@dr[[reduction.type]]@cell.embeddings)[2]
	sro <- FindClusters(sro, reduction.type=reduction.type,dims.use=1:(nb.pc),resolution=resolution, k.param=k,force.recalc=TRUE)
	name <- paste("res.", 0.8,sep="")
	tmptmp <- as.factor(sro@meta.data[[name]])
	sro@meta.data[[name]] <- factor(as.character(tmptmp), levels=as.character(0:(length(levels(tmptmp))-1)))
	sro@ident <- sro@meta.data[[name]]

	sro <- RunTSNE(sro, reduction.use=reduction.type, dims.use=1:(nb.pc),check_duplicates=FALSE,force.recalc=TRUE)
	if (do.diffusion){
		sro <- RunDiffusion(sro, dims.use=1:(nb.pc), reduction.use=reduction.type)
	}
return(sro)}

# assumes NormalizeSeurat was run first
DoThatMonocleRun <- function(input, vargenelist, doSizeFactor = T){
	library("monocle")
	scaling = T
	if (class(input) == "SingleCellExperiment"){
		if (length(vargenelist) == 1){
			d <- t(reducedDims(input)[[vargenelist]])
			rownames(d) <- paste("DR_", 1:nrow(d),sep="")
			geneNames <- rownames(d)
			norm_method = "none"
			scaling=F
		}else stop("gene subset to do")
	}else{

	if ((class(input) != "seurat")&(class(input) != "Seurat")) stop("Unknown class for input, expected a seurat object")
		library(Seurat)
	
		if (length(vargenelist) == 1){
			if ("dr" %in% names(attributes(sro))) drslot = "dr"
			else drslot = "reductions"
			d <- t(slot(sro,drslot)[[vargenelist]]@cell.embeddings)
			geneNames <- rownames(d)
			norm_method = "none"
			scaling=F
		}else{
			d <- input@data[which(rownames(input@data) %in% vargenelist), ]
			if (class(d) != "matrix"){
				d <- as.matrix(d)
			}
			d <- d[!duplicated(rownames(d)), ]
			colnames(d) <- rownames(input@meta.data)
			geneNames <- vargenelist
			rownames(d) <- vargenelist
			norm_method= "log"
		}
	}
	#pd <- data.frame(timepoint = cellLabels)
	#pd <- new("AnnotatedDataFrame", data=pd)
	fd <- new("AnnotatedDataFrame", data= data.frame(gene_short_name = geneNames))
	rownames(fd) <- geneNames
	print("Step 1")
	dCellData <- newCellDataSet(d, featureData = fd, expressionFamily = tobit())
	print("Step 2")
	dCellData <- setOrderingFilter(dCellData, which(geneNames %in% vargenelist))
	print("Step 3")
	if (doSizeFactor) dCellData <- estimateSizeFactors(dCellData)
	print("Step 4")
	dCellDataSet <- reduceDimension(dCellData, max_components = 2, pseudo_expr = 1, norm_method=norm_method,scaling=scaling)
	print("Step 5")
	dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
	print("Step 6")
	plot_cell_trajectory(dCellDataSet)
	print("DONE")
	if (class(input) == "SingleCellExperiment"){
		colData(input)[["monocle_state"]] <- phenoData(dCellDataSet)$State
		colData(input)[["monocle_time"]] <- phenoData(dCellDataSet)$Pseudotime
		for(i in colnames(colData(ssce))) pData(dCellDataSet)[[i]] <- colData(ssce)[[i]]
		tmp2 <- input
	}else{
		smeta <- data.frame(data=phenoData(dCellDataSet)$State)
		rownames(smeta) <- rownames(phenoData(dCellDataSet))
		colnames(smeta) <- c("monocle_state")
		tmp <- AddMetaData(input,metadata=smeta)
		smeta <- data.frame(data=phenoData(dCellDataSet)$Pseudotime)
       		rownames(smeta) <- rownames(phenoData(dCellDataSet))
		colnames(smeta) <- c("monocle_time")
	        tmp2 <- AddMetaData(tmp,metadata=smeta)
	}
return(list(sro=tmp2,monocle=dCellDataSet))}
PlotMonocle <- function(monocle){
	if ((class(monocle) == list)&('monocle' %in% names(monocle))) plot_cell_trajectory(monocle$monocle)
	else plot_cell_trajectory(monocle)
}
getPseudotimeOrdering <-function(input, statefilter =""){
	if (!class(input) == "CellDataSet") stop("Expected monocle CellDataSet class as input")
	return(intersect(order(phenoData(input)$State,phenoData(input)$Pseudotime), grep(statefilter,phenoData(input)$State)))
}
PrePlotProcess <- function(input, zero.as.na=F, transform="none", row.normalize="none",col.normalize="none", row.resample=0,col.resample=0, row.reorder=c(), col.reorder=c()){
	data <- as.matrix(input)
	if (!is.null(row.reorder)) {data <- data[row.reorder,]}
  if (!is.null(col.reorder)) {data <- data[,col.reorder]}
	if (zero.as.na){data[data ==0] <- NA}
	if (transform=="log"){
		data <- log(data)
	}
	if (col.normalize == "erf"){
		for(i in 1:ncol(input)){
			nf <- !is.na(data[,i])
			m <- mean(data[nf,i])
			v <- var(data[nf,i]) ^-0.5
			if (is.na(v)||(v == 0.0)||(is.infinite(v))){
				data[nf,i] <- 0
			}else{
				data[nf,i] <- 2*pnorm((data[nf,i] - m)*v)-1
			}
		}
	}
  if (row.normalize == "erf"){
    for(i in 1:nrow(input)){
      nf <- !is.na(data[i,])
      m <- mean(data[i,nf])
      v <- var(data[i,nf])  ^-0.5
      if (is.na(v)||(v == 0)||is.infinite(v)){
              data[i,nf] <- 0
      }else{
              data[i,nf] <- 2*pnorm((data[i,nf] - m)*v)-1
      }
    }
  }
	if ((col.resample != 0)&&(col.resample < ncol(data))){
		outdata <-matrix(0,nrow(data), col.resample)
		nlabel <- c()
		f <- 1;
		ind = c()
		for(i in 1:col.resample){
			nlabel <- c(nlabel, colnames(data)[f])
			e <- floor(i * ncol(data) / col.resample)
			if (e != f) outdata[,i] = rowMeans(data[, f:e ], na.rm=T)
			else outdata[,i] = data[, e]
			f <- e +1
		}
		colnames(outdata) <- nlabel
		rownames(outdata) <- rownames(data)
		data <- outdata
	}
	if ((row.resample != 0)&&(row.resample < nrow(data))){
                outdata <-matrix(0,row.resample,ncol(data))
                nlabel <- c()
		f <- 1;
                for(i in 1:row.resample){
			nlabel <- c(nlabel, rownames(data)[f])
                        e <- i * nrow(data) / row.resample
                        if (e !=f)outdata[i,] = rowMeans(data[ f:e,], na.rm=T)
			else outdata[i,] = data[e,]
                        f <- e +1
                }
		rownames(outdata) <- nlabel
		colnames(outdata) <- colnames(data)
                data <- outdata
        }
return(data)}
getMarkerInOrdering <- function(input, out.number = 32, col.reorder= c(), min.na=0){
	perminput <- t(input[,col.reorder])
	w <- perminput
	w[is.na(perminput)] <- 0
	w[!is.na(perminput)] <- 1
	if (min.na !=0){
		print(sum(colSums(is.na(perminput)) >= min.na) )
		perminput[,(colSums(is.na(perminput)) < min.na)] <- 0
		#print(perminput)
	}
	perminput[is.na(perminput)] <- 0
	return(rownames(input)[order(abs(lm.fit( perminput, 1:length(col.reorder) )$coefficients),decreasing=T,na.last=F)[1:out.number]] )
}

getFileExtension <- function(path){
	tmp <- basename(path)
	return(substr(tmp,regexpr('\\.',tmp),nchar(tmp)))
}

# predecated
makeAnnotation <- function(input, col.resample=0,col.reorder=c()){
  if (!class(input) == "CellDataSet") stop("Expected monocle CellDataSet class as input")

	if (col.resample==0) {
		if (!is.null(col.reorder)) filt = col.reorder
		else filt = 1:length(phenoData(input)$State)}
	else{
		#filt = floor(( 0:(col.resample-1) * (length(phenoData(input)$State)-1)) /(col.resample-1)) +1
		filt = (floor(0:(col.resample-1) * length(phenoData(input)$State) / col.resample))+1
		if (!is.null(col.reorder)){
			 filt <- col.reorder[filt]
		}
	}
	annot <- data.frame(Var1 = phenoData(input)$State[filt] )
	rownames(annot) <- colnames(input)[filt]
	colnames(annot)	<- c("State")

#	colorss <- c("#FF0000", "#AAAA00","#00C000","#0040FF","#FF00FF" )
 #       names(colorss) <- c("1","2","3","4","5")
       colorss <- rainbow(5)
        names(colorss) <- c("1","2","3","4","5")

	return(list(annotation_col=annot, annotation_colors =list(State =colorss)))
}

appendAnnotation <- function(input, annotation, annot.name , is.rowannotation=T, row.names=NULL, do.sort=F, labels= c()){
	library("gtools")
	#tmp <- annotation
	#tmp[annotation == ""] <- NA

	if (do.sort){
		if (is.na(as.numeric(unique(annotation)[1]))) annotation <- factor(annotation, levels=unique(annotation)[order(as.character(unique(annotation)))])
		else annotation <- factor(annotation, levels=unique(annotation)[order(as.numeric(unique(annotation)))])
	}
	

	newannot <- data.frame(Var1 = annotation)
	if (is.null(row.names)) rownames(newannot) <- rownames(annotation)
	else rownames(newannot) <- row.names
    colnames(newannot) <- c(annot.name)

	if (is.rowannotation){
		if ('annotation_row' %in% names(input)) input$annotation_row <- cbind(input$annotation_row, newannot)
		else input$annotation_row <- newannot
	}else{
		if ('annotation_col' %in% names(input)) input$annotation_col <- cbind(input$annotation_col, newannot)
		else input$annotation_col <- newannot
	}
	if (length(labels)==0){
		tmp <- unique(sort(annotation, na.last = T))
		#print(tmp)
		#tmp <- tmp[!is.na(tmp)]
		colorss <- myrainbow(length(tmp))
		#print(colorss)
		#print(length(tmp))
		#names(colorss) <- attributes(annotation)$levels[attributes(annotation)$levels != ""]

		if (class(annotation) == "factor"){
					if (do.sort){
						if (is.na(as.numeric(levels(annotation)[1]))) annotation <- factor(annotation, levels=levels(annotation)[order(as.character(levels(annotation)))])
						else annotation <- factor(annotation, levels=levels(annotation)[order(as.numeric(levels(annotation)))])
					}
				names(colorss) <- attributes(annotation)$levels
		}else{
			names(colorss) <- unique(annotation)
		}
		#print(colorss)
		input$annotation_colors[[annot.name]] <- colorss
	}else{
		colorss <- myrainbow(length(labels))
		if (class(annotation) == "factor"){
				if (do.sort){
					if (is.na(as.numeric(levels(annotation)[1]))) annotation <- factor(annotation, levels=levels(annotation)[order(as.character(levels(annotation)))])
					else annotation <- factor(annotation, levels=levels(annotation)[order(as.numeric(levels(annotation)))])
				}				
				names(colorss) <- attributes(annotation)$levels
		}
		names(colorss) <- labels

		#print(colorss)
		input$annotation_colors[[annot.name]] <- colorss

		if (is.rowannotation){
			input$annotation_row[is.na(match(input$annotation_row[,c(annot.name)], labels)),c(annot.name)] <- NA
		}else{
			input$annotation_col[is.na(match(input$annotation_col[,c(annot.name)], labels)),c(annot.name)] <- NA
		}

	}
	return(input)
}

makeDefaultAnnotations <- function(pair){
	annot_ct <- readAnnotationFile("~/neuron_top500_markers_full.tsv")
	annot_cy <- readAnnotationFile("~/cell_cycle_list_full.tsv")

	monoannot <- makeAnnotation(pair$monocle)
	monoannot <- appendAnnotation(monoannot, annot_cy[,"stage"], "stage",row.names=rownames(annot_cy))
	monoannot <- appendAnnotation(monoannot, annot_ct[,"type"], "type", row.names=rownames(annot_ct))
	monoannot <- appendAnnotation(monoannot, pair$sro@meta.data[,"res.0.8"], "Cluster",is.rowannotation=F )
return(monoannot)}

readAnnotationFile <- function(path){
	ext <- getFileExtension(path)
	if (ext == ".csv") data <- read.csv(file=path,row.names=1)
	else data <- read.csv(file=path,row.names=1,sep='\t')
return(data)}

makeDefaultAnnotationsSeurat <- function(sro, use.cycle=T,use.neuromarks=T,sample.color=c(),cell.labels=c("orig.ident"), do.sort.label=c(T)){
	annot <- list(annotation_col = data.frame(row.names=colnames(sro)))
	if (use.cycle){		
		annot_cy <- readAnnotationFile("/nfs/users/nfs_l/lh20/cell_cycle_list_full.tsv")
		annot <- appendAnnotation(annot, annot_cy[,"stage"], "stage",row.names=rownames(annot_cy))
	}
	if (use.neuromarks){
		annot_ct <- readAnnotationFile("/nfs/users/nfs_l/lh20/neuron_top500_markers_full.tsv")
		annot <- appendAnnotation(annot, annot_ct[,"type"], "type", row.names=rownames(annot_ct))
	}

	if (length(do.sort.label) != length(cell.labels)) do.sort.label <-rep(F,length(cell.labels))
	if (length(cell.labels) > 0) {
		for(i in 1:length(cell.labels)){
			annot <- appendAnnotation(annot, sro@meta.data[,cell.labels[i]], cell.labels[i],is.rowannotation=F, do.sort=do.sort.label[i])
			if (cell.labels[i] %in% names(sro@misc$meta.color)){
				print("trying out colors!")
				annot$annotation_colors[[cell.labels[i]]]  <- sro@misc$meta.color[[cell.labels[i]]]
			}
		}
	}else annot <- appendAnnotation(annot, sro@ident, "Cluster",is.rowannotation=F, do.sort=T)

	if (length(sample.color) != 0){
		names(sample.color) <- levels(sro@meta.data[,"orig.ident"])
		annot$annotation_colors[["orig.ident"]] <- sample.color
	}
return(annot)}


# outputs 1 tabdelimited file, containing top markers and all markers
exportSeuratAllMarkerStruct <- function(allmarkerstruct, tsvpath, top_pos = 100, top_neg = 10){
	cat("Top ", top_pos ," Significant Positive Gene Markers per cluster:\n", file = tsvpath, append=FALSE)
	for(i in 1:length(levels(allmarkerstruct[[6]]))){
		cat(paste("\t", levels(allmarkerstruct[[6]])[i], sep=''), file = tsvpath, append=TRUE)
		sub <- ((allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]] > 0))
		l <- sum(sub)
		if (l >= top_pos){
			if (i == 1) topdf <- data.frame(allmarkerstruct[[7]][sub][1:top_pos])
			else topdf <- cbind(topdf, data.frame(allmarkerstruct[[7]][sub][1:top_pos]))
		}else{
			if (i == 1) topdf <- data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_pos-l) ) )
			else topdf <- cbind(topdf, data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_pos-l) )))
		}
	}
	cat("\n", file = tsvpath, append=TRUE)
	write.table(topdf, tsvpath,sep = "\t",col.names=F, append = TRUE)
	cat("\nTop ", top_pos ," High Foldchange Positive Gene Markers per cluster:\n", file = tsvpath, append=TRUE)
	for(i in 1:length(levels(allmarkerstruct[[6]]))){
		cat(paste("\t", levels(allmarkerstruct[[6]])[i], sep=''), file = tsvpath, append=TRUE)
		sub <- ((allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]] > 0))
		sub <- sub[order(-allmarkerstruct[[2]])]
		l <- sum(sub)
		if (l >= top_pos){
			if (i == 1) topdf <- data.frame(allmarkerstruct[[7]][sub][1:top_pos])
			else topdf <- cbind(topdf, data.frame(allmarkerstruct[[7]][sub][1:top_pos]))
		}else{
			if (i == 1) topdf <- data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_pos-l) ) )
			else topdf <- cbind(topdf, data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_pos-l) )))
		}
	}
	cat("\n", file = tsvpath, append=TRUE)
	write.table(topdf, tsvpath,sep = "\t",col.names=F, append = TRUE)
	cat("\nTop", top_neg ," Negative Gene Markers per cluster:\n", file = tsvpath, append=T)
	for(i in 1:length(levels(allmarkerstruct[[6]]))){
		cat(paste("\t", levels(allmarkerstruct[[6]])[i], sep=''), file = tsvpath, append=TRUE)
		sub <- ((allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]] < 0))
		l <- sum(sub)
		if (l >= top_neg){
			if (i == 1) topdf <- data.frame(allmarkerstruct[[7]][sub][1:top_neg])
			else topdf <- cbind(topdf, data.frame(allmarkerstruct[[7]][sub][1:top_neg]))
		}else{
			if (i == 1) topdf <- data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_neg-l) ) )
			else topdf <- cbind(topdf, data.frame(c ( allmarkerstruct[[7]][sub][1:l], rep("",top_neg-l) )))
		}
	}

	cat("\n", file = tsvpath, append=TRUE)
	write.table(topdf, tsvpath,sep = "\t",col.names=F, append = TRUE)
	cat("\nAll Gene Positive Markers\n\t", file = tsvpath, append=TRUE)
	write.table(allmarkerstruct[(allmarkerstruct[[2]] > 0),], tsvpath,sep = "\t",append = TRUE)
	cat("\nAll Gene Negative Markers\n\t", file = tsvpath, append=TRUE)
	write.table(allmarkerstruct[(allmarkerstruct[[2]] > 0),], tsvpath,sep = "\t",append = TRUE)
}

listtopinAllMarkerstruct <- function(allmarkerstruct,  markperclass=10){
  genelist = c();
  for(i in 1:length(levels(allmarkerstruct[[6]]))){
	nbpos <- sum((allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]]@.Data > 0))
	if (nbpos < markperclass) {
		genelist <- c(genelist,allmarkerstruct[[7]][(allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]]@.Data > 0) ])
	}else{
		genelist <- c(genelist,allmarkerstruct[[7]][(allmarkerstruct[[6]]@.Data == i)&(allmarkerstruct[[2]]@.Data > 0) ][1:markperclass])
	}
  }
return(genelist)}


# needs "allmarkerstruct" can be a list, is so it is assumed to be a genelist
makeDefaultMarkerHeatMap <- function(sro, allmarkerstruct = c(), markperclass=10,fontsize_row = 6, use.cycle=T,use.neuromarks=T, use.map=F, sample.color=c(), order.label= c("res.0.8", "orig.ident"),  do.sort.label=c(T,F), overpresentation.class = c(),  overpresentation.faction = 0.5,useraw=F, cluster_rows=F, use.cells=c()){
	if ((class(sro) != "seurat")&(class(sro) != "Seurat")) stop("Expects Seurat Object");
	genelist = c();
	grlen = c()
	if (is.null(use.cells)) use.cells <- rep(T, colnames(sro))
	if (is.null(allmarkerstruct)){
		if ("allmarkers" %in% names(sro@misc)){
			allmarkerstruct <- sro@misc[["allmarkers"]]
		}else stop("list of marker or Structure outputted by 'Seurat::FindMarker' required")
	}
	if (class(allmarkerstruct) == "data.frame"){
		for(i in levels(allmarkerstruct[[6]])){
			nbpos <- sum((allmarkerstruct[[6]] == i)&(allmarkerstruct[[2]]@.Data > 0))
			if (nbpos < markperclass) {
				genelist <- c(genelist,allmarkerstruct[[7]][(allmarkerstruct[[6]] == i)&(allmarkerstruct[[2]]@.Data > 0) ])/
				grlen = c(grlen,rep(i, nbpos))
			}else{
				genelist <- c(genelist,allmarkerstruct[[7]][(allmarkerstruct[[6]] == i)&(allmarkerstruct[[2]]@.Data > 0) ][1:markperclass])
				grlen = c(grlen,rep(i, markperclass))
			}
		}
	}else if (class(allmarkerstruct) == "list"){
		genelist <- as.character(allmarkerstruct)
	}else if (class(allmarkerstruct) == "character"){
		if (length(allmarkerstruct) == 1) genelist <- parselist(allmarkerstruct)
		else genelist <- allmarkerstruct
	}else stop(parse("not legal class for input: ", class(allmarkerstruct)))
	print("heh")
	
	daorder <- switch(as.character(length(order.label)),
"0"= order(as.numeric(sro@ident@.Data[use.cells]),sro@ident),
"1"= order(sro@meta.data[use.cells,order.label[1]]),
"2"= order(sro@meta.data[use.cells,order.label[1]],sro@meta.data[use.cells,order.label[2]]),
"3"= order(sro@meta.data[use.cells,order.label[1]],sro@meta.data[use.cells,order.label[2]],sro@meta.data[use.cells,order.label[3]]),
"4"= order(sro@meta.data[use.cells,order.label[1]],sro@meta.data[use.cells,order.label[2]],sro@meta.data[use.cells,order.label[3]],sro@meta.data[use.cells,order.label[4]]))

	print("heh")

	if (useraw) {
		map <- !is.na(match(genelist, SeuratRawData(sro, "R")))
		grlen <- grlen[map]
		genelist <- genelist[map]
		print(range(daorder))
		print(sum(use.cells))
		print(dim(SeuratRawData(sro)))
		print(dim( SeuratRawData(sro)[genelist,use.cells]))
		sdata <- PrePlotProcess( SeuratRawData(sro)[genelist,use.cells], row.normalize="erf", col.resample=1024, col.reorder=daorder)
	}else{
		map <- !is.na(match(genelist, rownames(sro)))
		grlen <- grlen[map]
		genelist <- genelist[map]
		sdata <- PrePlotProcess(sro@data[genelist,use.cells], row.normalize="erf", col.resample=1024, col.reorder=daorder)
	}
	print("heh")
	annotation <- makeDefaultAnnotationsSeurat(sro, use.cycle=use.cycle,use.neuromarks=use.neuromarks,sample.color=sample.color, cell.labels=order.label, do.sort.label=do.sort.label)
	printf("gothere")
	annot <- matrix(NA,length(rownames(annotation$annotation_row)),1)
	rownames(annot) <- rownames(annotation$annotation_row)
	if (use.map){
		annotation <- appendAnnotation(annotation, sro@meta.data[,"scmap_MDS"], "Mapped",is.rowannotation=F)
	}

	if (length(grlen) != 0){
		if (length(genelist) > 0){
			for(i in 1:length(genelist)){
				if (genelist[i] %in% rownames(annot)){
					if (is.na(annot[genelist[i],1])) annot[genelist[i],1] <- grlen[i]
					else annot[genelist[i],1] <- "Multiple"
				}
			}
			annotation <- appendAnnotation(annotation, annot, "MarkCluster", row.names=rownames(annot), labels = levels(sro@ident))
		}else stop("nogenes!!")
	}

	PlotHeat(sdata,annotation=annotation,cluster_rows=cluster_rows,fontsize_row=fontsize_row)
return(genelist);}




makeTopMarkerList <- function(allmarkerstruct, nbmark = 10){
	#library("gProfileR")
	genelist <- c()
	tmp <- allmarkerstruct[[6]]
	for(i in levels(allmarkerstruct[[6]])){
		nbpos <- sum(tmp == i)
		if (nbpos < nbmark) genelist <- c(genelist, allmarkerstruct[[7]][tmp == i])
		else genelist <- c(genelist, allmarkerstruct[[7]][tmp == i][1:nbmark])
	}
	return(genelist)
}
makeDefaultHeat <- function(pair, genelist, cluster_rows=T){
    annot <- makeDefaultAnnotations(pair)
    order <- getPseudotimeOrdering(pair$monocle)
    data <- PrePlotProcess(pair$sro@data[genelist,], row.normalize="erf",col.resample=1000,col.reorder=order)
    library(pheatmap)
    pheatmap(data,color= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)
    ,cluster_rows=cluster_rows,cluster_cols=F	,annotation_names_row=F,annotation_names_col=F		,annotation_colors=annot$annotation_colors, annotation_row=annot$annotation_row, annotation=annot$annotation_col,show_rownames=T,show_colnames=F)

#PrePlotProcess <- function(input, zero.as.na=F, transform="none", row.normalize="none",col.normalize="none", row.resample=0,col.resample=0, row.reorder=c(), col.reorder=c()){
    return(annot)
}


PlotHeat <- function(data,annotation=list(),cluster_rows=T,cluster_cols=F,...){
	library(pheatmap)
	data[is.na(data)] <- 0
	pheatmap(data
		,color= colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)
		,cluster_rows=cluster_rows,cluster_cols=cluster_cols
		,annotation_names_row=F,annotation_names_col=T
		,show_rownames=T,show_colnames=F,
		,annotation_colors=annotation$annotation_colors
		,annotation_col=annotation$annotation_col
		,annotation_row=annotation$annotation_row, legend_breaks= c(1),fontsize = 6
		,...
		)
}

myFlFormat <- function(value){
	str1 <- format(value,digits=2,scientific=T)
	str2 <- format(value,digits=3)
	if (nchar(str1) <nchar(str2)) return(str1)
	else return(str2)
}

myCustomDatagrid <- function(table, rowselect, ctselect, batchselect, which= "DEseq", nbrows=0,decreasing=F, special.x=c()){
	input <- list();
	colselect <- paste(rep(ctselect,each =length(batchselect)),rep(batchselect,length(ctselect)),sep="_")
	if (which == "DEseq"){
		input$data <- table$DEseq$deseq.log2FC[rowselect,colselect]
		input$w <- table$DEseq$deseq.log10pvalue[rowselect,colselect]
		
	}else{
		input$data <- 1.0 - 2.0 / (1 + exp(table$Wilcox$wilcox.logitAuroc[rowselect,colselect]))
		input$w <- table$Wilcox$wilcox.log10pval[rowselect,colselect]
	}
	input$x <- ((table$log2FC_FAD_to_Empty[rowselect,colselect] > 0)|(table$log2FC_CTRL_to_Empty[rowselect,colselect] > 0))
	input$x <- input$x & ((table$Wilcox$dropoutPosClass[rowselect,colselect] > 0.1)|(table$Wilcox$dropoutNegClass[rowselect,colselect] > 0.1))
	
	input$y <- matrix(0, length( rowselect), length(colselect))

	if (!is.null(special.x)){
		toord <- input$data
		toord[!input$x] <- NA
		if (length(special.x) != 1) toord <- toord[,special.x]

		ord <- order(rowSums(toord,na.rm=T),decreasing=decreasing)
	}else ord <- order(rowSums(input$data),decreasing=decreasing)
	if ((nbrows != 0)&(nbrows< length(ord))) ord <- ord[1:nbrows]
	input$data <- input$data[ord,]
	input$w <- input$w[ord,]
	input$x <- input$x[ord,]
	input$y <- input$y[ord,]

return(plotDataGrid(data=input, transform=list(w= "log10pval"),bgcolor=c("#FFFFFF")))} # #008800

plotDataGrid <- function(data, wdata= c(), xdata = c(), ydata =c(), transform=c(), plot.attribs=c(),do.zero.center=T, bgcolor = "#BBBBBB"){
	library(ggplot2)
	if (class(data) != "list") {
		data <- list(data=data)
		if (is.null(wdata)) data$w <-matrix(1, dim(data$data)[1], dim(data$data)[2])
		else data$w <- dataw
		if (is.null(xdata)) data$x <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		else data$x <- datax
		if (is.null(ydata)) data$y <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		else data$y <- datay
	}else{
		if (!"data" %in% names(data)) stop("input list has 'data' as mandatory field")
		for( i in names(data)){
			if (!i %in% c("data", "x", "y", "w")) print(paste("The nknown field",i," is ignored, valid fields are \"data\", \"x\", \"y\", \"w\" only"))
		}
		
		if (!"w" %in% names(data)) data$w <- matrix(1, dim(data$data)[1], dim(data$data)[2])
		if (!"x" %in% names(data)) data$x <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		if (!"y" %in% names(data)) data$y <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
	}
 	dd <- dim(data$data)
	cliprect <- c(0,0,dd[2],dd[1])


	if (is.null(plot.attribs)) plot.attribs <- list(flags=c())
	if ("flags" %in% names(plot.attribs)) pflags <- plot.attribs[["flags"]]
	else pflags <- c()
	
	trformval <-c()
 	aurange <- range(as.vector(data$data),na.rm=T);
        if (do.zero.center){
                if (abs(aurange[1]) > abs(aurange[2])) aurange[2] <- -aurange[1]  
                else aurange[1] <- -aurange[2]
        }
	for(i in names(transform)){
		if (!i %in% c("data", "x", "y", "w")) stop("allowed transform fiels are \"data\", \"x\", \"y\", \"w\" only")
		if (transform[i] == "pval"){
			data[[i]] <- -0.8 + 0.09 / (0.05 + data[[i]])
			data[[i]][data[[i]] < 0.1] <- 0.1
		}else if (transform[i] == "log10pval"){
			data[[i]] <- 0.05 / (0.05 + exp(data[[i]]/log(10)))
		}else if (transform[i] == "lerfp1"){
			tmp <- as.matrix(log(data[[i]] + 1))
			m <- mean(tmp,na.rm=T)
			v <- var(as.vector(tmp),na.rm=T) ^-0.5
			if (is.na(v)){
				data[[i]] <- 0.5
			}else{
				data[[i]]  <- pnorm((tmp - m)*v)
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "colwise.erf"){
			trformval <- matrix(0, 5, dd[2])
			for(j in 1:dd[2]){
				tmp <- data[[i]][,j]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][,j] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][,j]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "rowwise.erf"){
			trformval <- matrix(0, 5, dd[1])
			for(j in 1:dd[2]){
				tmp <- data[[i]][j,]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][j,] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][j,]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "threshold"){
			flt <- (data[[i]] <= 0.05)
			flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 1;
			data[[i]][!flt] <- 0;
		}else if (transform[i] == "threshold2"){
			flt <- (data[[i]] <= 0.05)
			flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 0;
			data[[i]][!flt] <- 1;
		}

	}

	fgdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 8))
	bgdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 8))
  for(j in 1:dd[2]){
    for(i in 1:dd[1]){
	offset <- (i-1+ (j-1) * dd[1])*8
	fgdata$C[(offset+1):(offset+8)] <- rep(data$data[i,j],8)
	bgdata$C[(offset+1):(offset+4)] <- rep(0,4)
	bgdata$C[(offset+5):(offset+8)] <- rep(1,4)

	
	fgdata$I[(offset+1):(offset+8)] <- rep(offset/8,8)
	bgdata$I[(offset+1):(offset+4)] <- rep(offset/4,4)
	bgdata$I[(offset+5):(offset+8)] <- rep(1+offset/4,4)
	radata <- c(1 - sqrt(data$w[i,j]), 1 - sqrt(data$w[i,j]/2)) /2
	fgdata$X[offset+1] <- j - 1.0 + radata[1];		fgdata$Y[offset+1] <- i - 0.5;
 	fgdata$X[offset+2] <- j - 1.0 + radata[2];		fgdata$Y[offset+2] <- i - radata[2];
	fgdata$X[offset+3] <- j - 0.5;			fgdata$Y[offset+3] <- i - radata[1];
	fgdata$X[offset+4] <- j - radata[2];	fgdata$Y[offset+4] <- i - radata[2];
	fgdata$X[offset+5] <- j - radata[1];	fgdata$Y[offset+5] <- i - 0.5;
	fgdata$X[offset+6] <- j - radata[2]; 	fgdata$Y[offset+6] <- i - 1.0 + radata[2];
	fgdata$X[offset+7] <- j - 0.5;			fgdata$Y[offset+7] <- i -1.0 + radata[1];
	fgdata$X[offset+8] <- j - 1.0 + radata[2];		fgdata$Y[offset+8] <- i - 1.0 + radata[2];
	
	bgdata$X[offset+1] <- j -1.0;			bgdata$Y[offset+1] <- i - 1.0;
 	bgdata$X[offset+2] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+2] <- i - 1.0;
	bgdata$X[offset+3] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+3] <- i - 1.0 + data$y[i,j];
	bgdata$X[offset+4] <- j -1.0;			bgdata$Y[offset+4] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+5] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+5] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+6] <- j;		 	bgdata$Y[offset+6] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+7] <- j;			bgdata$Y[offset+7] <- i;
	bgdata$X[offset+8] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+8] <- i;
    }
  }
 


	p <- ggplot(data = fgdata,mapping=aes(fill= C, group = I, y=Y, x=X) )
	p <- p + theme(axis.text.x=element_text(angle=90,vjust=0.5))
	p <- p + scale_x_discrete(limits= (1:dd[2])-0.5, labels= colnames(data$data) )# + xlab(NULL)
	p <- p + scale_y_discrete(limits= (1:dd[1])-0.5, labels= rownames(data$data) )# + ylab(NULL) 
	p <- p + geom_polygon(data=bgdata, mapping=aes(group = I, y=Y, x=X), fill = rep(c(bgcolor), dd[1] * dd[2]*8))
	
	if (!is.null(trformval)){
		newbot = cliprect[2]  - (cliprect[4] / 10)
		cbdata <- data.frame(row.names = 1:164)
		for(i in 1:41){
			cbdata$X[i *4 -3] = 0; cbdata$X[i *4 -2] = 0; cbdata$X[i *4 -1] = cliprect[3]; cbdata$X[i *4] = cliprect[3];
			cbdata$Y[i *4 -3] = (newbot * i)/41; cbdata$Y[i *4 -2] = (newbot * (i-1))/41; cbdata$Y[i *4 -1] = (newbot * (i-1))/41; cbdata$Y[i *4] = (newbot * i)/41;
			cbdata$I[i *4 -3] = i; cbdata$I[i *4 -2] = i; cbdata$I[i *4 -1] = i; cbdata$I[i *4] = i;
			cbdata$C[i *4 -3] = (41-i)/40; cbdata$C[i *4 -2] = (41-i)/40; cbdata$C[i *4 -1] = (41-i)/40; cbdata$C[i *4] = (41-i)/40;
		}
		p <- p + geom_polygon(data=cbdata, mapping=aes(group = I, y=Y, x=X, fill=C))
		cliprect[2] = newbot
		cbdata <- data.frame(row.names = 1:(5*dd[2]))
		cbdata$C = 1:(5*dd[2])
		for(i in 1:dd[2]){
			cbdata$X[i *5 -4] = i-0.5; cbdata$X[i *5 -3] = i-0.5; cbdata$X[i *5 -2] = i-0.5; cbdata$X[i *5 -1] = i-0.5; cbdata$X[i *5] = i-0.5;
			cbdata$Y[i *5 -4] = newbot * 0.1; cbdata$Y[i *5 -3] = newbot * 0.3; cbdata$Y[i *5 -2] = newbot * 0.5; cbdata$Y[i *5 -1] = newbot * 0.7; cbdata$Y[i *5] = newbot * 0.9;
			cbdata$I[i *5 -4] = myFlFormat(trformval[1,i]); cbdata$I[i *5 -3] = myFlFormat(trformval[2,i]); cbdata$I[i *5 -2] = myFlFormat(trformval[3,i]); cbdata$I[i *5 -1] = myFlFormat(trformval[4,i]); cbdata$I[i *5] = myFlFormat(trformval[5,i]);
		}
		p <- p + geom_text(data=cbdata, mapping=aes(label = I, y=Y, x=X), color=rep(c("#FFFFFF"), dd[2]*5),size=rep(2,dd[2]*5))
	}


  	p <- p + geom_polygon()
	p <- p + theme(axis.text=element_text(color="#000000",face="bold",size=10))
	p <- p + scale_fill_gradientn(colours=c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"),limits=aurange)
	p <- p + scale_size(limits=c(0,1))
	p <- p + coord_cartesian(xlim=c(cliprect[1],cliprect[3]),ylim=c(cliprect[2],cliprect[4]),expand=F)
	return(changeStyle(p,plot.attribs = plot.attribs))
}

plotCollumnData <- function(data, transform=c(), plot.attribs=c()){
	library(ggplot2)
 	dd <- dim(data)

	if (is.null(plot.attribs)) plot.attribs <- list(flags=c())
	if ("flags" %in% names(plot.attribs)) pflags <- plot.attribs[["flags"]]
	else pflags <- c()
	
	for(i in names(transform)){
		if (transform[i] == "pval"){
			data[,i] <- -0.8 + 0.09 / (0.05 + data[,i])
			data[,i][data[,i] < 0.1] <- 0.1
		}else if (transform[i] == "log10pval"){
			data[,i] <- 0.05 / (0.05 + exp(data[,i]/log(10)))
		}else if (transform[i] == "lerfp1"){
			tmp <- as.matrix(log(data[,i] + 1))
			m <- mean(tmp,na.rm=T)
			v <- var(as.vector(tmp),na.rm=T) ^-0.5
			if (is.na(v)){
				data[,i] <- 0.5
			}else{
				data[,i]  <- pnorm((tmp - m)*v)
			}
		}

	}

	fgdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 8))
	bgdata <- data.frame(row.names = 1:(dd[2] * 160))

  for(j in 1:dd[2]){
    for(i in 1:dd[1]){
	offset <- (i-1+ (j-1) * dd[1])*8
	fgdata$C[(offset+1):(offset+8)] <- rep(data$data[i,j],8)

	
	fgdata$I[(offset+1):(offset+8)] <- rep(offset/8,8)
	radata <- c(1 - sqrt(data$w[i,j]), 1 - sqrt(data$w[i,j]/2)) /2
	fgdata$X[offset+1] <- i + radata[1];		fgdata$Y[offset+1] <- j + 0.5;
 	fgdata$X[offset+2] <- i + radata[2];		fgdata$Y[offset+2] <- j + 1.0 - radata[2];
	fgdata$X[offset+3] <- i + 0.5;			fgdata$Y[offset+3] <- j + 1.0 - radata[1];
	fgdata$X[offset+4] <- i + 1.0 - radata[2];	fgdata$Y[offset+4] <- j + 1.0 - radata[2];
	fgdata$X[offset+5] <- i + 1.0 - radata[1];	fgdata$Y[offset+5] <- j + 0.5;
	fgdata$X[offset+6] <- i + 1.0 - radata[2]; 	fgdata$Y[offset+6] <- j + radata[2];
	fgdata$X[offset+7] <- i + 0.5;			fgdata$Y[offset+7] <- j + radata[1];
	fgdata$X[offset+8] <- i + radata[2];		fgdata$Y[offset+8] <- j + radata[2];
    }
  }
 


  aurange <- range(as.vector(data));
	p <- ggplot(data = fgdata,mapping=aes(fill= C, group = I, y=Y, x=X) )
	p <- p + theme(axis.text.x=element_text(angle=90,vjust=0.5))
	p <- p + scale_x_discrete(limits= (1:dd[1])+0.5, labels= rownames(data$data) )# + xlab(NULL)
	p <- p + scale_y_discrete(limits= (1:dd[2])+0.5, labels= colnames(data$data) )# + ylab(NULL) 
	p <- p + geom_polygon(data=bgdata, mapping=aes(group = I, y=Y, x=X), fill = rep(c("#BBBBBB"), dd[1] * dd[2]*8))
  	p <- p + geom_polygon()
	p <- p + theme(axis.text=element_text(color="#000000",face="bold",size=10))
	p <- p + scale_fill_gradientn(colours=c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"),limits=aurange,position=none)
	p <- p + scale_size(limits=c(0,1))
	p <- p + coord_cartesian(xlim=c(1,dd[1]+1),ylim=c(1,dd[2]+1),expand=F)
	return(changeStyle(p,plot.attribs = plot.attribs))
}

plotLabels <- function(xdata, ydata, names = c(), color = c(), alpha = c(), filter = c(), xlim=c(), ylim=c(),plot.attribs=c(),label.size=4){
	if (is.null(plot.attribs)) plot.attribs <- list()
	if (!is.null(xlim)) plot.attribs$xlim = xlim
	if (!is.null(ylim)) plot.attribs$ylim = ylim
	if (is.null(filter)) filter <- rep(TRUE,length(xdata))

	if ("label.size" %in% names(plot.attribs)) label.size = plot.attribs$label.size
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()

	if (is.null(names)) names <- names(xdata);
	if (is.null(color)) color <- rep("#000000", length(xdata))
	if (is.null(alpha)) alpha <- rep(1, length(xdata))


	hehe <- cbind(xdata[filter],ydata[filter],names[filter], color[filter], alpha[filter])
	rownames(hehe) <- NULL
	colnames(hehe) <- c("X","Y","Label", "Color", "Alpha")
	hehe <- data.frame(hehe)
	hehe[,1] <- as.numeric(as.character(hehe[,1]))
	hehe[,2] <- as.numeric(as.character(hehe[,2]))
	hehe[,4] <- as.numeric(as.character(hehe[,4]))
	hehe[,5] <- as.numeric(as.character(hehe[,5]))

	

	library(ggplot2)
	p <- ggplot(hehe, aes(x=X,y=Y))
	if ("xlim" %in% names(plot.attribs)) p <- p + scale_x_continuous(limits = plot.attribs$xlim)
	if ("ylim" %in% names(plot.attribs)) p <- p + scale_y_continuous(limits = plot.attribs$ylim)

	colscale <- unique(color)
	names(colscale) <- colscale
#	p <- p + scale_color_manual(name="color",values=colscale)

	p <- p + geom_text(aes(label=Label,color=Color,alpha=Alpha),size=label.size)
return(changeStyle(p,plot.attribs))}


#' Plot Inferred Network
#'
#' @param networkoutput: list generated by 'Infern0_IdentifyNetwork'
#'
#' @export
Infern0plotNetwork <- function(out,vertex.size.base=0.5,vertex.size.incr=3,edge.width=10,edge.color="#000000FF", layout= c(), vertex.height=5,label.cex=0.75, incl_isolates=F){
	library(igraph);
	if (out$nbedge == 0) {
		      range <- c();
		      lone = 1:dim(out$covar)[1];
	}else{
		      range <- 3:(out$nbedge*2 +2);
		      lone = setdiff(1:dim(out$covar)[1], (unique(as.vector(t(out$EdgeListCoor))[3:(out$nbedge*2 +2)])+1));
	}
	if (incl_isolates){
		gr <- graph(edges=rownames(out$covar)[ as.vector(t(out$EdgeListCoor + 1))[range] ], isolates = rownames(out$covar)[lone], directed = F); 
	}else{
		gr <- graph(edges=rownames(out$covar)[ as.vector(t(out$EdgeListCoor + 1))[range] ], directed = F); 
	}	
	nbnode = length(V(gr))

	if (is.null(layout)) layout <- layout_nicely(gr)

	xrange <- c( min(layout[,1]), max(layout[,1]))
	yrange <- c( min(layout[,2]), max(layout[,2]))

	xfact <- xrange[2] - xrange[1];
	yfact <- yrange[2] - yrange[1];

	diff <- out$LL_incrs[2:(out$nbedge+1),1] - out$LL_incrs[1:(out$nbedge),1];
	diff <- ((diff - min(diff)) / (2.0 * max(diff) - min(diff))) + 0.5

	par(mar=c(1,1,1,1) + 0.1)
	V(gr)$label.cex <- 0.0001
	plot(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange,vertex.size= 0,edge.width=diff * edge.width,edge.color=edge.color);
	par(new=T)
	
	V(gr)$label.cex <- label.cex
	plot(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange, vertex.size= (nchar(names(V(gr)) )*vertex.size.incr + vertex.size.base) * xfact, vertex.shape="rectangle" ,vertex.size2=rep(vertex.height * yfact,nbnode),edge.width=0, label.font=4);
}

myFindMarkerGenes <- function(sro, metaname="",clusterlist= c()){
	if (metaname == "") curclass <- as.factor(sro@ident)
	else curclass <- as.factor(sro@meta.data[[metaname]])
	if (is.null(clusterlist)) todo <- levels(curclass)
	else{
		todo <- match(clusterlist,levels(curclass) );
		todo <- levels(curclass)[todo[!is.na(todo)]]
	}
	library(InferN0)
	for(i in 1:length(todo)){
		output <- InferN0ZeroCorrectedWilcox(sro@raw.data, curclass == todo[i], curclass != todo[i])
		

	}
return(output)
}



findGoEnrichments <- function(named.input, size=0, filter.list=c(), return.raw=F,do.plot.grid=T,plot.maxrow=10){
	library(gProfileR)
	if (is.null(filter.list)) filter.list = rep(TRUE,length(named.input))
	if (size == 0) size = length(named.input) / 20
	ord <- order(named.input[filter.list])
	
	fout <- list(low.genes= names(named.input[filter.list])[ord[1:size]], high.genes=names(named.input[filter.list])[ord[sum(filter.list):(sum(filter.list)-size+1)]])
#	printlist(fout$low.genes,append.comma=F)
#	print("hehehehehehe")
#	printlist(fout$high.genes,append.comma=F)
#	return(FALSE)
	outputL <- gprofiler(fout$low.genes,custom_bg=names(named.input))
	outputL <-outputL[order(outputL$p.value),]
	outputH <- gprofiler(fout$high.genes,custom_bg=names(named.input))
	outputH <-outputH[order(outputH$p.value),]
	if (return.raw) {
		fout[["low.output"]] <- outputL
		fout[["high.output"]] <- outputH
		return(fout);
	}
	poutput <- data.frame(outputL$term.name)
	colnames(poutput) <- c("Name")
	rownames(poutput) <- outputL$term.id
	poutput[["pvalue"]] <- outputL$p.value
	poutput[["enrich"]] <- outputL$overlap.size * length(named.input) / (outputL$term.size * length(fout$low.genes) )
	fout[["low.annot"]] <- poutput
	poutput <- data.frame(outputH$term.name)
	colnames(poutput) <- c("Name")
	rownames(poutput) <- outputH$term.id
	poutput[["pvalue"]] <- outputH$p.value
	poutput[["enrich"]] <- outputH$overlap.size * length(named.input) / (outputH$term.size * length(fout$high.genes) )
	fout[["high.annot"]] <- poutput
	if (do.plot.grid){
#		library(gridExtra)
#		library(gtable)
#		library(grid)
#		pt.low <- textGrob("Low")
#		tmptable <- gtable_add_rows(tableGrob(fout[["low.annot"]][1:plot.maxrow,]), heights=(grobHeight(pt.low)+unit(5,"mm"))*2)
#		p.low <-gtable_add_grob(tmptable,pt.low,1,1,1,ncol(tmptable))
#		pt.high <- textGrob("High")
#		tmptable <- gtable_add_rows(tableGrob(fout[["high.annot"]][1:plot.maxrow,]), heights=(grobHeight(pt.high)+unit(5,"mm"))*2)
#		p.high <- gtable_add_grob(tmptable,pt.high,1,1,1,ncol(tmptable))
#		grid.newpage()
#		#grid.draw(p.low)
#		print(grid.arrange(p.low,p.high,ncol=2))
		library(gridExtra)
		mytheme <- gridExtra::ttheme_default(
    		core = list(fg_params=list(cex = 0.8)),
    		colhead = list(fg_params=list(cex = 0.8)),
    		rowhead = list(fg_params=list(cex = 0.8)))
		if (nrow(fout[["low.annot"]]) < plot.maxrow) p.low <- tableGrob(fout[["low.annot"]], theme = mytheme)
		else p.low <- tableGrob(fout[["low.annot"]][1:plot.maxrow,], theme = mytheme)
		if (nrow(fout[["low.annot"]]) < plot.maxrow) p.high <- tableGrob(fout[["high.annot"]], theme = mytheme)
		else p.high <- tableGrob(fout[["high.annot"]][1:plot.maxrow,], theme = mytheme)
		print(grid.arrange(p.low,p.high,ncol=2, top="Enriched Annotations"))
	}
return(fout);
}


funkyZscoreTest <- function(scptr,Lnn, Lnp, Lpn, Lpp, cur,ord){
fout <- list()
fout[["outCp"]] <-  InferN0ZeroCorrectedWilcox2(scptr, intersect(c(Lnp, Lpp), cur), intersect(c(Lnn, Lpn), cur), ord, use.clusterID=T,do.quartile.average=T)
fout[["outCn"]] <- InferN0ZeroCorrectedWilcox2(scptr,   setdiff(c(Lnp, Lpp), cur),   setdiff(c(Lnn, Lpn), cur), ord, use.clusterID=T,do.quartile.average=T)
fout[["outN"]] <-  InferN0ZeroCorrectedWilcox2(scptr, intersect(Lnp, cur), intersect(Lnn, cur), ord, use.clusterID=T,do.quartile.average=T)
fout[["outP"]] <-  InferN0ZeroCorrectedWilcox2(scptr, intersect(Lpp, cur), intersect(Lpn, cur), ord, use.clusterID=T,do.quartile.average=T)

return(fout)
}

findNearestNeighbors <- function(sro){
	library(RANN);
	
}

fftRtsne <- function(X, 
		     dims=2, perplexity=30, theta=0.5,
		     check_duplicates=TRUE,
		     max_iter=1000,
		     fft_not_bh = TRUE,
		     ann_not_vptree = TRUE,
		     stop_early_exag_iter=250,
		     exaggeration_factor=12.0, no_momentum_during_exag=FALSE,
		     start_late_exag_iter=-1.0,late_exag_coeff=1.0,
             mom_switch_iter=250, momentum=.5, final_momentum=.8, learning_rate=200,
		     n_trees=50, search_k = -1,rand_seed=-1,
		     nterms=3, intervals_per_integer=1, min_num_intervals=50, 
		     K=-1, sigma=-30, initialization=NULL,
		     data_path=NULL, result_path=NULL,
		     load_affinities=NULL,
		     fast_tsne_path="/lustre/scratch117/cellgen/team218/lh20/FIt-SNE/bin", nthreads=0, perplexity_list = NULL, 
                     get_costs = FALSE, df = 1.0,... ) {
        version_number = '1.1.0'

	if (is.null(fast_tsne_path)) {
		if(.Platform$OS.type == "unix") {
			fast_tsne_path = sprintf('%s/bin/fast_tsne', FAST_TSNE_SCRIPT_DIR )
		} else {
			fast_tsne_path = sprintf('%s/bin/FItSNE.exe', FAST_TSNE_SCRIPT_DIR)
		}
	}

	if (is.null(data_path)) {
		data_path <- tempfile(pattern='fftRtsne_data_', fileext='.dat')
	}
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='fftRtsne_result_', fileext='.dat')
	}
	if (is.null(fast_tsne_path)) {
		fast_tsne_path <- system2('which', 'fast_tsne', stdout=TRUE)
	}
	fast_tsne_path <- normalizePath(fast_tsne_path)
	if (!file_test('-x', fast_tsne_path)) {
		stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
	}

	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	if (!is.numeric(theta) || (theta<0.0) || (theta>1.0) ) { stop("Incorrect theta.")}
	if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
	if (!is.matrix(X)) { stop("Input X is not a matrix")}
	if (!(max_iter>0)) { stop("Incorrect number of iterations.")}
	if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter<0) { stop("stop_early_exag_iter should be a positive integer")}
	if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
	if (!is.numeric(df)) { stop("df should be numeric")}
	if (!is.wholenumber(dims) || dims<=0) { stop("Incorrect dimensionality.")}
	if (search_k == -1) {
       if (perplexity>0) {
          search_k = n_trees*perplexity*3
       } else if (perplexity==0) {
          search_k = n_trees*max(perplexity_list)*3
       } else { 
          search_k = n_trees*K
       }
    }

	if (fft_not_bh){
	  nbody_algo = 2;
	}else{
	  nbody_algo = 1;
	}

	if (is.null(load_affinities)) {
		load_affinities = 0;
	} else {
		if (load_affinities == 'load') {
			load_affinities = 1;
		} else if (load_affinities == 'save') {
			load_affinities = 2;
		} else {
			load_affinities = 0;
		}
	}
	
	if (ann_not_vptree){
	  knn_algo = 1;
	}else{
	  knn_algo = 2;
	}
	tX = c(t(X))

	f <- file(data_path, "wb")
	n = nrow(X);
	D = ncol(X);
	writeBin(as.integer(n), f,size= 4)
	writeBin( as.integer(D),f,size= 4)
	writeBin( as.numeric(theta), f,size= 8) #theta
	writeBin( as.numeric(perplexity), f,size= 8) #theta

    if (perplexity == 0) {
    	writeBin( as.integer(length(perplexity_list)), f, size=4)
	    writeBin( perplexity_list, f) 
    }

	writeBin( as.integer(dims), f,size=4) #theta
	writeBin( as.integer(max_iter),f,size=4)
	writeBin( as.integer(stop_early_exag_iter),f,size=4)
	writeBin( as.integer(mom_switch_iter),f,size=4)
	writeBin( as.numeric(momentum),f,size=8)
	writeBin( as.numeric(final_momentum),f,size=8)
	writeBin( as.numeric(learning_rate),f,size=8)
	writeBin( as.integer(K),f,size=4) #K
	writeBin( as.numeric(sigma), f,size=8) #sigma
	writeBin( as.integer(nbody_algo), f,size=4)  #not barnes hut
	writeBin( as.integer(knn_algo), f,size=4) 
	writeBin( as.numeric(exaggeration_factor), f,size=8) #compexag
	writeBin( as.integer(no_momentum_during_exag), f,size=4) 
	writeBin( as.integer(n_trees), f,size=4) 
	writeBin( as.integer(search_k), f,size=4) 
	writeBin( as.integer(start_late_exag_iter), f,size=4) 
	writeBin( as.numeric(late_exag_coeff), f,size=8) 
	
	writeBin( as.integer(nterms), f,size=4) 
	writeBin( as.numeric(intervals_per_integer), f,size=8) 
	writeBin( as.integer(min_num_intervals), f,size=4) 
	tX = c(t(X))
	writeBin( tX, f) 
	writeBin( as.integer(rand_seed), f,size=4) 
        writeBin(as.numeric(df), f, size=8)
	writeBin( as.integer(load_affinities), f,size=4) 
	if (! is.null(initialization)){ writeBin( c(t(initialization)), f) }		
        print(df)
	close(f) 

	flag= system2(command=fast_tsne_path, args=c(version_number,data_path, result_path, nthreads));
	if (flag != 0) {
		stop('tsne call failed');
	}
	f <- file(result_path, "rb")
	n <- readBin(f, integer(), n=1, size=4);
	d <- readBin(f, integer(), n=1, size=4);
	Y <- readBin(f, numeric(), n=n*d);
        Y <- t(matrix(Y, nrow=d));
        if (get_costs ) {
            tmp <- readBin(f, integer(), n=1, size=4);
            costs <- readBin(f, numeric(), n=max_iter,size=8);
            Yout <- list( Y=Y, costs=costs);
        }else {
            Yout <- Y;
        }
        close(f)
        file.remove(data_path)
        file.remove(result_path)
        return(Yout)
}

#' Multinomial counts best significance test
#'
#' Initialize Infern0 scope from input data, which can be a matrix, a singlecelleperiment or Seurat object
#'
#' @return list of arrays with length matching the number of row in input. Each contain the Smallest pvalue for high count within a row. 
#' @param count data matrix
#' @export
funTagEval <- function(data){
	freq <- Matrix::colSums(data)
	freq <- freq / sum(freq);
	freq <- 1 - freq
	best = 0;
	best <- rep(0, nrow(data))
	pval <- rep(0, nrow(data))
	mlogpval <- rep(0, nrow(data))

	if (ncol(data) == 2){
	for(i in 1:nrow(data)){
		tsum <- sum(data[i,])
		for(j in 1:ncol(data)){
			tmp <- pbinom(tsum - data[i,j],tsum, freq[j])
			if ((j == 1)||(pval[i] > tmp)) {best[i]<-j; pval[i] <- tmp;}
		}
	}
	return(list(b = best, p=pval))
	}
	spval <- rep(0, nrow(data))
	sbest <- rep(0, nrow(data))
	mlogpval <- rep(0, nrow(data))
	smlogpval <- rep(0, nrow(data))

	for(i in 1:nrow(data)){
		tsum <- sum(data[i,])
		tmp <- pbinom(tsum - data[i,1],tsum, freq[1])
		pval[i] <- pbinom(tsum - data[i,2],tsum, freq[2])
		if (tmp < pval[i]) {
			spval[i] <- pval[i]	
			pval[i] <- tmp
			best[i] <- 1
		}else{
			spval[i] <- tmp
			best[i] <- 2
		}
		sbest[i] <- 3 - best[i]
		for(j in 3:ncol(data)){
			tmp <- pbinom(tsum - data[i,j],tsum, freq[j])
			if (spval[i] > tmp) {
				if (pval[i] > tmp){
					spval[i] <- pval[i]
					sbest[i] <- best[i]
					pval[i] <- tmp
					best[i] <- j
				}else{
					spval[i] <- tmp
					sbest[i] <- j
				}
			}
		}
		mlogpval[i] <- log10(pval[i]) / tsum
		smlogpval[i] <- log10(spval[i]) / tsum
	}
	tmp <- min(mlogpval,na.rm=T)
	mlogpval[is.na(mlogpval)] <- tmp
	tmp <- min(smlogpval,na.rm=T)
	smlogpval[is.na(smlogpval)] <- tmp

return(list(best = best, pvalue=pval,second.best=sbest,second.pvalue= spval, normalized.logpvalue=mlogpval, second.normalized.logpvalue=smlogpval))
}

cleanQuit <- function(){
	print("rm(list = setdiff(ls(), lsf.str())); gc(); q()")
}

testInfern0 <- function(scp =c(), genenames = c("MEG3", "LINGO1", "EMC4", "FCN1"), nb.correlated= 20, cell.filt =c(), cell.meta=list(celltype.broad = "Neuron_cortical")){
	library(InferN0)
	if (is.null(scp)) scp <- InferN0LoadFile("ev4.infr5.scp");
	cell.filt <- InferN0:::.celllistquery(scp, cell.meta, cell.filt)
	if (nb.correlated !=0){
		out <- InferN0FindCoexpressed(scp, gene.list=genenames, cell.list=scp$cell.names[cell.filt], nb.output=nb.correlated)
		genenames <- c(genenames , out$gene.list)
	}
return(InferN0ComputeCovariance(scp, genenames, method="Model", nb.crosseval.partition=8))}
runSoupX <- function(sro, celltype, batch, batch_folder_index= c(), cr_folders=c(), genes_to_regress = c(), soupQuantile = 0.15, tfidfMin = 0.2, filterclustername= "filtered", check_mapping_only=F){
	library(SoupX); library(Seurat)
    if (is.null(cr_folders)) cr_folders = paste("/lustre/scratch117/cellgen/team292/lh20/I-O-", c("1_13-N", "1_13", "11_12-N", "11_12", "2_3-N", "2_3", "4-N", "4", "5_8-N", "5_8", "6_9-N", "6_9", "7-N", "7", "10-N", "10") ,"_C-0/outs/raw_feature_bc_matrix", sep="")
    if (is.null(batch_folder_index)) batch_folder_index = c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,15,16,9,10,11,12,13,14,9,10,11,12) 
    #if (is.null(batch_folder_index)) batch_folder_index = c(1,1,2,2,15,16,3,3,4,4,5,5,6,6,7,8,9,9,10,10,11,11,12,12,13,14) 

    if (is.null(genes_to_regress)) genes_to_regress <- 33568:33759
    cellnames = c()
    for(i in 1:length(cr_folders)){
        curbatches <- levels(sro@meta.data[[batch]])[which(batch_folder_index == i)]
        print(paste("The following samples are associated to the output folder" , cr_folders[i]))
        print(curbatches)
	if (!check_mapping_only){
	emptymat <- Read10X(cr_folders[i])[["Antibody Capture"]]
	flt <- sro@meta.data[[batch]] %in% curbatches
	flt[as.character(sro@meta.data[[celltype]]) == filterclustername] <- F
	clust = paste( as.character(sro@meta.data[[batch]][flt]),  as.character(sro@meta.data[[celltype]][flt]))
	clust <- as.factor(clust)
	mDat = data.frame(clusters=clust@.Data, row.names=rownames(sro@meta.data)[flt])
	toc = sro@assays$RNA@counts[genes_to_regress, flt]
	rownames(emptymat) <- rownames(toc)
	print(paste("Processing soup in ", sum(flt)," cells and ", length(levels(clust)) ," clusters", sep=""))
	sc = SoupChannel(emptymat, toc = toc, metaData = mDat, calcSoupProfile=F)
	sc = estimateSoup(sc, soupRange=c(4,100))
	sc = autoEstCont(sc,soupQuantile= soupQuantile, tfidfMin = tfidfMin,forceAccept=TRUE )
	out = adjustCounts(sc)
	cellnames = c(cellnames, rownames(sro@meta.data)[flt])
	if (i == 1) output <- out
	else output = Matrix::cbind2(output,out)
}}
	if (check_mapping_only) return;
	fout = matrix(0, nrow(output), nrow(sro@meta.data))
	map = match(rownames(sro@meta.data), cellnames)
	print(paste("reordering available data for", ncol(output), "genes and", sum(!is.na(map)), "cells out of", length(map)))
	fout[,!is.na(map)] = as.matrix(output[,map[!is.na(map)]])
return(fout)}
#mat <- dasuperconcat(paste(prefix, "_C-0/outs/filtered_feature_bc_matrix/", sep = ""), prefix)
dasuperconcat <- function(pathlist, prefixlist){
	library(Seurat)
	bdata = c()
	ctgene = c()
	ctcount =c()
	mtcount =c()
	for( i in 1:length(pathlist)){
		tmp <- Read10X(pathlist[i])
		rownames(tmp[["Antibody Capture"]]) <- sub("^","Antibody_", rownames(tmp[["Antibody Capture"]]))
		ctcount = c(ctcount, Matrix::colSums(tmp[["Antibody Capture"]]))
		ctgene = c(ctgene, Matrix::colSums(tmp[["Antibody Capture"]] != 0))
		flt <- grep( "^GRCh38___MT-", rownames(tmp[["Gene Expression"]]))
		mtcount <-  c(mtcount, Matrix::colSums(tmp[["Gene Expression"]][flt,,drop=F]))
		tmp <- Matrix::rbind2(tmp[["Gene Expression"]], tmp[["Antibody Capture"]])
		colnames(tmp) <- sub("^", prefixlist[i], colnames(tmp))
		if (i == 1) out <- tmp
		else out <- Matrix::cbind2(out, tmp)
		bdata <- c(bdata, rep(prefix[i] , ncol(tmp)))
	}
	rownames(out) <- sub("GRCh38___","", rownames(out))
	out <- CreateSeuratObject(out);
	out@meta.data$nCiteseq_count <- ctcount
	out@meta.data$nCiteseq_features <- ctgene
	out@meta.data$MT_fraction <- mtcount / (out@meta.data$nCount_RNA - ctcount) 
	out@meta.data$Batch <- bdata


return(out);}

showInfernalMarkers <- function(output, selected.col, top = 50, specific.constraint=T, negative.tail = F){
	if (negative.tail){

	if (specific.constraint){
		flt = apply(output$Zscore, 1, function(r){which.min(r)}) == selected.col
	}else{
		flt = rep(T,nrow(output$Zscore))
	}
	flt <- flt & (!(grepl("^RPL",rownames(output$Zscore)) | grepl("^RPS",rownames(output$Zscore))))
	print(sum(flt))
	print(output$Zscore[flt,selected.col][order(output$Zscore[flt,selected.col], decreasing=F)[1:top]])


	}else{


	if (specific.constraint){
		flt = apply(output$Zscore, 1, function(r){which.max(r)}) == selected.col
	}else{
		flt = rep(T,nrow(output$Zscore))
	}
	flt <- flt & (!(grepl("^RPL",rownames(output$Zscore)) | grepl("^RPS",rownames(output$Zscore))))
	print(sum(flt))
	print(output$Zscore[flt,selected.col][order(output$Zscore[flt,selected.col], decreasing=T)[1:top]])
	}
}

showInfernalTFIDF <- function(outlist,zthr = 3.0){
	fout = data.frame(Gene = character(), Best =numeric(), Zscore = numeric(), TFIDF = numeric(), NBDE= numeric())
	sout = fout
	# nbDEgenes ( length of documents)
	docsize = matrix(0, length(outlist), ncol(outlist[[1]]$Zscore))
	
	nbcells = table(outlist[[1]]$MarkPartition)
	
	
	for(j in 1:nrow(docsize)){
		for (i in 1:ncol(docsize)){
		docsize[j,i] <- sum(outlist[[j]]$Zscore[,i] > zthr, rm.na=T)
		}
	}
	for( i in 1:nrow(outlist[[1]]$Zscore)){
		j = 1;
		tmp = outlist[[j]]$Zscore[i,]
		size = !is.na(tmp)
		z = rep(0, length(tmp))
		if (sum(!is.na(tmp)) != 0) {
			z[!is.na(tmp)] = tmp[!is.na(tmp)]
			nbpos = sum(tmp > zthr, rm.na=T)
		}else nbpos =0
		if (length(outlist) >1){for(j in 2:length(outlist)){
			tmp = outlist[[j]]$Zscore[i,]
			size = size + !is.na(tmp)
			nbpos = nbpos + sum(tmp > 3, rm.na=T)
			if (sum(!is.na(tmp)) != 0)  z[!is.na(tmp)] = z[!is.na(tmp)] + tmp[!is.na(tmp)] 
		}}
		if (sum(nbpos) > 0){
			for(j in 1:length(nbpos)){
				if (nbpos[j] != 0){
					if (nrow(sout) == 0)  sout <- data.frame(Gene = rownames(outlist[[1]]$Zscore)[i], Best =j, Zscore =  outlist[[1]]$Zscore[j], TFIDF = (outlist[[1]]$PosCoverage[i,j] / sum(outlist[[1]]$PosCoverage[,j])) * (log(sum(nbcells)) - log(sum(outlist[[1]]$PosCoverage[i,]))), NBDE = nbpos[j])
					else sout <- rbind(sout, c(rownames(outlist[[1]]$Zscore)[i], j, outlist[[1]]$Zscore[j], (outlist[[1]]$PosCoverage[i,j] / sum(outlist[[1]]$PosCoverage[,j])) * (log(sum(nbcells)) - log(sum(outlist[[1]]$PosCoverage[i,]))), nbpos[j] ) )
					if (nrow(sout) == 100) {
						print("merge")
						print(colnames(fout))
						print(colnames(sout))
						fout <- rbind(fout, sout)
						sout = data.frame(Gene = character(), Best =numeric(), Zscore = numeric(), TFIDF = numeric(), NBDE= numeric())
					}
				}
			}
		}
	}
	fout <- rbind(fout, sout)
return(fout)}


exportSeuratToFolder <- function(sro, pathprefix, what = c("mat", "meta", "rep")){
	library(Seurat)
	library(Matrix)
	if ("mat" %in% what) writeMM(sro@assays$RNA@counts, paste(pathprefix, "raw.mtx", sep=""))
	if ("meta" %in% what) write.csv(sro@meta.data , paste(pathprefix, "meta.csv",sep=""),quote=F)	
	if ("mat" %in% what) write.table(rownames(sro@assays$RNA@counts) , paste(pathprefix, "feature.csv",sep=""),quote=F,row.names=F,col.names=F)
	if ("rep" %in% what){
		for(i in names(sro@reductions)){
			write.csv(sro@reductions[[i]]@cell.embeddings , paste(pathprefix,"obsm_", i, ".csv",sep=""),quote=F)	
		}
	}
}

importSeuratFromFolder <- function(pathprefix, input_sro = c(), filter = c("meta", "data")){
	library(Seurat)
	library(Matrix)

	if ("data" %in% filter){
		mat <- readMM(paste(pathprefix, "raw.mtx",sep=""))
		rownames(mat) <-read.csv(paste(pathprefix, "feature.csv",sep=""),header=F)$V1
		meta <- read.csv(paste(pathprefix, "meta.csv",sep=""),row.names=1)	
		colnames(mat) <- rownames(meta)
		if (!is.null(input_sro)) stop("if data is loaded, a starting Seurat object cannot be provided")
		sro <- CreateSeuratObject(mat)
	}else{
		if (is.null(input_sro)) stop("if data is not loaded, a starting Seurat object is needed")
		sro = input_sro
	}
	if ("meta" %in% filter) {
		if  (!"data" %in% filter) meta <- read.csv(paste(pathprefix, "meta.csv",sep=""),row.names=1)	
		for( i in colnames(meta)) sro@meta.data[[i]] <- meta[[i]]
	}
return(sro)}

showCum <- function(arga = 1, argb = 3 , argc= 2/3, xval = (1:80) / 10, type = "d",fake=0, order = 16){
	d = rep(0,length(xval))
	if (type == "d"){
		for(i in 1:length(xval)) d[i] = callCumulant(arga, argb , argc , 0.3, xval[i],fake=fake,order=order)$density;
		if (arga == 2) return(plot(xval , d) +  lines(xval, dgamma(xval,  argb, 1.0 / argc)))
		else return(plot(xval , d) +  lines(xval,dnbinom(xval,3, 3/2)))
	}else if (type == "l"){
		for(i in 1:length(xval)) d[i] = callCumulant(arga, argb , argc , 0.3, xval[i],fake=fake,order=order)$LogQuantile;
		library(EQL)
		cumulants = gammaCumulants(argb, argc)

		K.deriv <- function(x) {return(argb * argc/(1 - argc*x))}

		if (arga == 2) {
#			return(cbind(xval,cumulants$mu.inv(xval),cumulants$mu.inv(xval),cumulants$mu.inv(xval)))

			return(cbind(xval,d,log(saddlepoint(xval,1,cumulants)$approx), log(1.0 - pgamma(xval,  argb, 1.0 / argc))))

#			return(plot(xval , d) +  lines(xval, cumulants$mu.inv(xval)))
#			return(plot(xval , d) +  lines(xval, cumulants$kappa2(cumulants$mu.inv(xval))))

#			return(plot(xval , d) +  lines(xval, log(1.0 - pgamma(xval,  argb, 1.0 / argc))))
		}else {
			return(plot(xval , d) +  lines(xval,log (1.0 - pnbinom(xval,3, 3/2))))
			
		}
	}else{
		for(i in 1:length(xval)) d[i] = callCumulant(arga, argb , argc , 0.3, xval[i],fake=fake,order=order)$Quantile;
		if (arga == 2) return(plot(xval , d) +  lines(xval, pgamma(xval,  argb, 1.0 / argc)))
		else return(plot(xval , d) +  lines(xval,dnbinom(xval,3, 3/2)))
	}
#	plot(xval, dgamma(xval,3, 3/2), pgamma(1:100/3,3/2))

return( plot(xval , d))}

runSlingShot <- function(sce, cellfilter, clusterLabels = "leiden_res2_scvi_batchv2_Stromal", subsample_group = "DonorID", max_cells=1000){
	library(slingshot, quietly=T)
	nflt <- rep(F,length(cellfilter))
	met <- colData(sce)
	for(i in levels(met[[subsample_group]])){
		print(i)
		ind <- which(met[[subsample_group]][cellfilter] == i)
		if (length(ind) > max_cells) ind <- sample(ind)[1:max_cells]
		nflt[cellfilter][ind]<-T
	}
	ssce <- sce[,nflt]
	print(table(colData(ssce)[[subsample_group]]))
	print(paste("Nb cells to fit trajectory", sum(nflt)))
	print(table(colData(ssce)[[clusterLabels]]))
	
	ssce <- slingshot(ssce, clusterLabels = clusterLabels, reducedDim = 'STR')
return(ssce)}

plotPointCloud <- function(xcoor, ycoor, color, sizes, rect = c(), plot.attribs=c()){
	tt <- data.frame(X= xcoor, Y=-ycoor, fill=color, size=sizes)
	p <- ggplot(tt, aes(X, Y, color=fill, fill=fill, size=size))
	p <- p + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + coord_fixed()
	return(changeStyle(p,plot.attribs)) 	
}

plotPointCloud2 <- function(coords, sizes, use.color = c(), rect = c(), plot.attribs=c(),nbshow = 3){
	if (ncol(sizes) < nbshow) nbshow = ncol(sizes)
	if (is.null(use.color)) use.color <- myrainbow(ncol(sizes))
	if (is.null(colnames(sizes)))  colnames(sizes) <- paste("Class", 1:ncol(sizes))
	tt <- data.frame(X= rep(coords[,1],nbshow), Y=rep(-coords[,2],nbshow), fill=rep(0, nbshow * nrow(coords)), size = rep(0.5, nbshow * nrow(coords)))
	for(i in 1:nrow(sizes)) {
		daorder = order(sizes[i,],decreasing=T)[1:nbshow]
		curs = 0;
		fact = 1.0 / sum(sizes[i,daorder])
		for(j in nbshow:1) {
			tt$fill[i + (j-1) * nrow(sizes)] = daorder[j]
			curs = curs + sizes[i,daorder[j]]
			tt$size[i + (j-1) * nrow(sizes)] = (curs * fact) ^ 2
		}
	}
	
	tt$fill <- factor(tt$fill, levels = 1:ncol(sizes))
	levels(tt$fill) <- colnames(sizes)
	p <- ggplot(tt, aes(X, Y, color=fill, fill=fill, size=size))
	p <- p + geom_point() + scale_color_manual(values=use.color,label=colnames(sizes),drop = FALSE) + scale_size_area(max_size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + coord_fixed()
	return(changeStyle(p,plot.attribs)) 	
}


#' Load the fitted exposures and normalise them
#'
#' Load's the exposure table, separates out the goodness-of-fit metrics and exposures and normalises the exposures to sum to 1 across each sample.
#'
#' @param tgt The relevant *_fitExposures.tsv file that results should be loaded from.  Can either be the exact file, or the base part to which _fitExposures.tsv will be added.
#' @return A list containing the normalised exposure table, raw exposure table, and a goodness-of-fit table
normaliseExposures = function(tgt){
  #Is this just the base?
  if(file.exists(paste0(tgt,'_fitExposures.tsv')))
    tgt = paste0(tgt,'_fitExposures.tsv')
  fit = read.table(tgt,sep='\t',header=TRUE)
  #Extract the goodness of fit rows
  gofNoms = c('pR2','fitCount','obsCount')
  gof = fit[gofNoms,]
  gof['log2(countRatio)',] = log2(unlist(gof['fitCount',]/gof['obsCount',]))
  #And the exposure table
  exposures = fit[-match(gofNoms,rownames(fit)),]
  #Normalise
  exposures = t(t(exposures)/colSums(exposures))
  #Done
  return(list(exposures=exposures,gof=gof,raw=fit[-match(gofNoms,rownames(fit)),]))
}

#' Plot normalised exposures
#' 
#' Takes the normalised exposure object produced by \code{normaliseExposures} and creates a basic heatmap to visualise the results.
#'
#' @param fit The normalised exposure object returned by \code{normaliseExposures}
#' @param exposureScale Range of exposures to show.
#' @param cluster_rows Should we cluster the rows.
#' @param ... Passed to Heatmap.
plotExposures = function(fit,exposureScale=c(0,0.5),cluster_rows=FALSE,use_raster = FALSE,show_column_names=F,extraannot = c()){
  library("ComplexHeatmap")
  #Colours for exposures
  hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
  #Colours for pR2 metric
  pR2Cols  = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  #Colours for log library sizes
  libCols = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
  #And library size ratio
  libRatCols = c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')
  #Create the bottom annotation
  gof = fit$gof
  gof['fitCount',] = log10(gof['fitCount',])
  gof['obsCount',] = log10(gof['obsCount',])
  rownames(gof) = gsub('(.*Count)','log10(\\1)',rownames(gof))
  #Decide on range for library size
  libRange = range(gof[grep('Count',rownames(gof)),])
  #Convert colours to ramps
  bCols = circlize::colorRamp2(seq(0,1,length.out=length(pR2Cols)),pR2Cols)
  lCols = circlize::colorRamp2(seq(libRange[1],libRange[2],length.out=length(libCols)),libCols)
  hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)
  lrCols = circlize::colorRamp2(seq(-1,1,length.out=length(libRatCols)),libRatCols)
  if (!is.null(extraannot)) {
	gof['tissue',] = extraannot@.Data
  	#dacol[["tissue"]] = mydoublerainbow(rep(1,length(levels(extraannot))))
  }
  botAnno = HeatmapAnnotation(df = t(gof),
                              annotation_name_side = 'left',
                              col = list(pR2 =bCols,
                                         `log10(fitCount)` = lCols,
                                         `log10(obsCount)` = lCols,
                                         `log2(countRatio)` = lrCols)
  )
  hm = Heatmap((fit$exposures),
               col=hmColObj,
               name='Exposures',
               bottom_annotation=botAnno,
               cluster_rows = cluster_rows,
	       show_column_names=show_column_names
  )
  return(hm)
}

genSimilyBulk = function(sce, ctannot = "Broad cell type", dnannot= "DonorID", exannot = ""){
	library(SingleCellExperiment)
	ctannot = colData(sce)[[ctannot]] 
	dnannot = colData(sce)[[dnannot]]
	dnindex = matrix(0,length(levels(dnannot)) * (length(levels(dnannot)) - 1) / 2, 2)
	concat = paste(as.character(ctannot@.Data), as.character(dnannot@.Data),sep= '_')

	dadepth <- colData(sce)[["log2p1_count"]]
	names(dadepth) <- 1:length(dadepth)
	scell_filter <- rep(T, dim(sce)[2])
	for(i in unique(concat)) {
		if (sum(concat[scell_filter] == i) > 480){
			daind <- as.integer(names(sort(dadepth[concat == i], decreasing=T)))
			scell_filter[daind[481:length(daind)]] <- F
		}
	}
	nbentries <- table(ctannot[scell_filter], dnannot[scell_filter])
	print(nbentries)


	k=1;
	for( i in 2:length(levels(dnannot))){
		for( j in 1:(i-1)){
			dnindex[k,] = c(i,j)
			k = k +1	
		}	
	}
	ctindex = matrix(0,length(levels(ctannot)) * (length(levels(ctannot)) - 1) / 2, 2)
	k=1;
	for( i in 2:length(levels(ctannot))){
		for( j in 1:(i-1)){
			ctindex[k,] = c(i,j)
			k = k +1	
		}	
	}
	
	fracbulkout= Matrix(0, dim(sce)[1], nrow(nbentries) * ncol(nbentries)*8)
	nbcell = rep(0,nrow(nbentries) * ncol(nbentries)*8)
	k=1
	for(i in 1:nrow(nbentries)){
		for(j in 1:ncol(nbentries)){
			whi <- as.integer(names(dadepth[(concat == paste(as.character(i), as.character(j),sep= '_')) & scell_filter]))
			das = sample(1:length(whi))
			limit = as.integer((0:8) * length(whi) / 8) + 1
			print(limit)
			for(l in 1:8){
				nbcell[k] = limit[l+1] - limit[l]
				if (nbcell[k] == 1) fracbulkout[,k] = assays(sce)[[1]][, whi[das[limit[l]:(limit[l+1]-1)]]]
				else if (nbcell[k] > 1) fracbulkout[,k] = Matrix::rowSums(assays(sce)[[1]][, whi[das[limit[l]:(limit[l+1]-1)]]])
				k <-k+1
			}
					
		}
	}
	print("Alive")
	k=1
	bulkout= Matrix(0, dim(sce)[1], nrow(dnindex) * nrow(ctindex))
	bulkindex = matrix(0, nrow(dnindex) * nrow(ctindex),4)
	bulkcounts = matrix(0, nrow(dnindex) * nrow(ctindex),4)
	for(i in 1:nrow(ctindex)){
		for(j in 1:nrow(dnindex)){
			bulkindex[k,] = c(dnindex[j,], ctindex[i,])
			curind = c( (ctindex[i,1]-1) * ncol(nbentries) + dnindex[j,1]-1, (ctindex[i,1]-1) * ncol(nbentries) + dnindex[j,2]-1, (ctindex[i,2]-1) * ncol(nbentries) + dnindex[j,1]-1, (ctindex[i,2]-1) * ncol(nbentries) + dnindex[j,2]-1 )  * 8
			das = c(sample(1:8)+curind[1],sample(1:8)+curind[2],sample(1:8)+curind[3],sample(1:8)+curind[4]) ;
			bulkout[,k] = fracbulkout[, das[1] ] + fracbulkout[, das[2] ] + fracbulkout[, das[3] ]+ fracbulkout[, das[4] ]+ fracbulkout[, das[9] ] + fracbulkout[, das[10] ] + fracbulkout[, das[11] ] + fracbulkout[, das[12] ] + fracbulkout[, das[17] ] + fracbulkout[, das[18] ] + fracbulkout[, das[19] ]+ fracbulkout[, das[20] ]+ fracbulkout[, das[25] ] + fracbulkout[, das[26] ] + fracbulkout[, das[27] ] + fracbulkout[, das[28] ]
			bulkcounts[k,] = c( nbcell[das[1]] + nbcell[das[2]] + nbcell[das[3]] + nbcell[das[4]], nbcell[das[9]] + nbcell[das[10]] + nbcell[das[11]] + nbcell[das[12]],nbcell[das[17]] + nbcell[das[18]] + nbcell[das[19]] + nbcell[das[20]],nbcell[das[25]] + nbcell[das[26]] + nbcell[das[27]] + nbcell[das[28]])


			k<-k+1
		}
	}

	sceout <- SingleCellExperiment(bulkout)
	colData(sceout)[["DonorA"]] = factor(bulkindex[,1], 1:length(levels(dnannot)))
	levels(colData(sceout)[["DonorA"]]) = levels(dnannot)
	colData(sceout)[["DonorB"]] = factor(bulkindex[,2], 1:length(levels(dnannot)))
	levels(colData(sceout)[["DonorB"]]) = levels(dnannot)
	colData(sceout)[["CtA"]] = factor(bulkindex[,3], 1:length(levels(ctannot)))
	levels(colData(sceout)[["CtA"]]) = levels(ctannot)
	colData(sceout)[["CtB"]] = factor(bulkindex[,4], 1:length(levels(ctannot)))
	levels(colData(sceout)[["CtB"]]) = levels(ctannot)
	colData(sceout)[["ctAdnA_cellcount"]] = bulkcounts[,1]
	colData(sceout)[["ctAdnB_cellcount"]] = bulkcounts[,2]
	colData(sceout)[["ctBdnA_cellcount"]] = bulkcounts[,3]
	colData(sceout)[["ctBdnB_cellcount"]] = bulkcounts[,4]	

	ctfrac = matrix(0, nrow(dnindex) * nrow(ctindex), length(levels(ctannot)))
	for(i in 1:length(levels(ctannot))){
		ctfrac[,i] = ((bulkindex[,3] == i) * (bulkcounts[,1] + bulkcounts[,2]) ) + ((bulkindex[,4] == i) * (bulkcounts[,3] + bulkcounts[,4]))
	}
	ctfrac = ctfrac / rowSums(ctfrac)
	
	dnfrac = matrix(0, nrow(dnindex) * nrow(ctindex), length(levels(dnannot)))
	for(i in 1:length(levels(dnannot))){
		dnfrac[,i] = ((bulkindex[,1] == i) * (bulkcounts[,1] + bulkcounts[,3]) ) + ((bulkindex[,2] == i) * (bulkcounts[,2] + bulkcounts[,4]))
	}
	dnfrac = dnfrac / rowSums(dnfrac)
	#reducedDims(sce) <- list(donor= )

	colnames(ctfrac) <- levels(ctannot)
	colnames(dnfrac) <- levels(dnannot)
	rownames(sceout) = rownames(sce)
	rowData(sceout) = rowData(sce)
	colnames(sceout) = paste(bulkindex[,1],bulkindex[,2],bulkindex[,3],bulkindex[,4] , sep="_")
	reducedDim(sceout, "celltype") <- ctfrac
	reducedDim(sceout, "donor") <- dnfrac
	
	
return(list(sce = sceout, data= bulkout, indexes = bulkindex, cellcounts = bulkcounts, ct = levels(ctannot), dn=levels(dnannot), dndn = dnfrac, ctct = ctfrac))}


testMusic = function(sce, bulk, sc_cellfilter = c(), genefilter = c(), sc_annot = "Broad cell type",  sc_dn_annot = "DonorID"){
	# assumes inputs are singlecellexperiments
	
	library(SingleCellExperiment)
	library(MuSiC)
	library(xbioc)
	dainput <- as.matrix(assays(sce[genefilter,sc_cellfilter])[[1]])
	rownames(dainput) <- rownames(sce)[genefilter]
	colnames(dainput) <- colnames(sce)[sc_cellfilter]
	scee = ExpressionSet(dainput)
	dainput <- as.matrix(assays(bulk[genefilter,])[[1]])
	rownames(dainput) <- rownames(bulk)[genefilter]
	colnames(dainput) <- colnames(bulk)
	bulke = ExpressionSet(dainput)

	scee@phenoData[["ct"]] <- droplevels(colData(sce)[[sc_annot]][sc_cellfilter])
	scee@phenoData[["dn"]] <- droplevels(colData(sce)[[sc_dn_annot]][sc_cellfilter])

	
	musout = music_prop(bulk.eset = bulke, scee, clusters = "ct", samples = "dn", verbose = T)
}

#library(InferN0); net <- InferN0IdentifyNetwork(readRDS("tmp.rds"), nb.thread=2)

prepThatEndoAtlas <- function(maxcell = 500){
	sce <- readRDS("/lustre/scratch117/cellgen/team292/lh20/endoaltas.sce.rds")
	scell_filter <- (colData(sce)[["Broad cell type"]] != 'Excluded')&!(colData(sce)$SampleID %in% c("4861STDY7387181", "4861STDY7387182", "4861STDY7387183", "4861STDY7771115" ,"4861STDY7771123"))


	concat <- paste(colData(sce)[["Broad cell type"]], colData(sce)[["DonorID"]])
	dadepth <- colData(sce)[["log2p1_count"]]
	dadepth[!scell_filter] <- 0
	names(dadepth) <- 1:length(dadepth)

	for(i in unique(concat[scell_filter])) {
		if (sum(concat[scell_filter] == i) > maxcell){
			daind <- as.integer(names(sort(dadepth[concat == i], decreasing=T)))
			scell_filter[daind[(maxcell+1):length(daind)]] <- F
		}
	}

return(list(sce=sce, scell_filter = scell_filter))}

sumRowsUsingAttribute <- function(data, attribute){
	if (length(attribute) != nrow(data)) stop("expects an attribute for every row, in the match order")

	outrownames <- unique(attribute)
	fout <- matrix(0,length(outrownames), ncol(data))
	rownames(fout) <- outrownames
	colnames(fout) <- colnames(data)
	oldnames = rep("",length(outrownames))
	for(i in 1:length(outrownames)){
		if ((i %% 100) == 0) print(outrownames[i]);
		fout[i,] = colSums(data[attribute == outrownames[i],,drop=F])
		oldnames[i] = paste(rownames(data)[attribute == outrownames[i]])
	}

return(list(data=fout, oldname=oldnames));}

prepThatGonadAtlas <- function(maxcell = 500){
	sce <- readRDS("/lustre/scratch117/cellgen/team292/lh20/gonadaltas.sce.rds")
	scell_filter <- (!(colData(sce)[["lineages_v1"]] %in% c('Doublets', 'lowQC', 'nan', 'Metanephros', 'Neural')))&(colData(sce)$cryopreserved == "No")


	concat <- paste(colData(sce)[["lineages_v1"]], colData(sce)[["sample"]])
	dadepth <- colData(sce)[["log2p1_count"]]
	dadepth[!scell_filter] <- 0
	names(dadepth) <- 1:length(dadepth)

	for(i in unique(concat[scell_filter])) {
		if (sum(concat[scell_filter] == i) > maxcell){
			daind <- as.integer(names(sort(dadepth[concat == i], decreasing=T)))
			scell_filter[daind[(maxcell+1):length(daind)]] <- F
		}
	}

return(list(sce=sce, scell_filter = scell_filter))}




