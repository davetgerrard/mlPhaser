

alleleDbToHaploTable <- function(alleleDb, haploMarkerNames=names(alleleDb), shareRatio=0.9)  {
	hlaAlleTableList <- list()

	fullAlleleList <- character()
	intersectAlleleList <- character()
	for(j in 1:length(haploMarkerNames))  {
		thisMarker <- haploMarkerNames[j]
		#thisMarker <- "HLA_A2"
		locusAlleleTable <- data.frame()
		for(i in 1:length(alleleDb[[thisMarker]]@variantMap))  {
			thisVariant <- alleleDb[[thisMarker]]@variantMap[i]
			alleleNames <- unlist(strsplit(thisVariant[[1]],"_", fixed=T))
			sequence <- names(thisVariant[1])
			tempTable <- data.frame(allele=alleleNames,sequence=sequence)
			locusAlleleTable <- rbind(locusAlleleTable,tempTable)
		}
		names(locusAlleleTable) <- c("allele", thisMarker) 
		hlaAlleTableList[[thisMarker]] <- locusAlleleTable
		fullAlleleList <- union(fullAlleleList ,locusAlleleTable$allele)
		if(j==1)  {
			intersectAlleleList <- as.character(locusAlleleTable$allele)
		} else {
			intersectAlleleList <- intersect(as.character(intersectAlleleList),as.character(locusAlleleTable$allele))
		}
		#print(j)
		#print(head(locusAlleleTable))
	}

	print(paste(length(intersectAlleleList), "of", length(fullAlleleList), "alleles shared"))
	ratio <- length(intersectAlleleList)/length(fullAlleleList)

	if(ratio < shareRatio)  {
		warning(paste("Proportion of shared allele names is less than", shareRatio))
	}

	for(i in 1:length(haploMarkerNames))  {
		thisMarker <- haploMarkerNames[i]
		if(i==1)  {
			hlaHaploTable <- hlaAlleTableList[[thisMarker]]
		} else {
			hlaHaploTable <- merge(hlaHaploTable, hlaAlleTableList[[thisMarker]], by="allele")
		}
	}

	## Need to convert haplotypes to haplotype tables with column for each marker
	## DONE: generalise this part.
	##hlaHaploTable <- merge(hlaAlleTableList[["HLA_C2"]], hlaAlleTableList[["HLA_C3"]], by="allele")
	##hlaHaploTable <- merge(hlaAlleTableList[["HLA_B2"]], hlaAlleTableList[["HLA_B3"]], by="allele")
	#hlaHaploTable <- merge(hlaAlleTableList[["HLA_A2"]], hlaAlleTableList[["HLA_A3"]], by="allele")
	rownames(hlaHaploTable) <- hlaHaploTable$allele
	hlaHaploTable <- subset(hlaHaploTable, select=-allele)
	return(hlaHaploTable)
}


calledGenotypesToGenoTable <- function(genotypeList,alleleDb,alleleNumber=2, markerNames=names(genotypeList)) {
	
	# Should put some checks in here that makers present in all lists.
	#alleleNumber <- 2
	# add alllele sequence to genotypeTables.
	for(thisMarker in markerNames)  {
		for(k in 1:alleleNumber) {
			genotypeList[[thisMarker]]@genotypeTable[,paste("seq",k,sep=".")] <-""
			for(i in 1:nrow(genotypeList[[thisMarker]]@genotypeTable) ) {
				seq <- names(which(alleleDb[[thisMarker]]@variantMap == genotypeList[[thisMarker]]@genotypeTable[i,paste("allele",k,sep=".")]))
				genotypeList[[thisMarker]]@genotypeTable[i,paste("seq",k,sep=".")] <- ifelse(length(seq) > 0, seq, NA)
			}
		}
	}

	#genotypeList[[thisMarker]]@genotypeTable

	genotypeList.refined <- data.frame()

	myTempGenotypes <- list()
	for(i in 1:length(markerNames))  {
		thisMarker <- markerNames[i]
		# sample id as rownames and one column per allele per marker (e.g. 2 markers/2 alleles = 4)
		myTempGenotypes[[thisMarker]] <- genotypeList[[thisMarker]]@genotypeTable

		myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")] <- myTempGenotypes[[thisMarker]][,paste("seq",1,sep=".")]
		myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")] <- myTempGenotypes[[thisMarker]][,paste("seq",2,sep=".")]
		myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")] <- ifelse(myTempGenotypes[[thisMarker]][,"status"] =="HOMOZYGOTE",myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")],myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")])
		# should add in checks if genotype status' are tooFewReads or complexVars.
		myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")] <- ifelse(myTempGenotypes[[thisMarker]][,"status"] =="complexVars",NA,myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")])
		myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")] <- ifelse(myTempGenotypes[[thisMarker]][,"status"] =="complexVars",NA,myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")])
		myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")] <- ifelse(myTempGenotypes[[thisMarker]][,"status"] =="tooFewReads",NA,myTempGenotypes[[thisMarker]][,paste(thisMarker,1,sep=".")])
		myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")] <- ifelse(myTempGenotypes[[thisMarker]][,"status"] =="tooFewReads",NA,myTempGenotypes[[thisMarker]][,paste(thisMarker,2,sep=".")])
		#rownames(myTempGenotypes[[thisMarker]]) <- myTempGenotypes[[thisMarker]][,"sample"]
		myTempGenotypes[[thisMarker]] <- subset(myTempGenotypes[[thisMarker]], select=c("sample", paste(thisMarker,1,sep="."),paste(thisMarker,2,sep=".")))
		if(i==1)  {
			genotypeList.refined <- myTempGenotypes[[thisMarker]]
		} else {
			genotypeList.refined <- merge(genotypeList.refined, myTempGenotypes[[thisMarker]], by="sample")
		}
	}

	#print(nrow(genotypeList.refined))
	#print(genotypeList.refined)
	# need to generalise this section
	#genotypeList.refined <- merge(myTempGenotypes[['HLA_A2']],myTempGenotypes[['HLA_A3']], by="sample")
	#genotypeList.refined <- merge(myTempGenotypes[['HLA_B2']],myTempGenotypes[['HLA_B3']], by="sample")
	#genotypeList.refined <- merge(myTempGenotypes[['HLA_C2']],myTempGenotypes[['HLA_C3']], by="sample")

	rownames(genotypeList.refined) <- genotypeList.refined[,'sample']
	genotypeList.refined <- subset(genotypeList.refined, select=-sample)
	return(genotypeList.refined)
}