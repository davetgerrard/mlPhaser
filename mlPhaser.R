

## creates genotypes from a table of hapotypes. Will use frequencies.
simGenoFromHaplo <- function(haploTable, haploFreqs, n=1, ploidy=2)  {
	
	## put in checks that the haplotypes in haploTable have freqs in haploFreqs


	#haplos <- rep(names(haploFreqs), haploFreqs)

	sampleHaplos <- sample(names(haploFreqs), n*ploidy, replace=T, prob=haploFreqs )

	genoTable <- data.frame()

	for(i in 1:ploidy)  {
		newCols <- haploTable[sampleHaplos[((i-1)*n)+1:n],]
		names(newCols) <- paste(names(newCols),i,sep=".")
		if(i==1) {
			genoTable <- newCols
		} else {
			genoTable <- cbind(genoTable, newCols)
		}
	}
		
	rownames(genoTable) <- 1:n
	#return(sampleHaplos)
	genoTable <- genoTable[,order(colnames(genoTable))]		# should columns be re-ordered?
	return(genoTable)
}



#unused!
getPloidyFromGenoTable <- function(genotype)  {
	
}


#####  this is based only on using known haplotypes. 
##  Each valid haplotype is found by intersecting across each locus.
### What about if there are no valid haplotypes?  NA?
### Account for novel haplotypes?
## Could have optional call to list all possible haplotypes. (and assign very low probabilities to novel forms).
## This function relies on grep of loci names. Not good idea when loci names are "A1", "A10" etc. 
listValidHaplotypes <- function(genotype,haplotypes, ploidy=2)  {
	# could generate haplotype list from fresh. Might be very long.	
	
	# TEST IF genotype and haplotypeList are list() or data.frame().

	if(class(haplotypes) == "data.frame")  {
		haploList <- tableHaploToList(haplotypes)
		haploTable <- haplotypes
	} else if(class(haplotypes) == "list") {
		haploList <- haplotypes
		haploTable <- listHaploToTable(haplotypes)
	} else {
		stop("bad format for haplotypes")
	}

	#lociNames <- c("A", "B", "C")	# TODO find these from genotype
	lociNames <- names(haploList[[1]])	# currently finds lociNames from first entry in haplotypeList

	if(class(genotype) == "data.frame")  {
		genoTable <- genotype
		genoList <- tableGenoToList(genotype, lociNames)
	} else if (class(genotype) == "list") {
		genoList <- genotype
		genoTable  <- listGenoToTable(genotype)
	} else {
		stop("bad format for genotype")
	}





	# OLD
	#validHaplos <- character()
	#for(i in 1:length(lociNames))  {
	#	thisLocus <- lociNames[i]
	#	## this is very poor way to find columns. 
	#	lociCols <- grep(thisLocus,colnames(genotype))
	##	genotypes <- genotype[lociCols]
	#	# rownames(haplotypes)[haplotypes[,"A"] == "a"]
	#	theseHaplos <- character()
	#	for(thisGenotype in genotypes) {
	#		theseHaplos <- c(theseHaplos, rownames(haplotypeList)[haplotypeList[,thisLocus] == thisGenotype])
	#	}
	#	if(i == 1)  {
	#		validHaplos <- theseHaplos
	#	} else  {
	#		validHaplos <- intersect(validHaplos,theseHaplos)
	#	}
	#	
	#}

	# NEW
	validHaplos <- character()
	for(i in 1:length(lociNames))  {
		thisLocus <- lociNames[i]
		genotypes <- genoList[[thisLocus]]
		theseHaplos <- character()
		for(thisGenotype in genotypes) {
			theseHaplos <- c(theseHaplos, rownames(haploTable)[haploTable[,thisLocus] == thisGenotype])
		}
		if(i == 1)  {
			validHaplos <- theseHaplos
		} else  {
			validHaplos <- intersect(validHaplos,theseHaplos)
		}
		
	}

	validHaploList  <- haploList[validHaplos]


	### need to turn list of validHaplotypes into consistent pairs/groups.
	# e.g. bb,bb,bb -> bbb/bbb homozygote
	#      ab,ab,ab -> aaa/bbb ; aab/bba ; abb/baa ; aba/bab etc. 
	# Return complex data. 	with vector of haplos and list of pairs, each a list.
	
	#validGroups <- list()
	# this bit needs to deal with ploidy > 2 . Recursive? 
	# With ploidy > 2, there are then multiple options after each first ploidy.
	#for(thisHaplo in validHaplos)  {
	#	thisRemGeno <- remGeno(haplotypeList[thisHaplo,], genotype)
	#}
	

	return(validHaploList)
}





# Tries to extract a single haplotype from a compound genotype and return, amongst other things, the remainder genotype.
# this version to test and return if remaining alleles.
remGeno <- function(haplotype, genotypeList)  {
	# haplotype should be a list, preferably named.
	# and entries as alleles at each locus.
	# genotype should be formatted as a list with loci as entries and present alleles under each entry.
	# can use tableGenoToList() to format.
	
	passTable <- logical(length=length(haplotype[[1]]))
	names(passTable) <- names(haplotype[[1]])
	remList <- list()
	for(thisLocus in names(haplotype[[1]])) {
		matchIndex <-  match(haplotype[[1]][[thisLocus]] , genotypeList[[thisLocus]])
		if(is.na(matchIndex))  {
			passTable[thisLocus] <- FALSE
		} else {
			passTable[thisLocus] <- TRUE
			remList[[thisLocus]] <- genotypeList[[thisLocus]][-which(genotypeList[[thisLocus]] == haplotype[[1]][[thisLocus]])[1]]
		}
	}
	if(sum(unlist(lapply(remList,length))) == 0) {
		remList <- NULL
	}

	return(list(haplotype=haplotype, passTable=passTable,remList = remList))
}

## test if a single haplotype is consistent with a genotype.
testHaploInGeno <- function(haplotype, genotypeList)  {
	# haplotype should be one line dataframe with columns as loci names
	# and entries as alleles at each locus.
	# genotype should be formatted as a list with loci as entries and present alleles under each entry.
	# can use tableGenoToList() to format.
	haploPresent <- FALSE
	testList <- logical()
	for(thisLocus in names(genotypeList)) {
		testList <- c(testList,which(genotypeList[[thisLocus]] == haplotype[,thisLocus]) > 0)
	}

	haploPresent <- all(testList)

	return(haploPresent)
}




# perhaps should be called extractHaploFromGeno
# return genoList with certain haplo taken out, error if haplo cannot be taken out.



## recursion
## requires access to a global list to store results:
### NOT THIS ANYMORE  groupStorageEnv <- new.env(parent = baseenv())
#assign("validHaploGroups", list(), globalenv())
## then call function
## then retrieve validGroups.
#get("validHaploGroups", globalenv())
recurseHaplos <- function(validHaplotypes, remGenotype, group) {
	#print("New recurse")
		#print("Rem geno, try more haplotypes")
		for(i in 1:length(validHaplotypes))  {
			thisHaplotype <- validHaplotypes[i]
			#print(names(validHaplotypes)[i] )
			remGeneList <- remGeno(thisHaplotype , remGenotype)
			#print(remGeneList[['passTable']])
			if(all(remGeneList[['passTable']])) {
				#print("Extracted haplotype")
				#myGroup <- c(group,remGeneList[['haplotype']])
				#myGroup <- c(group,list(thisHaplotype))
				myGroup <- c(group,thisHaplotype)
				if( is.null(remGeneList[['remList']]) )  {
					#print("Valid pair")
					#print(myGroup)
					#groupStore <- list()
					#for(g in 1:length(myGroup))  {
					#	groupStore[[g]] <- myGroup[[g]]
					#}
					#groupStoreList <- list(groupStore)
					groupStoreList <- list(myGroup)
					# store the valid group within a previously set up globally accesible list.	List of list of list.
					#assign("validHaploGroups", c(get("validHaploGroups", groupStorageEnv),groupStoreList), groupStorageEnv)
					assign("validHaploGroups", c(get("validHaploGroups", globalenv()),groupStoreList), globalenv())
					#print("breaking")
					break()	
				} else  {  # further genotypes to extract
					recurseHaplos(validHaplotypes,remGeneList[['remList']],myGroup)
				}
			}
		}
		#print("tried all haplos or broke")
}

##wrapper function to set up and control the recursive search for groups of haplotyps, each of which are consistent with the genotype in question.
getValidHaploGroups <- function(genotype, haplotypes)  {

	if(class(haplotypes) == "data.frame")  {
		haploList <- tableHaploToList(haplotypes)
		haploTable <- haplotypes
	} else if(class(haplotypes) == "list") {
		haploList <- haplotypes
		haploTable <- listHaploToTable(haplotypes)
	} else {
		stop("bad format for haplotypes")
	}

	if(class(genotype) == "data.frame")  {
		genoTable <- genotype
		genoList <- tableGenoToList(genotype, names(haploTable))
	} else if (class(genotype) == "list") {
		genoList <- genotype
		genoTable  <- listGenoToTable(genotype)
	} else {
		stop("bad format for genotype")
	}


	
	startGroup=list()
	#startRemGeno <- list('remList'=genoList )
	startRemGeno <- list()
	startRemGeno[['remList']] <- genoList

	#print(genoList)
	#print(genoTable)
	#print(haploList)
	#print(haploTable)

	#haploList<- haplotypes[listValidHaplotypes(thisGenotype, haplotypes )]
	#haploList.valid <- listValidHaplotypes(thisGenotype, haplotypes )
	#haploList.valid <- listValidHaplotypes(genoList[1], haploTable )
	haploList.valid <- listValidHaplotypes(genoTable, haploTable )

	assign("validHaploGroups", list(), globalenv() )
	# then call function
	# then retrieve validGroups.
	recurseHaplos(validHaplotypes=haploList.valid,startRemGeno[['remList']] ,
				group=startGroup )
	validHaploGroups <- get("validHaploGroups", envir=globalenv())
	#rm(groupStorageEnv)
	validHaploGroups.nonRedund <- reduceRedundantList(validHaploGroups )
	rm(validHaploGroups, envir=globalenv())
	return(validHaploGroups.nonRedund)
}


# creates a list of lists replacing haplotypes with their names from a full, named list of haplotypes.
# will use this function 'temporarily' to match returned haplotypes with the originals.
haploGroupsNamed <- function(haploListOfLists, fullHaploList)  {
	# assume this is a list of list of lists
	# TODO: check there is any data coming in (will not be if getValidHaploGroups() produced empty list).
	resultList <- list()
	if(length(haploListOfLists) < 1)  {
		warning("Empty list of valid groups")
	} else {
		for(i in 1:length(haploListOfLists)) {
			resultList[[i]] <- list()
			thisGroup <- haploListOfLists[[i]]
			for(j in 1:length(thisGroup)) {
				resultList[[i]][j] <- names(fullHaploList[ match(thisGroup[j], fullHaploList)])
			}
		}
	}
	return(resultList)

}

#namedGroups <- haploGroupsNamed(result,fullHaploList )

## first attempt at assigning probabilities/likelihoods to competing haplotype groups.
## prints and summed log-likelihood and reconstituted (exp()) likelihood.
printHaploProbs <- function(namedHaploGroups, haploFrequencies) {
	for(i in 1:length(namedHaploGroups)) {	
		thisGroup <- namedHaploGroups[[i]]
			sumLikelihood <- 0
			for(thisHaplo in thisGroup)  {
				thisProb <- haploFrequencies[thisHaplo]
				sumLikelihood <- sumLikelihood + log(thisProb)
				cat(thisHaplo)
				cat("/")
			}
			cat("\t")
			cat(sumLikelihood)
			cat("\t")
			cat(exp(sumLikelihood))
			cat("\n")
	}
}

#printHaploProbs(namedGroups,haploFreqs )


# simple function to print out list of list haplotypes, one per line
printHaploGroups <- function(haploListOfLists)  {
	# assume this is a list of list of lists
	for(i in 1:length(haploListOfLists)) {
		thisGroup <- haploListOfLists[[i]]
		for(j in 1:length(thisGroup)) {
			thisHaplo <- thisGroup[[j]]
				for(k in 1:length(thisHaplo))  {
					cat(thisHaplo[[k]])
					cat(".")
				}
			cat("/")
		}
		cat("\n")

	}

}

#printHaploGroups(result )



# and some way of scoring/ranking them based on relative probabilities.



## The recursive method of finding groups of consistent haplotypes does not differentiate,
## re-arranged versions of the same set. e.g. keeps aaa/bbb AND bbb/aaa.
##  This function removes that redundancy from the results.
## uses the length of intersect to determine if two lists contain all the same elements. 
reduceRedundantList <- function(startList)  {	
	listLength <- length(startList)
	if(listLength < 2) {
		return(startList)

	} else {
		removeIndex <- integer()

		for(i in 1:(listLength-1)) {
			if(!is.na(match(i,removeIndex))) {
				next()
			} else {
				for(j in (i+1):listLength)  {
					if(!is.na(match(j,removeIndex))) {	
						next()
					} else {
						if(length(intersect(startList[[i]],startList[[j]])) == length(startList[[j]])) {
							removeIndex <- c(removeIndex, j)
						}
					}
				}
			}		
		}	
		return(startList[-removeIndex])
	}
}


### NOT USED
phasedGenotypes <- function(genotypeTable, haplotypeList, haploFreqs, genoHaploMap, popOfOrigin, popFreqs, ordering) {
	
	# check table of genotypes.


	# order loci.

	# look for perfect homozygotes as these are easy to score.



	#
	validHaplos <- listValidHaplotypes()

}

# converts a table genotype with multi columns per locus to a list with one item per locus, each listing the alleles present.
tableGenoToList <- function(genoTable, locusNames)  {
	# need to check if only one row in genoTable
	
	genoList <- list()
	for(thisLocus in locusNames)  {
		locusColumns <- grep(thisLocus, names(genoTable))
		genoList[[thisLocus]] <- as.character(unlist(genoTable[,locusColumns]))
	}
	return(genoList)

}

# converts a table (with rownames) of haplotypes to a list of haplotypes.
tableHaploToList <- function(haploTable, locusNames=colnames(haploTable))  {
	haploList <- list()
	for(thisRow in 1:nrow(haploTable)) {
		thisList <- list()
		thisHaplo <- rownames(haploTable)[thisRow]
		for(thisLocus in locusNames)  {
			#locusColumns <- grep(thisLocus, names(haploTable))
			thisList[[thisLocus]] <- as.character(unlist(haploTable[thisRow,thisLocus]))
		}
		haploList[[thisHaplo]] <- thisList
	}
	return(haploList )
	
}

# tableHaploToList(haplotypes[listValidHaplotypes(thisGenotype, haplotypes ),])

listGenoToTable <- function(genoList) {
	exportTable <- data.frame()
	for(i in 1:length(genoList)) {
		thisGeno <- genoList[i]
		genoTable <- as.data.frame(thisGeno)
		thisRow <- data.frame(dummy=NA)		# require dummy column to use cbind with first row.
		for(j in 1:nrow(genoTable)) {
			partRow <- genoTable[j,]
			colnames(partRow) <- paste(colnames(genoTable),j,sep=".")
			#print(partRow)
			thisRow <- cbind(thisRow,partRow)
			#print(thisRow)
			
		}
		thisRow <- subset(thisRow, select=-dummy)	# need to remove dummy column again.
		exportTable <- rbind(exportTable,thisRow)
	}
	exportTable <- exportTable[,order(names(exportTable))]

	return(exportTable)
}

#listGenoToTable(list(tableGenoToList(thisGenotype, locusNames=c("A","B","C"))))


listHaploToTable <- function(haploList)  {
	exportTable <- data.frame()
	for(i in 1:length(haploList)) {
		thisHaplo  <- haploList[i]
		thisRow <- data.frame(thisHaplo[[1]])
		rownames(thisRow) <- names(thisHaplo)
		exportTable <- rbind(exportTable , thisRow)
	}
	return(exportTable)
}

#listHaploToTable(fullHaploList)