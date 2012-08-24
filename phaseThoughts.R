


###############CURRENT RUN EXAMPLES
source("C:/Users/Dave/HalfStarted/mlPhaser/mlPhaser.R")



haplotypes <- data.frame(	A= c("a","b","c","a","b","c","b"),
					B= c("a","b","c","b","c","a","a"),
					C= c("a","b","c","b","c","a","a") )
		
rownames(haplotypes) <- apply(haplotypes, 1, paste,sep="" , collapse="")
# test for uniqueness.

haploFreqs <- c(0.4, 0.3, 0.15, 0.07,0.05, 0.02, 0.01)
names(haploFreqs) <- rownames(haplotypes)

# have separate 

#### N.B. THIS IS DIFFERENT EVERY TIME.
my.genotypes <- simGenoFromHaplo(haploTable=haplotypes, haploFreqs=haploFreqs , 20) 

source("C:/Users/Dave/HalfStarted/mlPhaser/mlPhaser.R")
thisGenotype <- my.genotypes[1,]
thisGenotype <- my.genotypes[11,]		# was double heterozygote when developing (random result).
 listValidHaplotypes(thisGenotype, haplotypes )

test <- listValidHaplotypes(thisGenotype, haplotypes )


getValidHaploGroups(thisGenotype , haplotypes=haplotypes)

result <- getValidHaploGroups(thisGenotype , haplotypes=haplotypes)

length(result)	# how many haplotype combinations are there?
printHaploGroups(result )	# what are the combinations (pretty formatting).


## need to match haplotypes with haplofreqs. 

fullHaploList <- tableHaploToList(haplotypes)		# could use this within above functions.


namedGroups <- haploGroupsNamed(result,fullHaploList )
printHaploProbs(namedGroups,haploFreqs )


# test to see if can work through all genotypes.

for(i in 1:nrow(my.genotypes)) {
	thisGenotype <- my.genotypes[i,]
	print(thisGenotype)
	result <- getValidHaploGroups(thisGenotype , haplotypes=haplotypes)
	namedGroups <- haploGroupsNamed(result,fullHaploList )
	printHaploProbs(namedGroups,haploFreqs )
}


phaseReport <- function(genotypes,haplotypes,haploFreqs,outFormat="all")   {
	## do conversion for genotypes and haplotypes.
	genoTable <- genotypes
	haploTable <- haplotypes	


	phaseResults <- data.frame()
	for(i in 1:nrow(genoTable)) {
		thisGenotype <- genoTable[i,]
		genoName <- rownames(thisGenotype)
		result <- getValidHaploGroups(thisGenotype , haplotypes=haploTable)
		for(j in 1:length(result))  {
			
			#for(k in 1:length(names(result[[j]])))  {			
			

			#}	
			thisLH=NA	
			#thisLH <- getHaploGroupProb()
			thisRow <- data.frame(id=genoName,n.validGroups = length(result),haploGroup=paste(sort(names(result[[j]])), collapse="/"), likelihood=thisLH)
			phaseResults <- rbind(phaseResults, thisRow)
		}
			
	}
	
	phaseResults <- phaseResults[order(phaseResults$id,phaseResults$likelihood),]
	if(outFormat=="top")  {
		phaseResults <- phaseResults[match(unique(phaseResults$id),phaseResults$id),]
	}
	return(phaseResults )
}

test <- phaseReport(genotypes=my.genotypes,haplotypes)
test <- phaseReport(genotypes=my.genotypes,haplotypes, outFormat="top")

## TODO
# Clean up!
# Think of a package name: mlPhaser, haplo.phaser, phaser, 

# Decide on use of lists or dataframes for haplo and genos.
# 	would be helpful to keep track of allele names
# 	Would be nice to allow users to supply either data frames or lists (with parameter).
# 	Then convert to 'list' internally if necessary. (would that affect function defintiion?)
# 	TODO use lists and try to use list[x] rather than list[[x]] to pass list objects -> better preserves names.
#	Also use c() to combine create new lists.
# Remove hard-coded parts (e.g. loci_names)
# create tools to allow user to build haplotype tables/lists (e.g. from mlgt output)
# Package!

# How will it cope with HLA, where mutltiple haplotypes have same 'allele' at sub-locus and 
#	the allele is defined/named by it's sequence. Should be ok if translate known haplotypes 
#	into common name-set. Could use concatentated allele names, variant names or even sequence!

# general probability framework so that users can add in arbitrary probability lists or dataframes of values. 
# 	would be able to weight each term. What if conditional? e.g. rely more on Chinese frequency IF ethnicity == Chinese.
#	Sounds like a formula. 

#  Does ploidy >2 ability have other advantages. e.g. for data where loci are confused?
# 	perhaps not if can't distinguish haplotypes. 

# is the order of loci important? Not until recombination is involved (which it won't be for a while). 


# how to deal with missing data (NA)?


# how to deal with unknown haplotypes. e.g. when remGeno fails.  Perhaps test if genotype is solvable with know haplotypes before starting?
# However, some may appear solvable but not be e.g.   ac,ab,cb   when available haplos are abb, aac, cbb
# Might need to add in code to remgeno or recursion to allow for failure. 
#Currently should fail (gracefully) to find combinations and provide empty list.


## similar functions in existing packages
# haplo.stats
#	summaryGeno()
#	get.geno.pairs()


#what about X chromosome?  second allele as NULL?


## Using sequence. Should be able to get sequence directly within a custom callGenotypes function as it is held by the mlgtResult.


############################################################################
################ Attempt to process mlgt data.##############################
### Overview (because it seems long and complicated currently)

# define the haplotype (which loci or sub-loci)
# Obtain/build a knownAlleleDb
# build haplotype table: loci as columns, haplotype names as rownames, alleles or sequences as data.
# summarise haplotype table? (does it matter if not all unique? - currently permissive of duplicates on the presumptin that they will later be ranked on external frequency data).

# obtain haplotype frequency data.  
# TODO have multiple columns (or tables?) for different lines of evidence?
#		e.g. if want to link ethnicity of sample to frequencies.

# load mlgt library and set BLAST/MUSCLE (last bit if going to produce knownAlleleDb).
# Call genotypes on an mlgtResult object
# Get sequence from genotypes for each locus and allele as seq.1, seq.2 etc (or locusName.1, locusName.2 etc).
# Refine based on genotype status (HOMOZYGOTE/HETEROZYGOTE/tooFewReads/complexVars)
# Combine genotypes for all loci within the haplotype into a single table. Multiple columns per locus giving distinc allele data (homozygotes repeated).

# load haploPhaser code (library)
# call getValidHaploGroups() on each genotype and retain the results in useable format.
#	may be 0,1 or more valid haplotype group combinations.
# score/rank resulting groups based on haplotype frequencies (TODO)
# output the results. 	(TODO)



##  could use custom callGenotypes to provide sequence as allele name. NO, could need to be replacement of 
## the master function as callGenotypes.default only received the table of frequencies, not the full mlgtResult object.

library(mlgt)
Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")

setwd("C:/Users/Dave/HalfStarted/mlgt/testProject/testPhasing")
#data("mlgtResult", package="mlgt") 	# no pairs of loci.
#my.mlgt.Result

#load("C:/Users/Dave/HalfStarted/mlgt/testProject/cleanRun/thisRun.mlgtResult.Rdata")

load("C:/Users/Dave/NextGen/DNAseqLab/TestProject/TestRun_09_03_2012_mlgt0.14/my.mlgt.Result.RData")
my.mlgt.Result.HLA_09_03_2012

myMarkerList <- read.fasta("C:/Users/Dave/NextGen/DNAseqLab/HLA_29_02_2012/HLA_HR_Markers_20120320.fasta", as.string=T)	
#myPhasingMarkers <- myMarkerList[c("HLA_A2", "HLA_A3")]
myPhasingMarkers <- myMarkerList[c("HLA_B2", "HLA_B3")]
myPhasingMarkers <- myMarkerList[c("HLA_C2", "HLA_C3")]
# the default method to callGenotypes
#my.genoytpes <- callGenotypes(my.mlgt.Result)


## need to build genotype table with 2 genotypes per sample per marker.
## Control NAs etc. 
## TODO build in NA control into the phaser

## Need to convert haplotypes to haplotype tables with column for each marker



#markerImgtFileTable <- read.delim(system.file("marker.imgt.msf.list.tab", package="mlgt"),
 	#			header=T)

markerImgtFileTable <- data.frame(marker=c("HLA_A2","HLA_A3"),imgtAlignFile=c("A_nuc.msf","A_nuc.msf"))
markerImgtFileTable <- data.frame(marker=c("HLA_A2","HLA_A3","HLA_B2","HLA_B3","HLA_C2","HLA_C3"),
					imgtAlignFile=c("A_nuc.msf","A_nuc.msf","B_nuc.msf","B_nuc.msf","C_nuc.msf","C_nuc.msf"))

alignFilesSource <- 'ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/'
# select a folder to store the alignments in. Here using current working directory.
 alignFilesDir <- getwd()	
 ## Download the allele alignments and create a 'variantMap' object for each marker and store them all in a list.
knownAlleleDb <- list()
for(thisMarker in names(myPhasingMarkers)) {	
	fileName <-  markerImgtFileTable$imgtAlignFile[match(thisMarker, markerImgtFileTable$marker)]
	alleleAlignUrl  <- paste(alignFilesSource , fileName , sep="/")
	alleleAlignFile <- paste(alignFilesDir , fileName , sep="/")
	download.file(alleleAlignUrl,alleleAlignFile)
 	knownAlleleDb[[thisMarker]] <- createKnownAlleleList(thisMarker,
 		myPhasingMarkers[[thisMarker]][1], alleleAlignFile)
}

#save(knownAlleleDb, file="knownAlleleDb.testPhasing.HLA_B.Rdata")
#save(knownAlleleDb, file="knownAlleleDb.testPhasing.Rdata")

#load(

### reformat the knownAlleleDb into a haploTypeTable.

hlaAlleTableList <- list()

for(thisMarker in names(myPhasingMarkers))  {

	#thisMarker <- "HLA_A2"
	locusAlleleTable <- data.frame()
	for(i in 1:length(knownAlleleDb[[thisMarker]]@variantMap))  {
		thisVariant <- knownAlleleDb[[thisMarker]]@variantMap[i]
		alleleNames <- unlist(strsplit(thisVariant[[1]],"_", fixed=T))
		sequence <- names(thisVariant[1])
		tempTable <- data.frame(allele=alleleNames,sequence=sequence)
		locusAlleleTable <- rbind(locusAlleleTable,tempTable)
	}
	names(locusAlleleTable) <- c("allele", thisMarker) 
	hlaAlleTableList[[thisMarker]] <- locusAlleleTable
}

## Need to convert haplotypes to haplotype tables with column for each marker
## TODO: generalise this part.
hlaHaploTable <- merge(hlaAlleTableList[["HLA_C2"]], hlaAlleTableList[["HLA_C3"]], by="allele")
#hlaHaploTable <- merge(hlaAlleTableList[["HLA_B2"]], hlaAlleTableList[["HLA_B3"]], by="allele")
#hlaHaploTable <- merge(hlaAlleTableList[["HLA_A2"]], hlaAlleTableList[["HLA_A3"]], by="allele")
rownames(hlaHaploTable) <- hlaHaploTable$allele
hlaHaploTable <- subset(hlaHaploTable, select=-allele)
 str(hlaHaploTable)



## need to test if all haplotypes are unique (probably not if alleles defined
## by region not within markers selected).
### N.B. perhaps not. If want to use frequencies, then leave all combinations in. Even identical ones.
## this part was to test for unique haplotyeps. Might not be necessary
#hlaHaploTable$fullHaplo <- ""
#for(thisMarker in names(myPhasingMarkers))  {
#	hlaHaploTable$fullHaplo <- paste(hlaHaploTable$fullHaplo, hlaHaploTable[,thisMarker ])
#}
#hlaHaploTable$fullHaplo <- paste(hlaHaploTable[,names(myPhasingMarkers)])
#length(hlaHaploTable$fullHaplo)
#length(unique(hlaHaploTable$fullHaplo))




## need to callGenotypes, leaving behind sequence or haplotype name as genotype allele.
## Sequence will be compatible with haploTable method above. 
## Give NA if NA. What about homozygotes?  

my.genotypes.hla <- callGenotypes(my.mlgt.Result.HLA_09_03_2012, mapAlleles=TRUE,
		alleleDb=knownAlleleDb)

names(my.genotypes.hla)

my.genotypes.hla[thisMarker]
alleleNumber <- 2
# add alllele sequence to genotypeTables.
for(thisMarker in names(myPhasingMarkers))  {
	for(k in 1:alleleNumber) {
		my.genotypes.hla[[thisMarker]]@genotypeTable[,paste("seq",k,sep=".")] <-""
		for(i in 1:nrow(my.genotypes.hla[[thisMarker]]@genotypeTable) ) {
			seq <- names(which(knownAlleleDb[[thisMarker]]@variantMap == my.genotypes.hla[[thisMarker]]@genotypeTable[i,paste("allele",k,sep=".")]))
			my.genotypes.hla[[thisMarker]]@genotypeTable[i,paste("seq",k,sep=".")] <- ifelse(length(seq) > 0, seq, NA)
		}
	}
}

#my.genotypes.hla[[thisMarker]]@genotypeTable


myTempGenotypes <- list()
for(thisMarker in names(myPhasingMarkers))  {
	# sample id as rownames and one column per allele per marker (e.g. 2 markers/2 alleles = 4)
	myTempGenotypes[[thisMarker]] <- my.genotypes.hla[[thisMarker]]@genotypeTable
	
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
}
# need to generalise this section
#my.genotypes.hla.refined <- merge(myTempGenotypes[['HLA_A2']],myTempGenotypes[['HLA_A3']], by="sample")
#my.genotypes.hla.refined <- merge(myTempGenotypes[['HLA_B2']],myTempGenotypes[['HLA_B3']], by="sample")
my.genotypes.hla.refined <- merge(myTempGenotypes[['HLA_C2']],myTempGenotypes[['HLA_C3']], by="sample")

rownames(my.genotypes.hla.refined) <- my.genotypes.hla.refined[,'sample']
my.genotypes.hla.refined <- subset(my.genotypes.hla.refined, select=-sample)


## need to consolidate table to use proper sequence for homozygotes. 
## NA for non-homozygotes/hetero


## filter for complete.cases?



## should now be callable with:-
hlaHaploTable 
my.genotypes.hla.refined 

listValidHaplotypes(my.genotypes.hla.refined[2,] , hlaHaploTable )

getValidHaploGroups(my.genotypes.hla.refined[3,] , hlaHaploTable)
getValidHaploGroups(my.genotypes.hla.refined[5,] , hlaHaploTable)

result <- getValidHaploGroups(my.genotypes.hla.refined[1,], haplotypes=haplotypes)

for(i in 1:nrow(my.genotypes.hla.refined))  {
	print(my.genotypes.hla.refined[i,])
	result <- getValidHaploGroups(my.genotypes.hla.refined[i,] , hlaHaploTable)
	print(length(result))
	print(result)
}


# could have used either dna sequence or list of shared alleles as 'allele name'


 #my.mlgt.Result.HLA_09_03_2012@markerSampleList['HLA_A3']


#names(knownAlleleDb)

## next need to add in frequencies.
## Looked at alleleFrequencies.net. Haplotype there are often defined across mulitple loci. 
## Need to do some work on this.
## e.g. which populations. How many loci per haplotype. 
## What resolution to compare to. e.g.  A*01:01:01:01 or just A*01:01
## May depend on how many calls come back as ambiguous or novel. 
##
## A decision needs to be made as to how to class 'novel' calls. Are they novel or are they sequencing/amplification errors?


## also structure the results better.
## How many valid groups were there (usually 0 or 1).
## Valid group definition: a set of known haplotypes that fully explain the observed genotype.

# zero to multiple valid combinations per sample. 
#	return all valid combinations. 
# 	rank multiple valid combinations based on rank (likelihood).
# If mixture of valid and unknown haplotypes, then none will be given.








stopifnot()
################################# TESTING



### Test what happens when genotype appears to have valid haplotypes, but doesn't really



## testing for valid haplotypes -> may return incompatible pair


ac,bb,cb   when available haplos are abb, aac, cbb

failHaplos <- data.frame( A=c("a","a","c"), B=c("b","a","b"), C=c("b","c","b"), stringsAsFactors=F )
rownames(failHaplos) <- c("abb","aac","cbb")

failGeno <- data.frame(A.1="a",A.2="c",B.1="b",B.2="b",C.1="c",C.2="b", stringsAsFactors=F)		


listValidHaplotypes(failGeno, failHaplos)		# lists all three as valid haplos.
listValidHaplotypes(thisGenotype, failHaplos)
listValidHaplotypes(failGeno, haplotypes)

test <- remGeno(tableHaploToList(failHaplos[3,])[[1]],tableGenoToList(failGeno, locusNames=c("A","B","C")))	
# now try to remove a haplotype that is not in the genotype. Returns passTable all FALSE
remGeno(tableHaploToList(failHaplos[3,])[[1]],test$remList)	 ## some passTable false.

###

result.fail <- getValidHaploGroups(failGeno , haplotypes=failHaplos)		# produced empty list.

fullHaploList.fail <- tableHaploToList(failHaplos)	


namedGroups <- haploGroupsNamed(result.fail,fullHaploList.fail)
















############################### DEVELOPMENT





### how to give probabilities to phased haplotypes.



## what info to combine and how.

point genotypes 
locusA	allele1	allele2	# what if >2 alleles?
locusB	allele1	allele2

common haplotypes
A1:B1  50%
A2:B2	 35%
A1:B2  10%
A2:B1  5%

qualify by population?
haplotype freqs by populations
Certainty over population membership?


## how many loci to combine (depends on prior info). Need general framework and map.
## e.g. if want to combine Aex2, Aex3, Bex2 & Cex3, need map which gives order and frequencies of haplotypes.

### Provision for fully interanl assessment? e.g. using frequencies within data. 
## becomes almost impossible to allow recombination for rare variants. 


## HOw to combine statistically?


## homozygote/heterozygote uncertainty.


## how to encode/present the input

genotypeTable

genotype <-> haplotype map  (1-to-many)

number of possible haplotypes:-

L loci with n1.. nL genotypes.

e.g. Four loci with 3,6,3,5 alleles = product
genoCounts <- c(3,6,3,5)
haplo.poss.count <- prod(genoCounts)		# 270
genoCounts <- c(300,60,30,50)
haplo.poss.count <- prod(genoCounts)		# 2.7e+07



## Individuality of genotypes/haplotypes is somewhat relative. e.g. protein-coding differences are more important than synonymous changes.
## Depends on ends use.


## How to encode/present the results. 





test <- tableGenoToList(thisGenotype, locusNames=c("A","B","C"))



## genotypes based on 

haploFreqsTable <- data.frame(haplotype=c("aaa", "bbb", "ccc", "ddd", "aac", "dbb", "cca","bcd", "bca", "acd", "cdb"),freq=0.1)




## random genotypes
my.genotypes <- data.frame(	A.1=sample(c("a","b","c","d"),20, replace=T), 
					A.2=sample(c("a","b","c","d"),20, replace=T),
					B.1=sample(c("a","b","c","c","d"),20, replace=T),
					B.2=sample(c("a","b","c","c","d"),20, replace=T),
					C.1=sample(c("a","b","b","b","c","d"),20, replace=T),
					C.2=sample(c("a","b","b","b","c","d"),20, replace=T)
)




#my.phased <- phasedGenotypes(my.genotypes)




## dictate haplotype map
## each haplotype can have only one genoytpe per locus
## each genotype can map to multiple loci.

# so haplotype table should link genotypes across loci.
# rownames for names





build intersect list for each ploidy
set haplotype of first to focal haplotype.




#### TODO validHaplos is not enough. Need to also find groups of haplos consistent with genotypes.
## e.g. "aaa" && "bcc"  for
#  A.1 A.2 B.1 B.2 C.1 C.2
#1   a   b   a   c   a   c




thisGenotype <- my.genotypes[11,]
#haploList <- list(haplotypes[listValidHaplotypes(thisGenotype, haplotypes )[1],], haplotypes[listValidHaplotypes(thisGenotype, haplotypes )[2],])
haploList <- tableHaploToList(haplotypes[listValidHaplotypes(thisGenotype, haplotypes ),])

startGroup=list()
startRemGeno <- list()
startRemGeno[['remList']] <- tableGenoToList(thisGenotype, locusNames=c("A","B","C"))
recurseHaplos(validHaplotypes=haploList,startRemGeno[['remList']] ,
				group=startGroup )

# working but giving too much output. Might need to strip out levels?

### need special env to store valid haplo groups.

groupStorageEnv <- new.env(parent = baseenv())
assign("validHaploGroups", list(), groupStorageEnv)
# then call function
# then retrieve validGroups.
recurseHaplos(validHaplotypes=haploList,startRemGeno[['remList']] ,
				group=startGroup )
validHaploGroups <- get("validHaploGroups", groupStorageEnv)

## need to retrieve this result and remove the environment
#rm(groupStorageEnv)

## then need to test for duplicated haplotype groups.

validHaploGroups <- reduceRedundantList(validHaploGroups )
length(reduceRedundantList(validHaploGroups ))


# now need a function to do this all together. 

thisHaplotype <- haplotypes[listValidHaplotypes(thisGenotype, haplotypes )[1],]
remGeneList <- remGeno(thisHaplotype , tableGenoToList(thisGenotype, locusNames=c("A","B","C")))
thisHaplotype <- haplotypes[listValidHaplotypes(thisGenotype, haplotypes )[2],]
remGeneList <- remGeno(thisHaplotype , remGeneList[['remList']])

testHaploInGeno(thisHaplotype , remGeneList)
