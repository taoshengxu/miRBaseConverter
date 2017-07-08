#' miRBase version convert for miRNA Names
#'
#' This function converts a group of any species' miRNA names (including precursor and mature miRNA) to the specified miRBase version if the miRNAs have been defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param targetVersion A character value representing the target miRBase version corresponding the source miRNA names.
#' Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @param exact Logical value. If true, the result will be the most exactly matched result.
#'  If FALSE, the result will include all the possible matched miRNA name. If one miRNA can match multiple names. All the matched names are
#'   concatenated with "&".
#' @param verbose Logical value. If true, it will print the multiple matched miRNA Names and Accessions to the console.
#' @note Please note: Due to some miRNA names changing many times in history. Even if choose the third parameter "exact"=TRUE, it may still have some miRNAs that can't
#' match the unique name in the target version. In order to return the accurate result as possible, we also concatenate the multiple matched miRNA names with "&". This is the rare case but it happens sometimes.
#'
#' @return
#' A data frame with three columns. The number of rows equal to the input miRNA names. The three columns are defined as below:
#'\itemize{
#'  \item  \strong{OriginalName} : The original miRNA names (Column 1).\cr
#'  \item \strong{TargetName} : The converted miRBase names (in specified version) corresponding to the original miRNA names (Column 2).\cr
#'  \item \strong{Accession} : The corresponding miRBase Accessions (Column 3).
#'}
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result1=miRNAVersionConvert(miRNANames,targetVersion="v13",exact=TRUE,verbose=TRUE)
#' result2=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE,verbose=TRUE)
#' result3=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=FALSE,verbose=TRUE)
#'
#' miRNANames=c( "hsa-let-7c","hsa-miR-3190-3p","hsa-let-7c","hsa-miR-34b","hsa-miR-378",
#' "hsa-miR-499a-3p","hsa-miR-499a-5p","hsa-miR-500","hsa-miR-516a-5p","hsa-miR-550","hsa-miR-589")
#' result4=miRNAVersionConvert(miRNANames, targetVersion="v21", exact=TRUE, verbose=TRUE)
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
miRNAVersionConvert <-function(miRNANames,targetVersion="v21",exact=TRUE,verbose=TRUE)
{
  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space

  VMAP <- .getmiRNAs(targetVersion)[, c("ACC", "SYM")]

  uid = unique(as.vector(miRNANames))
  SYM_ID = match(uid, SYM)
  df <- data.frame(uid = uid, SYM = SYM_ID, stringsAsFactors=FALSE)
  df <- merge(df, ACC_SYM)[, c("uid", "ACC")]
  df <- unique( merge(df, VMAP, by="ACC") )

  target <- data.frame(
    OriginalName = df$uid,
    TargetName = SYM[df$SYM],
    Accession = ACC[df$ACC],
    stringsAsFactors = FALSE
  )
  if (exact) {
    idx <- (target$OriginalName == target$TargetName) | (!target$OriginalName %in% target$TargetName)
    target <- target[idx, , drop=FALSE]
  }
  ## collapse 1:many maps
  splitpaste <- function(x, f) {
    result <- vapply(split(x, f), paste, character(1), collapse="&")
    result[!nzchar(result)] <- NA
    result
  }
  f <- factor(target$OriginalName, levels=uid)
  target <- data.frame(
    OriginalName = uid,
    TargetName = splitpaste(target$TargetName, f),
    Accession = splitpaste(target$Accession, f),
    row.names=NULL, stringsAsFactors = FALSE
  )
  if(verbose)
  {
    multiindex=grep('&',target$TargetName)
    if(length(multiindex)>0)
    {
      output=target[multiindex,]
      colnames(output)<-c("original",paste("Version",targetVersion),"Accession")
      rownames(output)<-NULL
      message("********************************************\n")
      message("The multiple matched miRNAs are list below: \n")
      # message("                                            \n")
      print.data.frame(output)
    }
  }
  target[match(miRNANames, target$OriginalName),]
}

#' miRBase Accession to miRNA Name in specified version
#'
#' This function converts a group of any species' miRNA Accessions (including precursor and mature miRNA) to a specified miRBase version if the Accessions have been defined in miRBase.
#'
#' @param Accessions A character vector representing the miRNA Accessions needed to be convert.
#' @param targetVersion A character value representing the target miRBase version corresponding the Accessions.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of rows equal to the input miRNA names. The two columns are defined as below:
#'\itemize{
#' \item \strong{Accession} : The Accession of miRNAs (Column 1).\cr
#' \item  \strong{TargetName} : The converted miRBase names (in specified version) corresponding to the Accessions (Column 2).\cr
#' }
#' @examples
#' data(miRNATest)
#' Accessions=miRNATest$Accession
#' result1=miRNA_AccessionToName(Accessions,targetVersion="v13")
#' result2=miRNA_AccessionToName(Accessions,targetVersion="v21")
#'
#' @author
#'  Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
miRNA_AccessionToName<-function(Accessions,targetVersion="v21")
{
  Accessions=as.character(Accessions)
  Accessions=gsub(" ","",Accessions)##Remove the possible space

  VMAP <- .getmiRNAs(targetVersion)[, c("ACC", "SYM")]

  uid=unique(as.vector(Accessions))
  ACC_ID=match(uid,ACC)
  df=data.frame(uid = uid, ACC = ACC_ID, stringsAsFactors=FALSE)
  df <- unique(merge(df, VMAP, by="ACC"))
  target <- data.frame(Accession = df$uid,TargetName = SYM[df$SYM],stringsAsFactors = FALSE)
  target <- target[match(Accessions, target$Accession),]
  target$Accession=Accessions
  rownames(target)= NULL
  target
}


#' The miRBase miRNA names with specified version to Accessions
#'
#' This function converts a group of any species' miRNA name to the Accessions defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param version A character value representing the version corresponding the miRNANames.
#' Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of rows equal to the input miRNA names. The two columns are defined as below:
#'
#' \itemize{
#'  \item \strong{miRNAName_\{Version\}} : The input miRNA names (Column 1).\cr
#'  \item \strong{Accession} : The convert Accession(Column 2).\cr
#'  }
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' version=checkMiRNAVersion(miRNANames,verbose=TRUE)
#' result1=miRNA_NameToAccession(miRNANames,version=version)
#' result2=miRNA_AccessionToName(result1[,2],targetVersion="v21")
#' result3=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#'
#' @author
#'  Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
miRNA_NameToAccession<-function(miRNANames,version="v21")
{
  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space

  VMAP <- .getmiRNAs(version)[, c("ACC", "SYM")]

  uid=unique(as.vector(miRNANames))
  SYM_ID = match(uid, SYM)
  df <- data.frame(uid = uid, SYM = SYM_ID, stringsAsFactors=FALSE)
  df <- unique( merge(df, VMAP, by="SYM"))
  target <- data.frame(miRNANames = df$uid,Accession = ACC[df$ACC],stringsAsFactors = FALSE)
  target <- target[match(miRNANames, target$miRNANames),]
  target$miRNANames=miRNANames
  rownames(target)= NULL
  colnames(target)[1]=paste("miRNAName_",version,sep="")
  target
}

#' Get the miRNA sequences
#'
#' This function returns the miRNA sequences for a list of miRNAs.
#'
#' @param Accessions A character vector representing the miRNA Accessions in miRBase.
#' @param targetVersion A character value representing the target miRBase version corresponding the Accessions.
#'  Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of row equals to input miRNAs. The two columns are defined as below:
#'
#' \itemize{
#'  \item \strong{Accession} : The original miRNA (Column 1).\cr
#'  \item \strong{miRNASequence_\{targetVersion\}} : The return miRNA sequence (in specified version) corresponding to the input miRNAs (Column 2).\cr
#'  }
#' @examples
#' #####1, The input are miRNA Accessions
#' data(miRNATest)
#' Accessions=miRNATest$Accession
#' result1=getMiRNASequence(Accessions,targetVersion="v13")
#' result2=getMiRNASequence(Accessions,targetVersion="v21")
#'
#' #####2, The input are miRNA Names
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result3=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#' Accessions=result3[,3]
#' result4=getMiRNASequence(Accessions,targetVersion="v21")
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getMiRNASequence<-function(Accessions,targetVersion="v21")
{
  Accessions=as.character(Accessions)
  Accessions=gsub(" ","",Accessions)##Remove the possible space

  VMAP <- .getmiRNAs(targetVersion)[, c("ACC", "SEQ")]

  uid=unique(as.vector(Accessions))
  ACC_ID=match(uid,ACC)
  df=data.frame(uid = uid, ACC = ACC_ID, stringsAsFactors=FALSE)
  df <- unique(merge(df, VMAP, by="ACC"))
  target <- data.frame(Accession = df$uid,miRNASequence = SEQ[df$SEQ],stringsAsFactors = FALSE)
  target <- target[match(Accessions, target$Accession),]
  target$Accession=Accessions
  rownames(target)= NULL
  colnames(target)[2]=paste("miRNASequence_",targetVersion,sep="")
  target
}


#' Get the detailed information of a single specified miRNA in all miRBase versions.
#'
#' This function returns all available miRBase versions' information of a single specified miRNA.
#'
#' @param Accession A character representing the single Accession.
#'
#' @return
#'  A data frame (21X7) including all the history information (Precursor, Mature, Sequence) of the specified miRNA.
#'  Each row represents a miRBase version.
#' @examples
#' #####1,The input is a miRNA Name
#' miRNAName="hsa-miR-26b-5p"
#' result1=miRNA_NameToAccession(miRNAName,version="v21")
#' Accession=result1[,2]
#' result2=getMiRNAHistory(Accession)
#'
#' #####2,The input is miRNA Accession
#' Accession="MIMAT0000765"
#' result3=getMiRNAHistory(Accession)
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getMiRNAHistory<-function(Accession)
{
  Accession=as.character(Accession)
  Accession=gsub(" ","",Accession)##Remove the possible space

  if(length(Accession)>1)
    stop("This function is only for a single miRNA query, Please query miRNAs one by one")

  ACC_ID=match(Accession,ACC)
  if(!is.na(ACC_ID))
  {
    target=data.frame(matrix(vector(),length(VER), 7,dimnames=list(c(), c("Precursor","PrecursorSequence","Mature1","Mature1Sequence","Mature2","Mature2Sequence","Status"))),stringsAsFactors=FALSE)
    rownames(target)<-VER
    for(i in 1:length(VER))
    {
      miRNAs=miRNA_data[[i]]
      ind1=which(miRNAs[,c(1,5,8)] ==ACC_ID , arr.ind = TRUE)[1]
      if(!is.na(ind1))
      {
        target[i,]=miRNAs[ind1,][c(2,4,6,7,9,10,3)]
      }
    }
    target$Precursor=SYM[target$Precursor]
    target$PrecursorSequence =SEQ[target$PrecursorSequence]
    target$Mature1=SYM[target$Mature1]
    target$Mature1Sequence=SEQ[target$Mature1Sequence]
    target$Mature2=SYM[target$Mature2]
    target$Mature2Sequence=SEQ[target$Mature2Sequence]
    target$Status=STA[target$Status]
    target
  }
  else
  {
    message("There is no records about this miRNA.")
    NULL
  }
}


#' check the miRNA Version in miRBase
#'
#' This function checks the most possible miRBase version for a list of miRNA names.
#' @importFrom stats na.omit
#' @param miRNANames A character vector representing the miRNA names.
#' @param verbose Logical value. If true, the detail version information is printed in the console for user reference.
#' @return
#'  A single character value or a character vector represent the most possible miRBase version for the list of miRNA names.
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' version=checkMiRNAVersion(miRNANames,verbose=TRUE)
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
checkMiRNAVersion<-function(miRNANames,verbose=TRUE)
{
  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space

  uid=unique(as.vector(miRNANames))
  SYM_ID=match(uid,SYM)
  SYM_ID=na.omit(SYM_ID)
  #result=data.frame(matrix(vector(),length(VER), 3,dimnames=list(c(), c("Version","Proportion","Recommend") )),stringsAsFactors=FALSE)
  #result$Version<-VER
  result=data.frame("Version"=VER,"Proportion"=NA,"Recommend"="",stringsAsFactors=FALSE)

  for(i in 1:length(VER))
  {
    allSYM=miRNA_data[[i]]
    allSYM=c(allSYM[,2],allSYM[,6],allSYM[,9])
    allSYM=na.omit(allSYM)
    num =length(intersect(allSYM, SYM_ID))
    result$Proportion[i]=round((num/length(uid))*100,2)
  }
  ind=which.max(result$Proportion)
  if(verbose)
  {
    result[result$Proportion == max(result$Proportion), "Recommend"] <- " ***BEST Matched***"
    result$Proportion=paste(result$Proportion,"%",sep="")
    print.data.frame(result)
  }
  VER[ind]
}

#' Open the miRBase webpages of the specified miRNAs
#'
#' This function redirects the miRBase webpage of the specified miRNAs
#' @param Accessions A character vector representing the miRNA Accessions in miRBase.
#' We restict the number of queried miRNAs each time. The maximum number of the input miRNAs is 15.
#'
#' @examples
#' #### 1. A step-loop
#' Accession1="MI0000447"
#' goTo_miRBase(Accession1)
#'
#' #### 2. A mature miRNA
#' Accession2="MIMAT0026477"
#' goTo_miRBase(Accession2)
#'
#' #### 3. A list of miRNAs
#' Accession3=miRNATest$Accession[1:10]
#' goTo_miRBase(Accession3)
#'
#' @return
#' No values
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
goTo_miRBase<-function(Accessions)
{
  Accessions=unique(Accessions)
  if(length(Accessions)>15)
  {
    message("The number of the queried miRNAs is out of 15. Please reduce the input miRNA number")
  }
  else
  {

    Accessions=as.character(Accessions)
    Accessions=gsub(" ","",Accessions)##Remove the possible space
    alive=checkMiRNAAlive(Accessions)
    for(i in alive)
    {
      if(grepl("MAT",i))
        URL=paste0("http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=",i)
      else
        URL=paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",i)
      utils::browseURL(URL)
    }
  }
}

#' Get all miRBase version information
#'
#' This fuction return a reference for all miRBase versions' information including Verson name,
#' Release date, miRNA number and Status.
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @return
#' The detailed version information is printed in the console for user reference.
#' @examples
#' getAllVersionInfo()
#' @export
#'
getAllVersionInfo<-function()
{
  print.data.frame(VER_INFO)
}

#' Get all species of miRNAs embodied in miRBase repository
#'
#' This fuction return a reference for all species of miRNAs including the
#' abbreviation and full name.
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @return
#' A dataframe
#' A data frame with two columns. The two columns are defined as below:
#' \itemize{
#'  \item \strong{Species} \cr
#'  \item \strong{FullName} \cr
#'  }
#' @examples
#' allSpecies=getAllSpecies()
#' @export
#'
getAllSpecies<-function()
{
  SPE1=SPE
  colnames(SPE1)<-c("Species","FullName")
  SPE1
}

#' Get the full miRNAs information table of the specified miRBase version
#'
#' This function returns the full miRNAs information table of the specified miRBase version
#'
#' @param version A character value representing the specified miRBase version for
#' retrieval. Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @param species A character value representing the abbreviation of species.  Users can apply
#' the function \strong{getAllSpecies()} to get the available  abbreviation of species. If species is set
#' to \strong{"all"}, the miRNAs of all species will return.
#' @examples
#' miRNA_Tab=getMiRNATable(version="v21",species="hsa")
#' @return
#' A data frame
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getMiRNATable<-function(version="v21",species="all")
{
  ver_index=match(tolower(version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  miRNA_Tab <- miRNA_data[[ver_index]]
  if(species=="all")
  {
    miRNA_Tab=miRNA_Tab
  }
  else
  {
    spe_index=match(species,SPE$SPE)
    if(is.na(spe_index))
      stop("It is a wrong species abbreviation, Please check it")
    miRNA_Tab=miRNA_Tab[which(spe_index==miRNA_Tab$Species),]
  }
  miRNA_Tab$Precursor_Acc=ACC[miRNA_Tab$Precursor_Acc]
  miRNA_Tab$Precursor=SYM[miRNA_Tab$Precursor]
  miRNA_Tab$Status=STA[miRNA_Tab$Status]
  miRNA_Tab$Precursor_Seq=SEQ[miRNA_Tab$Precursor_Seq]
  miRNA_Tab$Mature1_Acc=ACC[miRNA_Tab$Mature1_Acc]
  miRNA_Tab$Mature1=SYM[miRNA_Tab$Mature1]
  miRNA_Tab$Mature1_Seq=SEQ[miRNA_Tab$Mature1_Seq]
  miRNA_Tab$Mature2_Acc=ACC[miRNA_Tab$Mature2_Acc]
  miRNA_Tab$Mature2=SYM[miRNA_Tab$Mature2]
  miRNA_Tab$Mature2_Seq=SEQ[miRNA_Tab$Mature2_Seq]

  miRNA_Tab[,c(1:10)]
}

#' Check the miRNA status(Alive or Dead)
#'
#' This function checks the miRNA status (Alive or Dead) in the latest miRBase verison.
#'
#' @param Accessions A character vector representing the miRNA Accessions in miRBase.
#' @param verbose Logical value. If true, the dead miRNAs will be printed the console.
#' @examples
#' data(miRNATest)
#' ## The input is miRNA Accessions
#' Accessions=miRNATest$Accession
#' alive_Accession1=checkMiRNAAlive(Accessions)
#'
#' ##The input is miRNA names
#' miRNANames=miRNATest$miRNA_Name
#' version=checkMiRNAVersion(miRNANames,verbose = TRUE)
#' result=miRNA_NameToAccession(miRNANames,version = version)
#' Accessions=result$Accession
#' alive_Accession2=checkMiRNAAlive(Accessions)
#' @return
#' A  character vector of Accessions for all alive miRNAs. The names of the return vector are the position indexes in the input Accessions.
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
checkMiRNAAlive<-function(Accessions, verbose=TRUE)
{
  Accessions=as.character(Accessions)
  Accessions=gsub(" ","",Accessions)##Remove the possible space
  Accession=match(Accessions,ACC)

  all_accession=.getmiRNAs()[,1]
  index = is.element(Accession,all_accession)
  alive=Accessions[index]
  names(alive)=which(index)
  if(length(which(!index)>0) && verbose)
  {
    dead=na.omit(unique(Accessions[!index]))
    dead=data.frame("Dead.miRNAs"=dead)
    message("The dead miRNAs are in below:")
    print.data.frame(dead)
  }
  alive
}

#' Get all miRNAs in the specified miRBase version
#'
#' This function gets all miRNAs in the specified miRBase version.
#'
#' @param version A character value representing the specified miRBase version for
#' retrieval. Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @param species A character value representing the abbreviation of species.  Users can apply
#' the \strong{getAllSpecies()} function to get the available  abbreviation of species. If species is set
#' to \strong{"all"}, the miRNAs of all species will return.
#' @param type A character value representing the miRNA type for retrieval.
#' \itemize{
#' \item \strong{"precursor"}  \cr
#' \item \strong{"mature"} \cr
#' \item \strong{"all"} : precursor and mature \cr
#' }
#' @examples
#' miRNAs=getAllMiRNAs(version="v21", type="all", species="hsa")
#' @return
#' A data frame with three columns. The three columns are defined as below:
#' \itemize{
#'  \item \strong{Accession} \cr
#'  \item \strong{Name} \cr
#'  \item \strong{Sequence} \cr
#'  }
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getAllMiRNAs<-function(version="v21",type="all",species="all")
{
  ver_index=match(tolower(version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  MiRNAs <- miRNA_data[[ver_index]]
  if(species=="all")
  {
    MiRNAs=MiRNAs
  }
  else
  {
    spe_index=match(species,SPE$SPE)
    if(is.na(spe_index))
      stop("It is a wrong species abbreviation, Please check it")
    MiRNAs=MiRNAs[which(spe_index==MiRNAs$Species),]
  }
  if(type=="all")
  {
    Precursor=MiRNAs[,c(1,2,4)]
    Mature1=MiRNAs[,c(5,6,7)]
    Mature2=MiRNAs[,c(8,9,10)]
    colnames(Precursor)=c("Accession","Name","Sequence")
    colnames(Mature1)=c("Accession","Name","Sequence")
    colnames(Mature2)=c("Accession","Name","Sequence")
    MiRNAs=rbind.data.frame(Precursor,Mature1,Mature2)
  }
  else if(type=="precursor")
  {
    MiRNAs=MiRNAs[,c(1,2,4)]
    colnames(MiRNAs)=c("Accession","Name","Sequence")
  }
  else if(type=="mature")
  {
    Mature1=MiRNAs[,c(5,6,7)]
    Mature2=MiRNAs[,c(8,9,10)]
    colnames(Mature1)=c("Accession","Name","Sequence")
    colnames(Mature2)=c("Accession","Name","Sequence")
    MiRNAs=rbind(Mature1,Mature2)
  }
  else
    stop("It is a wrong miRNA type parameter, Please check it")

  #MiRNAs=na.omit(MiRNAs)
  ind <- apply(MiRNAs, 1, function(x) all(is.na(x)))
  MiRNAs=MiRNAs[-ind,]
  MiRNAs=unique(MiRNAs)

  MiRNAs$Accession=ACC[MiRNAs$Accession]
  MiRNAs$Name=SYM[MiRNAs$Name]
  MiRNAs$Sequence=SEQ[MiRNAs$Sequence]

  MiRNAs=MiRNAs[order(MiRNAs$Name),]
  rownames(MiRNAs)=NULL

  MiRNAs
}

#' Convert the precursors to the corresponding mature miRNAs
#'
#' This function converts the precursors to the corresponding mature miRNAs in the specified miRBase version.
#' @param miRNANames A character vector representing the miRNA names.
#' @param version The default is \strong{NULL} representing the most possible latest version of the input miRNA Names will be checked automatically.
#' Otherwise, a character value representing the version corresponding the input miRNA Names.
#' Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21".\cr
#' @return
#' A data frame(nx3). The number of rows equal to the input miRNA names. The three columns are defined as below:
#' \itemize{
#'  \item \strong{OriginalName} : The input miRNA Names.\cr
#'  \item \strong{Mature1} : The corresponding mature miRNAs (always "-5p") .\cr
#'  \item \strong{Mature2} : The corresponding mature miRNAs (always "-3p") .\cr
#'  }
#'
#' @note
#'  If the input miRNA Names mingle mature miRNA names, the mature miRNA names will match to themselves in the output.
#'
#' @examples
#' miRNANames=c("pma-mir-100a","sko-mir-92a","hsa-mir-6131","mtr-MIR2655i",
#' "mmu-mir-153","mtr-MIR2592am","mml-mir-1239","xtr-mir-128-2","oan-mir-100",
#' "mmu-mir-378b","hsa-miR-508-5p","mmu-miR-434-3p")
#' result=miRNA_PrecursorToMature(miRNANames)
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
miRNA_PrecursorToMature<-function(miRNANames,version=NULL)
{
  if(is.null(version))
  {
    message("miRNA version check information:")
    c_version=rev(checkMiRNAVersion(miRNANames,verbose = FALSE))[1]
  }
  else
    c_version=version

  ver_index=match(tolower(c_version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  VMAP <-miRNA_data[[ver_index]][,c(2,6,9)]
  VMAP[,1]=SYM[VMAP[,1]]
  VMAP[,2]=SYM[VMAP[,2]]
  VMAP[,3]=SYM[VMAP[,3]]

  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space
  uid = unique(as.vector(miRNANames))

  ind=apply(VMAP,2,function(x){match(uid,x)})
  ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),2]
  ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),3]

  target <- data.frame(
    OriginalName = uid,
    Mature1 = VMAP[ind[,1],2],
    Mature2 = VMAP[ind[,1],3],
    row.names=NULL, stringsAsFactors = FALSE)
  target=target[match(miRNANames, target$OriginalName),]
}

#' Convert the mature miRNAs to the corresponding precursors
#'
#' This function converts the mature miRNAs to the corresponding precursors in the specified miRBase version.
#' @param miRNANames A character vector representing the miRNA names.
#' @param version The default is \strong{NULL} representing the most possible latest version of the input miRNA Names will be checked automatically.
#' Otherwise, a character value representing the version corresponding the input miRNA Names.
#' Users can apply the function \strong{getAllVersionInfo()} to get the available miRNA version names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1","v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21".\cr
#' @return
#' A data frame(nx2). The number of rows equal to the input miRNA Names. The two columns are defined as below:
#' \itemize{
#'  \item \strong{OriginalName} : The input miRNA Names.\cr
#'  \item \strong{Precursor} : The corresponding precursors of the mature miRNAs.\cr
#'  }
#'
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result=miRNA_MatureToPrecursor(miRNANames)
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export


miRNA_MatureToPrecursor<-function(miRNANames,version=NULL)
{
  if(is.null(version))
  {
    message("miRNA version check information:")
    c_version=rev(checkMiRNAVersion(miRNANames,verbose = FALSE))[1]
  }

  else
    c_version=version

  ver_index=match(tolower(c_version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  VMAP <-miRNA_data[[ver_index]][,c(2,6,9)]
  VMAP[,1]=SYM[VMAP[,1]]
  VMAP[,2]=SYM[VMAP[,2]]
  VMAP[,3]=SYM[VMAP[,3]]

  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space
  uid = unique(as.vector(miRNANames))

  ind=apply(VMAP,2,function(x){match(uid,x)})
  ind[which(is.na(ind[,2])),2]=ind[which(is.na(ind[,2])),3]
  ind[which(is.na(ind[,2])),2]=ind[which(is.na(ind[,2])),1]


  target <- data.frame(
    OriginalName = uid,
    Precursor = VMAP[ind[,2],1],
    row.names=NULL, stringsAsFactors = FALSE)
  target=target[match(miRNANames, target$OriginalName),]
}

#' Check the miRNA family
#'
#' This function checks the miRNA family for a list of miRNA Names.
#'
#' @param Accessions A character vector representing the miRNA Accessions in miRBase.
#' @examples
#' data(miRNATest)
#' ## The input is miRNA Accessions
#' Accessions=miRNATest$Accession
#' Family_Info1=checkMiRNAFamliy(Accessions)
#'
#' ##The input is miRNA names
#' miRNANames=miRNATest$miRNA_Name
#' version=checkMiRNAVersion(miRNANames,verbose = TRUE)
#' result=miRNA_NameToAccession(miRNANames,version=version)
#' Accessions=result$Accession
#' Family_Info2=checkMiRNAFamliy(Accessions)
#'
#' @return
#'  A data frame with four columns. The number of rows equal to the input Accessions. The four columns are defined as below:
#' \itemize{
#'  \item  \strong{Accession} : The input miRNA accessions.\cr
#'  \item \strong{miRNAName_v21} : The miRNA names (version 21) corresponding to the Accession.\cr
#'  \item \strong{FamliyAccession} : The accession of the family .\cr
#'  \item \strong{Family} : The family name.\cr
#' }
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
checkMiRNAFamliy<-function(Accessions)
{
  Accessions=as.character(Accessions)
  Accessions=gsub(" ","",Accessions)##Remove the possible space

  uid=unique(as.vector(Accessions))
  uid=na.omit(uid)
  miRNANames=miRNA_AccessionToName(uid)
  ACC_ID=match(uid,ACC)

  VMAP <-miRNA_data[[length(VER)]][,c(1,5,8,12)]
  ind=apply(VMAP,2,function(x){match(ACC_ID,x)})
  ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),2]
  ind[which(is.na(ind[,1])),1]=ind[which(is.na(ind[,1])),3]

  xx=FAM[VMAP$FAM[ind[,1]],]
  target=data.frame(miRNANames,xx,stringsAsFactors=FALSE)
  target=target[match(Accessions,target$Accession),]
  rownames(target)=NULL
  colnames(target)=c("Accession","miRNAName_v21","FamilyAccession","Family")
  target
}

#' Open the miRNA family webpages of the specified miRNAs
#'
#' This function redirects the miRBase miRNA family webpages of the specified miRNA families
#'
#' @param FamilyAccessions A character vector representing the miRNA family Accessions in miRBase.
#' We restict the queried number of miRNA family each time. The maximum number of the input miRNA families is 15.
#' @param verbose Logical value. If true, the invalid miRNA Family will be printed the console.
#'
#' @examples
#' data(miRNATest)
#' Accessions=miRNATest$Accession
#' Family_Info=checkMiRNAFamliy(Accessions)
#' FamilyAccessions=Family_Info$FamilyAccession[1:15]
#' goTo_miRNAFamily(FamilyAccessions)
#' @return
#' No values
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
goTo_miRNAFamily<-function(FamilyAccessions, verbose = TRUE)
{
  FamilyAccessions=unique(FamilyAccessions)
  FamilyAccessions=na.omit(FamilyAccessions)
  if(length(FamilyAccessions)>15)
  {
    message("The number of the queried miRNA familys is out of 15. Please reduce the input miRNA family number")
  }
  else
  {
    ind=match(FamilyAccessions,FAM$FAM_ACC)
    ind1=which(is.na(ind))
    ind2=which(!is.na(ind))
    dead=FamilyAccessions[ind1]
    alive=FamilyAccessions[ind2]
    if(length(dead)>0 && verbose)
    {
      dead=data.frame("Invalid.miRNAFamily"=dead)
      message("The invalid Family Accessions are in below: ")
      print.data.frame(dead)
    }
    for(i in alive)
    {
      URL=paste0("http://www.mirbase.org/cgi-bin/mirna_summary.pl?fam=",i)
      utils::browseURL(URL)
    }
  }
}

###internal function
.getmiRNAs<-function(version="v21")
{
  ver_index=match(tolower(version),VER)
  if (is.na(ver_index))
    stop("It is a wrong version name, Please check it")
  MiRNAs <- as.matrix(miRNA_data[[ver_index]])
  MiRNAs=rbind(MiRNAs[,c(1,2,4)],MiRNAs[,c(5,6,7)],MiRNAs[,c(8,9,10)])
  colnames(MiRNAs)=c("ACC","SYM","SEQ")
  #MiRNAs=na.omit(MiRNAs)
  ##check the rows with all NA
  ind <- apply(MiRNAs, 1, function(x) all(is.na(x)))
  MiRNAs=MiRNAs[-ind,]

  MiRNAs=unique(MiRNAs)
  MiRNAs=as.data.frame(MiRNAs)
  MiRNAs
}
