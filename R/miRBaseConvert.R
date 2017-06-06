#' miRBase name version convert
#'
#' This function converts a group of any species' miRNA names (including precursor and mature miRNA) to the specified miRBase version if the miRNAs have been defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param targetVersion A character value representing the target miRBase version corresponding the source miRNA names.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @param exact logical value. If true, the result will be the most exactly matched result.
#'  If FALSE, the result will include all the possible matched miRNA name. If one miRNA can match multiple names. All the matched names are
#'   concatenated with "&".
#' @param verbose logical value. If true, it will print the multiple matched miRNA names and Accession IDs to the console.
#' @note Please note: Due to some miRNA names changing many times in history. Even if choose the third parameter "exact"=TRUE, it may still have some miRNAs that can't
#' match the unique name in the target version. In order to return the accurate result as possible, we also concatenate the multiple matched miRNA names with "&". This is the rare case but it happens sometimes.
#'
#' @return
#' A data frame with three columns. The three columns are defined as below:
#'\itemize{
#'  \item  \strong{OriginalName} : The original miRNA names (Column 1).\cr
#'  \item \strong{TargetName} : The converted miRBase names (in specified version) corresponding to the original miRNA names (Column 2).\cr
#'  \item \strong{AccessionID} : The corresponding miRBase accession IDs (Column 3).
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

  ver_index=match(tolower(targetVersion),VER[,2])
  if (is.na(ver_index))
    stop("It is a wrong target version, Please check it")
  VMAP <- miRNA_data[[ver_index]][, c("ACC", "SYM")]

  uid = unique(as.vector(miRNANames))
  SYM_ID = match(uid, SYM$SYM)
  df <- data.frame(uid = uid, SYM = SYM_ID, stringsAsFactors=FALSE)
  df <- merge(df, ACC_SYM)[, c("uid", "ACC")]
  df <- unique( merge(df, VMAP, by="ACC") )

  target <- data.frame(
    OriginalName = df$uid,
    TargetName = SYM$SYM[df$SYM],
    AccessionID = ACC$ACC[df$ACC],
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
    AccessionID = splitpaste(target$AccessionID, f),
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

#' miRBase Accession ID to miRNA Name in specified version
#'
#' This function converts a group of any species' miRNA Accession IDs (including precursor and mature miRNA) to a specified miRBase version if the Accession IDs have been defined in miRBase.
#'
#' @param AccessionIDs A character vector representing the Accession IDs needed to be convert.
#' @param targetVersion A character value representing the target miRBase version corresponding the Accession IDs.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of rows equal to the input miRNA names. The two columns are defined as below:
#'\itemize{
#' \item \strong{AccessionID} : The Accession ID of miRNAs (Column 1).\cr
#' \item  \strong{TargetName} : The converted miRBase names (in specified version) corresponding to the Accession IDs (Column 2).\cr
#' }
#' @examples
#' data(miRNATest)
#' AccessionIDs=miRNATest$AccessionID
#' result1=miRNA_AccessionToName(AccessionIDs,targetVersion="v13")
#' result2=miRNA_AccessionToName(AccessionIDs,targetVersion="v21")
#'
#' @author
#'  Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
miRNA_AccessionToName<-function(AccessionIDs,targetVersion="v21")
{
  AccessionIDs=as.character(AccessionIDs)
  AccessionIDs=gsub(" ","",AccessionIDs)##Remove the possible space

  ver_index=match(tolower(targetVersion),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")
  VMAP <- miRNA_data[[ver_index]][, c("ACC", "SYM")]

  uid=unique(as.vector(AccessionIDs))
  ACC_ID=match(uid,ACC$ACC)
  df=data.frame(uid = uid, ACC = ACC_ID, stringsAsFactors=FALSE)
  df <- unique(merge(df, VMAP, by="ACC"))
  target <- data.frame(AccessionID = df$uid,TargetName = SYM$SYM[df$SYM],stringsAsFactors = FALSE)
  target <- target[match(AccessionIDs, target$AccessionID),]
  target$AccessionID=AccessionIDs
  rownames(target)= NULL
  target
}


#' The miRBase miRNA names with specified version to Accession IDs
#'
#' This function converts a group of any species' miRNA name to the Accession IDs defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param Version A character value representing the version corresponding the miRNANames.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of rows equal to the input miRNA names. The two columns are defined as below:
#'
#' \itemize{
#'  \item \strong{miRNAName_\{Version\}} : The input miRNA names (Column 1).\cr
#'  \item \strong{AccessionID} : The convert Accession ID (Column 2).\cr
#'  }
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' checkMiRNAVersion(miRNANames)
#' result1=miRNA_NameToAccession(miRNANames,Version="v18")
#' result2=miRNA_AccessionToName(result1[,2],targetVersion="v21")
#' result3=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#'
#' @author
#'  Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
miRNA_NameToAccession<-function(miRNANames,Version="v21")
{
  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space

  ver_index=match(tolower(Version),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrongn version, Please check it")
  VMAP <- miRNA_data[[ver_index]][, c("ACC", "SYM")]

  uid=unique(as.vector(miRNANames))
  SYM_ID = match(uid, SYM$SYM)
  df <- data.frame(uid = uid, SYM = SYM_ID, stringsAsFactors=FALSE)
  df <- unique( merge(df, VMAP, by="SYM"))
  target <- data.frame(miRNANames = df$uid,AccessionID = ACC$ACC[df$ACC],stringsAsFactors = FALSE)
  target <- target[match(miRNANames, target$miRNANames),]
  target$miRNANames=miRNANames
  rownames(target)= NULL
  colnames(target)[1]=paste("miRNAName_",Version,sep="")
  target
}

#' Obtain the miRNA sequences
#'
#' This function returns the miRNA sequences for a list of miRNAs.
#'
#' @param AccessionIDs A character vector representing the miRNA Accession IDs in miRBase.
#' @param targetVersion A character value representing the target miRBase version corresponding the Accession IDs.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data frame. The number of row equals to input miRNAs. The two columns are defined as below:
#'
#' \itemize{
#'  \item \strong{AccessionID} : The original miRNA (Column 1).\cr
#'  \item \strong{miRNASequence_\{targetVersion\}} : The return miRNA sequence (in specified version) corresponding to the input miRNAs (Column 2).\cr
#'  }
#' @examples
#' #####1, The input are miRNA Names
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result1=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#' AccessionIDs=result1[,3]
#' result2=getMiRNASequence(AccessionIDs,targetVersion="v21")
#'
#' #####2, The input are miRNA Accession IDs
#' data(miRNATest)
#' AccessionIDs=miRNATest$AccessionID
#' result3=getMiRNASequence(AccessionIDs,targetVersion="v13")
#' result4=getMiRNASequence(AccessionIDs,targetVersion="v21")
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getMiRNASequence<-function(AccessionIDs,targetVersion="v21")
{
  AccessionIDs=as.character(AccessionIDs)
  AccessionIDs=gsub(" ","",AccessionIDs)##Remove the possible space

  ver_index=match(tolower(targetVersion),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")
  VMAP <- miRNA_data[[ver_index]][, c("ACC", "SEQ")]

  uid=unique(as.vector(AccessionIDs))
  ACC_ID=match(uid,ACC$ACC)
  df=data.frame(uid = uid, ACC = ACC_ID, stringsAsFactors=FALSE)
  df <- unique(merge(df, VMAP, by="ACC"))
  target <- data.frame(AccessionID = df$uid,miRNASequence = SEQ$SEQ[df$SEQ],stringsAsFactors = FALSE)
  target <- target[match(AccessionIDs, target$AccessionID),]
  target$AccessionID=AccessionIDs
  rownames(target)= NULL
  colnames(target)[2]=paste("miRNASequence_",targetVersion,sep="")
  target
}


#' Obtain all miRBase versions' information for a single specified miRNA.
#'
#' This function returns all miRBase versions' information of a single specified miRNA.
#'
#' @param AccessionID A character representing the single AccessionID.
#'
#' @return
#'  A 21X7 data frame including all the history information (Precursor, Mature, Sequence) of the specified miRNA.
#'  Each row represents a miRBase version.
#' @examples
#' #####1,The input is a miRNA Name
#' name="hsa-miR-26b-5p"
#' result1=miRNAVersionConvert(name,targetVersion="v21",exact=TRUE)
#' AccessionID=result1[,3]
#' result2=getMiRNAHistory(AccessionID)
#'
#' #####2,The input is miRNA Accession Id
#' AccessionID="MIMAT0000765"
#' result3=getMiRNAHistory(AccessionID)
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
getMiRNAHistory<-function(AccessionID)
{
  AccessionID=as.character(AccessionID)
  AccessionID=gsub(" ","",AccessionID)##Remove the possible space

  if(length(AccessionID)>1)
    stop("This function is only for a single miRNA query, Please query miRNAs one by one")

  ACC_ID=match(AccessionID,ACC[,2])
  if(!is.na(ACC_ID))
  {
    target=data.frame(matrix(vector(),nrow(VER), 7,dimnames=list(c(), c("Version","Precursor","PrecursorSequence","Mature1","Mature1Sequence","Mature2","Mature2Sequence") )),stringsAsFactors=FALSE)
    target$Version<-VER[,2]
    for(i in 1:nrow(VER))
    {
      Precursor=NA
      seq0=NA
      Mature1=NA
      seq2=NA
      Mature2=NA
      seq3=NA

      ind1=match(ACC_ID,miRNA_data[[i]]$ACC)
      if(!is.na(ind1))
      {
        ind2=which(miRNA_data[[i]]$NUM==miRNA_data[[i]]$NUM[ind1])
        for(j in ind2)
        {
          SYM_ID=match(miRNA_data[[i]]$SYM[j],SYM[,1])
          SEQ_ID=match(miRNA_data[[i]]$SEQ[j],SEQ[,1])

          if(miRNA_data[[i]]$TYPE[j]==1)
          {
            Precursor=SYM[SYM_ID,2]
            seq0=SEQ[SEQ_ID,2]
          }

          if(miRNA_data[[i]]$TYPE[j]==2)
          {
            Mature1=SYM[SYM_ID,2]
            seq1=SEQ[SEQ_ID,2]
          }

          if(miRNA_data[[i]]$TYPE[j]==3)
          {
            Mature2=SYM[SYM_ID,2]
            seq2=SEQ[SEQ_ID,2]
          }
        }
      }
      target$Precursor[i]=Precursor
      target$PrecursorSequence[i]=seq0
      target$Mature1[i]=Mature1
      target$Mature1Sequence[i]=seq1
      target$Mature2[i]=Mature2
      target$Mature2Sequence[i]=seq2
    }
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
#' This function checks miRBase Version for a list of miRNA names
#' @importFrom stats na.omit
#' @param miRNANames A character vector representing the miRNA names.
#' @return
#'  The version information is printed in the console.
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' checkMiRNAVersion(miRNANames)
#'
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
#'
checkMiRNAVersion<-function(miRNANames)
{
  miRNANames=as.character(miRNANames)
  miRNANames=gsub(" ","",miRNANames)##Remove the possible space

  uid=unique(as.vector(miRNANames))
  SYM_ID=match(uid,SYM[,2])
  SYM_ID=na.omit(SYM_ID)
  result=data.frame(matrix(vector(),nrow(VER), 3,dimnames=list(c(), c("Version","Proportion","Recommend") )),stringsAsFactors=FALSE)
  result$Version<-VER[,2]

  for(i in 1:nrow(VER))
  {
    num =length(intersect(miRNA_data[[i]]$SYM, SYM_ID))
    result$Proportion[i]=round((num/length(uid))*100,2)
  }
  result$Recommend=""
  result[result$Proportion == max(result$Proportion), "Recommend"] <- " ***BEST Matched***"
  result$Proportion=paste(result$Proportion,"%",sep="")
  print.data.frame(result)
}

#' Open the miRBase webpage of the specified miRNA
#'
#' This function locates the miRBase webpage of the specified miRNA
#' @param AccessionID A character representing a single miRNA Accession ID.
#'
#' @examples
#' #### 1. A step-loop
#' AccessionID1="MI0000447"
#' goToMiRBase(AccessionID1)
#'
#' #### 2. A mature miRNA
#' AccessionID2="MIMAT0026477"
#' goToMiRBase(AccessionID2)
#' @return
#' No values
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
goToMiRBase<-function(AccessionID)
{
  AccessionID=as.character(AccessionID)
  AccessionID=gsub(" ","",AccessionID)##Remove the possible space
  if(is.na(match(AccessionID,ACC$ACC)))
  {
    stop("The input Accession ID is wrong. Please check")
  }
  else
  {
    if(grepl("MAT",AccessionID))
      URL=paste0("http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=",AccessionID)
    else
      URL=paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",AccessionID)
  }
  utils::browseURL(URL)
}
