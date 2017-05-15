#' miRBase name version convert
#'
#' This function converts a group of any species' miRNA names (including precursor and mature miRNA) to a specified miRBase version if the miRNAs have been defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param targetVersion A character value representing the target miRBase version name corresponding the source miRNA name.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @param exact logical value. If true, the result will be the most exactly matched result.
#'  If FALSE, the result will include all the possible matched miRNA name. If one miRNA can match multiple names. All the names will be
#'   concatenated with "&".
#' @note Please note: Due to some miRNAs' name changing many times in history. Even if choose the third parameter "exact"=TRUE, maybe it may still have some miRNAs that cannot
#' match the unique name in target version. In order to return the accurate result as possible, We also concatenate the multiple matched miRNA name with "&". This is rare case but it happens sometimes.
#'
#' @return
#' A nx3 data matrix. The number of row equals to input miRNA name. The three columns are defined as below:
#'\itemize{
#'  \item  \strong{OrignalName} : The original miRNA name (Column 1).\cr
#'  \item \strong{TargetName} : The converted miRBase name (in specified version) corresponding to the original miRNA name (Column 2).\cr
#'  \item \strong{AccessionID} : The corresponding miRBase accession ID (Column 3).
#'}
#' @examples
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result1=miRNAVersionConvert(miRNANames,targetVersion="v13",exact=TRUE)
#' result2=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#' result3=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=FALSE)
#' @author
#' Xu, Taosheng \email{taosheng.x@@gmail.com}
#' @export
miRNAVersionConvert<-function(miRNANames,targetVersion="v21",exact=TRUE)
{
  #version=paste("miRNA_",tolower(targetVersion),sep="")
  ver_index=match(tolower(targetVersion),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")

  temp=unique(as.vector(miRNANames))
  SYM_ID=match(temp,SYM[,2])

  temp_target=data.frame("orignal"=temp,"SYM"=NA,"ACC"=NA,stringsAsFactors=FALSE)
  for(i in 1:length(SYM_ID))
  {
    if(!is.na(SYM_ID[i]))
    {
      ind1=which(SYM_ID[i]==ACC_SYM[,2])
      ACC_ID=ACC_SYM[ind1,1]
      ind2=match(ACC_ID,miRNA_data[[ver_index]]$ACC)
      ACC_ID=ACC_ID[!is.na(ind2)]
      ind2=ind2[!is.na(ind2)]

      ind3= miRNA_data[[ver_index]]$SYM[ind2]
      ind3=unique(ind3)
      ###############################
      if(exact)
      {
        ind4=match(temp[i],SYM[ind3,2])
        if(!is.na(ind4))
        {
          mm=SYM[ind3,2][ind4]
          nn=ACC[ACC_ID,2][ind4]
        }
        else
        {
          mm=SYM[ind3,2]
          nn=ACC[ACC_ID,2]
        }
      }
      else
      {
        mm=SYM[ind3,2]
        nn=ACC[ACC_ID,2]
      }
      ###################################
      if(length(mm)>0)
        temp_target$SYM[i]=paste(mm, collapse = '&')
      if(length(nn)>0)
        temp_target$ACC[i]=paste(nn, collapse = '&')
    }
  }
  ####################
  multiindex=grep('&',temp_target[,2])
  if(length(multiindex)>0)
  {
    aa=temp_target[multiindex,]
    colnames(aa)<-c("orignal",paste("Version",targetVersion),"accession")
    rownames(aa)<-NULL
    message("********************************************\n")
    message("The multiple matched miRNAs are list below: \n")
    message("                                            \n")
    print.data.frame(aa)
  }
  #####################

  target= cbind("OrignalName"=miRNANames,"TargetName"=NA,"AccessionID"=NA)
  for(i in 1:length(temp))
  {
    ind5=which(miRNANames==temp[i])
    #if(!is.na(ind5))
    target[ind5,2]=temp_target[i,2]
    target[ind5,3]=temp_target[i,3]
  }
  colnames(target)[2]=paste("Version",targetVersion)
  target
}

#' miRBase Accession ID to miRNA Name in specified version
#'
#' This function converts a group of any species' miRNA Accession ID (including precursor and mature miRNA) to a specified miRBase version name if the Accession ID have been defined in miRBase.
#'
#' @param AccessionIDs A character vector representing the Accession ID needed to be convert.
#' @param targetVersion A character value representing the target miRBase version name corresponding the Accession ID.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr

#' @return
#' A nx2 data matrix. The number of row equals to input miRNA name. The two columns are defined as below:
#'\itemize{
#' \item \strong{AccessionID} : The Accession ID of miRNAs (Column 1).\cr
#' \item  \strong{TargetName} : The converted miRBase name (in specified version) corresponding to the Accession ID (Column 2).\cr
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
  ver_index=match(tolower(targetVersion),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")
  #AccessionIDs=miRNATest$AccessionID
  temp=unique(as.vector(AccessionIDs))

  temp_target=data.frame("ACC"=temp,"SYM"=NA,stringsAsFactors=FALSE)

  ACC_ID=match(temp,ACC[,2])
  for(i in 1:length(ACC_ID))
  {
    if(!is.na(ACC_ID[i]))
    {
      ind1=match(ACC_ID[i],miRNA_data[[ver_index]]$ACC)
      if(!is.na(ind1))
      {
        SYM_ID=miRNA_data[[ver_index]]$SYM[ind1]
        temp_target[i,2]=SYM[SYM_ID,2]
      }
    }

  }
  target= cbind("AccessionID"=AccessionIDs,"TargetName"=NA)
  for(i in 1:length(temp))
  {
    ind2=which(AccessionIDs==temp[i])
    target[ind2,2]=temp_target[i,2]
  }
  colnames(target)[2]=paste("Version",targetVersion)
  target
}


#' The miRBase miRNA Name with specified version to Accession ID
#'
#' This function converts a group of any species' miRNA name to the Accession ID defined in miRBase.
#'
#' @param miRNANames A character vector representing the source miRNA names needed to be convert.
#' @param Version A character value representing the version name corresponding the miRNANames.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr
#' @return
#' A nx2 data matrix. The number of row equals to input miRNA name. The two columns are defined as below:
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
  ver_index=match(tolower(Version),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")
  #miRNANames=miRNATest$miRNA_Name
  temp=unique(as.vector(miRNANames))

  SYM_ID=match(temp,SYM[,2])
  ind1=match(SYM_ID,miRNA_data[[ver_index]]$SYM)
  ind2=miRNA_data[[ver_index]]$ACC[ind1]
  ACC_ID=ACC[ind2,2]
  temp_target=data.frame("orignal"=temp,"ACC"=ACC_ID,stringsAsFactors=FALSE)


  target= cbind("miRNAName"=miRNANames,"AccessionID"=NA)
  for(i in 1:length(temp))
  {
    ind3=which(miRNANames==temp[i])
    target[ind3,2]=temp_target[i,2]
  }
  colnames(target)[1]=paste("miRNAName_",Version,sep="")
  target
}

#' Obtain the miRNA sequence
#'
#' This function return the miRNA sequence based on the miRNA name or Accession ID.
#'
#' @param AccessionIDs A character vector representing the miRNA Accession ID in miRBase.
#' @param targetVersion A character value representing the target miRBase version name corresponding the Accession ID.
#' The optional values are in below:\cr
#' "v6","v7_1","v8","v8_1","v8_2","v9","v9_1",\cr
#' "v9_2","v10","v10_1","v11","v12","v13","v14",\cr
#' "v15","v16","v17","v18","v19","v20","v21"\cr

#' @return
#' A nx2 data matrix. The number of row equals to input miRNA. The two columns are defined as below:
#'
#' \itemize{
#'  \item \strong{AccessionID} : The original miRNA (Column 1).\cr
#'  \item \strong{miRNASequence_\{targetVersion\}} : The return miRNA sequence (in specified version) corresponding to the input miRNAs (Column 2).\cr
#'  }
#' @examples
#' #####1,The input aer miRNA Names
#' data(miRNATest)
#' miRNANames=miRNATest$miRNA_Name
#' result1=miRNAVersionConvert(miRNANames,targetVersion="v21",exact=TRUE)
#' AccessionIDs=result1[,3]
#' result2=getMiRNASequence(AccessionIDs,targetVersion="v21")
#'

#' #####2,The input are miRNA Accession IDs
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
  ver_index=match(tolower(targetVersion),VER[,2])
  if(is.na(ver_index))
    stop("It is a wrong target version, Please check it")

  temp=unique(as.vector(AccessionIDs))
  temp_target=data.frame("ACC"=temp,"SEQ"=NA,stringsAsFactors=FALSE)

  ACC_ID=match(temp,ACC[,2])

  for(i in 1:length(ACC_ID))
  {
    if(!is.na(ACC_ID[i]))
    {
      ind1=match(ACC_ID[i],miRNA_data[[ver_index]]$ACC)
      if(!is.na(ind1))
      {
        sEQ_ID=miRNA_data[[ver_index]]$SEQ[ind1]
        temp_target[i,2]=SEQ[sEQ_ID,2]
      }
    }
  }
  target= cbind("AccessionID"=AccessionIDs,"miRNASequence"=NA)
  for(i in 1:length(temp))
  {
    ind2=which(AccessionIDs==temp[i])
    target[ind2,2]=temp_target[i,2]
  }
  colnames(target)[2]=paste("miRNASequence_",targetVersion,sep="")
  target
}


#' Obtain all versions information for a single specified miRNA.
#'
#' This function return all versions information of a single specified miRNA.
#'
#' @param AccessionID A character representing the single AccessionID.

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

  if(length(AccessionID)>1)
    stop("This function is only for a single miRNA query, Please query miRNAs one by one")

  ACC_ID=match(AccessionID,ACC[,2])
  if(!is.na(ACC_ID))
  {
    target=data.frame(matrix(vector(),nrow(VER), 7,dimnames=list(c(), c("Version","Precursor","PrecursorSequence","Mautur1","Mautur1Sequence","Mautur2","Mautur2Sequence") )),stringsAsFactors=FALSE)
    target$Version<-VER[,2]
    for(i in 1:nrow(VER))
    {
      Precursor=NA
      seq0=NA
      Mautur1=NA
      seq2=NA
      Mautur2=NA
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
            Mautur1=SYM[SYM_ID,2]
            seq1=SEQ[SEQ_ID,2]
          }

          if(miRNA_data[[i]]$TYPE[j]==3)
          {
            Mautur2=SYM[SYM_ID,2]
            seq2=SEQ[SEQ_ID,2]
          }
        }
      }
      target$Precursor[i]=Precursor
      target$PrecursorSequence[i]=seq0
      target$Mautur1[i]=Mautur1
      target$Mautur1Sequence[i]=seq1
      target$Mautur2[i]=Mautur2
      target$Mautur2Sequence[i]=seq2
    }
    target
  }
  else
  {
    message("There is no records about this miRNA.")
    NULL
  }
}


#' check miRNA Version in miRBase
#'
#' This function check miRBase Version for a list of miRNA names
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
  temp=unique(as.vector(miRNANames))
  SYM_ID=match(temp,SYM[,2])
  SYM_ID=na.omit(SYM_ID)
  result=data.frame(matrix(vector(),nrow(VER), 3,dimnames=list(c(), c("Version","Proportion","Recommend") )),stringsAsFactors=FALSE)
  result$Version<-VER[,2]

  for(i in 1:nrow(VER))
  {
    num =length(intersect(miRNA_data[[i]]$SYM, SYM_ID))
    result$Proportion[i]=round((num/length(temp))*100,2)
  }
  result$Recommend=""
  result[result$Proportion == max(result$Proportion), "Recommend"] <- "*****"
  result$Proportion=paste(result$Proportion,"%",sep="")
  print.data.frame(result)
}




