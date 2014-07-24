
TopList <- function(object,topnum=10,sorted.by="RValue")  {
    ####  Example usage:
    ####  toplist(rvs,topnum=15,sorted.by="MLE")

    if(sorted.by=="RValue")  {
        ans <- object$main[1:topnum,]
    }
    else if(sorted.by=="PostMean") {
        object$main <- object$main[order(object$main$PM.rank),]
        ans <- object$main[1:topnum,]
    }
    else if(sorted.by=="MLE") {
        object$main <- object$main[order(object$main$MLE.rank),]
        ans <- object$main[1:topnum,]
    }
    else if(sorted.by=="PVal") {
        object$main <- object$main[order(object$main$PVal.rank),]
        ans <- object$main[1:topnum,]
    }
    return(ans)
}
