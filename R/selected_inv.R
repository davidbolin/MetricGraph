#' Selected Inverse Calculation
#' 
#' @param Q A sparse matrix in dgCMatrix format
#' @return A numeric vector containing the selected inverse
#' @export
selected_inv <- function(Q) {
    if (!inherits(Q, "sparseMatrix")) {
        stop("Q must be a sparse matrix (e.g., dgCMatrix, dsCMatrix, or dgTMatrix)")
    }
    if (!inherits(Q, "dgCMatrix")) {
        Q <- as(Q, "dgCMatrix")  
    }
    
    result <- selected_inv_cpp(Q)
    return(result)
}