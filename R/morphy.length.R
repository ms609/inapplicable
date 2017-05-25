#' @title Count a tree length
#'
#' @description Count a tree length using the Morphy Lib implementation
#'
#' @param tree a \code{phylo} object.
#' @param matrix a discrete morphological matrix
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

morphy.length <- function(tree, matrix) {
    # tree <- read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));") ; stop("DEBUG")
    # ## Data
    # matrix <- matrix(data = c("2", "3", "-", "-", "1", "?", "?", "-", "-", "0", "3", "2"), ncol = 1) ; stop("DEBUG")

    ## Check matrix and replace get dimensions

    ## morphy variables
    rawmatrix <- paste(apply(matrix, 2, paste, collapse = ""), ";", sep = "")
    ntax <- ape::Ntip(tree)
    nnode <- ape::Nnode(tree)
    nchar <- ncol(matrix)
    
    morphy_object <- mpl_new_Morphy()
    catch_error <- mpl_init_Morphy(ntax, nchar, morphy_object)
    catch_error <- mpl_attach_rawdata(rawmatrix, morphy_object)
    catch_error <- mpl_set_num_internal_nodes(ntax+1, morphy_object)
    catch_error <- mpl_apply_tipdata(morphy_object)
    

    tips_ancestors <- c(12, 12, 13, 14, 15, 16, 21, 20, 19, 18, 17, 17) ; #tips_ancestors <- tips_ancestors +1
    ancestors      <- c(13, 14, 15, 16, 22, 18, 19, 20, 21, 22, 23)     ; ancestors <- ancestors +1
    nodes          <- c(12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23) ; #nodes <- nodes +1
 left_descendants  <- c(0,  12, 13, 14, 15, 10,  9,  8,  7,  6, 16)     ; #left_descendants <- left_descendants +1
right_descendants  <- c(1,   2,  3,  4,  5, 11, 17, 18, 19, 20, 21)     ; #right_descendants <- right_descendants +1
    

    ## Initialise length
    length <- 0

    ## First up pass
    for (i in 1:(ntax-1)) {
        length <- length + mpl_first_down_recon(nodes[i], left_descendants[i], right_descendants[i], morphy_object)
    }
    
    ## Update root
    catch_error <- mpl_update_lower_root(23, 22, morphy_object)
    
    ## First up pass
    for (i in (ntax-1):1) {
        length <- length + mpl_first_up_recon(nodes[i], left_descendants[i], right_descendants[i], ancestors[i], morphy_object)
    }
    
    ## Update tips
    for (i in 1:ntax) {
        catch_error <- mpl_update_tip(i, tips_ancestors[i], morphy_object)
    }

    ## Second down pass
    for (i in 1:(ntax-1)) {
        length <- length + mpl_second_down_recon(nodes[i], left_descendants[i], right_descendants[i], morphy_object)
    }
    paste("length after 2nd down:", length)
    
    ## Second up pass
    for (i in (ntax-1):1) {
        length <- length + mpl_second_up_recon(nodes[i], left_descendants[i], right_descendants[i], ancestors[i], morphy_object)
    }
    paste("length after 2nd up:", length)
    
    mpl_delete_Morphy(morphy_object)


    # ## First downpass
    # sapply_first_down <- function(i, nodes, left_descendants, right_descendants, morphy_object) {
    #     return(mpl_first_down_recon((nodes[i]), (left_descendants[i]), (right_descendants[i]), morphy_object))
    # }
    # increment <- sum(sapply(1:nnode, sapply_first_down, nodes, left_descendants, right_descendants, morphy_object))
    # length <- length + increment #TG: note that this is pretty useless (no length increment possible here). Just left it for consistency.

    # ## Update the root
    # catch_error <- mpl_update_lower_root(ntax+nnode, ntax, morphy_object) # mpl_update_lower_root(23, 22, m1);

    # ## First uppass
    # sapply_first_up <- function(i, nodes, left_descendants, right_descendants, ancestors, morphy_object) {
    #     return(mpl_first_up_recon((nodes[i]), (left_descendants[i]), (right_descendants[i]), (ancestors[i]), morphy_object))
    # }
    # increment <- sum(sapply(nnode:1, sapply_first_up, nodes, left_descendants, right_descendants, ancestors, morphy_object))
    # length <- length + increment #TG: note that this is pretty useless (no length increment possible here). Just left it for consistency.

    # ## Update the tips
    # sapply_update_tip <- function(i, tips_ancestors, morphy_object) {
    #     return(mpl_update_tip(i, (tips_ancestors[i]), morphy_object))
    # }
    # catch_error <- sapply(1:ntax, sapply_update_tip, tips_ancestors, morphy_object)

    # ## Second downpass
    # sapply_second_down <- function(i, nodes, left_descendants, right_descendants, morphy_object) {
    #     return(mpl_second_down_recon((nodes[i]), (left_descendants[i]), (right_descendants[i]), morphy_object))
    # }
    # increment <- sum(sapply(1:nnode, sapply_second_down, nodes, left_descendants, right_descendants, morphy_object))
    # print(increment)
    # length <- length + increment

    # ## Second uppass
    # sapply_second_up <- function(i, nodes, left_descendants, right_descendants, ancestors, morphy_object) { 
    #     return(mpl_second_up_recon((nodes[i]), (left_descendants[i]), (right_descendants[i]), (ancestors[i]), morphy_object)) #TG: must call the dummy root here
    # }
    # increment <- sum(sapply(nnode:1, sapply_second_up, nodes, left_descendants, right_descendants, ancestors, morphy_object))
    # print(increment)
    # length <- length + increment
 
    return(length)
}