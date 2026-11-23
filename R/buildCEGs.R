#' buildCEGs: A Tool for Constructing, Manipulating and Analysing Chain Event Graph Models
#'
#' @name buildCEGs
#' @importFrom magrittr %>%
NULL

#' Save a Plot as a High-Resolution PNG
#'
#' Saves a plot of a tree or Chain Event Graph as a PNG file
#' in the working directory.
#'
#' Allows the user to specify the dimensions for the PNG file.
#' The file is saved directly in the working directory.
#'
#' @param plot A tree or CEG plot
#' @param filename A filename for the PNG file without the .png extension
#' @param width Width of the output in pixels. Defaults to 4000.
#' @param height Height of the output in pixels. Defaults to 8000.
#' @return A .png file in the working directory
#' @export
PlotToPNG <- function(plot,filename,width=4000,height=8000){
  svg <- DiagrammeRsvg::export_svg(plot)
  rsvg::rsvg_png(charToRaw(svg),
           paste0(filename,".png"),
           width = width,
           height = height)
}


#' Generate a Table Specifying an Event Tree
#'
#' Generates a table documenting each node in the event tree and information about these nodes to uniquely specify an event tree.
#'
#' The table specifies the ID, the label, the parent, the incoming edge label, the level and the incoming edge probability for each node.
#' It also features a placeholder for future stage and position information for creating staged trees and CEGs.
#' Where a node has no parent, NA values are used to indicate this for the parent, edge and edge probability columns.
#' By default, this function generates a full event tree including every configuration of states across the inputted variables.
#'
#' @param variables A list of named variables with their possible states
#'
#' @return A table uniquely describing the event tree
#' @export
#'
#' @examples
#' variables <- list(
#'      Bite       =   c("No","Yes"),
#'      Antidote   =   c("No","Yes"),
#'      Die        =   c("No","Yes")
#' )
#'
#' EventTreeTable(variables)
#'
EventTreeTable <- function(variables){
  nodes <- data.frame(
    id = 0,
    label = "s0",
    parent = NA,
    edge = NA,
    level = 0,
    stringsAsFactors = FALSE
  )
  nextID <- 1
  recurse <- function(parentID, level){
    if(level > length(variables)) return()
    for(state in variables[[level]]){
      thisID <- nextID
      nodes <<- rbind(
        nodes,
        data.frame(
          id = thisID,
          label = paste0("s", thisID),
          parent = parentID,
          edge = state,
          level = level,
          stringsAsFactors = FALSE
        )
      )
      nextID <<- nextID + 1
      recurse(thisID, level + 1)
    }
  }
  recurse(0,1)
  # Breadth-first ordering
  order_index <- order(nodes$level, nodes$id)

  # Old IDs in breadth-first order
  old_ids <- nodes$id[order_index]

  # New IDs
  new_ids <- seq_along(old_ids) - 1

  # Mapping: old ID -> new ID
  id_map <- stats::setNames(new_ids, old_ids)

  # Renumber ids
  nodes$id <- id_map[as.character(nodes$id)]

  # Renumber parents
  new_parent <- id_map[as.character(nodes$parent)]
  new_parent[is.na(nodes$parent)] <- NA
  nodes$parent <- as.integer(new_parent)

  # Update labels
  nodes$label <- paste0("s", nodes$id)

  # Sort by the new ids
  nodes <- nodes[order(nodes$id), ]
  row.names(nodes) <- NULL
  nodes$edgeprob <- NA
  nodes$stage <- NA
  nodes$position <- NA
  nodes$edgeprob <- as.numeric(nodes$edgeprob)
  nodes$id <- as.integer(nodes$id)
  nodes$level <- as.integer(nodes$level)
  nodes$stage <- as.integer(nodes$stage)
  nodes$position <- as.integer(nodes$position)
  nodes
}


#' Get the Descendants of a Set of Nodes
#'
#' Returns a vector of descendants of an inputted node or set of nodes, including the nodes inputted.
#'
#' The descendants of a node are the nodes that can be reached from the specified node following the edges of the tree.
#' The function finds all such nodes for a set of inputted nodes.
#'
#' @param infoTable A table uniquely specifying an event tree
#' @param nodes A node ID or vector of node IDs
#'
#' @return A vector of node IDs that are descendants of at least one inputted node ID
#' @export
#'
#' @seealso [ancestors()] to return the set of ancestors
#' @examples
#' variables <- list(
#'      Bite       =   c("No","Yes"),
#'      Antidote   =   c("No","Yes"),
#'      Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' # to get the descendants of a node
#'
#' descendants(tree,2)
#' descendants(tree,c(2))
#'
#' # to get descendants across a vector of nodes
#'
#' descendants(tree,c(2,3))
#'
#'
descendants <- function(infoTable, nodes){
  if(all(nodes %in% infoTable$id)==FALSE){
    stop("There are node IDs that don't exist!")
  }
  children <- infoTable$id[infoTable$parent %in% nodes]
  if(length(children) == 0){
    return(nodes)
  }
  desc <- nodes
  for(child in children){
    desc <- unique(c(desc, descendants(infoTable, child)))
  }
  sort(desc)
}


#' Get the Ancestors of a Set of Nodes
#'
#' Returns a vector of ancestors of an inputted node or set of nodes, including the nodes inputted.
#'
#' The ancestors of a node are the nodes from which the specified node can be reached following the edges of the tree.
#' The function finds all such nodes for a set of inputted nodes.
#'
#' @inheritParams descendants
#'
#' @return A vector of node IDs that are ancestors of at least one inputted node ID
#'
#' @export
#'
#' @seealso [descendants()] to return the set of descendants
#' @examples
#' variables <- list(
#'      Bite       =   c("No","Yes"),
#'      Antidote   =   c("No","Yes"),
#'      Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' # to get the ancestors of a node
#'
#' ancestors(tree,11)
#' ancestors(tree,c(11))
#'
#' # to get descendants across a vector of nodes
#'
#' ancestors(tree,c(1,11))
#'
ancestors <- function(infoTable, nodes){
  if(all(nodes %in% infoTable$id)==FALSE){
    stop("There are node IDs that don't exist!")
  }
  anc <- nodes
  current <- nodes
  while (length(current) > 0) {
    rows <- match(current, infoTable$id)
    parents <- unique(infoTable$parent[rows])
    parents <- parents[!is.na(parents)]
    current <- setdiff(parents, anc)
    anc <- unique(c(anc, current))
  }
  sort(anc)
}


#' Add or Update Edge Probabilities
#'
#' Input a vector of incoming edge probabilities for a vector of child nodes.
#'
#' The vector of edge probabilities must correspond to the vector of node IDs in the order they are specified.
#' A separate function is required to check validitity of edge probabilities - see `checkEdgeProbs`.
#'
#' @param infoTable A table uniquely specifying an event tree
#' @param nodes A node ID or vector of node IDs
#' @param edgeprobs A vector of edge probabilities
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table specifying the same event tree with added or updated edge probabilities
#'
#' @export
#'
#' @seealso
#' \itemize{
#'     \item [addUniformProbs()] to add uniform probability vectors across the edges in specified florets
#'     \item [checkEdgeProbs()] to check the validity of inputted edge probabilities
#' }
#' @examples
#' variables <- list(
#'      Bite       =   c("No","Yes"),
#'      Antidote   =   c("No","Yes"),
#'      Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' # suppose the floret rooted at (i.e. the edges coming out of) the root node has probability
#' # vector (0.75,0.25), and the floret at nodes 1 and 2 both have probability vector (0.8,0.2)
#'
#' nodes <- tree$id[which(tree$parent %in% c(0,1,2))]
#'
#' tree <- addEdgeProbs(tree,nodes,c(0.75,0.25,0.8,0.2,0.8,0.2))
#'
#' # if we know the child node IDs directly, we can input them directly
#'
#' tree <- addEdgeProbs(tree,7:14,c(0.95,0.05,0.95,0.05,0.9,0.1,0.95,0.05))
#'
#' # entries can be updated and corrected
#'
#' tree <- addEdgeProbs(tree,c(11,12),c(0.1,0.9))
addEdgeProbs <- function(infoTable,nodes,edgeprobs,deleteStages = TRUE){
  infoTable$edgeprob[match(nodes,infoTable$id)] <- edgeprobs
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  infoTable
}


#' Plot an Event Tree
#'
#' Plots the event tree specified by the tree information table outputted by `EventTreeTable` (and other functions).
#'
#' Produce a plot of the event tree which can then be saved locally using `PlotToPNG`, as well as dot representation which may be used by other graphical manipulation functions.
#' Features several optional graphical visualisation arguments including variable labels and path probability labels.
#' The leafprobs argument can simply be a vector of path probabilities (ensuring the correct ordering), or a table matching each path probability to each leaf node, as per the output of `getPathProbs`.
#'
#' @param infoTable A table uniquely specifying an event tree
#' @param varnames Optional vector of variable names corresponding to the levels of the tree
#' @param leafprobs Optional vector of path probabilities associated with the set of leaf nodes
#' @param padding Optional parameter to modify boundary space around the tree
#' @param leafnodesep Optional parameter to modify vertical separation between nodes
#' @param varsep Optional parameter to modify horizontal separation between nodes
#' @param nodewidth Optional parameter specifying node width
#' @param nodefont Optional parameter specifying node font size
#' @param edgefont Optional parameter specifying font size of edge labels
#' @param probfont Optional parameter specifying font size of labelled path probabilities
#' @param varfont Optional parameter specifying font size of variable name labels
#'
#' @return A plot of the inputted event tree alongside its dot representation
#'
#' @export
#'
#'
#' @examples
#' variables <- list(
#'      Bite       =   c("No","Yes"),
#'      Antidote   =   c("No","Yes"),
#'      Die        =   c("No","Yes")
#' )
#'
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#' tree <- addEdgeProbs(tree,1:14,c(0.75,0.25,0.8,0.2,0.8,0.2,0.95,0.05,
#'                                  0.95,0.05,0.1,0.9,0.95,0.05))
#' pathprobs <- getPathProbs(tree)
#'
#' plotEventTree(tree,varnames,pathprobs,padding=1)
plotEventTree <- function(
    infoTable,
    varnames = NULL,
    leafprobs = NULL,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  dot <- sprintf("
    digraph EventTree {
    graph [
      rankdir=LR,
      nodesep=%s,
      ranksep=%s,
      pad=%s
    ]
    node [
    shape=circle,
    fixedsize=true,
    width=%s,
    fontsize=%s
    ]
    edge [fontsize=%s]
    ", leafnodesep, varsep, padding,nodewidth,nodefont,edgefont)
  if(is.na(infoTable[1,]$parent)) infoTable[1,]$label<-"s0"
  for(i in 2:nrow(infoTable)) {

    parent <- infoTable$parent[i]
    parentLabel <- infoTable$label[infoTable$id == parent]
    childLabel <- infoTable$label[i]

    parentlevel <- infoTable$level[infoTable$id == parent]
    childlevel <- infoTable$level[i]

    leveldiff <- childlevel - parentlevel

    dot <- paste0(
      dot,
      sprintf(
        '"%s" -> "%s" [label="%s (%s)", minlen=%d];\n',
        parentLabel,
        childLabel,
        infoTable$edge[i],
        as.character(round(infoTable$edgeprob[i],2)),
        leveldiff
      )
    )
  }

  for (lev in sort(unique(infoTable$level))) {

    rows <- which(infoTable$level == lev)

    dot <- paste0(dot, "{rank=same; ")

    for (r in rows) {
      dot <- paste0(dot, sprintf('"%s"; ', infoTable$label[r]))
    }

    dot <- paste0(dot, "}\n")
  }

  if (!is.null(varnames)) {

    for (lev in seq_along(varnames)) {

      rows <- which(infoTable$level == lev - 1)

      pad <- ceiling(nchar(varnames[lev]) / 3)

      nbsp <- paste(rep("&nbsp;", pad), collapse = "")

      dot <- paste0(
        dot,
        sprintf(
          '"var%d" [label=<%s%s<BR/><BR/><BR/>>, shape=plaintext, fontsize=%s, fontname="Computer Modern Roman"];\n',
          lev,
          varnames[lev],
          nbsp,
          varfont
        )
      )

      dot <- paste0(dot, "{rank=same; ")
      dot <- paste0(dot, sprintf('"var%d"; ', lev))

      for (r in rows) {
        dot <- paste0(dot, sprintf('"%s"; ', infoTable$label[r]))
      }

      dot <- paste0(dot, "}\n")
    }
  }

  if (!is.null(leafprobs)) {

    if(!is.null(dim(leafprobs))) leafprobs <- leafprobs[,2]

    leaves <- infoTable$id[!(infoTable$id %in% infoTable$parent)]

    for(i in seq_along(leaves)) {

      leaf_id <- leaves[i]
      leaf_label <- infoTable$label[infoTable$id == leaf_id]

      prob <- leafprobs[i]
      prob_node <- paste0(leaf_label, "_prob")

      dot <- paste0(
        dot,
        sprintf('"%s" [label="%.4f", shape=plaintext, fontsize=%s];\n',
                prob_node, prob,probfont)
      )

      dot <- paste0(
        dot,
        sprintf('"%s" -> "%s" [style=invis, arrowhead=none, weight=100];\n',
                leaf_label, prob_node)
      )
    }
  }
  dot <- paste0(dot, "\n}")
  list(plot=DiagrammeR::grViz(dot),dot=dot)
}

#' Input Uniform Edge Probabilities for Specified Parents
#'
#' Adds a uniform probability vector across the outgoing edges of each inputted node
#'
#' This can also be described as inputting a uniform edge probability vector for each floret based at one of the inputted nodes.
#' Unlike many functions in this package, the inputted nodes here refer to parent nodes and their outgoing edges rather than child nodes and their incoming edges.
#'
#' @inheritParams addEdgeProbs
#'
#' @return A table specifying the same event tree updated with uniform probability vectors across outgoing edges for specified nodes.
#' @export
#'
#' @seealso
#' \itemize{
#'    \item [addEdgeProbs()] for specifying non-uniform probability vectors
#'    \item [checkEdgeProbs()] to check the validity of inputted edge probabilities
#' }
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' # to add uniform probability vectors to only some florets:
#' tree <- addUniformProbs(tree,c(0,1,2))
#'
#' # to add uniform probability vectors to all florets, all node IDs can be inputted,
#' # though a warning will appear that leaf nodes have been ignored:
#' tree <- addUniformProbs(tree,tree$id)
#'
#' # to get the same result without the warning:
#' tree <- addUniformProbs(tree,tree$id[which(tree$level<max(tree$level))])
#'
addUniformProbs <- function(infoTable,nodes,deleteStages=TRUE){
  if(max(infoTable$level) %in% infoTable$level[which(infoTable$id %in% nodes)]){
    warning("You have inputted leaf nodes which will be ignored. Inputted nodes should have children.")
  }
  if(!all(nodes %in% infoTable$id)){
    warning("There are node IDs that don't exist!")
  }
  nodes <- intersect(nodes,infoTable$id)
  for(i in nodes){
    childrows <- which(infoTable$parent==i)
    infoTable$edgeprob[childrows] <- 1/length(childrows)
  }
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  infoTable
}

#' Check the Validity of Edge Probabilities
#'
#' Checks whether the specified edge probabilities satisfy the sum-to-one condition across each floret in the tree.
#'
#' A floret refers to a node and the set of its outgoing edges. The probabilities associated with these edges must sum to one.
#' This function checks this condition, specifying valid and invalid parent nodes - i.e. the nodes whose florets are probabilistically valid/invalid.
#'
#' @inheritParams plotEventTree
#'
#' @return A vector of nodes whose outgoing edges are valid, and for those whose outgoing edges are invalid, alongside a message highlighting whether there are invalid probabilities or not.
#' @export
#'
#' @seealso
#' \itemize{
#'    \item [addEdgeProbs()] to add edge probabilities
#'    \item [addUniformProbs()] to add uniform probability vectors to specified florets
#' }
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' tree <- addEdgeProbs(tree,1:20,c(0.75,0.25,0.6,0.2,0.2,0.5,0.2,0.2,0.95,0.05,
#'                                  0.95,0.05,0.95,0.05,0.1,0.9,0.5,0.5,0.95,0.05))
#' checkEdgeProbs(tree)
#'
#' # This gives us a warning stating there are invalid probabilities,
#' # and tells us that the floret from node 2 is invalid. We can then fix the mistake.
#'
#' tree <- addEdgeProbs(tree,6,0.6)
#' checkEdgeProbs(tree)
checkEdgeProbs <- function(infoTable){
  parents <- unique(infoTable$parent[-1])
  valid <- numeric()
  invalid <- numeric()
  for(i in parents){
    if(abs(sum(infoTable$edgeprob[infoTable$parent == i], na.rm = TRUE) - 1) < 1e-8){
      valid <- c(valid,i)
    } else {
      invalid <- c(invalid,i)
    }
  }
  if(length(invalid) > 0) message("CAUTION - there are invalid edge probabilities. If working with a pruned tree, check that the invalid edge probabilities are not present.")
  if(length(invalid)==0) message("These edge probabilities are valid. Take care to nevertheless verify these probabilities.")
  list(valid=valid,invalid=invalid)
}


#' Relabel an Edge
#'
#' Gives a new label to the edge leading to the specified child node
#'
#' @inheritParams addEdgeProbs
#' @param child The node ID that the target edge leads to
#' @param newname The new label for the edge
#'
#' @return A table specifying the same event tree updated with the new edge label.
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#'
#' tree <- EventTreeTable(variables)
#'
#' tree <- renameEdge(tree,7,"Antidote 1")
#' tree <- renameEdge(tree,8,"Antidote 2")
#'
#'
renameEdge <- function(infoTable,child,newname,deleteStages=TRUE){   # this is currently only for children with only one incoming edge i.e. fine for a tree
  if(length(which(infoTable$id==child))>1) stop("This function can only edit edge labels of children with one incoming edge!")
  if(!(child %in% infoTable$id)) stop("This child doesn't exist!")
  if(child == 0) stop("The root node has no incoming edge!")
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  infoTable[infoTable$id==child,]$edge <- newname
  infoTable
}



#' Remove a Node by Connecting its Parent to its Children
#'
#' Bypasses a given node by deleting the node and its incident edges while adding an edge from the node's parent to each of its children.
#'
#' Used to represent context-dependent conditional independence statements in which a variable subsequent to a particular node is irrelevant for the modelling of later variables in the tree.
#' For example, if, given \eqn{X_0=0} and \eqn{X_1=1}, \eqn{X_3} is independent of \eqn{X_2}, we can bypass the node situated at the end of the partial path \eqn{\{X_0=0,X_1=1\}}.
#' Bypass multiple nodes simulataneously using the function `bypassNodes`.
#' @inheritParams plotEventTree
#' @param node A node ID
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#'
 bypassNode <- function(
    infoTable,
    node,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  if(!(node %in% infoTable$id)) stop("This node ID doesn't exist!")
  if(length(which(infoTable$id==node))>1) stop("This function can only edit edge labels of children with one incoming edge!")
  if(infoTable$level[which(infoTable$id==node)]==0) stop("You cannot bypass the root node!")

  row <- which(infoTable$id==node)
  bypassrow <- which(infoTable$parent==node)

  if(length(bypassrow)==0) stop("You cannot bypass a leaf node!")
  if(length(unique(infoTable$level[bypassrow])) > 1) stop("You cannot bypass a node whose children are in different levels because only uninformative nodes can be bypassed! A node whose non-zero-porbability outcomes lead to different levels must be informative. If this node is uninformative, first modify its subtree structure to reflect this.")

  excdesc <- descendants(infoTable,infoTable$id[bypassrow[-1]])
  infoTable <- infoTable[-which(infoTable$id %in% excdesc),]
  parent <- infoTable$parent[row]
  edgename <- infoTable$edge[row]
  edgeprob <- infoTable$edgeprob[row]
  bypassrow <- bypassrow[1]
  infoTable$edge[bypassrow] <- edgename
  infoTable$parent[bypassrow] <- parent
  infoTable$edgeprob[bypassrow] <- edgeprob
  infoTable[c(row,bypassrow),] <- infoTable[c(bypassrow,row),]
  infoTable <- infoTable[-bypassrow,]
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  plot <- plotEventTree(infoTable, varnames, leafprobs,
                        padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                        nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                        probfont = probfont, varfont = varfont)

  list(plot = plot$plot, dot = plot$dot, table = infoTable)
}


#' Remove Nodes by Connecting Each Parent to Each of its Children
#'
#' Bypasses a node or each of a set of nodes by deleting each node and its incident edges while adding an edge from each node's parent to each of its children.
#'
#' Used to represent context-dependent conditional independence statements in which a variable subsequent to a particular node is irrelevant for the modelling of later variables in the tree.
#' For example, if, given \eqn{X_0=0} and \eqn{X_1=1}, \eqn{X_3} is independent of \eqn{X_2}, we can bypass the node situated at the end of the partial path \eqn{\{X_0=0,X_1=1\}}.
#' Implements the hidden function `bypassNode` on a set of node IDs.
#' @inheritParams plotEventTree
#' @param nodes A node ID or vector of node IDs
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' # The antidote is irrelevant conditional upon Bite = "No", so we reflect this
#' # by bypassing the child node for Bite = "No".
#'
#' # We can identify the relevant node(s) to bypass using `getNodes`:
#' child <- getNodes(tree,varnames,"Bite","No")
#'
#' model <- bypassNodes(tree,child,varnames)
#' tree <- model$table
#' plotEventTree(tree,varnames)
#'
#' # Nodes can be passed directly into the function too:
#' model <- bypassNodes(tree,c(1,8),varnames)
#'
#' # If `infoTable` is a staged tree, you can preserve the node stage
#' # information during this operation:
#' model <- bypassNodes(tree,c(1,8),varnames,deleteStages=FALSE)
bypassNodes <- function(
    infoTable,
    nodes,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  failed_nodes <- c()
  failed_nodes <- c(nodes[infoTable$level[match(nodes,infoTable$id)]==max(infoTable$level)],nodes[nodes %in% infoTable$id[is.na(infoTable$parent)]])
  nodes <- setdiff(nodes,failed_nodes) # remove any leaves
  nodes <- nodes[order(infoTable$level[match(nodes, infoTable$id)], decreasing = TRUE)]

  for(i in nodes){
    model <- tryCatch(
      bypassNode(infoTable, i, deleteStages=deleteStages),
      error = function(e){
        failed_nodes <<- c(failed_nodes, i)
        warning(paste0("Node ", i, " could not be bypassed: ", conditionMessage(e)), call.=FALSE)
        NULL
      }
    )
    if(!is.null(model)) infoTable <- model$table
  }

  plot <- plotEventTree(infoTable, varnames, leafprobs,
                        padding=padding, leafnodesep=leafnodesep, varsep=varsep,
                        nodewidth=nodewidth, nodefont=nodefont, edgefont=edgefont,
                        probfont=probfont, varfont=varfont)

  result <- list(plot=plot$plot, dot=plot$dot, table=infoTable)
  if(length(failed_nodes)>0) result$failed_nodes <- failed_nodes
  result
}


#' Remove a set of Nodes and their Descendants From the Tree
#'
#' Prune a tree by deleting a set of nodes and the subtrees rooted at these nodes
#'
#' Used, for example, when an edge necessarily has zero probability thus all paths including this edge can be deleted.
#' In this case, the node that the zero-probability edge arrives at would be inputted.
#' The descendants of a node are the nodes that can be reached from the specified node following the edges of the tree.
#'
#' @inheritParams plotEventTree
#' @param cutnodes A node ID or vector of node IDs
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#' @export
#'
#' @examples
#'  variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' pruned <- getNodes(tree,varnames,"Antidote","Weak Antidote")
#' model <- PrunedTree(tree,pruned,varnames)
#' model
PrunedTree <- function(
    infoTable,
    cutnodes = NULL,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  if(0 %in% cutnodes) stop("You can't prune the root node!")
  if (length(cutnodes)==0) {
    prunedtable <- infoTable
  } else {
    if(deleteStages){
      infoTable$stage <- NULL
      infoTable$position <- NULL
    }
    prunedtable <- infoTable[-which(infoTable$id %in% descendants(infoTable,cutnodes)),]
  }
  plot <- plotEventTree(prunedtable, varnames, leafprobs,
                        padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                        nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                        probfont = probfont, varfont = varfont)

  list(plot = plot$plot, dot = plot$dot, table = prunedtable)
}

#' Generate the Subtree(s) Rooted at the Inputted Node(s)
#'
#' Generates a plot and information table for the event tree formed by the inputted nodes and their descendants.
#'
#' This may generate disconnected subtrees if inputting multiple nodes whose descendant sets do not overlap. However, this function is designed to handle such cases.
#' The descendants of a node are the nodes that can be reached from the specified node following the edges of the tree.
#'
#' @inheritParams plotEventTree
#' @param rootnodes A node ID or vector of node IDs
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the generated subtree and a plot of the generated subtree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' # To focus on just the paths downstream of Bite = "Yes":
#' root <- getNodes(tree,varnames,"Bite","Yes")
#' model <- SubTree(tree,root,varnames,padding=1)
#' model
SubTree <- function(
    infoTable,
    rootnodes = NULL,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  if(any(duplicated(unlist(lapply(rootnodes, function(x) descendants(infoTable, x)))))){
    warning("You have selected root nodes that are descendants of other root nodes! These nodes will be ignored.")
  }
  desc <- unique(unlist(lapply(rootnodes, function(x) descendants(infoTable, x))))

  subtree <- infoTable[which(infoTable$id %in% desc), ]
  if(length(varnames)==max(subtree$level)) varnames <- varnames[unique(subtree$level)[-1]]
  subtree$level <- subtree$level - min(subtree$level)
  if(deleteStages){
    subtree$stage <- NULL
    subtree$position <- NULL
  }
  plot <- plotEventTree(subtree,varnames,
                        padding=padding,leafnodesep=leafnodesep,varsep = varsep,
                        nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                        probfont = probfont, varfont = varfont)
  list(plot = plot$plot, dot = plot$dot, table = subtree)
}

#' Generate an Ancestral Tree
#'
#' Generates a plot and information table for the event tree formed by the inputted nodes and their ancestors.
#'
#' This may generate disconnected subtrees if inputting multiple nodes whose descendant sets do not overlap. However, this function is designed to handle such cases.
#' The ancestors of a node are the nodes from which the specified node can be reached following the edges of the tree.
#'
#' @inheritParams plotEventTree
#' @param nodes A node ID or vector of node IDs
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the generated subtree and a plot of the generated subtree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' # To focus on just the paths upstream of Die = "Yes":
#' anc <- getNodes(tree,varnames,"Die","Yes")
#' model <- AncestralTree(tree,anc,varnames,padding=1)
#' model
AncestralTree <- function(
    infoTable,
    nodes,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){
  if(identical(0,nodes)) stop("You are not allowed to input just the root node here! Anyway, why would you even want to?")
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  inc <- ancestors(infoTable,nodes)
  exc <- setdiff(infoTable$id,inc)
  PrunedTree(infoTable,exc,varnames,leafprobs = NULL,padding=padding,leafnodesep=leafnodesep,varsep = varsep,
             nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,probfont = probfont, varfont = varfont)
}

#' Connect a Set of Nodes Direct to a Set of Leaf Nodes
#'
#' Bypasses all intermediate nodes between inputted nodes and leaf nodes. Assigns one leaf node to each inputted nodes, and deletes intermediate nodes and edges.
#'
#' Bypassing edges do not get assigned a label automatically. Use function `renameEdge` to assign such edges their labels if desired.
#'
#' @inheritParams plotEventTree
#' @param nodes A node ID or vector of node IDs
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#' @param edges Optional character vector assigning labels to the new edges resulting from the operation
#'
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#' @export
#'
#' @examples
#'  variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#'
#' tree <- EventTreeTable(variables)
#'
#' # If death can only occur through a bite, and the antidote is truly perfect:
#' model <- collapseSubTree(tree,c(1,8),rep("Survival",2),varnames = varnames)
#' model
collapseSubTree <- function(
    infoTable,
    nodes,
    edges = NULL,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){

  nodes <- setdiff(nodes,c(0,infoTable$id[which(infoTable$level==max(infoTable$level))]))
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  if(length(nodes)==0){
    collapsetree <- infoTable
  } else{
  newrows <- list()
  for(i in seq_along(nodes)){

    node <- nodes[i]

    # Identify descendants of node (excluding node itself)
    desc <- descendants(infoTable, node)[-1]

    # Find a leaf among the descendants
    desclevels <- infoTable$level[match(desc, infoTable$id)]
    leafid <- desc[which(desclevels == max(desclevels))][1]

    # Create the bypass edge
    newrow <- infoTable[match(leafid, infoTable$id), , drop = FALSE]

    newrow$id <- leafid
    newrow$label <- paste0("s", leafid)
    newrow$parent <- node
    newrow$edge <- if (is.null(edges)) NA else edges[i]
    newrow$level <- max(desclevels)
    newrow$edgeprob <- 1

    newrows[[i]] <- newrow
  }

  all_desc <- unique(unlist(
    lapply(nodes, function(node) descendants(infoTable, node)[-1])
  ))

  # Remove collapsed subtrees
  tree2 <- infoTable[!(infoTable$id %in% all_desc), ]

  # Add bypass edges
  newrows <- do.call(rbind, newrows)

  collapsetree <- rbind(tree2, newrows)

  # Sort tree
  collapsetree <- collapsetree[
    order(collapsetree$level, collapsetree$id),
  ]

  rownames(collapsetree) <- NULL
  }
  plot <- plotEventTree(collapsetree,varnames, leafprobs,
                        padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                        nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                        probfont = probfont, varfont = varfont)
  list(plot = plot$plot, dot = plot$dot, table = collapsetree)
}

#' Generate a Tree Only Containing Paths that Contain the inputted nodes
#'
#' Generates the subtree formed by the root-to-leaf paths that pass through the selected nodes.
#'
#' A root-to-leaf path is any ordered set of nodes that can be traversed on the tree following its directed edges, originating at the root node and terminating at any leaf node.
#'
#' @inheritParams collapseSubTree
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' # To focus on the two antidote options:
#' nodes <- getNodes(tree,varnames,"Antidote",c("Weak Antidote","Perfect Antidote"))
#' model <- PathsThroughNodes(tree,nodes,varnames)
#' model
PathsThroughNodes <- function(
    infoTable,
    nodes,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){
  if(deleteStages){
    infoTable$stage <- NULL
    infoTable$position <- NULL
  }
  inc <- unique(c(ancestors(infoTable,nodes),descendants(infoTable,nodes)))
  exc <- setdiff(infoTable$id,inc)
  PrunedTree(infoTable,exc,varnames,leafprobs = NULL,padding=padding,leafnodesep=leafnodesep,varsep = varsep,
             nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,probfont = probfont, varfont = varfont)
}

#' Add Edges to Complete Partial Florets Within a Subtree
#'
#' Checks whether any nodes have an incomplete set of outgoing edges and, if any edges are missing, reintroduces them for the sake of checking sum-to-one probability conditions.
#'
#' A floret is a node together with its set of outgoing edges.
#' When manipulating a tree, some non-zero-probability edges may be pruned.
#' This function reintroduces these edges, but does not add subsequent edges to these, in order to validate the probability vector assigned to each node in the tree.
#'
#' @inheritParams plotEventTree
#' @param fulltable A table describing the complete event tree from which the subtree or modified tree has been derived
#'
#' @return A table describing the modified tree and a plot of the modified tree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' nodes <- getNodes(tree,varnames,"Antidote",c("Weak Antidote","Perfect Antidote"))
#' model <- PathsThroughNodes(tree,nodes,varnames)
#'
#' modtree <- model$table
#' CompletedFlorets(tree,modtree,varnames)
CompletedFlorets <- function(
    fulltable,
    infoTable,
    varnames = NULL,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30){
  ptParents <- c(NA,infoTable$id[which(infoTable$id %in% infoTable$parent)]) # which nodes in prunedtable are parents
  inc <- fulltable$id[which(fulltable$parent %in% ptParents)] # which children do these parents have in fulltable, inc. the prunedtable root
  exc <- setdiff(fulltable$id,inc)
  PrunedTree(fulltable,exc,varnames,leafprobs = NULL,padding=padding,leafnodesep=leafnodesep,varsep = varsep,
             nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,probfont = probfont, varfont = varfont)
}

#' Get the Node IDs Corresponding to an Outcome of a Variable
#'
#' Finds the IDs of the lowest-level nodes for which the specified variable takes the specified state, or one of the specified states.
#' That is, these nodes and their descendants lie along paths in which the inputted variable takes the inputted state(s) as its outcome.
#' The subtrees rooted at these nodes describe the potential unfoldings conditional upon the specified variable having this state/these states.
#'
#' Used to help identify node IDs for modifying a tree
#'
#' @inheritParams plotEventTree
#' @param variable The name of modelled variable
#' @param state A state or vector of states of the specified variable
#'
#' @return The node IDs whose subtrees correspond to the paths that necessarily contain the specified variable taking (one of) the specified state(s)
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#' tree <- EventTreeTable(variables)
#'
#' getNodes(tree,varnames,"Die","No")
#'
#' getNodes(tree,varnames,"Antidote",c("Weak Antidote","Perfect Antidote"))
#'
#' getNodes(tree,names(variables),names(variables)[2],variables[[2]][1])
#'
getNodes <- function(infoTable,varnames,variable,state){
  varlevel <- which(varnames==variable)
  subtable <- infoTable[infoTable$level==varlevel,]
  subtable$id[subtable$edge%in%state]
}

#' Look Up Nodes in an Event Tree
#'
#' Returns the nodes downstream of a specified combination of variable states.
#'
#' @param infoTable A table describing the event tree.
#' @param varnames A vector of variable names in level order.
#' @param states A named vector of selected states, with NA indicating any state.
#'
#' @return A vector of node IDs.
#' @export
TreeLookup <- function(infoTable, varnames, states) {

  if (is.null(names(states)) || any(!names(states) %in% varnames)) {
    stop("states must be a named vector using names in varnames.")
  }

  # Retain only specified states and identify their variable levels
  specified <- states[!is.na(states)]
  names(specified) <- as.character(match(names(specified), varnames))

  # Find the root node
  root <- infoTable$id[infoTable$level == 0][1]

  # Find all descendants of a node, including the node itself
  get_descendants <- function(node) {

    children <- infoTable$id[
      !is.na(infoTable$parent) & infoTable$parent == node
    ]

    if (length(children) == 0) return(node)

    c(node, unlist(lapply(children, get_descendants)))
  }

  # If no states are specified, return every node downstream of the root
  if(length(specified) == 0){
    result <- get_descendants(root)[-1]
    return(list(all = result, nextlevel = result))
  }

  # Traverse the tree, retaining paths compatible with the specified states
  traverse <- function(node, remaining){

    children <- infoTable[
      !is.na(infoTable$parent) & infoTable$parent == node,
    ]

    if(nrow(children) == 0){
      return(list(all = integer(0), nextlevel = integer(0)))
    }

    found <- integer(0)
    nextlevel <- integer(0)

    parent_level <- infoTable$level[match(node, infoTable$id)]
    edge_variable <- parent_level + 1

    for(j in seq_len(nrow(children))){

      child <- children[j, ]
      next_remaining <- remaining

      # Check whether this edge represents a specified variable
      variable <- as.character(edge_variable)

      if(variable %in% names(remaining)){

        if(as.character(child$edge) != remaining[variable]) next

        next_remaining <- remaining[names(remaining) != variable]
      }

      # Reject paths that bypass a variable with a specified state
      skipped <- as.numeric(names(next_remaining))

      if(any(skipped > edge_variable & skipped <= child$level)) next

      # Once all specified states have been matched, collect downstream nodes
      if(length(next_remaining) == 0){
        nextlevel <- c(nextlevel, child$id)
        found <- c(found, get_descendants(child$id))
      } else {
          deeper <- traverse(child$id, next_remaining)
          found <- c(found, deeper$all)
          nextlevel <- c(nextlevel, deeper$nextlevel)
        }
      }

    list(
      all = unique(found),
      nextlevel = unique(nextlevel)
    )
  }

  result <- traverse(root, specified)

  list(
    all = infoTable$id[infoTable$id %in% result$all],
    nextlevel = infoTable$id[infoTable$id %in% result$nextlevel]
  )
}


#' Reassign Node IDs After Tree Modification
#'
#' Reassigns the integers \eqn{1:N} as the node IDs in the tree.
#' After performing tree operations, the original node IDs are preserved, even if it results in integers being skipped.
#' This function reverses this effect.
#'
#' @inheritParams plotEventTree
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @return A table describing the renumbered tree and a plot of the renumbered tree, alongside its dot representation
#' @export
#'
#' @examples
#' variables <- list(
#'   Bite       =   c("No","Yes"),
#'   Antidote   =   c("None","Weak Antidote","Perfect Antidote"),
#'   Die        =   c("No","Yes")
#' )
#' varnames <- names(variables)
#'
#' tree <- EventTreeTable(variables)
#'
#' model <- collapseSubTree(tree,c(1,8),rep("Survival",2))
#' model <- PrunedTree(model$table,7,varnames)
#' model <- reassignID(model$table,varnames)
#' model
reassignID <- function(
    infoTable,
    varnames = NULL,
    leafprobs = NULL,
    deleteStages = TRUE,
    padding = 0.5,
    leafnodesep = 0.15,
    varsep = 2,
    nodewidth = 0.5,
    nodefont = 18,
    edgefont = 20,
    probfont = 18,
    varfont = 30) {
  idmap <- stats::setNames(seq_len(nrow(infoTable)) - 1, sort(infoTable$id))
  infoTable$id <- unname(idmap[as.character(infoTable$id)])
  infoTable$parent <- unname(idmap[as.character(infoTable$parent)])
  infoTable$label <- paste0("s", infoTable$id)
  plot <- plotEventTree(infoTable, varnames, leafprobs,
                        padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                        nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                        probfont = probfont, varfont = varfont)

  list(plot = plot$plot, dot = plot$dot, table = infoTable)
}

#' Plot an Existing Tree with an Added Stage Colouring
#'
#' Plots an existing tree from its dot representation but adds a stage formed of inputted nodes shown graphically by colouring the node with the inputted colour.
#'
#' The function outputs the plot and the new modified dot information representing the colouring. The inputted dot representation may already include coloured stages.
#'
#' @param dot The dot representation of a plot for a tree
#' @param nodes A node ID or vector of node IDs
#' @param colour A colour
#' @param shape A character shape input, defaults to "ellipse"
#'
#' @return A plot of the tree with the additional stage colouring and its dot representation
#' @export
addTreeStage <- function(dot, nodes, colour, shape = "ellipse"){

  if(!inherits(nodes, "character")) nodes <- paste0("s", nodes)

  additions <- ""

  for (node in nodes) {
    additions <- paste0(
      additions,
      sprintf(
        '"%s" [style=filled, fillcolor="%s", shape="%s"];\n',
        node,
        colour,
        shape
      )
    )
  }

  dot <- sub("\\}\\s*$", paste0(additions, "}"), dot)

  list(plot=DiagrammeR::grViz(dot),dot=dot)
}

#' Plot an Existing Tree with Added Stage Colourings
#'
#' Plots an existing tree from its dot representation but adds stages formed of inputted nodes shown graphically by colouring the nodes with their corresponding colour.
#'
#' The function outputs the plot and the new modified dot information representing the colouring. The inputted dot representation may already include coloured stages.
#'
#' @param dot The dot representation of a plot for a tree
#' @param nodelist A list of vectors of node IDs
#' @param colourvec An vector of colours
#'
#' @return A plot of the tree with the additional stage colouring and its dot representation
#' @export
addTreeStageList <- function(
    dot,
    nodelist,
    colourvec = c(
      "skyblue", "firebrick3", "plum", "chartreuse2",
      "darkorange", "grey67", "peachpuff", "hotpink",
      "gold", "violetred2", "goldenrod1", "olivedrab4",
      "steelblue2", "ivory2", "darkorchid2", "khaki3",
      "aquamarine3", "tan2"
    )){

  shapes <- c("ellipse", "rectangle", "hexagon", "doublecircle", "doubleoctagon")

  for(i in seq_along(nodelist)){

    colour_index <- ((i - 1) %% length(colourvec)) + 1
    shape_index <- ((i - 1) %/% length(colourvec)) + 1

    shape <- shapes[((shape_index - 1) %% length(shapes)) + 1]

    newtree <- addTreeStage(
      dot,
      nodelist[[i]],
      colourvec[colour_index],
      shape
    )

    dot <- newtree$dot
  }

  list(plot=DiagrammeR::grViz(dot),dot=dot)
}

#' Get Stages from Edge Probabilities
#'
#' Assigns each node its stage based on its vector of outgoing edge probabilities and the stages of its children.
#'
#' Two nodes are in the same stage if they have identical probability vectors over their florets.
#'
#' @inheritParams TreeLookup
#'
#' @return A table uniquely describing the staged tree
#' @export
getStages <- function(infoTable){
  stagedtree <- infoTable
  stagedtree$stage <- 0
  counter <- 0
  for(j in (max(stagedtree$level)-1):0){
    floretlist <- list()
    nodes <- stagedtree$id[which(stagedtree$level==j)]
    for(i in nodes){
    children <- sort(stagedtree$id[which(stagedtree$parent==i)])
    #childstages <- stagedtree$stage[which(stagedtree$id%in%children)]
    childedges <- stagedtree$edge[which(stagedtree$id%in%children)]
    childprobs <- stagedtree$edgeprob[which(stagedtree$id%in%children)]
    floretlist[[as.character(i)]] <- c(childedges,childprobs) # childstages was in this vector before - TEST TO CHECK WHAT IT SHOULD BE
    }
    stagegroups <- unname(split(names(floretlist),sapply(floretlist, paste, collapse = ",")))
    for(k in seq_along(stagegroups)){
      if(length(stagegroups[[k]])>1) stagedtree$stage[which(stagedtree$id %in% stagegroups[[k]])] <- min(stagedtree$stage)-1
    }
  }
  stagedtree$stage <- -stagedtree$stage
  posStage <- stagedtree$stage>0
  stagedtree$stage[posStage] <- max(stagedtree$stage[posStage])-stagedtree$stage[posStage]+1
  stagedtree$stage[stagedtree$level==max(stagedtree$level)] <- Inf
  stagedtree
}


#' Plot a Staged Tree from the Staged Tree Table
#'
#' Plots a staged tree using the stage information in the table describing the staged tree.
#'
#' If the table does not have stage information, this simply plots the event tree.
#'
#' @inheritParams plotEventTree
#' @param colourvec A vector of colours
#' @param deleteStages If the table specifies a staged tree, whether to remove the stage information. Defaults to TRUE.
#'
#' @export
plotStagedTree <- function(infoTable,
                           varnames = NULL,
                           leafprobs = NULL,
                           deleteStages = TRUE,
                           padding = 0.5,
                           leafnodesep = 0.15,
                           varsep = 2,
                           nodewidth = 0.5,
                           nodefont = 18,
                           edgefont = 20,
                           probfont = 18,
                           varfont = 30,
                           colourvec = c(
                             "skyblue", "firebrick3", "plum", "chartreuse2",
                             "darkorange", "grey67", "peachpuff", "hotpink",
                             "gold", "violetred2", "goldenrod1", "olivedrab4",
                             "steelblue2", "ivory2", "darkorchid2", "khaki3",
                             "aquamarine3", "tan2"
                           )){

  if(is.null(infoTable$stage) || all(infoTable$stage==0) || all(is.na(infoTable$stage))){
    plotEventTree(infoTable, varnames, leafprobs,
                  padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                  nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                  probfont = probfont, varfont = varfont)
  }
  else{
    stagedmodel <- plotEventTree(infoTable, varnames, leafprobs,
                                 padding = padding, leafnodesep = leafnodesep, varsep = varsep,
                                 nodewidth = nodewidth, nodefont = nodefont, edgefont = edgefont,
                                 probfont = probfont, varfont = varfont)
    nonleaf_nodes <- infoTable$id[
      infoTable$id %in% unlist(infoTable$parent)
    ]

    staged_nodes <- infoTable$id[
      infoTable$stage > 0 & infoTable$id %in% nonleaf_nodes
    ]

    stagelist <- split(
      staged_nodes,
      infoTable$stage[match(staged_nodes, infoTable$id)]
    )
    addTreeStageList(stagedmodel$dot,stagelist,colourvec)
  }
}


#' Get Position Information from Node Stages
#'
#' Uses node stage information to infer the position of each node.
#'
#' Two nodes are in the same position if they are in the same stage and have topologically and probabilistically identical subtrees.
#' That is, their subtrees have an identical structure and the stages found in equivalent positions across the two structures are the same.
#'
#' @param stagedTable A table uniquely specifying the staged tree
#'
#' @return A table uniquely specifying the staged tree including each node's position
#' @export
PositionsFromStages <- function(stagedTable){

  stagedTable$position <- NA  # create a vector indicating a node's CEG position label
  stagedTable$position[which(stagedTable$level==max(stagedTable$level))] <- nrow(stagedTable)   # used to indicate that all leaves have equivalent subtrees (trivially) though these aren't proper positions

  levels <- (max(stagedTable$level)-1):0
  for(i in levels){
    nodes <- stagedTable$id[which(stagedTable$level==i)]
    levelparentinfo <- list()   #PARENTS THAT ARE ONE LEVEL UP ONLY!
    for(j in nodes){
      parentstage <- stagedTable$stage[which(stagedTable$id==j)]
      if(parentstage==0 || is.na(parentstage)){
        levelparentinfo[[as.character(j)]] <- paste0("NA_", j)}  # unique stage => unique position
      else{
        rows <- which(stagedTable$parent==j)
        # get parent stage (must be in same stage for same position) and child position
        levelparentinfo[[as.character(j)]] <- c(parentstage,stagedTable$position[rows])    # position => same subtree structure inc future positions => same stage and same next positions (in order!) (stage already checks edge labels, edge probs)
      }
    }   # levelparent info gives the stage of the parent and then a vector of subsequent positions it can lead to
    positiongroups <- split(
      names(levelparentinfo),
      sapply(levelparentinfo, paste, collapse = ",")
    )
    positiongroups <- lapply(positiongroups, as.numeric)
    positiongroups <- positiongroups[order(sapply(positiongroups, min),decreasing=TRUE)]
    for(k in seq_along(positiongroups)){
      stagedTable$position[stagedTable$id%in%positiongroups[[k]]] <- min(stagedTable$position,na.rm=TRUE)-1  ### THIS ASSUMES A ROW+1 NUMBERING ON THE NODES
    }
  }
  stagedTable$position <- stagedTable$position-min(stagedTable$position)
  stagedTable$position[which(stagedTable$level==max(stagedTable$level))] <- Inf
  stagedTable
}

#' Get Position Information from Event Tree Table
#'
#' Uses edge probability information to infer the stage and position of each node.
#'
#' Two nodes are in the same position if they are in the same stage and have topologically and probabilistically identical subtrees.
#' That is, their subtrees have an identical structure and the stages found in equivalent positions across the two structures are the same.
#'
#' @param infoTable A table uniquely specifying the event tree
#'
#' @return A table uniquely specifying the staged tree including each node's position
#' @export
PositionsFromEventTree <- function(infoTable){
  if(length(checkEdgeProbs(infoTable)$invalid)>0) stop("There are invalid edge probabilities!")
  PositionsFromStages(getStages(infoTable))
}

#' Generate Information About the Nodes of the CEG Formed by a Staged Tree
#'
#' Produces a table uniquely describing the nodes of the CEG.
#'
#' The outputted table documents the node ID, the node label, its level in the CEG and the stage ID it represents for each node in the CEG.
#' Where a stage comprising one node is represented by a node, its stage is noted as NA.
#'
#' @param positionTable A table uniquely specifying a staged tree including the position of each node
#'
#' @return A table documenting the nodes of the generated CEG, including ID, label, level and stage information
#' @export
getCEGnodes <- function(positionTable){
  positions <- unique(positionTable$position)
  CEGnodes <- data.frame(
    id = positions,
    label = paste0("w", positions),
    level = sapply(positions, function(p){max(positionTable$level[positionTable$position == p])}),
    stage = sapply(positions, function(p) unique(positionTable$stage[positionTable$position == p])),
    parents = I(lapply(positions, function(p) {
      nodes <- positionTable$id[positionTable$position == p]
      parent_ids <- unique(positionTable$parent[positionTable$id %in% nodes])
      unique(positionTable$position[positionTable$id %in% parent_ids])
    }))
  )
  CEGnodes
}

#' Generate Information About the Edges of the CEG Formed by a Staged Tree
#'
#' Produces a table uniquely describing the edges of the CEG.
#'
#' The outputted table documents the edge origin, the edge destination, the edge label and the edge probability for each edge in the CEG.
#'
#' @param positionTable A table uniquely specifying a staged tree including the position of each node
#'
#' @return A table documenting the edges of the generated CEG, including origin, destination, label and edge probability information
#' @export
getCEGedges <- function(positionTable){

  positionTable$position_from <- as.numeric(positionTable$position[match(positionTable$parent, positionTable$id)])
  positionTable <- positionTable[order(positionTable$position_from, positionTable$position,na.last=FALSE), ]
  CEGedges <- unique(positionTable[,c("position_from","position","edge","edgeprob")])
  CEGedges <- CEGedges[CEGedges$position!=0,]
  rownames(CEGedges) <- NULL
  names(CEGedges) <- c("pos_from","pos_to","edge","edgeprob")
  CEGedges$edge_id <- seq_len(nrow(CEGedges))
  CEGedges
}

#' Obtain Information About the Nodes and Edges of a CEG Formed by a Staged Tree
#'
#' Produces a table uniquely describing the nodes of the CEG, and a table uniquely describing the edges of the CEG.
#'
#' This simultaneously applies the functions `getCEGnodes` and `getCEGedges`. The two tables uniquely specify a CEG.
#'
#' @param positiontable A table uniquely specifying a staged tree including the position of each node
#'
#' @return A table specifying the nodes of the generated CEG, and a table specifying the edges of the generated CEG.
#' @export
getCEGtables <- function(positiontable){
  list(
    CEGnodes = getCEGnodes(positiontable),
    CEGedges = getCEGedges(positiontable)
  )
}

#' Plot a Chain Event Graph
#'
#' Plots the CEG specified by the table outputs of `getCEGnodes` and `getCEGedges` (or `getCEGtables`).
#'
#' @inheritParams plotEventTree
#' @param CEGnodes A table specifying the nodes of the CEG (see `getCEGnodes`)
#' @param CEGedges A table specifying the edges of the CEG (see `getCEGedges`)
#' @param AllPosColoured A logical value indicating whether to colour nodes whose position comprises only one node from the staged tree
#' @param colourvec A vector of colours
#' @param splines A logical value indicating how edges are drawn in the CEG plot
#'
#' @return A plot of the CEG alongside its dot representation
#' @export
plotCEG <- function(
    CEGnodes,
    CEGedges,
    varnames = NULL,
    AllPosColoured = FALSE,
    colourvec = c(
      "skyblue", "firebrick3", "plum", "chartreuse2",
      "darkorange", "grey67", "peachpuff", "hotpink",
      "gold", "violetred2", "goldenrod1", "olivedrab4",
      "steelblue2", "ivory2", "darkorchid2", "khaki3",
      "aquamarine3", "tan2"
    ),
    padding = 0.5,
    leafnodesep = 0.5,
    varsep = 1.5,
    nodewidth = 0.75,
    nodefont = 32,
    edgefont = 16,
    varfont = 30,
    splines = TRUE
){
  dot <- sprintf("
  digraph CEG {

  graph [
    rankdir=LR,
    nodesep=%s,
    ranksep=%s,
    pad=%s,
    splines=%s    # 'true','line','polyline'
  ]

  node [
    shape=circle,
    fixedsize=true,
    width=%s,
    fontsize=%s
  ]

  edge [
    fontsize=%s
  ]
  ",leafnodesep, varsep, padding, splines, nodewidth, nodefont,edgefont)

  CEGnodes$stage <- as.numeric(CEGnodes$stage)

  if(AllPosColoured){
    CEGnodes$uniqueStage <- CEGnodes$stage
    CEGnodes$uniqueStage[CEGnodes$uniqueStage==0] <- seq(max(CEGnodes$stage[is.finite(CEGnodes$stage)]) + 1,length.out = sum(CEGnodes$stage==0))
    if(max(CEGnodes$uniqueStage[is.finite(CEGnodes$uniqueStage)])>length(colourvec)) stop("There are not enough colours in your colour vector!")
  } else if(nrow(CEGnodes)>length(colourvec)) stop("There are not enough colours in your colour vector!")

  # Add nodes
  for(i in seq_len(nrow(CEGnodes))){

    stage <- CEGnodes$stage[i]

    if(CEGnodes$label[i] == "winf"){
      node_label <- "<<I>w</I><SUB><FONT POINT-SIZE='28'>&#8734;</FONT></SUB>>"
      html_label <- TRUE
      dot <- paste0(
        dot,
        sprintf(
          '"%s" [label=%s, shape=plaintext];\n',
          as.character(CEGnodes$id[i]),
          node_label
        )
      )
    } else{
      node_label <- CEGnodes$label[i]
      html_label <- FALSE
    }

    if(stage>0 && CEGnodes$label[i]!="winf"){
      if(html_label){
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label=%s, style=filled, fillcolor="%s"];\n',
            as.character(CEGnodes$id[i]),
            node_label,
            colourvec[as.numeric(stage)]
          )
        )
      } else{
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label="%s", style=filled, fillcolor="%s"];\n',
            as.character(CEGnodes$id[i]),
            node_label,
            colourvec[as.numeric(stage)]
          )
        )
      }
    }

    if(stage==0 && !AllPosColoured){
      if(html_label){                     # don't need this now as html_label should only be for leaf node, but kept in case of label changes
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label=%s];\n',
            as.character(CEGnodes$id[i]),
            node_label
          )
        )
      } else{
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label="%s"];\n',
            as.character(CEGnodes$id[i]),
            node_label
          )
        )
      }
    }

    if(stage==0 && AllPosColoured){
      uniquestage <- as.integer(CEGnodes$uniqueStage[i])
      if(html_label){
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label=%s, style=filled, fillcolor="%s"];\n',
            as.character(CEGnodes$id[i]),
            node_label,
            colourvec[uniquestage]
          )
        )
      } else{
        dot <- paste0(
          dot,
          sprintf(
            '"%s" [label="%s", style=filled, fillcolor="%s"];\n',
            as.character(CEGnodes$id[i]),
            node_label,
            colourvec[uniquestage]
          )
        )
      }

    }
  }

  # Add variable labels above each level
  if(!is.null(varnames)){
    for(lev in seq_along(varnames)){
      rows <- which(CEGnodes$level == lev - 1)
      dot <- paste0(
        dot,
        sprintf(
          '"var%d" [label=<%s>, shape=plaintext, fontsize=%s, fontname="Computer Modern Roman"];\n',
          lev,
          varnames[lev],
          varfont
        )
      )

      # Put variable label and nodes at this level on same rank
      dot <- paste0(dot, "{rank=same; ")

      dot <- paste0(dot, sprintf('"var%d"; ', lev))

      for(r in rows){
        dot <- paste0(
          dot,
          sprintf('"%s"; ', as.character(CEGnodes$id[r]))
        )
      }
      dot <- paste0(dot, "}\n")
    }
  }

  # Add edges
  pad <- paste(rep("&nbsp;", 3), collapse = "")
  for(i in seq_len(nrow(CEGedges))){

    from <- as.character(CEGedges$pos_from[i])
    to   <- as.character(CEGedges$pos_to[i])

    dot <- paste0(
      dot,
      sprintf(
        '"%s" -> "%s" [label=<%s (%s)%s>];\n',
        from,
        to,
        CEGedges$edge[i],
        as.character(CEGedges$edgeprob[i]),
        pad
      )
    )

  }
  dot <- paste0(dot, "\n}")
  list(
    CEG = DiagrammeR::grViz(dot),
    dot = dot
  )

}


#' Create an Interactive Plot of the CEG
#'
#' Creates an interactive plot of the generated CEG.
#'
#' @inheritParams plotCEG
#'
#' @param levelSeparation Optional parameter specifying spacing between variable levels
#' @param nodeSpacing Optional parameter specifying vertical spacing between nodes within a level
#' @param height Optional character parameter specifying the height of the interactive plot widget, in the form `"850px"'
#' @param sinklabels Optional character vector to plot multiple sink nodes, defaults to NULL
#' @param edgesink Optional vector indicating which edges connect to which sink, defaults to NULL
#' @param nodeprobs Optional Boolean input that makes the probability of arrival for a node appear when the node is hovered over, defaults to TRUE
#' @param returnnodes Optional Boolean input that makes the function output the CEG node information, defaults to FALSE
#' @param widthbyprob Optional Boolean input that scales the edge width by the probability of its traversal within its floret, defaults to FALSE
#' @param positions Optional named list of co-ordinates that forces node positions in the plot. Primarily designed for app functionality. Defaults to NULL.
#' @param view_position Optional named list of the co-ordinates of the plot's centre. Primarily designed for app functionality. Defaults to NULL.
#'
#' @return An interactive plot of the generated CEG.
#' @export
#'
plotCEG_vis <- function(
    CEGnodes,
    CEGedges,
    varnames = NULL,
    colourvec = c(
      "skyblue", "firebrick3", "plum", "chartreuse2",
      "darkorange", "grey67", "peachpuff", "hotpink",
      "gold", "violetred2", "goldenrod1", "olivedrab4",
      "steelblue2", "ivory2", "darkorchid2", "khaki3",
      "aquamarine3", "tan2"
    ),
    nodeprobs = TRUE,
    sinklabels = NULL,
    edgesink = NULL,
    returnnodes = FALSE,
    widthbyprob = FALSE,
    positions = NULL,
    view_position = NULL,
    nodefont = 24,
    edgefont = 16,
    varfont = 22,
    levelSeparation = 180,
    nodeSpacing = 100,
    height = "850px"
) {

  # --------------------------------------------------
  # 1. Prepare input tables
  # --------------------------------------------------

  nodes <- CEGnodes
  nodes$id <- as.character(nodes$id)
  nodes$level <- as.numeric(nodes$level)
  nodes$stage <- suppressWarnings(as.numeric(nodes$stage))

  CEGedges$pos_from <- as.character(CEGedges$pos_from)
  CEGedges$pos_to <- as.character(CEGedges$pos_to)

  # --------------------------------------------------
  # 2. Separate sink nodes, if requested
  # --------------------------------------------------

  separate_sinks <- !is.null(sinklabels) && length(sinklabels) > 0

  if (separate_sinks) {

    sinklabels <- trimws(as.character(sinklabels))

    if (anyNA(sinklabels) ||
        any(sinklabels == "") ||
        anyDuplicated(sinklabels)) {
      stop("Sink labels must be non-empty and unique.")
    }

    sink_row <- which(nodes$id == "inf" | nodes$id == Inf)

    if (length(sink_row) != 1) {
      stop("CEGnodes must contain exactly one sink with id 'inf'.")
    }

    sinkedges <- which(CEGedges$pos_to == "inf" | CEGedges$pos_to == Inf)

    if (is.null(edgesink)) {
      edgesink <- rep(1L, nrow(CEGedges))
    }

    edgesink_numeric <- suppressWarnings(as.numeric(edgesink))

    if(length(edgesink_numeric) != nrow(CEGedges) ||
        anyNA(edgesink_numeric) ||
        any(!is.finite(edgesink_numeric)) ||
        any(edgesink_numeric != floor(edgesink_numeric)) ||
        any(edgesink_numeric < 0 |
            edgesink_numeric > length(sinklabels))) {
      stop("Invalid sink assignment for one or more CEG edges.")
    }

    sink_level <- nodes$level[sink_row]
    new_ids <- paste0("sink_", seq_along(sinklabels))

    if (any(new_ids %in% nodes$id)) {
      stop("A generated sink ID already exists in CEGnodes.")
    }

    # Create sink rows preserving the original node columns.
    sink_nodes <- nodes[
      rep(sink_row, length(sinklabels)),
      ,
      drop = FALSE
    ]

    sink_nodes$id <- new_ids
    sink_nodes$label <- sinklabels
    sink_nodes$level <- sink_level
    sink_nodes$stage <- NA_real_

    # Replace the original sink in the local node table.
    remaining_nodes <- nodes[-sink_row, , drop = FALSE]

    for (nm in setdiff(names(remaining_nodes), names(sink_nodes))) {
      sink_nodes[[nm]] <- NA
    }

    for (nm in setdiff(names(sink_nodes), names(remaining_nodes))) {
      remaining_nodes[[nm]] <- NA
    }

    sink_nodes <- sink_nodes[, names(remaining_nodes), drop = FALSE]

    nodes <- rbind(remaining_nodes, sink_nodes)

    # Redirect the incoming sink edges in the local edge table.
    CEGedges$pos_to[sinkedges] <-
      new_ids[edgesink_numeric[sinkedges]]
  }

  # --------------------------------------------------
  # 3. Calculate node probabilities
  # --------------------------------------------------

  if (nodeprobs) {

    nodes$nodeprob <- 0
    nodes$nodeprob[nodes$level == 0] <- 1

    for (i in seq_len(max(nodes$level))) {

      levelnodes <- nodes$id[nodes$level == i]

      for (j in levelnodes) {

        parentrows <- CEGedges[
          CEGedges$pos_to == j,
          ,
          drop = FALSE
        ]

        if (nrow(parentrows) > 0) {

          parent_probs <- nodes$nodeprob[
            match(parentrows$pos_from, nodes$id)
          ]

          nodes$nodeprob[match(j, nodes$id)] <-
            sum(parent_probs * parentrows$edgeprob)

        } else {

          nodes$nodeprob[match(j, nodes$id)] <- 0

        }
      }
    }
  }

  # --------------------------------------------------
  # 4. Prepare node labels and colours
  # --------------------------------------------------

  node_labels <- as.character(nodes$label)
  node_labels[node_labels == "winf"] <- "w\u221e"
  node_labels[node_labels == "wInf"] <- "w\u221e"

  node_colours <- rep("white", nrow(nodes))
  node_shapes <- rep("circle", nrow(nodes))

  coloured <- !is.na(nodes$stage) &
    nodes$stage > 0 &
    !nodes$id %in% c("Inf", "inf") &
    !grepl("^sink_", nodes$id)

  # Identify distinct stages in ascending order.
  stage_ids <- sort(unique(nodes$stage[coloured]))

  # Assign colours and shapes by stage.
  shapes <- c(
    "ellipse", "box", "circle","database"
  )

  if (length(stage_ids) > 0) {

    colourvec <- grDevices::rgb(
      t(grDevices::col2rgb(colourvec)),
      maxColorValue = 255
    )

    for (i in seq_along(stage_ids)) {

      stage_nodes <- coloured & nodes$stage == stage_ids[i]

      colour_index <- ((i - 1) %% length(colourvec)) + 1
      shape_index <- ((i - 1) %/% length(colourvec)) %% length(shapes) + 1

      node_colours[stage_nodes] <- colourvec[colour_index]
      node_shapes[stage_nodes] <- shapes[shape_index]
    }
  }

  vis_nodes <- data.frame(
    id = nodes$id,
    label = node_labels,
    title = if (nodeprobs) {
      paste0(
        "Node probability: ",
        sprintf("%.3f", nodes$nodeprob)
      )
    } else {
      rep(NA_character_, nrow(nodes))
    },
    level = nodes$level,
    shape = node_shapes,
    color.background = node_colours,
    color.border = "black",
    size = 22,
    font.size = nodefont,
    stringsAsFactors = FALSE
  )

  if (!is.null(positions)) {
    matched <- match(vis_nodes$id, names(positions))
    has_position <- !is.na(matched)

    vis_nodes$x <- NA_real_
    vis_nodes$y <- NA_real_

    vis_nodes$x[has_position] <- vapply(
      positions[matched[has_position]],
      function(pos) pos$x,
      numeric(1)
    )

    vis_nodes$y[has_position] <- vapply(
      positions[matched[has_position]],
      function(pos) pos$y,
      numeric(1)
    )
  }

  # --------------------------------------------------
  # 5. Prepare edges
  # --------------------------------------------------

  edges <- CEGedges
  edges$from <- as.character(edges$pos_from)
  edges$to <- as.character(edges$pos_to)
  edges$prob <- as.character(edges$edgeprob)

  edges$edge_label <- ifelse(
    is.na(edges$edge),
    "",
    as.character(edges$edge)
  )

  # Distinguish edges sharing the same endpoints.
  pair <- paste(edges$from, edges$to, sep = "\r")

  edge_number <- stats::ave(
    seq_along(pair),
    pair,
    FUN = seq_along
  )

  edge_count <- stats::ave(
    seq_along(pair),
    pair,
    FUN = length
  )

  curve_type <- ifelse(
    edge_number %% 2 == 1,
    "curvedCW",
    "curvedCCW"
  )

  roundness <- pmin(
    0.15 + 0.20 * ceiling(edge_number / 2),
    0.85
  )

  vis_edges <- data.frame(
    id = edges$edge_id,
    from = edges$from,
    to = edges$to,
    label = paste0(edges$edge_label, " (", formatC(as.numeric(edges$prob), format = "f", digits = 3), ")"),
    title = paste0("Edge ID: ", edges$edge_id),
    font.size = edgefont,
    stringsAsFactors = FALSE
  )

  vis_edges$smooth <- I(lapply(seq_len(nrow(edges)), function(i) {
    list(
      enabled = edge_count[i] > 1,
      type = curve_type[i],
      roundness = roundness[i]
    )
  }))

  if(widthbyprob){
    vis_edges$value <- 0.5 + 3 * as.numeric(edges$prob)
  }

  # --------------------------------------------------
  # 6. Add variable labels
  # --------------------------------------------------

  if(!is.null(varnames) && length(varnames) > 0){

    var_levels <- seq_along(varnames) - 1

    # Only retain variable labels for levels containing nodes.
    existing_levels <- unique(vis_nodes$level)
    keep <- var_levels %in% existing_levels

    if (any(keep)) {

      var_nodes <- data.frame(
        id = paste0("var_", which(keep)),
        label = as.character(varnames[keep]),
        title = NA_character_,
        level = var_levels[keep],
        shape = "text",
        color.background = "rgba(0,0,0,0)",
        color.border = "rgba(0,0,0,0)",
        size = 1,
        font.size = varfont,
        stringsAsFactors = FALSE
      )

      if ("x" %in% names(vis_nodes)) {
        var_nodes$x <- NA_real_
        var_nodes$y <- NA_real_
      }

      vis_nodes <- rbind(vis_nodes, var_nodes)
    }
  }

  # --------------------------------------------------
  # 7. Render interactive graph
  # --------------------------------------------------

  p <- visNetwork::visNetwork(
    vis_nodes,
    vis_edges,
    width = "100%",
    height = height
  ) %>%
    visNetwork::visNodes(
      font = list(
        vadjust = 0,
        align = "center"
      ),
    ) %>%
    visNetwork::visEdges(
      arrows = "to",
      color = list(color = "grey40"),
      scaling = list(min = 0.5, max = 3.5),
      selectionWidth = 4,
      labelHighlightBold = TRUE
    )

  p <- p %>%
    visNetwork::visHierarchicalLayout(
      direction = "LR",
      levelSeparation = levelSeparation,
      nodeSpacing = nodeSpacing,
      sortMethod = "directed",
      shakeTowards = "roots",
      blockShifting = FALSE,
      edgeMinimization = FALSE,
      parentCentralization = FALSE
    )

  saved_positions_js <- if (!is.null(positions)) {
    jsonlite::toJSON(positions, auto_unbox = TRUE)
  } else {
    "null"
  }

  saved_view_position_js <- if (!is.null(view_position)) {
    jsonlite::toJSON(view_position, auto_unbox = TRUE)
  } else {
    "null"
  }

  p <- p %>%
    visNetwork::visInteraction(
      dragNodes = TRUE,
      dragView = TRUE,
      zoomView = TRUE,
      navigationButtons = TRUE,
      selectable = TRUE,
      multiselect = TRUE
    ) %>%
    visNetwork::visEvents(
      stabilizationIterationsDone = paste0("
    function () {
      var savedPositions = ", saved_positions_js, ";
      var savedViewPosition = ", saved_view_position_js, ";
      var network = this;
      var positions = network.getPositions();
          var allNodes = network.body.data.nodes.get();

          var cegNodes = allNodes.filter(function(node) {
            return !String(node.id).startsWith('var_') &&
                   !String(node.id).startsWith('sink_');
          });

          var levelPositions = {};

          cegNodes.forEach(function(node) {
            var pos = positions[node.id];
            if (!pos) return;

            var level = Number(node.level);

            if (!levelPositions[level]) {
              levelPositions[level] = { x: [], y: [] };
            }

            levelPositions[level].x.push(pos.x);
            levelPositions[level].y.push(pos.y);
          });

          var highestY = Infinity;

          Object.keys(levelPositions).forEach(function(key) {
            var ys = levelPositions[key].y;
            var minLevelY = Math.min.apply(null, ys);

            if (minLevelY < highestY) {
              highestY = minLevelY;
            }
          });

          var maxX = -Infinity;
          var minY = Infinity;
          var maxY = -Infinity;

          cegNodes.forEach(function(node) {
            var pos = positions[node.id];
            if (!pos) return;

            maxX = Math.max(maxX, pos.x);
            minY = Math.min(minY, pos.y);
            maxY = Math.max(maxY, pos.y);
          });

          var sinkNodes = allNodes.filter(function(node) {
            return String(node.id).startsWith('sink_');
          });

          if (sinkNodes.length > 0 &&
              isFinite(maxX) &&
              isFinite(minY) &&
              isFinite(maxY)) {

            var sinkX = positions[sinkNodes[0].id].x;
            var centreY = (minY + maxY) / 2;
            var sinkGap = 70;

            var defaultSink = sinkNodes.find(function(node) {
              return String(node.id) === 'sink_1';
            });

            var otherSinks = sinkNodes.filter(function(node) {
              return String(node.id) !== 'sink_1';
            });

            otherSinks.sort(function(a, b) {
              return Number(String(a.id).replace('sink_', '')) -
                     Number(String(b.id).replace('sink_', ''));
            });

            var medianIndex = Math.floor(sinkNodes.length / 2);
            var orderedSinks = otherSinks.slice();

            if (defaultSink) {
              orderedSinks.splice(medianIndex, 0, defaultSink);
            }

            orderedSinks.forEach(function(node, index) {
              network.moveNode(
                node.id,
                sinkX,
                centreY +
                  (index - (orderedSinks.length - 1) / 2) * sinkGap
              );
            });
          }

          var variableNodes = allNodes.filter(function(node) {
            return String(node.id).startsWith('var_');
          });

          variableNodes.forEach(function(node) {
            var index = Number(
              String(node.id).replace('var_', '')
            ) - 1;

            var level = levelPositions[index];

            if (!level ||
                level.x.length === 0 ||
                !isFinite(highestY)) {
              return;
            }

            var meanX = level.x.reduce(function(a, b) {
              return a + b;
            }, 0) / level.x.length;

            network.moveNode(node.id, meanX + ", levelSeparation / 2.5, ", highestY - 75);
          });

         if (savedPositions) {
            Object.keys(savedPositions).forEach(function(id) {
              var pos = savedPositions[id];
              var currentPos = positions[id];

              if(network.body.nodes[id] && pos && isFinite(pos.y) && !String(id).startsWith('var_')){

                if (currentPos && isFinite(currentPos.x)) {
                  network.moveNode(id, currentPos.x, pos.y);
                } else {
                  network.moveNode(id, pos.x, pos.y);
                }
              }
            });
          }

          if (savedViewPosition) {
            network.moveTo({
              position: savedViewPosition,
              animation: false
            });
          }

          network.setOptions({
            physics: { enabled: false }
          });
          Shiny.setInputValue(
              'ceg_node_positions',
              network.getPositions(),
              {priority: 'event'}
            );
        console.log('CEG positions captured');
        }
      "),

      selectEdge = "
        function(properties) {
          var network = this;
          var selected = (properties.edges || []).slice();

          network.__cegSelectedEdges = selected;

          Shiny.setInputValue(
            'ceg_selected_edges',
            selected,
            {priority: 'event'}
          );
        }
      ",

      click = "
        function(properties) {
          var network = this;
          var nodeId = properties.nodes && properties.nodes.length
            ? String(properties.nodes[0])
            : '';

          if (!nodeId.startsWith('sink_')) {
            return;
          }

          var selected = (
            network.__cegSelectedEdges || []
          ).slice();

          if (selected.length === 0) {
            Shiny.setInputValue(
              'ceg_sink_assignment',
              {sink: nodeId, edges: []},
              {priority: 'event'}
            );
            return;
          }

          Shiny.setInputValue(
            'ceg_sink_assignment',
            {sink: nodeId, edges: selected},
            {priority: 'event'}
          );

          network.__cegSelectedEdges = [];
          network.unselectAll();
        }
      ",

      dragEnd = "
        function(properties) {
          var network = this;

          if (properties.nodes && properties.nodes.length > 0) {
            Shiny.setInputValue(
              'ceg_node_positions',
              network.getPositions(),
              {priority: 'event'}
            );
          } else {
            Shiny.setInputValue(
              'ceg_view_position',
              network.getViewPosition(),
              {priority: 'event'}
            );
          }
        }
      "
    )

  if(returnnodes) return(list(plot = p,nodes = nodes))

  p
}

#' Look Up CEG Edge IDs from Variable Outcomes
#'
#' Identifies the CEG edges corresponding to specified outcomes of event-tree
#' variables.
#'
#' @param CEGnodes A data frame containing the CEG node information.
#' @param CEGedges A data frame containing the CEG edge information.
#' @param varnames A character vector giving the ordered event-tree variable names.
#' @param states A named list containing the selected outcomes for each variable.
#'
#' @return A vector of CEG edge IDs corresponding to the selected outcomes.
#'
#' @export
CEGedgeLookup <- function(CEGnodes, CEGedges, varnames, states) {

  edge_levels <- CEGnodes$level[
    match(CEGedges$pos_from, CEGnodes$id)
  ] + 1

  selected_edges <- unlist(
    lapply(seq_along(varnames), function(i) {
      if (is.null(states[[varnames[i]]])) {
        return(NULL)
      }

      CEGedges$edge_id[
        edge_levels == i &
          CEGedges$edge %in% states[[varnames[i]]]
      ]
    })
  )

  unique(selected_edges)
}

#' Get the Descendants of a Set of Nodes in the CEG
#'
#' Returns a vector of descendants of an inputted node or set of nodes, including the nodes inputted.
#'
#' The descendants of a node are the nodes that can be reached from the specified node following the edges of the CEG.
#' The function finds all such nodes for a set of inputted nodes.
#'
#' @param CEGinfo A table specifying the nodes of the CEG (see `getCEGnodes`), or the edges of the CEG (see `getCEGedges'), or the list of both (see `getCEGtables')
#' @param nodes A node ID or vector of node IDs
#'
#' @return A vector of node IDs that are descendants of at least one inputted node ID in the CEG
#'
#' @export
CEGdescendants <- function(CEGinfo, nodes){
  if(is.null(CEGinfo)) stop("Please input the CEG node or edge information table.")
  if(is.null(nodes)) stop("Please input at least one node whose descendants are to be returned.")

  if(inherits(CEGinfo, "list")) CEGinfo <- CEGinfo[[1]]
  if("parents" %in% names(CEGinfo)){
    CEGnodes <- CEGinfo
    if(all(nodes %in% CEGnodes$id)==FALSE){
      stop("There are node IDs that don't exist!")
    }

    children <- CEGnodes$id[sapply(CEGnodes$parents, function(p) any(p %in% nodes))]

    if(length(children) == 0){
      return(nodes)
    }

    desc <- nodes
    for(child in children){
      desc <- unique(c(desc, CEGdescendants(CEGinfo, child)))
    }
  } else if("pos_from" %in% names(CEGinfo)){
    CEGedges <- CEGinfo
    if(all(nodes %in% unique(c(CEGedges$pos_from, CEGedges$pos_to)))==FALSE){
      stop("There are node IDs that don't exist!")
    }

    children <- unique(CEGedges$pos_to[CEGedges$pos_from %in% nodes])

    if(length(children) == 0){
      return(nodes)
    }

    desc <- nodes
    for(child in children){
      desc <- unique(c(desc, CEGdescendants(CEGinfo, child)))
    }
  } else stop("The format of the inputted CEG information is not compatible with this function. Ensure you input CEGnodes, CEGedges or CEGtables as outputted from one of the respective CEG functions.")


  sort(desc)
}

#' Generate All Paths for a CEG
#'
#' Finds all root-to-sink paths (or paths with manually inputted start and sink node) for an inputted CEG.
#'
#' A path is an ordered set of nodes or edges from a root node to a terminal node following the edges of the graph. This function produces a list
#' of these paths for an inputted CEG by node and then by edge. A CEG uses a multigraph structure, enabling multiple edges to exist between a pair of
#' nodes. This means that an ordered set of nodes does not necessarily uniquely specify a path through the CEG. That is to say that a unique path
#' through the CEG must often be specified by the edges traversed. However, for some functions, such as finding descendants, paths specified by
#' node is sufficient. The list of paths for a CEG by node features duplicate entries when a `node path' corresponds to multiple `edge paths'.
#' The outputted list of node paths corresponds directly to the outputted list of edge paths, i.e. the node path in a given position in the node path
#' list is the node path obtained when following the edge path in the corresponding position of the edge path list.
#'
#' @param CEGedges A table specifying the edges of the CEG (see `getCEGedges`)
#' @param start A node ID
#' @param sink A node ID downstream of the start node
#'
#' @return A list of paths in the CEG from the start node to the specified sink node by nodes traversed and then by edges traversed
#' @export
CEGpaths <- function(CEGedges, start = 0, sink = Inf){
  if(inherits(CEGedges, "list")){
    CEGedges <- CEGedges[[which(sapply(CEGedges, function(x) "edge" %in% names(x)))]]
  } else if("parents" %in% names(CEGedges)){
    warning("Only CEG node information was supplied. Edge paths cannot be calculated and the outputted paths may correspond to multiple paths in the graph where there are parallel edges between nodes.")

    pathsFrom <- function(node){

      if(node == sink){
        return(list(node))
      }

      children <- CEGedges$id[sapply(CEGedges$parents, function(p) node %in% p)]

      if(length(children) == 0){
        return(list())
      }

      paths <- lapply(children, function(child){
        continuations <- pathsFrom(child)

        lapply(continuations, function(path){
          c(node, path)
        })
      })

      unlist(paths, recursive = FALSE)
    }

    paths <- pathsFrom(start)

    return(list(
      nodes = paths,
      edges = NULL
    ))

  } else if(!"pos_from" %in% names(CEGedges)){
    stop("The input must be a CEGnodes table, CEGedges table or CEGtables list.")
  }

  pathsFrom <- function(node){

    if(node == sink){
      return(list(list(nodes = node, edges = character(0))))
    }

    children <- which(CEGedges$pos_from == node)

    if(length(children) == 0){
      return(list())
    }

    paths <- lapply(children, function(i){
      continuations <- pathsFrom(CEGedges$pos_to[i])

      lapply(continuations, function(path){
        list(
          nodes = c(node, path$nodes),
          edges = c(CEGedges$edge[i], path$edges)
        )
      })
    })

    unlist(paths, recursive = FALSE)
  }

  paths <- pathsFrom(start)

  list(
    nodes = lapply(paths, `[[`, "nodes"),
    edges = lapply(paths, `[[`, "edges")
  )
}

#' Calculate Probability of Arriving at a Node
#'
#' Calculates the probability that a root-to-leaf path will pass through a given node.
#'
#' This is the product of the edge probabilities along the path between the root and the specified node.
#'
#' @inheritParams bypassNode
#'
#' @return Probability that a root-to-leaf path passes through a given node.
#' @export
getArrivalProb <- function(infoTable,node){
  prod(infoTable$edgeprob[match(ancestors(infoTable, node)[-1], infoTable$id)])
}

#' Calculate Root-to-leaf Path Probabilities
#'
#' Calculates the probability of occurrence for each root-to-leaf path.
#'
#' For a path, this is the product of the edge probabilities along the path.
#'
#' @inheritParams plotEventTree
#'
#' @return A table showing the leaf node for each path and its corresponding path probability
#' @export
getPathProbs <- function(infoTable){
  leaves <- infoTable$id[!(infoTable$id %in% infoTable$parent)]
  leaves <- leaves[order(leaves)]
  pathprobs <- numeric()
  for(i in leaves){
    pathprobs <- c(pathprobs,getArrivalProb(infoTable,i))
  }
  cbind(leaf=leaves,pathprob=pathprobs)
}

#' Calculate the Probability of Arrival for a Set of Leaf Nodes
#'
#' Calculates the probability that a root-to-leaf path will arrive at one of the specified leaf nodes.
#'
#' This is calculated by taking a sum over the product of edge probabilities along each path that arrives at one of the leaf nodes specified.
#'
#' @inheritParams plotEventTree
#' @param leaves A leaf node ID or vector of leaf node IDs
#'
#' @return The probability that a root-to-leaf path arriving at one of the inputted leaf nodes will be taken
#' @export
getPathSetProbs <- function(infoTable,leaves){
  if(!(all(leaves %in% infoTable$id[!(infoTable$id %in% infoTable$parent)]))){
    warning("There are node IDs that don't correspond to leaves in the tree! These nodes have been ignored on behalf of your own peril.")
  }
  pathprobs <- getPathProbs(infoTable)
  sum(pathprobs[pathprobs[,1] %in% leaves,2])
}


#' Add Evidence to the CEG
#'
#' Updates the CEG structure and edge probabilities through the CEG propagation algorithm using inputted evidence.
#'
#' Positive evidence is evidence in favour of a node or edge. Where a level contains at least one positive node, any root-to-sink
#' path must pass through one of these positive nodes. Where a floret (i.e. the edges emanating from a node) contains at least one
#' positive edge, any root-to-sink path arriving at the parent node of that floret must pass through one of the positive emanating
#' edges of the floret. Negative evidence is evidence against a node or edge. No root-to-sink path may pass through a negative node
#' or edge. Positive evidence enables a root-to-sink path to pass through, whereas negative evidence prevents any root-to-sink path
#' from passing through. There is a correspondence between positive and negative evidence; if there are four positions (nodes) in a
#' level and two of them are specified as positive, this is equivalent to specifying the other two as negative.
#' If the inputted evidence is negative, the inputted nodes and/or edges are removed from the CEG, as well as any other nodes
#' or edges that must pass through at least one negative node or edge as inputted. If the inputted evidence is positive, any
#' node in the same level as at least one positive node but itself is not specified as a positive node will be deleted. Any
#' edge in the same floret as that of at least one positive edge (i.e. sharing a parent node) but itself is not specified as
#' a positive edge will be deleted. This forms the subCEG structure corresponding to the evidence. A propagation algorithm
#' (Thwaites et al., 2008) is implemented on this subCEG structure in order to update the edge probabilities in line with this
#' inputted evidence. The updated CEG information is outputted, including both updated structural and probabilistic information.
#'
#' @param CEGnodes A table specifying the nodes of the CEG (see `getCEGnodes`)
#' @param CEGedges A table specifying the edges of the CEG (see `getCEGedges`)
#' @param positive A Boolean input identifying whether the evidence is positive or negative
#' @param nodes A node ID or vector of node IDs
#' @param edges An edge ID or vector of edge IDs
#'
#' @return A table specifying the nodes of the updated CEG, and a table specifying the edges of the updated CEG.
#'
#' @export
CEGevidence <- function(CEGnodes,CEGedges,positive=TRUE,nodes=NULL,edges=NULL){

  if(!all(nodes %in% CEGnodes$id)) stop("You have inputted node IDs that do not exist!")
  if(!all(edges %in% CEGedges$edge_id)) stop("You have inputted edge IDs that do not exist!")
  CEGedges$pos_from <- as.numeric(as.character(CEGedges$pos_from))
  CEGedges$pos_to <- as.numeric(as.character(CEGedges$pos_to))
  paths <- CEGpaths(CEGedges)$nodes

  if(!is.null(nodes)){
    nodes <- sort(nodes)
    levels <- CEGnodes$level[which(CEGnodes$id%in%nodes)]
    if(positive==TRUE){
    inc <- nodes
    exc <- setdiff(CEGnodes$id[which(CEGnodes$level%in%levels)],nodes)
    } else{
      exc <- nodes
      inc <- setdiff(CEGnodes$id[which(CEGnodes$level%in%levels)],nodes)
    }
  incpaths <- paths[sapply(paths, function(p) !any(p %in% exc))]
  incnodes <- sort(unique(unlist(incpaths)))
  if(length(incnodes)==0) stop("The inputted evidence is invalid. There are no root-to-leaf paths that adhere to the positively evidenced nodes. All paths must pass through at least one negatively evidenced node or edge.")

  newnodes <- CEGnodes[which(CEGnodes$id%in%incnodes),]
  newedges <- CEGedges[which(CEGedges$pos_from %in% incnodes & CEGedges$pos_to %in% incnodes),]
  newnodes$parents <- I(lapply(newnodes$id, function(id) {unique(newedges$pos_from[newedges$pos_to == id])}))

  } else{
    newnodes <- CEGnodes
    newedges <- CEGedges
  }

  if(!is.null(edges)){
    florets <- unique(newedges$pos_from[which(newedges$edge_id %in% edges)])  # relevant florets for edge evidence given they are not rooted at negative nodes (in which case the edge evidence is irrelevant)
    if(positive==TRUE){
      inc <- edges[which(edges %in% newedges$edge_id)]            # if a floret has positive edges, they must be taken
      exc <- setdiff(newedges$edge_id[which(newedges$pos_from %in% florets)],edges)
    } else{
      exc <- edges[which(edges %in% newedges$edge_id)]           # if a floret has negative edges, they must not be taken
      inc <- setdiff(newedges$edge_id[which(newedges$pos_from %in% florets)],edges)
    }
    newedges <- newedges[!(newedges$edge_id %in% exc),]
  }

  incpaths <- CEGpaths(newedges)$nodes
  incnodes <- sort(unique(unlist(incpaths)))
  if(length(incnodes)==0) stop("The inputted evidence is invalid. There are no root-to-leaf paths that adhere to the positively evidenced nodes. All paths must pass through at least one negatively evidenced node or edge.")

  newnodes <- newnodes[newnodes$id %in% incnodes,]
  newedges <- newedges[newedges$pos_from %in% incnodes & newedges$pos_to %in% incnodes,]
  if(length(CEGpaths(newedges)$nodes)==0) stop("The inputted evidence is invalid. There are no root-to-leaf paths that adhere to the positively evidenced nodes. All paths must pass through at least one negatively evidenced node or edge.")
  newnodes[newnodes=="wInf"] <- "winf"
  newnodes$parents <- I(lapply(newnodes$id, function(id) {unique(newedges$pos_from[newedges$pos_to == id])}))
  evidence <- incpaths

  newnodes$emphasis <- rep(1,nrow(newnodes))
  newnodes$actuated <- rep(0,nrow(newnodes))
  newedges$potential <- newedges$edgeprob

  # backwards collection step
      # actuate sink node(s) to 1
  newnodes$actuated[sapply(newnodes$id, function(id) !any(sapply(newnodes$parents, function(p) id %in% p)))] <- 1

      # loop until all nodes are actuated
  while(any(newnodes$actuated==0)){
      # find the nodes all of whose children are actuated but themselves aren't
  nextnodes <- newnodes$id[sapply(newnodes$id, function(id){
    children <- which(sapply(newnodes$parents, function(p) id %in% p)); length(children) > 0 && all(newnodes$actuated[children] == 1)
    })]
  nextnodes <- newnodes$id[which(newnodes$id %in% nextnodes & newnodes$actuated==0)]
  if(length(nextnodes)==0) stop("Unable to actuate all nodes. Check that there are no disconnected components of your CEG.")
      # calculate edge potentials
  newedges$potential[newedges$pos_from %in% nextnodes] <- newedges$edgeprob[newedges$pos_from %in% nextnodes] * newnodes$emphasis[match(newedges$pos_to[newedges$pos_from %in% nextnodes], newnodes$id)]

      # calculate new emphases for newly actuated edges
  newnodes$emphasis[which(newnodes$id %in% nextnodes)] <-  sapply(nextnodes, function(id) sum(newedges$potential[newedges$pos_from == id]))
  newnodes$actuated[which(newnodes$id %in% nextnodes)] <- 1
  }

  # forward distribution step
  for(i in seq_len(nrow(newedges))){
    newedges$edgeprob[i] <- newedges$potential[i]/newnodes$emphasis[match(newedges$pos_from[i],newnodes$id)]
  }

  newnodes <- newnodes[, !names(newnodes) %in% c("emphasis", "actuated")]
  newedges <- newedges[, !names(newedges) %in% "potential"]

  list(CEGnodes = newnodes, CEGedges = newedges)
}
