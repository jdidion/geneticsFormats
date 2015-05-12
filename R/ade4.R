read.tree.ade4 <- function(newick.file) {
    return(newick2phylog(readLines(newick.file)))
}

# Convert a phylog object to a binary tree, which is simply a list of lists, where each list
# has at most two elements.
as.binary.tree <- function(node, edges=NA, leaf.names=NULL) {
    if (is.null(node)) {
        return(NULL)
    }
    else if (inherits(node, "phylog")) {
        return(as.binary.tree(node$parts$Root, node$parts, leaf.names))
    }
    
    l <- list()
    for (child in node) {
        name <- leaf.names[[child]]
        if (is.null(name)) {
            name <- child
        }
        if (child %in% names(edges)) {
            l[[name]] <- as.binary.tree(edges[[child]], edges, leaf.names)
        }
        else {
            l[[name]] <- NA
        }
    }
    return(l)
}

# Find the nearest common ancestor of all specified leaf nodes.
find.ancestor <- function(p, leaves) {
    path <- reduce(p$paths[leaves], intersect)
    return(path[length(path)])
}

# Returns the subset of nodes that are decended from root.
find.descendants <- function(p, root) {
    cur <- root
    rows <- NULL
    while (length(cur) > 0) {
        rows <- c(rows, cur)
        kids <- unlist(p$parts[cur])
        cur <- kids[substr(kids, 1, 1) == 'I']
    }
    return(p$parts[rows])
}

# Returns information on the minimal subtree of nodes containing all leaf nodes
# in leaves. Return value is a list: root=name of root node, tree=subset of nodes matrix that
# forms the subtree, and cousins=leaf nodes within the subtree that are not in leaves.
find.minimal.subtree <- function(p, leaves) {
    root <- find.ancestor(p, leaves)
    tree <- find.descendants(p, root)
    temp <- unlist(tree)
    cous <- setdiff(temp[substr(temp, 1, 1) == 'X'], leaves)
    return(list(root=root, tree=tree, cousins=cous))
}

# Return TRUE if leaves are the only leaf nodes in the subtree consisting
# of all descendants of the nearest common ancestor of leaves.
is.minimal.subtree <- function(p, leaves) {
    root <- find.ancestor(p, leaves)
    desc <- unlist(find.descendants(p, root))
    return(setequal(leaves, desc[substr(desc, 1, 1) == 'X']))
}
