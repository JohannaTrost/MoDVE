#' Find dead branches in a shoot dataset
#'
#' Compares which trees appear in the first time step but are absent from the
#' next one and returns branches belonging to such trees.
#'
#' @param shoot_dt a data frame containing information on shoots, as produced by MoF3D.
#' @param shoot_dt_next the corresponding dataset for the next time step.
#'
#' @returns a vector containing the shoot IDs of dead branches
#' @export
#'
find_dead_segments <- function(shoot_dt, shoot_dt_next) {

  #check_shoot_dt(shoot_dt)
  #check_shoot_dt(shoot_dt_next)

  # Get all branch segments that die during time step
  dead_segments <- shoot_dt$shootID[!is.element(shoot_dt$treeID, shoot_dt_next$treeID)]

  return(dead_segments)
}

#' Find dead trees in a trunk dataset
#'
#' Compares which trees appear in the first time step but are absent from the
#' next one and returns which trunks belong to such trees.
#'
#' @param trunk_dt a data frame containing information on trunks, as produced
#' by MoF3D
#' @param trunk_dt_next the corresponding dataset for the next time step.
#'
#' @returns a vector containing the trunk IDs of dead trees
#' @export
#'
find_dead_trees <- function(trunk_dt, trunk_dt_next) {

  #check_trunk_dt(trunk_dt)
  #check_trunk_dt(trunk_dt_next)

  # Get all branch segments that die during time step
  DeadSegments <- trunk_dt$treeID[!is.element(trunk_dt$treeID, trunk_dt_next$treeID)]

  return(DeadSegments)
}

