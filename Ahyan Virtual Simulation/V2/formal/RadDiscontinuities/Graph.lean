import RadDiscontinuities.Operators

namespace RadDiscontinuities

/-!
Graph abstraction for cell topology.

Cell removal is modeled as deletion from a present-cell predicate plus deletion
of every incident edge. This is the graph-theoretic counterpart of physically
removing a RAD cell or disabling a module in another programmable mechanical
architecture.
-/

structure CellGraph (Cell : Type u) where
  present : Cell -> Bool
  adjacent : Cell -> Cell -> Bool

namespace CellGraph

def removeCell [DecidableEq Cell]
    (cell : Cell) (graph : CellGraph Cell) : CellGraph Cell where
  present := fun other =>
    if other = cell then false else graph.present other
  adjacent := fun first second =>
    if first = cell then false
    else if second = cell then false
    else graph.adjacent first second

theorem removed_cell_not_present [DecidableEq Cell]
    (cell : Cell) (graph : CellGraph Cell) :
    (removeCell cell graph).present cell = false := by
  simp [removeCell]

theorem removal_preserves_other_presence [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell)
    (hdifferent : Not (other = cell)) :
    (removeCell cell graph).present other = graph.present other := by
  simp [removeCell, hdifferent]

theorem removal_deletes_incident_edges_left [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    (removeCell cell graph).adjacent cell other = false := by
  simp [removeCell]

theorem removal_deletes_incident_edges_right [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    (removeCell cell graph).adjacent other cell = false := by
  by_cases hsame : other = cell
  case pos =>
    simp [removeCell, hsame]
  case neg =>
    simp [removeCell, hsame]

def oneStepReachable (graph : CellGraph Cell) (source target : Cell) : Bool :=
  graph.present source && graph.present target && graph.adjacent source target

theorem removal_deletes_one_step_from_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    oneStepReachable (removeCell cell graph) cell other = false := by
  simp [oneStepReachable, removeCell]

theorem removal_deletes_one_step_to_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    oneStepReachable (removeCell cell graph) other cell = false := by
  by_cases hsame : other = cell
  case pos =>
    simp [oneStepReachable, removeCell, hsame]
  case neg =>
    simp [oneStepReachable, removeCell, hsame]

theorem removal_preserves_nonincident_one_step [DecidableEq Cell]
    (removed source target : Cell) (graph : CellGraph Cell)
    (hsource : Not (source = removed))
    (htarget : Not (target = removed)) :
    oneStepReachable (removeCell removed graph) source target =
      oneStepReachable graph source target := by
  simp [oneStepReachable, removeCell, hsource, htarget]

inductive PositivePath (graph : CellGraph Cell) : Cell -> Cell -> Prop where
  | edge {source target : Cell} :
      oneStepReachable graph source target = true ->
      PositivePath graph source target
  | tail {source middle target : Cell} :
      oneStepReachable graph source middle = true ->
      PositivePath graph middle target ->
      PositivePath graph source target

theorem positive_path_source_present
    {graph : CellGraph Cell} {source target : Cell}
    (path : PositivePath graph source target) :
    graph.present source = true := by
  cases path with
  | edge hstep =>
      simp [oneStepReachable] at hstep
      exact hstep.1.1
  | tail hstep _ =>
      simp [oneStepReachable] at hstep
      exact hstep.1.1

theorem positive_path_target_present
    {graph : CellGraph Cell} {source target : Cell}
    (path : PositivePath graph source target) :
    graph.present target = true := by
  induction path with
  | edge hstep =>
      simp [oneStepReachable] at hstep
      exact hstep.1.2
  | tail _ _ ih =>
      exact ih

theorem removal_deletes_positive_path_from_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    Not (PositivePath (removeCell cell graph) cell other) := by
  intro path
  have hpresent := positive_path_source_present path
  simp [removeCell] at hpresent

theorem removal_deletes_positive_path_to_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    Not (PositivePath (removeCell cell graph) other cell) := by
  intro path
  have hpresent := positive_path_target_present path
  simp [removeCell] at hpresent

inductive Reachable (graph : CellGraph Cell) : Cell -> Cell -> Prop where
  | refl {cell : Cell} :
      graph.present cell = true ->
      Reachable graph cell cell
  | positive {source target : Cell} :
      PositivePath graph source target ->
      Reachable graph source target

theorem reachable_source_present
    {graph : CellGraph Cell} {source target : Cell}
    (path : Reachable graph source target) :
    graph.present source = true := by
  cases path with
  | refl hpresent =>
      exact hpresent
  | positive positivePath =>
      exact positive_path_source_present positivePath

theorem reachable_target_present
    {graph : CellGraph Cell} {source target : Cell}
    (path : Reachable graph source target) :
    graph.present target = true := by
  cases path with
  | refl hpresent =>
      exact hpresent
  | positive positivePath =>
      exact positive_path_target_present positivePath

theorem removal_deletes_reachable_from_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    Not (Reachable (removeCell cell graph) cell other) := by
  intro path
  have hpresent := reachable_source_present path
  simp [removeCell] at hpresent

theorem removal_deletes_reachable_to_removed [DecidableEq Cell]
    (cell other : Cell) (graph : CellGraph Cell) :
    Not (Reachable (removeCell cell graph) other cell) := by
  intro path
  have hpresent := reachable_target_present path
  simp [removeCell] at hpresent

structure CellConstraintMap (Cell : Type u) (Constraint : Type v) where
  touches : Constraint -> Cell -> Bool

def removeCellConstraints [DecidableEq Cell]
    (cell : Cell)
    (map : CellConstraintMap Cell Constraint)
    (mode : Mode Constraint) : Mode Constraint :=
  fun constraint =>
    if map.touches constraint cell then false else mode constraint

theorem removed_cell_constraints_inactive [DecidableEq Cell]
    (cell : Cell)
    (map : CellConstraintMap Cell Constraint)
    (mode : Mode Constraint)
    (constraint : Constraint)
    (htouches : map.touches constraint cell = true) :
    removeCellConstraints cell map mode constraint = false := by
  simp [removeCellConstraints, htouches]

theorem removed_cell_constraints_clear_group_operator
    [DecidableEq Cell] [DecidableEq Constraint]
    (cell : Cell)
    (map : CellConstraintMap Cell Constraint)
    (constraints : List Constraint)
    (enabledValue : Bool)
    (mode : Mode Constraint)
    (constraint : Constraint)
    (htouches : map.touches constraint cell = true) :
    removeCellConstraints cell map
        ((LocalModeOperator.groupOperator constraints enabledValue).apply mode)
        constraint = false := by
  simp [removeCellConstraints, htouches]

theorem removed_cell_constraints_not_active_after_group_operator
    [DecidableEq Cell] [DecidableEq Constraint]
    (cell : Cell)
    (map : CellConstraintMap Cell Constraint)
    (constraints : List Constraint)
    (enabledValue : Bool)
    (mode : Mode Constraint)
    (constraint : Constraint)
    (htouches : map.touches constraint cell = true) :
    Not
      (active
        (removeCellConstraints cell map
          ((LocalModeOperator.groupOperator constraints enabledValue).apply mode))
        constraint) := by
  intro hactive
  simp [active, removeCellConstraints, htouches] at hactive

end CellGraph

end RadDiscontinuities
