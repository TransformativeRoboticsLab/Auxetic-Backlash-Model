import RadDiscontinuities.EventAlgebra
import RadDiscontinuities.Graph
import RadDiscontinuities.Mechanics

namespace RadDiscontinuities

/-!
Finite examples for the first RAD abstraction layer.

These examples do not claim fabrication-accurate RAD geometry. They provide
small compiled witnesses that one-cell, two-cell, and small-lattice abstractions
fit the mode/operator/graph vocabulary.
-/

def oneCell : Fin 1 := { val := 0, isLt := by decide }
def twoCellLeft : Fin 2 := { val := 0, isLt := by decide }
def twoCellRight : Fin 2 := { val := 1, isLt := by decide }
def smallLatticeFirst : Fin 4 := { val := 0, isLt := by decide }

def emptyOneCellMode : FiniteMode 1 := fun _ => false

theorem one_cell_lock_idempotent_example :
    MechanicalSystem.lock oneCell
      (MechanicalSystem.lock oneCell emptyOneCellMode) =
        MechanicalSystem.lock oneCell emptyOneCellMode :=
  MechanicalSystem.lock_idempotent oneCell emptyOneCellMode

def twoCellGraph : CellGraph (Fin 2) where
  present := fun _ => true
  adjacent := fun first second =>
    decide
      ((first = twoCellLeft /\ second = twoCellRight) \/
        (first = twoCellRight /\ second = twoCellLeft))

theorem two_cell_removal_deletes_left_edge :
    (CellGraph.removeCell twoCellLeft twoCellGraph).adjacent
      twoCellLeft twoCellRight = false :=
  CellGraph.removal_deletes_incident_edges_left
    twoCellLeft twoCellRight twoCellGraph

def smallLatticeMode : FiniteMode 4 := fun _ => true

theorem small_lattice_unlock_idempotent_example :
    MechanicalSystem.unlock smallLatticeFirst
      (MechanicalSystem.unlock smallLatticeFirst smallLatticeMode) =
        MechanicalSystem.unlock smallLatticeFirst smallLatticeMode :=
  MechanicalSystem.unlock_idempotent smallLatticeFirst smallLatticeMode

end RadDiscontinuities
