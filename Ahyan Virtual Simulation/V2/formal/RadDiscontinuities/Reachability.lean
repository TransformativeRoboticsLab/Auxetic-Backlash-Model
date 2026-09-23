import RadDiscontinuities.Operators

namespace RadDiscontinuities

variable {Constraint : Type u}

/-!
Reachability in a finite mode graph.

The numerical simulator can search actuation/lock sequences. This Lean layer
only records the discrete mode-transition skeleton: a mode is reachable when a
finite sequence of indexed operators maps the initial mode to it.
-/

inductive Reachable {Index : Type v}
    (step : Index -> Mode Constraint -> Mode Constraint) :
    Mode Constraint -> Mode Constraint -> Prop where
  | refl (mode : Mode Constraint) : Reachable step mode mode
  | tail {start current : Mode Constraint} (index : Index)
      (path : Reachable step start current) :
      Reachable step start (step index current)

theorem reachable_trans {Index : Type v}
    {step : Index -> Mode Constraint -> Mode Constraint}
    {start middle finish : Mode Constraint}
    (first : Reachable step start middle)
    (second : Reachable step middle finish) :
    Reachable step start finish := by
  induction second with
  | refl =>
      exact first
  | tail index path ih =>
      exact Reachable.tail index ih

structure ModeGraph (Constraint : Type u) where
  vertices : List (Mode Constraint)
  edges : List (Prod (Mode Constraint) (Mode Constraint))

namespace ModeGraph

def Vertex (graph : ModeGraph Constraint) (mode : Mode Constraint) : Prop :=
  Membership.mem graph.vertices mode

def Edge (graph : ModeGraph Constraint)
    (edge : Prod (Mode Constraint) (Mode Constraint)) : Prop :=
  Membership.mem graph.edges edge

theorem finite_mode_graph_vertices (graph : ModeGraph Constraint) :
    Exists
      (fun vertices : List (Mode Constraint) =>
        forall mode, graph.Vertex mode -> Membership.mem vertices mode) := by
  exact Exists.intro graph.vertices (by intro mode hmode; exact hmode)

theorem finite_mode_graph_edges (graph : ModeGraph Constraint) :
    Exists
      (fun edges : List (Prod (Mode Constraint) (Mode Constraint)) =>
        forall edge, graph.Edge edge -> Membership.mem edges edge) := by
  exact Exists.intro graph.edges (by intro edge hedge; exact hedge)

end ModeGraph

end RadDiscontinuities
