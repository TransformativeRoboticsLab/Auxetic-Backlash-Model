import RadDiscontinuities.Operators

namespace RadDiscontinuities

/-!
Algebra of mode-level events.

This file treats a programmable discontinuity as a mode-to-mode map. It proves
identity, composition, associativity, and finite sequence application facts.
These are algebraic statements about event order, independent of any specific
RAD geometry or solver.
-/

abbrev ModeOperator (Constraint : Type u) := Mode Constraint -> Mode Constraint

def identityOperator {Constraint : Type u} : ModeOperator Constraint :=
  fun mode => mode

def composeOperators {Constraint : Type u}
    (second first : ModeOperator Constraint) : ModeOperator Constraint :=
  fun mode => second (first mode)

theorem identity_compose_left {Constraint : Type u}
    (operator : ModeOperator Constraint) :
    composeOperators identityOperator operator = operator := by
  rfl

theorem identity_compose_right {Constraint : Type u}
    (operator : ModeOperator Constraint) :
    composeOperators operator identityOperator = operator := by
  rfl

theorem composeOperators_assoc {Constraint : Type u}
    (third second first : ModeOperator Constraint) :
    composeOperators third (composeOperators second first) =
      composeOperators (composeOperators third second) first := by
  rfl

def applyEventSequence {Constraint : Type u} :
    List (ModeOperator Constraint) -> ModeOperator Constraint
  | [] => identityOperator
  | operator :: rest => composeOperators (applyEventSequence rest) operator

theorem applyEventSequence_nil {Constraint : Type u} :
    applyEventSequence ([] : List (ModeOperator Constraint)) =
      identityOperator := by
  rfl

theorem applyEventSequence_cons {Constraint : Type u}
    (operator : ModeOperator Constraint)
    (rest : List (ModeOperator Constraint)) :
    applyEventSequence (operator :: rest) =
      composeOperators (applyEventSequence rest) operator := by
  rfl

theorem applyEventSequence_append {Constraint : Type u}
    (left right : List (ModeOperator Constraint)) :
    applyEventSequence (left ++ right) =
      composeOperators (applyEventSequence right) (applyEventSequence left) := by
  induction left with
  | nil =>
      rfl
  | cons operator rest ih =>
      calc
        applyEventSequence ((operator :: rest) ++ right)
            = composeOperators (applyEventSequence (rest ++ right)) operator := by
              rfl
        _ = composeOperators
              (composeOperators (applyEventSequence right) (applyEventSequence rest))
              operator := by
              rw [ih]
        _ = composeOperators
              (applyEventSequence right)
              (composeOperators (applyEventSequence rest) operator) := by
              exact
                (composeOperators_assoc
                  (applyEventSequence right)
                  (applyEventSequence rest)
                  operator).symm
        _ = composeOperators
              (applyEventSequence right)
              (applyEventSequence (operator :: rest)) := by
              rfl

end RadDiscontinuities
