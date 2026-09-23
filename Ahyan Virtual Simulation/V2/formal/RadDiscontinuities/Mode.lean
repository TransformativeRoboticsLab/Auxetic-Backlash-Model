import Std

namespace RadDiscontinuities

/-!
Finite active-set mechanics.

This file formalizes the first, purely structural layer of programmable
discontinuities. A mode records whether each available constraint is active.
The physics of each constraint is intentionally abstract: a constraint is
represented only by the predicate saying whether a state satisfies it.
-/

abbrev Mode (Constraint : Type u) := Constraint -> Bool

abbrev FiniteMode (numberOfConstraints : Nat) := Mode (Fin numberOfConstraints)

def active {Constraint : Type u} (mode : Mode Constraint)
    (constraint : Constraint) : Prop :=
  mode constraint = true

def supportedBy {Constraint : Type u}
    (constraints : List Constraint) (mode : Mode Constraint) : Prop :=
  forall constraint, active mode constraint -> Membership.mem constraints constraint

theorem finite_active_constraints {Constraint : Type u}
    (constraints : List Constraint) (mode : Mode Constraint)
    (hsupported : supportedBy constraints mode) :
    Exists
      (fun finiteConstraints : List Constraint =>
        forall constraint,
          active mode constraint -> Membership.mem finiteConstraints constraint) := by
  exact Exists.intro constraints hsupported

structure MechanicalSystem (State : Type u) (Constraint : Type v) where
  satisfies : Constraint -> State -> Prop

namespace MechanicalSystem

variable {State : Type u} {Constraint : Type v}

def PredSubset (left right : State -> Prop) : Prop :=
  forall state, left state -> right state

def admissible (system : MechanicalSystem State Constraint)
    (mode : Mode Constraint) : State -> Prop :=
  fun state =>
    forall constraint, active mode constraint -> system.satisfies constraint state

def lock [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) : Mode Constraint :=
  fun other => if other = constraint then true else mode other

def unlock [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) : Mode Constraint :=
  fun other => if other = constraint then false else mode other

theorem lock_idempotent [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) :
    lock constraint (lock constraint mode) = lock constraint mode := by
  funext other
  by_cases h : other = constraint
  case pos =>
    simp [lock, h]
  case neg =>
    simp [lock, h]

theorem unlock_idempotent [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) :
    unlock constraint (unlock constraint mode) = unlock constraint mode := by
  funext other
  by_cases h : other = constraint
  case pos =>
    simp [unlock, h]
  case neg =>
    simp [unlock, h]

theorem admissible_subset_of_mode_subset
    (system : MechanicalSystem State Constraint)
    {smaller larger : Mode Constraint}
    (hsubset :
      forall constraint,
        active smaller constraint -> active larger constraint) :
    PredSubset (admissible system larger) (admissible system smaller) := by
  intro state hstate constraint hconstraint
  exact hstate constraint (hsubset constraint hconstraint)

theorem lock_shrinks_admissible
    [DecidableEq Constraint]
    (system : MechanicalSystem State Constraint)
    (constraint : Constraint) (mode : Mode Constraint) :
    PredSubset
      (admissible system (lock constraint mode))
      (admissible system mode) := by
  intro state hstate other hother
  apply hstate other
  by_cases h : other = constraint
  case pos =>
    subst other
    simp [active, lock]
  case neg =>
    simpa [active, lock, h] using hother

theorem unlock_enlarges_admissible
    [DecidableEq Constraint]
    (system : MechanicalSystem State Constraint)
    (constraint : Constraint) (mode : Mode Constraint) :
    PredSubset
      (admissible system mode)
      (admissible system (unlock constraint mode)) := by
  intro state hstate other hother
  apply hstate other
  by_cases h : other = constraint
  case pos =>
    simp [active, unlock, h] at hother
  case neg =>
    simpa [active, unlock, h] using hother

theorem adding_constraints_shrinks_admissible
    [DecidableEq Constraint]
    (system : MechanicalSystem State Constraint)
    (constraint : Constraint) (mode : Mode Constraint) :
    PredSubset
      (admissible system (lock constraint mode))
      (admissible system mode) :=
  lock_shrinks_admissible system constraint mode

theorem removing_constraints_enlarges_admissible
    [DecidableEq Constraint]
    (system : MechanicalSystem State Constraint)
    (constraint : Constraint) (mode : Mode Constraint) :
    PredSubset
      (admissible system mode)
      (admissible system (unlock constraint mode)) :=
  unlock_enlarges_admissible system constraint mode

end MechanicalSystem

end RadDiscontinuities
