import RadDiscontinuities.Mode

namespace RadDiscontinuities

variable {Constraint : Type u}

/-!
Mode-level discontinuity operators.

`LocalModeOperator` is a deliberately conservative first abstraction. It
overwrites active/inactive status only on a declared support and leaves every
constraint outside that support unchanged. This is strong enough to prove the
locality and commutation facts needed for the first RAD formal layer without
claiming a full mechanics theorem.
-/

structure LocalModeOperator (Constraint : Type u) where
  support : Constraint -> Bool
  enabled : Constraint -> Bool

namespace LocalModeOperator

def supportContainedIn (constraints : List Constraint)
    (operator : LocalModeOperator Constraint) : Prop :=
  forall constraint,
    operator.support constraint = true -> Membership.mem constraints constraint

theorem finite_operator_support (constraints : List Constraint)
    (operator : LocalModeOperator Constraint)
    (hsupported : supportContainedIn constraints operator) :
    Exists
      (fun finiteSupport : List Constraint =>
        forall constraint,
          operator.support constraint = true ->
            Membership.mem finiteSupport constraint) := by
  exact Exists.intro constraints hsupported

def apply (operator : LocalModeOperator Constraint)
    (mode : Mode Constraint) : Mode Constraint :=
  fun constraint =>
    if operator.support constraint then operator.enabled constraint
    else mode constraint

theorem apply_outside_support (operator : LocalModeOperator Constraint)
    (mode : Mode Constraint) (constraint : Constraint)
    (houtside : operator.support constraint = false) :
    operator.apply mode constraint = mode constraint := by
  simp [apply, houtside]

theorem apply_inside_support (operator : LocalModeOperator Constraint)
    (mode : Mode Constraint) (constraint : Constraint)
    (hinside : operator.support constraint = true) :
    operator.apply mode constraint = operator.enabled constraint := by
  simp [apply, hinside]

theorem disjoint_support_operators_commute
    (left right : LocalModeOperator Constraint) (mode : Mode Constraint)
    (hdisjoint :
      forall constraint,
        left.support constraint = true -> right.support constraint = false) :
    left.apply (right.apply mode) = right.apply (left.apply mode) := by
  funext constraint
  cases hleft : left.support constraint <;> cases hright : right.support constraint
  case false.false =>
    simp [apply, hleft, hright]
  case false.true =>
    simp [apply, hleft, hright]
  case true.false =>
    simp [apply, hleft, hright]
  case true.true =>
    have impossible := hdisjoint constraint hleft
    rw [hright] at impossible
    contradiction

def lockOperator [DecidableEq Constraint]
    (constraint : Constraint) : LocalModeOperator Constraint where
  support := fun other => decide (other = constraint)
  enabled := fun _ => true

def unlockOperator [DecidableEq Constraint]
    (constraint : Constraint) : LocalModeOperator Constraint where
  support := fun other => decide (other = constraint)
  enabled := fun _ => false

def groupOperator [DecidableEq Constraint]
    (constraints : List Constraint) (enabledValue : Bool) :
    LocalModeOperator Constraint where
  support := fun other => decide (Membership.mem constraints other)
  enabled := fun _ => enabledValue

def groupLockOperator [DecidableEq Constraint]
    (constraints : List Constraint) : LocalModeOperator Constraint :=
  groupOperator constraints true

def groupUnlockOperator [DecidableEq Constraint]
    (constraints : List Constraint) : LocalModeOperator Constraint :=
  groupOperator constraints false

theorem groupOperator_apply_inside [DecidableEq Constraint]
    (constraints : List Constraint) (enabledValue : Bool)
    (mode : Mode Constraint) (constraint : Constraint)
    (hinside : Membership.mem constraints constraint) :
    (groupOperator constraints enabledValue).apply mode constraint = enabledValue := by
  simp [groupOperator, apply, hinside]

theorem groupOperator_apply_outside [DecidableEq Constraint]
    (constraints : List Constraint) (enabledValue : Bool)
    (mode : Mode Constraint) (constraint : Constraint)
    (houtside : Not (Membership.mem constraints constraint)) :
    (groupOperator constraints enabledValue).apply mode constraint = mode constraint := by
  simp [groupOperator, apply, houtside]

theorem groupOperator_support_inside [DecidableEq Constraint]
    (constraints : List Constraint) (enabledValue : Bool)
    (constraint : Constraint)
    (hinside : Membership.mem constraints constraint) :
    (groupOperator constraints enabledValue).support constraint = true := by
  simp [groupOperator, hinside]

theorem groupOperator_support_outside [DecidableEq Constraint]
    (constraints : List Constraint) (enabledValue : Bool)
    (constraint : Constraint)
    (houtside : Not (Membership.mem constraints constraint)) :
    (groupOperator constraints enabledValue).support constraint = false := by
  simp [groupOperator, houtside]

theorem groupOperator_empty_apply [DecidableEq Constraint]
    (enabledValue : Bool) (mode : Mode Constraint) :
    (groupOperator ([] : List Constraint) enabledValue).apply mode = mode := by
  funext constraint
  simp [groupOperator, apply]

theorem groupOperator_append_support_left [DecidableEq Constraint]
    (left right : List Constraint) (enabledValue : Bool)
    (constraint : Constraint)
    (hleft : Membership.mem left constraint) :
    (groupOperator (left ++ right) enabledValue).support constraint = true := by
  have happend : Membership.mem (left ++ right) constraint := by
    exact List.mem_append_left right hleft
  simpa [groupOperator, List.mem_append] using happend

theorem groupOperator_append_support_right [DecidableEq Constraint]
    (left right : List Constraint) (enabledValue : Bool)
    (constraint : Constraint)
    (hright : Membership.mem right constraint) :
    (groupOperator (left ++ right) enabledValue).support constraint = true := by
  have happend : Membership.mem (left ++ right) constraint := by
    exact List.mem_append_right left hright
  simpa [groupOperator, List.mem_append] using happend

theorem groupOperators_with_disjoint_lists_commute [DecidableEq Constraint]
    (left right : List Constraint)
    (leftEnabled rightEnabled : Bool)
    (mode : Mode Constraint)
    (hdisjoint :
      forall constraint,
        Membership.mem left constraint ->
          Not (Membership.mem right constraint)) :
    (groupOperator left leftEnabled).apply
        ((groupOperator right rightEnabled).apply mode) =
      (groupOperator right rightEnabled).apply
        ((groupOperator left leftEnabled).apply mode) := by
  apply disjoint_support_operators_commute
  intro constraint hleftSupport
  have hleftMem : Membership.mem left constraint := by
    simpa [groupOperator] using hleftSupport
  have hrightNot := hdisjoint constraint hleftMem
  by_cases hright : Membership.mem right constraint
  case pos =>
    exact False.elim (hrightNot hright)
  case neg =>
    simp [groupOperator, hright]

theorem lockOperator_apply [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) :
    (lockOperator constraint).apply mode =
      MechanicalSystem.lock constraint mode := by
  funext other
  by_cases h : other = constraint
  case pos =>
    simp [lockOperator, apply, MechanicalSystem.lock, h]
  case neg =>
    simp [lockOperator, apply, MechanicalSystem.lock, h]

theorem unlockOperator_apply [DecidableEq Constraint]
    (constraint : Constraint) (mode : Mode Constraint) :
    (unlockOperator constraint).apply mode =
      MechanicalSystem.unlock constraint mode := by
  funext other
  by_cases h : other = constraint
  case pos =>
    simp [unlockOperator, apply, MechanicalSystem.unlock, h]
  case neg =>
    simp [unlockOperator, apply, MechanicalSystem.unlock, h]

end LocalModeOperator

end RadDiscontinuities
