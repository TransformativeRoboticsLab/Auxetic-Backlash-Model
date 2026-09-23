import Std

namespace RadDiscontinuities

/-!
Backlash dead-zone map.

This is the scalar law shared by the RAD preprint, the Python simulator, and the
browser implementation:

  max(0, x - b) + min(x + b, 0)

The first Lean target is deliberately local: prove exactly the dead-zone and
outside-gap branch behavior before formalizing lattice propagation. This
self-contained scaffold proves the ordered integer version so the package can
typecheck without mathlib; a real-valued version can be added later once the
mathlib cache is available.
-/

noncomputable def backlash (b x : Int) : Int :=
  max 0 (x - b) + min (x + b) 0

theorem backlash_zero_inside {b x : Int}
    (hupper : x <= b) (hlower : -b <= x) :
    backlash b x = 0 := by
  unfold backlash
  have hmax : x - b <= 0 := by omega
  have hmin : 0 <= x + b := by omega
  rw [Int.max_eq_left hmax, Int.min_eq_right hmin]
  omega

theorem backlash_positive_branch {b x : Int}
    (hb : 0 <= b) (houtside : b <= x) :
    backlash b x = x - b := by
  unfold backlash
  have hmax : 0 <= x - b := by omega
  have hmin : 0 <= x + b := by omega
  rw [Int.max_eq_right hmax, Int.min_eq_right hmin]
  omega

theorem backlash_negative_branch {b x : Int}
    (hb : 0 <= b) (houtside : x <= -b) :
    backlash b x = x + b := by
  unfold backlash
  have hmax : x - b <= 0 := by omega
  have hmin : x + b <= 0 := by omega
  rw [Int.max_eq_left hmax, Int.min_eq_left hmin]
  omega

structure RealValuedBacklashTargets (Scalar : Type u)
    [OfNat Scalar 0] [LE Scalar] [Add Scalar] [Sub Scalar] [Neg Scalar] where
  backlash : Scalar -> Scalar -> Scalar
  zero_inside :
    forall b x, x <= b -> -b <= x -> backlash b x = 0
  positive_branch :
    forall b x, 0 <= b -> b <= x -> backlash b x = x - b
  negative_branch :
    forall b x, 0 <= b -> x <= -b -> backlash b x = x + b

end RadDiscontinuities
