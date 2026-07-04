# AeroFlex — Code Review: Physical/Mathematical Consistency and Performance

Scope: full read of `main/` (Airplane.m, Element.m, Node.m, Engine.m, RigidFuselage.m,
dinamicaflex.m, aeroforceandmoment.m, Peter.m, getexpG.m, getJhbp.m, OMEGA.m,
matrixcross.m, atmosfera.m, flutter_speed.m, StructuralModes.m, create_uniform_structure.m,
plotaest3d.m, changedatarate.m) plus the examples that exercise them.

Severity legend: 🔴 likely wrong results · 🟠 inconsistency / silent approximation · 🟡 robustness · 🔵 performance only

---

## 1. Physical / mathematical inconsistencies

### 1.1 🔴 Structural Jacobians: airplane-level matrices frozen while member-level ones update
`Airplane.m`, `updateStrJac` (lines 102–139).

`ap.Jhep` and `ap.Jpep` are reassigned on **every** call, but `ap.Jthetaep`, `ap.Jhb`,
`ap.Jpb`, `ap.Jthetab` are only assigned **once** (guarded by `if isempty(ap.Jthetaep)`,
line 133). Meanwhile the member-level copies (`memberJhb`, `memberJpb`, `memberJthetaep`, …)
*are* refreshed every call and are the ones used by `aeroforceandmoment.m` to compute nodal
velocities.

Consequences when the structure deforms (`updateStrJac = 1` during trim, or `= 2` always):

1. **Virtual-work inconsistency.** The aerodynamic velocities are evaluated with the
   *deformed* Jacobians (`memberJpb`, `memberJthetaep`), but the generalized forces are
   projected with a mix of deformed (`ap.Jpep`) and *initial-geometry* (`ap.Jthetaep`,
   `ap.Jpb`, `ap.Jthetab`, `ap.Jhb`) Jacobians in `dinamicaflex.m` lines 105–106:
   ```matlab
   RAEROF = ap.Jpep'*ap.B*FAERO + ... + ap.Jthetaep'*ap.B*MAERO + ap.Jhep'*ap.N*GRAVITY;
   RAEROB = ap.Jpb'*ap.B*FAERO  + ... + ap.Jthetab'*ap.B*MAERO + ap.Jhb'*ap.N*GRAVITY + ...
   ```
   The distributed moment `MAERO` on a deformed wing is therefore projected onto strain
   coordinates through the *undeformed* `Jthetaep`. For the large deflections this toolbox
   targets (tip deflections of 25–50 % span), the torsion generalized force — the driver of
   both static aeroelastic twist and flutter — is projected with the wrong geometry.
2. **Inertia coupling inconsistency.** `MBB`, `MFB`, `MBF` (lines 45–50 of
   `dinamicaflex.m`) use the frozen `ap.Jhb`, while `MFF` uses the updated `ap.Jhep`.
3. **Wasted computation** (see §3.4): the expensive member-level recomputation is done
   every call and then discarded at the `ap` level.

**Fix.** Either update everything consistently:
```matlab
ap.Jhep = Jhep;  ap.Jpep = Jpep;
ap.Jthetaep = Jthetaep;  ap.Jhb = Jhb;  ap.Jpb = Jpb;  ap.Jthetab = Jthetab;
```
or, if the freeze is an intentional mean-axes-style approximation, freeze the member-level
matrices too (and skip their recomputation), so the same kinematics are used on both sides
of the virtual-work statement. Document whichever choice is made. This deserves a
regression test: static tip torsion of a deformed wing vs. a fine-mesh reference.

### 1.2 🔴 Fuselage gyroscopic moment is identically zero
`RigidFuselage.m`, `getCRB` (lines 29–32):
```matlab
CRB = [rigidfus.m*omega,             rigidfus.m*omega*rigidfus.ttpcm;
       rigidfus.m*rigidfus.tpcm*omega, rigidfus.I*omega];
```
`omega` is the skew matrix ω̃ (built in `dinamicaflex.m` line 59), and `CRB` multiplies
`beta = [v; ω]`. The rotational block therefore produces `I*ω̃*ω = I*(ω×ω) = 0`. The
correct Euler gyroscopic term is `ω×(Iω) = ω̃*I*ω`, i.e. the block must be `omega*rigidfus.I`
(ω̃ on the *left*):
```matlab
CRB = [rigidfus.m*omega,               rigidfus.m*omega*rigidfus.ttpcm;
       rigidfus.m*rigidfus.tpcm*omega, omega*rigidfus.I];
```
As written, any airplane whose fuselage inertia is non-spherical loses its gyroscopic
pitch/yaw/roll coupling (affects Dutch roll, spin, and any maneuvering simulation with
`fus.I ≠ 0`). The translational blocks are correct (`m ω̃ v` Coriolis force and
`m ω̃ (ω × c)`-type centripetal terms check out with the `matrixcross` sign convention).

### 1.3 🟠 Acknowledged missing rate-of-Jacobian terms — should be documented as model limitations
`dinamicaflex.m` lines 46, 58, 67 and `aeroforceandmoment.m` line 12:

- `CFF = ap.CG` — comment says it should include `Jhep'*Me*Jhepp` (Jacobian time
  derivative); `Element.setJhepp` exists but is commented out of `update()`.
- `CBF = ... * 0` — explicitly zeroed with the comment *"na realidade: = Jhb'*Me*Jhepp"*.
- `CBB` — comment: *"faltam termos aqui ainda!!!!! ver eq. 2.31 do Shearer, 2.33 do Weihua"*.
- Node accelerations `ppp = Jpep*strainpp + Jpb*betap` omit `J̇pep*strainp` and `J̇pb*beta`
  convective terms.
- `getJthetaep` (Element.m line 251) assumes the body frame does not rotate
  (*"VER SHEARER 183, D.39"*).

These are quadratic-rate terms: they vanish at trim and are small for slow structural
motion, but they enter the damping matrix, so they can bias flutter margins and large-
amplitude simulations. Recommendation: (a) list them in README/CLAUDE.md as known model
omissions; (b) since `setJhepp` and `getJhepep` already exist, wire them in behind a
`softPARAMS` flag and quantify the effect on example 2's flutter speed. The finite-
difference `getJhepep` is expensive (8 `getexpG` calls per element per node); an analytic
`d(dexpGdEps)/dε` from the same closed-form expansion used in `getexpG.m` would make it
affordable.

### 1.4 🔴 `trimairplanefull` / `equilibracorpo('full')`: steady-turn kinematics are wrong and the function cannot return
`Airplane.m` lines 326–354 and 396–414.

1. **Copy-paste bug in the turn-rate kinematics.** Both places set
   ```matlab
   beta(6) = psid*sin(phi)*cos(theta);   % should involve cos(phi)
   beta(4) = psid*sin(phi)*cos(theta);
   ```
   With the internal convention `Q = beta(4)`, `R = -beta(6)` (see §1.6), the standard
   steady-turn relations `q = ψ̇ sinφ cosθ`, `r = ψ̇ cosφ cosθ`, `p = −ψ̇ sinθ` give:
   ```matlab
   beta(4) =  psid*sin(phi)*cos(theta);        % q  (correct as-is)
   beta(6) = -psid*cos(phi)*cos(theta);        % -r (currently duplicates beta(4))
   re      = -beta(6);                          % = r
   beta(5) = -re*tan(theta)/cos(phi);           % p — formula is then consistent
   ```
   The existing `beta(5)` formula is only correct once `beta(6)` is fixed. As written, any
   trim with bank angle produces inconsistent body rates. (Both cases happen to agree when
   `phi = 0`, which is why level trim still works.)
2. **Gravity uses a stale attitude.** In the `'full'` case, `vec(1)` (θ) and `vec(6)` (φ)
   are trim unknowns but `kinetic(1)`/`kinetic(2)` are never updated inside
   `equilibracorpo`, so `dinamicaflex` computes the gravity vector with the initial guess
   (θ = 0.0537, φ = 0) forever. The `'longitudinal'` case does update `kinetic(1)`; the
   `'full'` case must do the same:
   ```matlab
   kinetic(1) = theta;  kinetic(2) = phi;
   ```
3. **Unassigned outputs.** The function signature is
   `function [vecEQ strainEQ] = trimairplanefull(...)` but only the local `vec` is
   assigned and `strainEQ` never is — calling it with any output errors out. It also never
   solves for the structural equilibrium (no `equilibraestrutura` stage). It should either
   be completed (mirror the longitudinal path) or clearly marked as WIP.
4. Hard-coded magic numbers (`aoa_trim = 0.0537`, `tracao = 6.4560`,
   `deltaflap = [-16.9703 0 0]`) are model-specific initial guesses baked into library
   code — move them to arguments with defaults.

### 1.5 🟠 Pinned DOFs handled by post-hoc zeroing of accelerations
`dinamicaflex.m` lines 110, 114, 120–122: for partially free rigid-body DOF sets, the
explicit branch solves the **fully coupled** `[MFF MFB; MBF MBB]` system and then multiplies
`bp` by `softPARAMS.vecFREEDEG` afterwards. That is not the constrained dynamics: the
flexible accelerations `xp` were obtained assuming the body accelerates freely, and then
the body acceleration is overwritten. The correct approach for prescribed (pinned) DOFs is
to partition: impose `β̇_pinned = 0`, solve
```
[MFF  MFB_free ] [ε̈      ]   [RHS_F]
[MBF_free MBB_ff] [β̇_free ] = [RHS_B,free]
```
i.e. delete the pinned rows/columns before the solve (keep them only for reaction-force
recovery). The implicit branch is closer to consistent (the residual forces `betap = 0` for
pinned DOFs), but the same partitioning would make both branches agree. Note the all-pinned
(`[0 0 0 0 0 0]`) and all-free cases are unaffected — which is what all shipped examples
use, so this has probably never been exercised.

Also fragile: `if softPARAMS.vecFREEDEG == [0 0 0 0 0 0]` (line 116) relies on
element-wise `==` inside `if`; use `if ~any(softPARAMS.vecFREEDEG)`.

### 1.6 🟠 State-ordering convention is undocumented and the example comments are wrong
From `dinamicaflex.m` lines 128–136 the actual convention is:

| index | `beta(i)` | meaning |
|-------|-----------|---------|
| 1 | V  | lateral velocity |
| 2 | U  | forward velocity (body y!) |
| 3 | −W | vertical |
| 4 | Q  | pitch rate |
| 5 | P  | roll rate |
| 6 | −R | yaw rate |

i.e. the body frame is (x spanwise/right, y forward, z up) and the rate ordering is
`[q p -r]`. Every example says `vecFREEDEG ... [u v w p q r]`, which mismatches on *all six*
entries. Anyone pinning "u" with `vecFREEDEG = [0 1 1 1 1 1]` actually pins lateral
velocity. Fix the comments in all examples and document the convention in one place
(README + a comment at the top of `dinamicaflex.m`). The gravity components (lines
100–103), Euler-rate equations (131–133) and `Hp` (137) are mutually consistent with this
convention — I verified them against the standard flat-Earth equations — so this is purely
a documentation/UX hazard, but a dangerous one.

### 1.7 🟠 Sectional inertia input convention is hybrid (off-diagonal sign)
`Node.m`, `calculateM` (lines 49–61) and `calculateMrb` (62–79).

The rotational block of the strain-based mass matrix needs the **second-moment tensor**
J = ∫ r rᵀ dm. The code builds its diagonal correctly from the inertia tensor,
`Jxx = ½(Iyy+Izz−Ixx)`, but copies the off-diagonals directly: `Jxy = +Ixy`. If `Inertia`
is the standard inertia *tensor* (off-diagonals = −∫xy dm), the correct conversion is
`Jxy = −Ixy`. As written, the code implicitly expects the user to supply **products of
inertia ∫xy dm with positive sign** in the off-diagonal slots while supplying tensor
moments on the diagonal — a hybrid convention nobody will guess. All shipped examples use
diagonal inertias, so it has never mattered, but the PAZY/experimental-wing work with tip
masses could hit it.

Additionally, `calculateMrb` applies parallel-axis terms (`+ m*rx^2`, `+ m*rx*ry`, …) —
so `rigid.I` is inertia **about the unit's CG** — whereas `calculateM` applies none, so
`Node.Inertia` is per-unit-length inertia **about the beam reference axis** even though a
CG offset `mcg` is also given. Both may be intended, but the two different conventions in
the same file need explicit doc comments (and ideally an off-diagonal sign fix +
convention note).

### 1.8 🔴 `trimairplane` iteration passes `deltaflap` with the wrong shape
`Airplane.m`: first structural solve (line 280) passes `[deltaflap 0 0]`, but inside the
`isITER` refinement loop (line 301) it passes the scalar `deltaflap`. With `ndelta = 1`
this silently works; with any surface indexed 2 or 3 it errors (`deltaflap(aero.ndelta)`
out of bounds) — and worse, the two calls trim against *different-sized* control vectors.
Fix line 301 to `[deltaflap 0 0]` (or, better, keep `deltaflap` a 3-vector throughout).

### 1.9 🟠 Strip-aerodynamics model: robustness and internal consistency
`aeroforceandmoment.m`, `getforceandmoment`:

- **`alphat = atan(-dotz/doty)`** (line 78) and the repeated `.../doty` terms: division by
  the local chordwise velocity is singular as `doty → 0` and picks the wrong branch for
  reversed local flow. Use `atan2(-dotz, doty)` and write the circulatory lift as
  `clalpha*rho*b*doty*(-dotz + (b/2-d)*dotalpha - lambda0)` (multiplying through by `doty`)
  so nothing is divided by a velocity. Near-hover/large-deformation cases (which this
  toolbox invites) currently blow up quietly.
- **Drag contains only `cd0`** — there is no induced/lift-dependent drag anywhere. For a
  free-flying trim, thrust is then equilibrated against parasite drag only, so trim thrust
  and the phugoid/spiral damping are underestimated. At minimum document it; better, add
  `cd = cd0 + cl^2/(pi*AR*e)` per strip or a leading-edge-suction-based tilt of the force
  vector.
- **Unsteady moment omits the inflow term.** In `case 3` the lift includes `−lambda0/doty`
  but the moment `M` doesn't include the corresponding `−0.5*doty*lambda0` circulatory
  term that appears in the (commented-out) "AED FULL DO WEIHUA" line 123. Since the moment
  is transferred to the elastic axis with the *full* lift (including λ₀) at line 127, the
  quarter-chord moment and the transfer term are consistent only if the circulatory moment
  about the ¾-chord-based AC truly vanishes — that holds for the flat-plate Theodorsen
  model, but then the same argument applies to the `(b/2−d)*dotalpha` term. Recommend
  re-deriving cases 1–3 from one reference (Peters/Hodges) and writing the chosen
  form in a comment with the equation number, so the mix of `clalpha`-based circulatory
  lift and π-based apparent-mass terms is auditable.
- **Wind is hard-coded along body-y** (`Vwind = [0;Vwindy;0]`, line 9) and is *added* to
  the structural velocity. Fine for the pinned-wing use case; misleading for free flight
  (a "wind" and a body speed `U` are equivalent only for a steady, uniform wind). Accepting
  a 3-vector would cost nothing and enable gust studies.
- `norm(Vwind)` on line 132 is `abs(doty)` of a scalar that was just renamed — harmless,
  but rename to `Uloc` for clarity.

### 1.10 🟠 Multiple engines all inherit the first engine's parameters
`Engine.m`, `getFPROP`: the parameter unpack (lines 25–34) happens **outside** the
`for i = 1:size(prop,2)` loop, so `Fmax, V0, rho0, nrho, nv, alphaf, betaf` come from
`prop(1)` and are applied to every engine. Two engines with different ratings/orientations
silently get identical forces. Move the unpack inside the loop
(`prop(i).params...`). Also: `CBb = eye(3)*Cwb` is a no-op multiply, and the follower-load
option exists only as a commented block — worth either implementing (the deformed direction
`h(4:6)…` is available) or deleting.

### 1.11 🟡 `flutter_speed` returns eigen-data that may belong to the stable bound
`flutter_speed.m`: after the bisection loop, `unstable_eig_value/vec` come from the **last
midpoint evaluated**, which is `Vwind_initial` (the stable side) whenever the final test was
stable, while `unstable_speed = Vwind_initial` is the *stable* bracket edge. So the
returned "unstable eigenvalue" can have negative real part. Also, if the initial interval
is already within `tol`, `eig_val/ind_max` are undefined (runtime error). Fix: after
convergence, evaluate once at `Vwind_final` (guaranteed unstable side) and return that
speed/eigenpair, or return both brackets. Bisection also assumes a single stable→unstable
crossing — fine, but state it in the header next to the existing caveats.

### 1.12 🟡 Spanwise coordinate normalized by a hard-coded 20
`create_uniform_structure.m` lines 11–13: `s = ((i-1)*ds)/20` — the 20 is a leftover wing
length. `Node.s` currently only feeds the commented-out tip-loss correction `fcorr` in
`aeroforceandmoment.m`, but the moment someone re-enables that correction on a wing that
isn't 20 m, the tip-loss factor is wrong. Use `.../Length` (available in the function).

### 1.13 🟡 Gravity constant inconsistencies
`atmosfera.m` uses `g0 = 9.81` internally; examples set `softPARAMS.g = 9.8`; the ISA
standard is 9.80665. The 0.1 % difference is irrelevant individually but means density and
weight are computed with *different* g. Pick 9.80665 in both places (or derive
`atmosfera`'s exponent from `softPARAMS.g`).

---

## 2. Numerical robustness

### 2.1 🔴 Linearization step size `delta = 1e-9` is far too small
`Airplane.m` line 229. For central differences in double precision the optimal relative
step is ~`eps^(1/3) ≈ 6e-6`; with `δ = 1e-9`, subtractive cancellation leaves only ~5–6
significant digits *at best*, and much less for states of magnitude ≫1 (velocities of
tens of m/s, altitude ~2×10⁴ m: the altitude column is differentiated with a *relative*
step of 5×10⁻¹⁴ — essentially noise). Flutter speeds from `flutter_speed` inherit this
noise through `max(real(eig(A)))`. Fix:
```matlab
delta_i = 1e-6 * max(1, abs(x(i)));   % per-state scaled step
```
Better: dinamicaflex is built from `expm`-free analytic operations, `atan`, and linear
algebra — it is complex-step safe or nearly so (replace `atan(-dotz/doty)` with a
complex-safe form and `norm(Vwind)` with `sqrt(Vwind^2)`). Complex-step differentiation
(`x + 1i*1e-100*e_i`) gives machine-precision Jacobians with the same call count.

### 2.2 🟡 fsolve tolerances below machine precision + magic residual scaling
`Airplane.m` lines 268, 308: `TolFun = 1e-15…1e-20`, `TolX = 1e-15…1e-17` are below what
double precision can deliver; fsolve just burns its iteration budget. The `10000*`
multiplier on the body-trim residual (lines 395, 413) exists only to fight these
tolerances. Use `TolFun/TolX ≈ 1e-10` with *scaled* residuals (divide force residuals by
`m·g`, moments by `m·g·c`), and drop the 10000 factor.

### 2.3 🟡 `getexpG` small-curvature switch
`getexpG.m` line 19: the series/closed-form switch at `λ² > 1e-10` (λ·ds ≈ 3×10⁻⁵ for
ds ≈ 1) is safe for the function value, but the closed-form coefficients
`a = ds/λ² − sin(λ ds)/λ³` suffer growing cancellation just above the threshold; and in the
series branch the derivative terms `dadEps, dbdEps` are set to exactly zero, which is only
the true limit *at* λ = 0. Use the next series terms instead:
```matlab
a = ds^3/6 - ds^5*lambda2/120;    b = ds^2/2 - ds^4*lambda2/24;
dadEps = -ds^5/120 * 2*[0;kx;ky;kz]-block;   % d a / d k = (da/dλ²)·2k
dbdEps = -ds^4/24  * 2*[0;kx;ky;kz]-block;
```
and raise the switch to `λ·ds < 1e-2` so both branches are accurate at the crossover.
(Also: the comment at line 47 recording that the theses' `dKEpsdkz` appears to have a sign
error is valuable — promote it from a code comment to documentation.)

### 2.4 🟡 `StructuralModes` can emit complex frequencies
`StructuralModes.m` line 27: `sqrt(autoval)` — near-rigid-body or ill-conditioned modes can
have tiny negative eigenvalues, producing complex "frequencies" that later plotting/sorting
handles oddly. Use `sqrt(max(autoval,0))` (or report a warning), and since `KFF`, `MFF` are
symmetric, call `eig(KFF, MFF, 'chol')` for guaranteed real spectra and better speed.

---

## 3. Computational cost

### 3.1 🔵🔴 Silent array growth: three functions allocate the wrong size
These all "work" only because MATLAB grows arrays on out-of-bounds assignment — repeatedly
reallocating and copying inside the ODE right-hand side:

| function | allocated | actually filled |
|---|---|---|
| `OMEGA.m` line 2 | `zeros(4n)` → 4n×4n | 12n×12n |
| `Element.getJpep` line 192 | `zeros(m*3, m*4)` | 9m×4m |
| `Element.getmemberJthetab` line 184 | `zeros(m*3, 6)` | 9m×6 |

`OMEGA` is the worst: it is called on **every** dynamic RHS evaluation
(`dinamicaflex.m` line 61) and grows a dense (12·3·NUMele)² matrix block-by-block. For 20
elements that is a 720×720 dense matrix rebuilt with ~240 growing assignments per time
step. Fix the allocations, and better, don't build OMEGA at all:
```matlab
% Hhb = OMEGA(omega,n)*Jhb  ==  apply 3x3 omega to each 3-row block of Jhb:
Hhb = reshape(omega * reshape(ap.Jhb.', 6*?, []) ...)   % or simply:
Hhb = kron(speye(4*nodesnumber), sparse(omega)) * ap.Jhb;
```
The `kron(speye, sparse)` one-liner is both correct-by-construction and ~O(n) instead of
O(n²) memory.

### 3.2 🔵 Everything is dense; the structure is block-sparse
`Me` (36n×36n), `B` (9n×9n), `KG`, `CG` (4n×4n) are block-diagonal; `Jhep` (36n×4n) is
block-lower-triangular; `Jhb` is dense-ish but tall-thin. Storing them dense makes every
triple product `Jhep'*Me*Jhep` (computed on **every** non-equilibrium RHS call,
`dinamicaflex.m` lines 45–67) cost O(n³)-ish dense flops. Two cheap wins:

1. Build the constant matrices sparse once in the constructor:
   `ap.Me = sparse(blkdiag(...))` (or collect blocks in a cell and call
   `blkdiag(cellblocks{:})` once — the current incremental `blkdiag` in a loop is itself
   O(n²) in copies).
2. Compute `MeJhep = ap.Me*ap.Jhep` once per Jacobian update and reuse it for
   `MFF = Jhep'*MeJhep`, `MFB/MBF`, `CFB` — currently `Me*Jhep`-type products are formed
   4–6 times per call.

### 3.3 🔵 Cache the generalized mass matrices (and their factorizations) when Jacobians are frozen
With `softPARAMS.updateStrJac ∈ {0,1}` the Jacobians are constant throughout a simulation,
so `MFF, MFB, MBF, MBB, CFF` are constant — yet `dinamicaflex.m` recomputes all of them on
every RHS evaluation, and re-factorizes (`MFF\…`, `MBB\…`, or the coupled backslash) each
time. Cache them (e.g. store in `ap` and invalidate inside `updateStrJac`), and store a
Cholesky factor of `MFF`/the coupled matrix: `ode15s/ode15i` typically make 10³–10⁵ RHS
calls, so this alone can be an order-of-magnitude speedup for the simulation path. Only
`CFB/CBB` (ω- and ḣ-dependent) genuinely change per call.

### 3.4 🔵 `updateStrJac` computes expensive Jacobians and throws them away
Because of the `isempty` guard (§1.1), on every call after the first, `getJthetaep`
(the function whose own comment says it takes ~0.1 s for a 10-element beam),
`getmemberJhb`, `getmemberJpb`, `getmemberJthetab` are fully recomputed only to be
discarded at the airplane level (the member-level assignment survives). Whatever the
resolution of §1.1, this waste should go: either use the results or don't compute them.

### 3.5 🔵 `getJthetaep` is an O(m²·4) scalar-loop hot spot
`Element.m` lines 237–276. Improvements, in order of effort:
1. The `ii` loop runs the full range and skips `i < ii`; make it `for ii = 1:i`.
2. The three dot products per (i, ii, j) reuse the same `h` vectors; hoist
   `W = [h(4:6) h(7:9) h(10:12)]` per node out of the (ii, j) loops.
3. Vectorize over j (all 4 strain components at once): each `dtdep*` is a 1×4 row from
   `Jhep(rows, cols)'*hvec`; assemble columns with one 3×4 matrix product per node instead
   of a scalar loop.
4. `Jthetab` (line 182) is a constant pattern — compute once at construction, not per
   update.

### 3.6 🔵 Trim: Jacobians are refreshed inside every residual evaluation, and the outer Picard loop duplicates that work
With `updateStrJac = 1`, `dinamicaflex.m` line 37 calls `updateStrJac(ap)` on **every**
`equilibraestrutura` residual evaluation — and fsolve, given no Jacobian, evaluates the
residual 4·NUMele+1 times per iteration just for finite differences. For 20 elements that
is ~80 full structural-Jacobian rebuilds per fsolve iteration, each including the O(m²)
`getJthetaep`. Meanwhile the `isITER` outer loop *already* implements the freeze-solve-
update Picard strategy. Recommendations:
1. During the fsolve solve, freeze the Jacobians (update only in the outer `isITER` loop);
   this changes the iterate path but not the converged equilibrium, since the residual at
   convergence is evaluated at the converged geometry after the final update.
2. Give fsolve the sparsity pattern (`JacobPattern`: `KFF` is block-diagonal, the aero
   coupling is member-local) or an analytic Jacobian `−KFF + ∂RAEROF/∂ε`.
3. Replace the 20 000-iteration budgets with sane limits once the tolerances (§2.2) are
   fixed.

### 3.7 🔵 `linearize` does finite differences for things it knows analytically
`Airplane.m` lines 233–243: `Mlin` columns for the `strainp` and `beta` blocks are obtained
by FD, but they are exactly `[MFF MFB; ...]`-type blocks already formed inside
`dinamicaflex`; similarly the first-block rows of `Alin` (`ε̇ = εp`) are identity. Even
keeping the FD design, you can skip the ~n columns/rows that are known, halving the cost;
each FD column currently costs two full RHS evaluations, each with a complete `update()`
(all element exponentials) — for the flutter bisection (≈10 `linearize` calls) this
dominates the run time.

### 3.8 🔵 Growing arrays in the aero assembly
`aeroforceandmoment.m`: `F`, `M`, `lambdap` are grown by indexed assignment without
preallocation, and `dinamicaflex.m` lines 70–82 grow `FAERO/MAERO/FLAMBDA` by
concatenation per member. Preallocate (`F = zeros(9*membersize,1)`, etc.). Also
`FLAMBDA = [FLAMBDA; []]` on line 80 silently produces a size mismatch if a *later* member
has aero states while an earlier one doesn't (lambda ordering assumes aero members come
first) — preallocating by `ap.membNAEDtotal` offsets fixes both.

### 3.9 🔵 Minor arithmetic reuse
- `Peter.m` constructor: `inv(obj.A)` is formed four times (lines 67–70); compute
  `Ainv = inv(obj.A)` once (or `obj.A\[B1 B2 ...]`).
- `Element.setexponentials`: `exp2G = expG*expG` (the exponentials commute) would save one
  of the two `getexpG` calls per element per update — but note the *derivative*
  `dexp2GdEps` is not simply related, so keep the analytic call for the derivative or use
  the product rule `d(expG²) = dexpG·expG + expG·dexpG` (block-wise).
- `Airplane.update`: `sum(ap.membSIZES(1:(i-1)))` recomputed per member — precompute
  cumulative offsets in the constructor (`[0 cumsum(...)]`).

---

## 4. Code health (quick list)

- **`avifile`/`addframe` were removed from MATLAB (R2014b+)** — `Airplane.airplanemovie`'s
  AVI branch is dead on any modern MATLAB. Port to `VideoWriter`. Also
  `elseif format == 'gif'` errors for any format string not 3 chars long — use
  `strcmpi(format,'gif')`.
- Dead/experimental code: `equilibratudo`, `getnumJhep`, `getnumexpG`, the `Jhepp`
  machinery (`setJhepp`, `product3dXvec`), commented `fcorr` tip-loss lines. Either wire
  them behind flags (see §1.3) or delete — they currently read like live physics.
- `global softPARAMS` everywhere makes results depend on hidden state; passing a params
  struct (or storing it in `Airplane`) would make runs reproducible and parallelizable
  (`parfor` over flutter speeds currently impossible without global broadcast hacks).
- Shadowed built-ins: `diff` (flutter_speed), `beta` (MATLAB's beta function) — cosmetic
  but confusing.
- Comment encoding is mangled (Latin-1 "Fl�vio", "n�") — re-save sources as UTF-8.
- `trimairplane` prints `theta/deltaflap/tracao` by leaving off semicolons — use
  `fprintf` guarded by a verbosity flag.
- No automated tests: `tests/exact.m` provides exact Euler–Bernoulli/torsion frequencies
  (I checked the torsion formula — the `/2/2` with the π absorbed into Hz conversion is
  correct) but nothing runs it automatically. A tiny script asserting
  `StructuralModes` frequencies within 1 % of `exact()` and example 2's flutter speed
  within tolerance would protect all the fixes above.

---

## 5. Suggested priority

1. **§1.1 Jacobian freeze inconsistency** — decide intended semantics, make member/airplane
   levels consistent (correctness of every deformed-configuration result).
2. **§1.2 `getCRB` gyroscopic term** — one-line physics fix.
3. **§1.4 + §1.8** — trim bugs (`beta(6)` copy-paste, stale `kinetic`, `deltaflap` shape,
   unassigned outputs).
4. **§2.1 + §2.2** — FD step and tolerances (accuracy of every linearization/flutter
   result).
5. **§3.1–§3.3** — allocation fixes, sparse/cached mass matrices (largest speedups,
   zero physics risk).
6. **§3.5–§3.7** — trim/linearize cost reduction.
7. **§1.3 / §1.9** — document (or implement) the omitted damping terms and aero-model
   caveats; add induced drag.
8. **§4** — regression tests, then the cleanups.
