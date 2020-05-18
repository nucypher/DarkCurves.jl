# Version history


## v0.2.0 (18 May 2020)

* CHANGED: renamed `curve_modulus()` to `curve_coordinate_modulus()`.
* CHANGED: renamed `EndomorphismType4`, `endomorphism()`, `curve_endomorphism_beta()`, `curve_endomorphism_lambda()` to `EndomorphismType4Curve`, `curve_endomorphism_type_4()`, `curve_endomorphism_type_4_beta()`, `curve_endomorphism_type_4_lambda()`, respectively.
* CHANGED: `EllipticCurvePoint` is now parametrized only by the curve type, but not by the coordinate type (to support custom curve point types).
* CHANGED: renamed `curve_coeff_a()` and `curve_coeff_b()` to `curve_weierstrass_coeff_a()` and `curve_weierstrass_coeff_b()`.
* CHANGED: not exporting point types or coefficient functions.
* ADDED: `rand()` for curve points.
* ADDED: a single endomorphism-based multiplication method.
* ADDED: wNAF-based batch mul without endomorphism.
* ADDED: `WeierstrassCurve` abstract type subtyping `EllipticCurve`.
* ADDED: `curve_point_type()` and `curve_scalar_type()` as main API entry points.
* ADDED: `StandardEllipticCurvePoint` abstract type for built-in `AffinePoint`, `JacobianPoint` and `ChudnovskyPoint`.


## v0.1.1 (6 March 2020)

* ADDED: exportng `EllipticCurve` and `EllipticCurvePoint`.
* ADDED: made typed curve metadata functions generated.
* ADDED: endomorphism type 4 support.
* ADDED: wNAF-based `batch_mul()`.


## v0.1.0 (3 March 2020)

Initial version.
