# BrownianMotion 0.4.0

## Breaking changes

* All exported functions now use snake_case names, so that no function name can
  be mistaken for an S3 method (e.g. `create.bm` previously read as a `create()`
  method for class `bm`). The old names have been removed:

  | Old name            | New name            |
  |---------------------|---------------------|
  | `create.bm()`       | `create_bm()`       |
  | `add.points()`      | `add_points()`      |
  | `add.segment()`     | `add_segment()`     |
  | `add.labels()`      | `add_labels()`      |
  | `delete.labels()`   | `delete_labels()`   |
  | `delete.skeleton()` | `delete_skeleton()` |
  | `concat.bm()`       | `concat_bm()`       |
  | `first.passage()`   | `first_passage()`   |

  `sim()`, `layers()`, `refine()` and `plot()` are unchanged.

## Bug fixes

* Multi-dimensional Brownian motion with a non-diagonal covariance matrix now
  has the requested covariance. Previously the transform used the transpose
  of the Cholesky factor, so paths had covariance `U %*% t(U)` instead of
  `cov = t(U) %*% U`.
* Simulating or refining within a localised bridge layer that ends at a jump
  (e.g. at a `concat_bm()` join) now draws from the correct distribution.
  Previously the right end of the auxiliary path was taken after the jump
  rather than before it, giving draws outside the layer bounds.
* `add_segment()` now moves the bounds of localised bridge layers by the
  correct amount.
* `delete_skeleton(type = "label")` now deletes only user labels; previously
  it also deleted the layers in the interval.
* `create_bm()` now rejects non-numeric and non-finite `t`, `W_t` and `cov`
  with a clear error, rather than failing later.
* Times and values in a path no longer pick up stray names (e.g. `q`) from
  simulation or from named times passed to `add_points()`. Names on times
  are still used as labels.
* Adding localised layers before existing localised bridge intervals now
  preserves the existing intervals' transformed layer bounds.
* Simulating or refining within localised layers requested via
  `layers(type = "localised")` no longer fails with "missing value where
  TRUE/FALSE needed".
* Indexing a `BrownianMotion` object with exactly one layer now works both at
  and off skeleton points.
* Conditional simulation within intersection layers no longer fails due to a
  missing internal helper.
* Updated ggplot2 usage for compatibility with current ggplot2 releases, and
  added missing namespace imports.

## Internal

* Internal helper functions renamed to avoid clashing with S3 dispatch
  (e.g. `sim.condbessel_()` is now `sim_condbessel_()`, and `update.lyr_()`,
  which R treated as a `stats::update()` method, is now `update_lyr_()`).

# BrownianMotion 0.3.0

* New `add.labels()` and `delete.labels()` for managing user labels on existing
  path times; `delete.skeleton()` gains the ability to delete labels.
* Tutorial vignette completed.

# BrownianMotion 0.2.0

* Multi-dimensional Brownian motion with arbitrary covariance matrix, via
  `create.bm(dim = , cov = )`, with access to the transformed path through `$Z`.
* Unified layer construction through `layers()` (Bessel, intersection and
  localised layers).
* New `add.points()`, `add.segment()` and `concat.bm()`.
* `delete.skeleton()` for removing path observations and layers.
* Labelling of path times.
* `[` indexing of Brownian motion objects by time, index, time range and label.
* Multi-dimensional pairwise plotting (requires patchwork; `sf` suggested for
  union layer plots).

# BrownianMotion 0.1.0

* Initial version: unconditional simulation, Brownian bridges, first passage
  simulation and Bessel layers.
