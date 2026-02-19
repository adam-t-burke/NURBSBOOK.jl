# Surfaces (Chapters 3–4)

## B-spline Surfaces (A3.5–A3.8)

```@docs
surface_point(::BSplineSurface, ::Real, ::Real)
surface_derivatives(::BSplineSurface, ::Real, ::Real, ::Int)
```

## NURBS Surfaces (A4.3–A4.4)

```@docs
surface_point(::NURBSSurface, ::Real, ::Real)
surface_derivatives(::NURBSSurface, ::Real, ::Real, ::Int)
```

## Derivative Control Points

```@docs
surface_deriv_control_points
```

## Surface Construction (Chapters 8, 10)

```@docs
make_bilinear_surface
make_ruled_surface
make_cylinder
make_revolved_surface
make_swung_surface
make_skinned_surface
make_swept_surface
make_gordon_surface
make_coons_surface
```
