# FaradayCages

## faradaycages.m
This has the most updated functionality, including:
- Visualization for magnitude of gradients(current bug: imagesc flipping image, behavior when wires are too big)
- Multiple point charges: change the value of zs to change/add/remove charges
- Different shapes, made with parametric equations
  -- Circle
  -- Heart
  -- Diamond
  -- Rose Curve
  -- Hypotrochoid(see https://www.desmos.com/calculator/r8bsgvsszo for how to edit the values to make a closed curve)

## faradaycages_multiplepointcharges.m
Is what it sounds like. This uses all the original code from the paper, but edits some functions to allow for additional point charges.

## faradaynoedits.m
The original code from the paper

## filteringquiver.m
This file tests how best to filter the gradient vectors when they are close to the point charge so that the visualization has a reasonable threshold.
