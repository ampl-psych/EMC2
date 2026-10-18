#!/usr/bin/env python3
"""The ONNX -> EMC2 card converter moved to the nleflow package.

    nleflow convert <model.onnx> --preset ddm_uniform_st [--verify 300]

(nleflow/ in the NLE repository, https://github.com/andrewheathcote/nle; the
steps are in nleflow/docs/lan-in.md there.) This file used to hold the converter; it
is kept as a pointer so that a call to it says where to go.
"""
import sys

sys.exit(__doc__ + "\nThis script no longer converts anything; use `nleflow convert`.")
