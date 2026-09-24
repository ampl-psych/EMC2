#!/usr/bin/env python3
"""Convert an MLP likelihood network in ONNX format into an EMC2 model card.

EMC2 does not run ONNX. Published likelihood approximation networks (LANs,
e.g. the HSSM / LANfactory models on huggingface.co) are plain multilayer
perceptrons; the .onnx file is only the container for their weights. This
script reads the container and writes a JSON model card that
`register_nn_model()` (kind "mlp_joint") evaluates with EMC2's own batched MLP
evaluator, in the compiled likelihood. No onnxruntime, no new system
dependency at run time.

What it does NOT know: an ONNX graph holds weights and an input width, not
what the inputs mean, on what scale the network was trained, or over what
region. You must say (or use a --preset whose values were checked against the
network's source):

  --params      the parameter inputs, in the network's order
  --box         the training box per parameter, e.g. v=-3:3,a=0.3:2.5
  --ll-floor    the floor of the log-density labels the network was trained
                on (rows outside the box get this value)
  --response-values  the network's input for responses 1 and 2 (EMC2's R = 1
                is the lower boundary; HSSM codes lower/upper as -1/+1); write
                a leading minus with an equals sign: --response-values=-1,1

The network input is assumed to be [params..., rt, response] with no
standardisation (the LAN convention; override with --layout if not). Only
tanh hidden layers are supported (the evaluator refuses anything else, and so
does this script) and one linear output, the log density.

    onnx_to_card.py ddm_uniform_st.onnx out.json --preset ddm_uniform_st --verify 2000

Requires `onnx` and numpy; `--verify` also needs onnxruntime and compares the
card's forward pass (numpy, from the JSON just written) with onnxruntime on
random inputs.

NOTE: converting does not validate. A LAN's calibration is unknown until it
has been checked (parameter recovery, simulation-based calibration); see
?register_nn_model.
"""
import argparse
import hashlib
import json
import sys

import numpy as np

# Values checked against the network's source. `ddm_uniform_st`: HSSM's
# "ddm-st-lans" DDM with uniform non-decision-time variability (huggingface.co/
# Eitanm/ddm-st-lans :: ddm_uniform_st.onnx). Input [v, a, z, t, st, rt,
# choice], choice in {+1 (upper), -1 (lower)}. EMC2 parameterisation
# (verified against the DDM density, tests/testthat/test-nn-mlp.R):
#   a_EMC2 = 2 a, Z = z, t0 = t - st, st0 = 2 st, sv = 0, SZ = 0, s = 1.
PRESETS = {
    "ddm_uniform_st": dict(
        params=["v", "a", "z", "t", "st"],
        box={"v": (-3.0, 3.0), "a": (0.3, 2.5), "z": (0.3, 0.7),
             "t": (0.25, 2.25), "st": (0.001, 0.25)},
        ll_floor=-16.11809565095832,      # log(1e-7)
        response_values=(-1.0, 1.0),
        notes={"source_model": "huggingface.co/Eitanm/ddm-st-lans :: ddm_uniform_st.onnx",
               "encoding": "input = [v,a,z,t,st,rt,choice] with choice in {+1,-1}; no standardisation",
               "map_to_emc2_ddm": "a_EMC2 = 2*a_lan; Z = z; t0 = t - st; st0 = 2*st; sv = 0; SZ = 0; s = 1",
               "calibration": "unknown; needs its own SBC cell before inference"}),
}

SUPPORTED_ACT = {"Tanh": "tanh"}


def load_chain(path):
    """The graph as [(W (n_in x n_out), b)] and the hidden activation name."""
    import onnx
    from onnx import numpy_helper

    g = onnx.load(path).graph
    init = {t.name: numpy_helper.to_array(t).astype(np.float64) for t in g.initializer}
    for n in g.node:                                  # Constant nodes are weights too
        if n.op_type == "Constant":
            for a in n.attribute:
                if a.name == "value":
                    init[n.output[0]] = numpy_helper.to_array(a.t).astype(np.float64)
    if len(g.input) - len([i for i in g.input if i.name in init]) != 1:
        sys.exit("expected exactly one network input")
    cur = [i.name for i in g.input if i.name not in init][0]
    layers, acts, pend = [], [], None

    def const(name):
        if name not in init:
            sys.exit(f"tensor '{name}' is computed, not a stored weight; not a plain MLP")
        return init[name]

    for n in g.node:
        if n.op_type == "Constant":
            continue
        if n.op_type == "MatMul" and n.input[0] == cur and pend is None:
            W = const(n.input[1])
            if W.ndim != 2:
                sys.exit("MatMul weight is not a matrix")
            pend = [W, None]
            cur = n.output[0]
        elif n.op_type == "Add" and pend is not None and pend[1] is None and cur in n.input:
            other = [x for x in n.input if x != cur]
            if len(other) != 1:
                sys.exit("unsupported Add")
            pend[1] = const(other[0]).reshape(-1)
            layers.append(tuple(pend))
            pend, cur = None, n.output[0]
        elif n.op_type == "Gemm" and n.input[0] == cur:
            at = {a.name: a for a in n.attribute}
            alpha = at["alpha"].f if "alpha" in at else 1.0
            beta = at["beta"].f if "beta" in at else 1.0
            if "transA" in at and at["transA"].i:
                sys.exit("Gemm with transA is not supported")
            W = const(n.input[1])
            if "transB" in at and at["transB"].i:
                W = W.T
            b = const(n.input[2]).reshape(-1) if len(n.input) > 2 else np.zeros(W.shape[1])
            layers.append((alpha * W, beta * b))
            cur = n.output[0]
        elif n.op_type in SUPPORTED_ACT and n.input[0] == cur:
            acts.append(SUPPORTED_ACT[n.op_type])
            cur = n.output[0]
        elif n.op_type == "Identity" and n.input[0] == cur:
            cur = n.output[0]
        else:
            sys.exit(f"unsupported operation '{n.op_type}' (inputs {list(n.input)}): "
                     "this converter handles MatMul/Add or Gemm layers with tanh activations only "
                     "(an activation other than tanh needs support in the EMC2 evaluator first)")
    if pend is not None:
        sys.exit("a MatMul without its bias Add")
    if not layers:
        sys.exit("no layers found")
    if cur != g.output[0].name:
        sys.exit("the layer chain does not end at the graph output")
    if len(acts) != len(layers) - 1:
        sys.exit(f"{len(layers)} layers need {len(layers) - 1} hidden activations; found {len(acts)} "
                 "(the output layer must be linear)")
    if len(set(acts)) > 1:
        sys.exit("mixed hidden activations are not supported")
    for i in range(1, len(layers)):
        if layers[i][0].shape[0] != layers[i - 1][0].shape[1]:
            sys.exit(f"layer {i + 1} input width does not match layer {i} output width")
    for W, b in layers:
        if b.shape[0] != W.shape[1]:
            sys.exit("a bias does not match its layer's output width")
    if layers[-1][0].shape[1] != 1:
        sys.exit("the network must have a single output (the log density)")
    return layers, (acts[0] if acts else "tanh")


def forward(layers, act, X):
    """numpy forward pass; X is n x n_in."""
    H = X
    for i, (W, b) in enumerate(layers):
        H = H @ W + b
        if i < len(layers) - 1:
            H = np.tanh(H)
    return H[:, 0]


def parse_box(text):
    box = {}
    for item in text.split(","):
        name, rng = item.split("=")
        lo, hi = rng.split(":")
        box[name.strip()] = (float(lo), float(hi))
    return box


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("onnx")
    ap.add_argument("out", help="model card to write (.json)")
    ap.add_argument("--preset", choices=sorted(PRESETS))
    ap.add_argument("--params", help="comma-separated parameter inputs, in the network's order")
    ap.add_argument("--box", help="training box, e.g. v=-3:3,a=0.3:2.5")
    ap.add_argument("--ll-floor", type=float, help="log-density floor of the training labels")
    ap.add_argument("--response-values", help="network inputs for responses 1 and 2, e.g. -1,1")
    ap.add_argument("--layout", help="input layout if not [params..., rt, R]; comma-separated of "
                    "parameter names and rt / log_rt / R")
    ap.add_argument("--transforms", help="per parameter scale the network was trained on "
                    "(identity|log|probit), e.g. a=log; default identity")
    ap.add_argument("--verify", type=int, default=0, metavar="N",
                    help="compare the card's forward pass with onnxruntime on N random rows")
    a = ap.parse_args()

    p = PRESETS[a.preset] if a.preset else {}
    params = a.params.split(",") if a.params else p.get("params")
    box = parse_box(a.box) if a.box else p.get("box")
    floor = a.ll_floor if a.ll_floor is not None else p.get("ll_floor")
    resp = tuple(float(x) for x in a.response_values.split(",")) if a.response_values else p.get("response_values")
    for what, val in (("--params", params), ("--box", box), ("--ll-floor", floor), ("--response-values", resp)):
        if val is None:
            sys.exit(f"{what} is required (or use --preset): an ONNX graph does not record it")
    if set(box) != set(params) or len(resp) != 2:
        sys.exit("--box must give a range for every parameter (and no others); --response-values needs two values")
    for k, (lo, hi) in box.items():
        if not lo < hi:
            sys.exit(f"empty box for {k}")
    tf = {k: "identity" for k in params}
    if a.transforms:
        for item in a.transforms.split(","):
            k, v = item.split("=")
            if k not in tf or v not in ("identity", "log", "probit"):
                sys.exit(f"bad --transforms entry '{item}'")
            tf[k] = v
    layout = a.layout.split(",") if a.layout else list(params) + ["rt", "R"]

    layers, act = load_chain(a.onnx)
    n_in = layers[0][0].shape[0]
    if len(layout) != n_in:
        sys.exit(f"the network has {n_in} inputs but the input layout has {len(layout)} entries "
                 f"({layout}); give --params/--layout")
    if not set(params) <= set(layout) or len(set(layout)) != len(layout) or \
            not set(layout) - set(params) <= {"rt", "log_rt", "R"}:
        sys.exit("layout must contain every parameter once, plus only rt / log_rt / R")

    with open(a.onnx, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    natural_lo = [box[k][0] for k in params]
    natural_hi = [box[k][1] for k in params]
    card = {
        "kind": "mlp_joint",
        "source": a.onnx.split("/")[-1],
        "source_sha256": sha,
        "context_names": params,
        "context_transforms": tf,
        "bounds_natural": {"lower": natural_lo, "upper": natural_hi},
        "input_layout": layout,
        "response_values": list(resp),
        "ll_floor_log": floor,
        "mlp": {"activation": act, "use_norm": False,
                "layers": [{"W": W.tolist(), "b": b.tolist()} for W, b in layers]},
        "converted_with": "inst/scripts/onnx_to_card.py",
    }
    card.update(p.get("notes", {}))
    with open(a.out, "w") as fh:
        json.dump(card, fh)
    dims = [n_in] + [W.shape[1] for W, _ in layers]
    print(f"wrote {a.out}: {'-'.join(map(str, dims))} {act} MLP, inputs {layout}, "
          f"{len(params)} parameters, ll floor {floor}, source sha256 {sha[:16]}")

    if a.verify:
        try:
            import onnxruntime as ort
        except ImportError:
            print("--verify: onnxruntime is not installed; skipped")
            return
        rng = np.random.default_rng(1)
        with open(a.out) as fh:
            back = json.load(fh)
        cl = [(np.array(l["W"]), np.array(l["b"])) for l in back["mlp"]["layers"]]
        n = a.verify
        X = np.zeros((n, n_in))
        for j, name in enumerate(layout):
            if name in box:
                X[:, j] = rng.uniform(*box[name], n)
            elif name == "rt":
                X[:, j] = rng.uniform(0.01, 5.0, n)
            elif name == "log_rt":
                X[:, j] = np.log(rng.uniform(0.01, 5.0, n))
            else:
                X[:, j] = rng.choice(resp, n)
        sess = ort.InferenceSession(a.onnx)
        name = sess.get_inputs()[0].name
        ref = np.array([sess.run(None, {name: X[i].astype(np.float32)})[0].reshape(-1)[0]
                        for i in range(n)], dtype=np.float64)
        mine = forward(cl, act, X)
        print(f"--verify: card forward pass vs onnxruntime on {n} random rows: "
              f"max |diff| = {np.max(np.abs(mine - ref)):.3g} (float32 vs float64 arithmetic)")


if __name__ == "__main__":
    main()
