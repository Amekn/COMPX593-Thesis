#!/usr/bin/env python3
"""Summarise weight drift against a pretrained baseline.

The bonito export command writes each parameter as a TorchScript ``.tensor`` archive.
This utility compares fine-tuned model tensors against the pretrained sup baseline and writes compact CSV tables for plotting.
"""

from __future__ import annotations

import argparse
import gc
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import torch

# Define a consistent order for models
MODEL_ORDER = [
    "Alpha",
    "AlphaQCH",
    "AlphaQCL",
    "Gamma",
    "GammaQCH",
    "GammaQCL",
    "Beta",
    "BetaQCH",
    "BetaQCL",
]

# Define a consistent order for neural network components
COMPONENT_ORDER = {
    "conv": 0,
    "transformer_self_attn": 1,
    "transformer_ffn": 2,
    "transformer_norm_deepnorm": 3,
    "upsample": 4,
    "crf": 5,
    "other": 99,
}

# Define easier to read labels for components
COMPONENT_LABELS = {
    "conv": "Convolution front-end",
    "transformer_self_attn": "Transformer self-attention",
    "transformer_ffn": "Transformer FFN",
    "transformer_norm_deepnorm": "Transformer norm/deepnorm",
    "upsample": "Upsample",
    "crf": "CRF head",
    "other": "Other",
}

HISTOGRAM_EDGES = np.concatenate(([0.0], np.logspace(-8, 8, 129), [np.inf]))

# Percentiles to compute for reporting 50th, 95th, 99th, and 99.9th percentiles of abs_delta and abs_dalta_z
QUANTILES = (0.50, 0.95, 0.99, 0.999)

# Thresholds for counting extreme "tail" values of abs_delta_z
TAIL_THRESHOLDS = (1.0, 2.0, 3.0, 5.0)

# Small constant to avoid division by zero in relative L2 and cosine similarity calculations
EPS = 1e-12


# immutable data container that holds the metadata for a single model configuration
@dataclass(frozen=True)
class ModelSpec:
    # name of the model
    name: str
    # filesystem path to the model's directory
    dorado_dir: Path
    # family of the model
    family: str
    # variant of the model (base,QCL, QCH)
    variant: str
    # base_model is the same as the name of the model
    base_model: str
    # seed_index = 1
    seed_index: int | None = None
    # seed_label = "s01" for example
    seed_label: str | None = None

    # a deterministic, user-defined sorting key for models, by using the MODEL_ORDER list.
    @property
    def sort_order(self) -> int:
        try:
            return MODEL_ORDER.index(self.base_model)
        except ValueError:
            return 999


# return a set of selected models; otherwise a empty set
def split_csv_arg(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def parse_seeded_name(name: str) -> tuple[str, int | None, str | None]:
    match = re.match(r"^(?P<base_model>.+)_s(?P<seed_index>\d{2})$", name)
    if not match:
        return name, None, None
    seed_index = int(match.group("seed_index"))
    return match.group("base_model"), seed_index, f"s{seed_index:02d}"


def infer_family_variant(base_model: str) -> tuple[str, str]:
    family_match = re.match(r"^(Alpha|Beta|Gamma)", base_model)
    family = family_match.group(1) if family_match else base_model
    if base_model.endswith("QCL"):
        variant = "QCL"
    elif base_model.endswith("QCH"):
        variant = "QCH"
    else:
        variant = "base"
    return family, variant


# discovers and constructs model specification from a metadata CSV file
def discover_model_specs(
    models_root: Path,
    metadata_path: Path,
    *,
    selected_models: set[str] | None = None,
) -> list[ModelSpec]:
    # read the metadata into a dataframe.
    metadata = pd.read_csv(metadata_path)
    # each model need to have a 'name' in model.csv
    if "name" not in metadata.columns:
        raise SystemExit(f"Model metadata is missing a name column: {metadata_path}")
    # if the model is pretained sup baseline skip it
    if "is_sup_baseline" in metadata.columns:
        metadata = metadata[
            ~metadata["is_sup_baseline"].fillna(False).astype(bool)
        ].copy()
    # create a list of model spec for each model
    specs: list[ModelSpec] = []
    # if no selected models provided use a empty set
    selected_models = selected_models or set()
    # for each model row
    for row in metadata.to_dict(orient="records"):
        # retrieve the name of the model
        name = str(row["name"])
        # assume selected_models exist (not empty), and name is not in selected models
        if selected_models and name not in selected_models:
            continue
        # otherwise locate the dorado export of the model
        dorado_dir = models_root / name / f"{name}_dorado"
        # a optional override that allows the CSV to specify a custom dorado directory when the default convention doesn't match nicely.
        if "dorado_model_dir" in row and pd.notna(row["dorado_model_dir"]):
            candidate = Path(str(row["dorado_model_dir"]))
            if candidate.exists():
                dorado_dir = candidate
        # retrieve the family of the model (Alpha, Beta, Gamma)
        family = str(row.get("family") or infer_family_variant(name)[0])
        # retrieve the variant of the model (base, QCL, QCH)
        variant = str(row.get("variant") or infer_family_variant(name)[1])
        # create a ModelSpec data container for the metadata and save to specs
        specs.append(
            ModelSpec(
                name=name,
                dorado_dir=dorado_dir.resolve(),
                family=family,
                variant=variant,
                base_model=name,
            )
        )

    # if there model selected for
    if selected_models:
        # actual models found
        discovered = {spec.name for spec in specs}
        # models that weren't found
        missing = sorted(selected_models - discovered)
        # raise a exception given there are missing models
        if missing:
            raise SystemExit(
                f"Requested models were not found in {metadata_path}: {', '.join(missing)}"
            )
    # return the sorted specs using the MODEL_ORDER defined.
    return sorted(specs, key=lambda spec: (spec.sort_order, spec.name))

# parse the model spec for seeded model just like for default models
def discover_seeded_specs(
    seeded_models_root: Path, *, selected_base_models: set[str] | None = None
) -> list[ModelSpec]:
    specs: list[ModelSpec] = []
    selected_base_models = selected_base_models or set()
    for dorado_dir in sorted(seeded_models_root.glob("*/*_dorado")):
        run_name = dorado_dir.parent.name
        base_model, seed_index, seed_label = parse_seeded_name(run_name)
        if (
            selected_base_models
            and base_model not in selected_base_models
            and run_name not in selected_base_models
        ):
            continue
        family, variant = infer_family_variant(base_model)
        specs.append(
            ModelSpec(
                name=run_name,
                dorado_dir=dorado_dir.resolve(),
                family=family,
                variant=variant,
                base_model=base_model,
                seed_index=seed_index,
                seed_label=seed_label,
            )
        )
    return sorted(
        specs, key=lambda spec: (spec.sort_order, spec.seed_index or 0, spec.name)
    )

# Prases a tensor filename into structured metadata by decoding the hierarchical pytorch module naming convention.
# It returns a flat file name like transformer_encoder.3.self_attn.q_proj.tensor into a disctionary with component, layers, and submodule information
def parse_tensor_name(filename: str) -> dict[str, object]:
    stem = filename.removesuffix(".tensor")
    component = "other"
    layer: int | None = None
    submodule = stem
    submodule_group = "Other"

    # if this ia convolution layer
    if stem.startswith("conv."):
        match = re.match(r"^conv\.(?P<layer>\d+)\.(?P<tail>.+)$", stem)
        component = "conv"
        submodule_group = "Conv"
        if match:
            # retrieve the layer index
            layer = int(match.group("layer"))
            # retrieve submodule
            submodule = match.group("tail").replace(".", "_")
    # else if this is a transformer layer
    elif stem.startswith("transformer_encoder."):
        match = re.match(r"^transformer_encoder\.(?P<layer>\d+)\.(?P<tail>.+)$", stem)
        if match:
            # retrieve layer index
            layer = int(match.group("layer"))
            # retrieve the tail
            tail = match.group("tail")
            # retrieve the submodule
            submodule = tail.replace(".", "_")
            # if this is a self-attention layer
            if tail.startswith("self_attn."):
                component = "transformer_self_attn"
                submodule_group = "Self-attention"
            # if this is a feed forward layer
            elif tail.startswith("ff."):
                component = "transformer_ffn"
                submodule_group = "Feed-forward"
            # if this is a normalisation layer or deepnorm layer
            elif tail.startswith("norm") or tail.startswith("deepnorm"):
                component = "transformer_norm_deepnorm"
                submodule_group = "Norm/deepnorm"
            # if this is anything else, just call it other
            else:
                component = "other"
                submodule_group = "Transformer other"
    # if this is a upsampling layer
    elif stem.startswith("upsample."):
        component = "upsample"
        submodule = stem.removeprefix("upsample.").replace(".", "_")
        submodule_group = "Upsample"
    # if this is a CRF head
    elif stem.startswith("crf."):
        component = "crf"
        submodule = stem.removeprefix("crf.").replace(".", "_")
        submodule_group = "CRF"

    # return the entry containing information about this layer.
    return {
        "tensor": filename,
        "tensor_stem": stem,
        "component": component,
        # for grouping and display purpose in analysis section
        "component_label": COMPONENT_LABELS.get(component, component),
        "component_order": COMPONENT_ORDER.get(component, 99),
        "layer": layer,
        "submodule": submodule,
        "submodule_group": submodule_group,
    }

# load the .tensor file into an actual pytorch Tensor object
def load_dorado_tensor(path: Path) -> torch.Tensor:
    # load the .tensor file into system memory (i.e., CPU memory, not GPU)
    # the .tensor file are actually serialized TorchScript modules
    module = torch.jit.load(str(path), map_location="cpu")
    # A torchscipt module stores its tensors as paramters/buffers. state_dict() return an OrderedDict mapping
    # paramter names to their tensor values
    state = module.state_dict()
    # the file should contain exactly one tensor
    if len(state) != 1:
        raise ValueError(f"Expected one tensor in {path}, found {len(state)}")
    # grab the single tensor from the dict, disconnect it from any computational graph, ensure it's on CPU, then cast to float32
    tensor = next(iter(state.values())).detach().cpu().float()
    # Free the torchscript module from memory
    del module
    # return the cleaned-up tensor.
    return tensor

# ensure the models being compared have layers 1:1 matching against the pretained sup model
def validate_tensor_set(model_dir: Path, baseline_names: set[str]) -> list[str]:
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Model Dorado directory does not exist: {model_dir}")
    names = {path.name for path in model_dir.glob("*.tensor")}
    missing = sorted(baseline_names - names)
    extra = sorted(names - baseline_names)
    if missing or extra:
        message = [f"Tensor set mismatch for {model_dir}"]
        if missing:
            message.append(
                f"missing: {', '.join(missing[:10])}"
                + (" ..." if len(missing) > 10 else "")
            )
        if extra:
            message.append(
                f"extra: {', '.join(extra[:10])}" + (" ..." if len(extra) > 10 else "")
            )
        raise ValueError("; ".join(message))
    return sorted(names)

# handle two edge cases
def safe_std(tensor: torch.Tensor) -> float:
    # 1. zero-division protection
    # std is undefined for a single value, calling torch.std() on a 1-element tensor returns nan
    if tensor.numel() <= 1:
        return 0.0
    # 2. population std, not sample std
    # The tensor represents the complete set of weight values, not a sample drawn from a larger distribution.
    return float(torch.std(tensor, unbiased=False).item())

# compute quantile stastistic for a tensor and returns them as a named dictionary
def quantile_values(values: torch.Tensor) -> dict[str, float]:
    # if a tensor has no element, return NaN for all quantile bins
    if values.numel() == 0:
        return {f"q{int(q * 1000):03d}": float("nan") for q in QUANTILES}
    # calculate the values at each percentile defined in the QUANTILES constant
    quantiles = torch.tensor(QUANTILES, dtype=torch.float32)
    result = torch.quantile(values.float(), quantiles)
    # return the quantile indx and the value
    return {
        f"q{int(q * 1000):03d}": float(value.item())
        for q, value in zip(QUANTILES, result, strict=True)
    }

# accumulates running toals into a shared dictionary.
# called once per tensor, and all tensors belonging to the same aggregation group
# feed into the same aggregate dictionary
def add_to_aggregate(
    aggregate: dict[str, object],
    *,
    tuned: torch.Tensor,
    baseline: torch.Tensor,
    delta: torch.Tensor,
    abs_z: torch.Tensor,
) -> None:
    aggregate["numel"] = int(aggregate.get("numel", 0)) + int(delta.numel())
    aggregate["delta_l1_sum"] = float(aggregate.get("delta_l1_sum", 0.0)) + float(
        delta.abs().sum().item()
    )
    aggregate["delta_sq_sum"] = float(aggregate.get("delta_sq_sum", 0.0)) + float(
        torch.dot(delta, delta).item()
    )
    aggregate["baseline_sq_sum"] = float(aggregate.get("baseline_sq_sum", 0.0)) + float(
        torch.dot(baseline, baseline).item()
    )
    aggregate["tuned_sq_sum"] = float(aggregate.get("tuned_sq_sum", 0.0)) + float(
        torch.dot(tuned, tuned).item()
    )
    aggregate["dot_sum"] = float(aggregate.get("dot_sum", 0.0)) + float(
        torch.dot(tuned, baseline).item()
    )
    aggregate["sign_flip_count"] = int(aggregate.get("sign_flip_count", 0)) + int(
        ((tuned * baseline) < 0).sum().item()
    )
    for threshold in TAIL_THRESHOLDS:
        key = f"tail_abs_delta_z_gt_{str(threshold).replace('.', 'p')}_count"
        aggregate[key] = int(aggregate.get(key, 0)) + int(
            (abs_z > threshold).sum().item()
        )
    counts, _ = np.histogram(abs_z.numpy(), bins=HISTOGRAM_EDGES)
    aggregate["hist_counts"] = (
        np.asarray(
            aggregate.get(
                "hist_counts", np.zeros(len(HISTOGRAM_EDGES) - 1, dtype=np.int64)
            )
        )
        + counts
    )

# converts accumulated sums into a final per-group static record.
def aggregate_record(
    # component, layer or model name etc
    base: dict[str, object],
    # the running total dictionary filled with accumulated record
    aggregate: dict[str, object],
    *,
    # total delta energy for the entire model, needed for energy share %
    model_delta_sq_sum: float | None = None,
) -> dict[str, object]:
    # retrieve the values from dictionary
    numel = int(aggregate.get("numel", 0))
    delta_sq = float(aggregate.get("delta_sq_sum", 0.0))
    baseline_sq = float(aggregate.get("baseline_sq_sum", 0.0))
    tuned_sq = float(aggregate.get("tuned_sq_sum", 0.0))
    dot = float(aggregate.get("dot_sum", 0.0))
    record = dict(base)
    record.update(
        {
            "numel": numel,
            "delta_l1_sum": float(aggregate.get("delta_l1_sum", 0.0)),
            # euclidean norm of the total weight change
            "delta_l2_norm": math.sqrt(delta_sq),
            # magnitude of sup weights
            "baseline_l2_norm": math.sqrt(baseline_sq),
            # mangitude of fine-tuned model weights
            "tuned_l2_norm": math.sqrt(tuned_sq),
            # how large the change is relative to the baseline
            "relative_l2": math.sqrt(delta_sq) / (math.sqrt(baseline_sq) + EPS),
            # directional alignment between tuned and baseline (1.0 = same direction)
            "cosine_similarity": dot
            / ((math.sqrt(tuned_sq) * math.sqrt(baseline_sq)) + EPS),
            # propotino of weights that changes sign
            "sign_flip_fraction": int(aggregate.get("sign_flip_count", 0)) / numel
            if numel
            else float("nan"),
        }
    )

    for threshold in TAIL_THRESHOLDS:
        suffix = str(threshold).replace(".", "p")
        count_key = f"tail_abs_delta_z_gt_{suffix}_count"
        count = int(aggregate.get(count_key, 0))
        # row count of elements exceeding threshold X
        record[f"tail_abs_delta_z_gt_{suffix}_count"] = count
        # fraction of that row count
        record[f"tail_abs_delta_z_gt_{suffix}_fraction"] = (
            count / numel if numel else float("nan")
        )
    hist_counts = np.asarray(
        aggregate.get("hist_counts", np.zeros(len(HISTOGRAM_EDGES) - 1, dtype=np.int64))
    )

    # estimate the 99.9th percentile of the z-score distirbution from histogram bin counts
    record["p999_abs_delta_z_approx"] = histogram_quantile(hist_counts, 0.999)
    # shows what fraction of the model's total weight cahgen is concentrated in this component/layer
    if model_delta_sq_sum is not None:
        record["delta_energy_share"] = (
            delta_sq / model_delta_sq_sum if model_delta_sq_sum else float("nan")
        )
    return record

# estimates a quantile value from histogram bin counts using linear interpolation within the bin
def histogram_quantile(counts: np.ndarray, quantile: float) -> float:
    total = int(counts.sum())
    if total == 0:
        return float("nan")
    # the element rank at the desired quantile
    target = total * quantile
    cumulative = np.cumsum(counts)
    # finds which bin containing that rank
    index = int(np.searchsorted(cumulative, target, side="left"))
    # prevent overflow
    index = min(index, len(counts) - 1)
    left = HISTOGRAM_EDGES[index]
    right = HISTOGRAM_EDGES[index + 1]
    previous = cumulative[index - 1] if index > 0 else 0
    count = counts[index]
    if count <= 0:
        return float(left)
    # linear interpolation, computes how far into the bin the quantile falls, then interpolrate between left and right
    fraction = (target - previous) / count
    if not np.isfinite(right):
        return float(left)
    return float(left + (right - left) * fraction)

# converts accumulated histogram bin counts into tabular row record
# one row per bin for each CSV output
def histogram_records(
    spec: ModelSpec, aggregate_key: dict[str, object], aggregate: dict[str, object]
) -> list[dict[str, object]]:
    counts = np.asarray(
        aggregate.get("hist_counts", np.zeros(len(HISTOGRAM_EDGES) - 1, dtype=np.int64))
    )
    total = int(counts.sum())
    cumulative = np.cumsum(counts)
    records = []
    for idx, count in enumerate(counts):
        left = HISTOGRAM_EDGES[idx]
        right = HISTOGRAM_EDGES[idx + 1]
        if np.isfinite(right):
            mid = math.sqrt(max(left, EPS) * right) if left > 0 else right / 2
        else:
            mid = left
        records.append(
            {
                "name": spec.name,
                "display_name": spec.name,
                "base_model": spec.base_model,
                "family": spec.family,
                "variant": spec.variant,
                "seed_index": spec.seed_index,
                "seed_label": spec.seed_label,
                "sort_order": spec.sort_order,
                **aggregate_key,
                "bin_index": idx,
                "bin_left": left,
                "bin_right": right,
                "bin_mid": mid,
                "count": int(count),
                "cumulative_count": int(cumulative[idx]),
                "cumulative_fraction": float(cumulative[idx] / total)
                if total
                else float("nan"),
                "total_count": total,
            }
        )
    return records

# pretained sup baseline ModelSpec record
def base_model_record(spec: ModelSpec) -> dict[str, object]:
    return {
        "name": spec.name,
        "display_name": spec.name,
        "base_model": spec.base_model,
        "family": spec.family,
        "variant": spec.variant,
        "seed_index": spec.seed_index,
        "seed_label": spec.seed_label,
        "sort_order": spec.sort_order,
    }

# core analysis and reporting function
# performs per-tensor statistical comparison between each tuned model and a baseline model
# then writes the results into six CSV files.
def write_outputs(
    specs: Sequence[ModelSpec],
    baseline_dorado_dir: Path,
    out_dir: Path,
    *,
    z_eps: float,
) -> dict[str, Path]:
    if not specs:
        raise SystemExit("No models selected for weight summarisation")
    if not baseline_dorado_dir.is_dir():
        raise SystemExit(
            f"Baseline Dorado directory does not exist: {baseline_dorado_dir}"
        )

    # 1. Baseline loading
    # loads all .tensor files from the pretained sup baseline directory into memory, these serve as the
    # reference point for all comparisons
    baseline_paths = sorted(baseline_dorado_dir.glob("*.tensor"))
    baseline_names = {path.name for path in baseline_paths}
    if not baseline_names:
        raise SystemExit(
            f"No .tensor files found in baseline Dorado directory: {baseline_dorado_dir}"
        )

    print(f"Loading {len(baseline_names)} baseline tensors from {baseline_dorado_dir}")
    # in-memory lookup table of all baseline tensor key by the filename (e.g., encoder.layer_0.weight.tensor)
    baseline_tensors = {path.name: load_dorado_tensor(path) for path in baseline_paths}

    tensor_rows: list[dict[str, object]] = []
    aggregate_map: dict[tuple[object, ...], dict[str, object]] = {}
    histogram_rows: list[dict[str, object]] = []

    # for each model required
    for spec in specs:
        # validate the layers/tensor in each fine-tuned model against pretained sup baseline
        names = validate_tensor_set(spec.dorado_dir, baseline_names)
        # number of tensor layers
        print(f"Summarising {spec.name}: {len(names)} tensors")
        # for each tensor file name
        for tensor_name in names:
            # parse the information encapsulated by the tensor layer
            metadata = parse_tensor_name(tensor_name)
            # load the actual raw tensor data for the layer
            tuned_raw = load_dorado_tensor(spec.dorado_dir / tensor_name)
            # the baseline raw tensor data for the same lyaer
            baseline_raw = baseline_tensors[tensor_name]
            # both layer must have the same dimension given that there are not architecture change
            # during fine-tuning
            if tuple(tuned_raw.shape) != tuple(baseline_raw.shape):
                raise ValueError(
                    f"Shape mismatch for {spec.name}/{tensor_name}: "
                    f"{tuple(tuned_raw.shape)} vs baseline {tuple(baseline_raw.shape)}"
                )
            # create human-readable shape string for the CSV output
            shape_text = "x".join(str(dim) for dim in tuned_raw.shape) or "scalar"
            # collpase any multidimentioanl tensor into a single vector, necessary for all subsequent statistics (means, standard deviations, dot products, element-wise delta, etc)
            tuned = tuned_raw.reshape(-1)
            baseline = baseline_raw.reshape(-1)

            # element-wise delta across flattened paramters
            delta = tuned - baseline
            # absolute delt across falttened paramters
            abs_delta = delta.abs()
            # internal spread of the baseline weights within this tensor, a large baseline_std means
            # the baseline has a wide range of weight values; small one means the weights are clustered tightly.
            baseline_std = safe_std(baseline)
            # change expressed in units of the sup baselins's own standard deviation
            # "Is this change large relative to the natural variation within this tensor"
            # abs_z ~ 0 (weights barely moved)
            # abs_z ~ 1 (weights shifted by rouhgly one standard deviation of the baseline)
            # abs_z > 1 (weights shifted dramatically)
            abs_z = abs_delta / (baseline_std + z_eps)
            # sqaured L2 norms (sum of squared elements) of each vector.
            # represent the total magnitude of each vector
            # TThese are intermediate values used later to compute:
                #  delta_l2_norm - overall distance betweenm tuned and baseline
                #  relative_l2 - change as fraction of baseline magnitude
                #  cosine_similarity - derived from the dot product below
            # delta_sq total amount of change across all weights
            delta_sq = float(torch.dot(delta, delta).item())
            # baseline_sq total energy of the pretained weights
            baseline_sq = float(torch.dot(baseline, baseline).item())
            # tuned_sq total energy of the fine-tuned weights
            tuned_sq = float(torch.dot(tuned, tuned).item())

            # alignment measure between tuned and baseline vector, used to compute cosine similarity
            dot = float(torch.dot(tuned, baseline).item())
            # total number of paramters in the tensor, used as a denominator for averageds and to guard against empty tensors
            numel = int(delta.numel())
            # distritbution summary of the absolute weight changes, not just the mean/std, so skewed distribution or outliers can be spotted.
            abs_delta_quantiles = quantile_values(abs_delta)
            abs_z_quantiles = quantile_values(abs_z)

            # create a row to capture the comparison information across this layer.
            row = {
                **base_model_record(spec),
                **metadata,
                "shape": shape_text,
                "numel": numel,
                "baseline_mean": float(baseline.mean().item()),
                "baseline_std": baseline_std,
                "tuned_mean": float(tuned.mean().item()),
                "tuned_std": safe_std(tuned),
                "delta_mean": float(delta.mean().item()),
                "delta_std": safe_std(delta),
                "delta_l1_sum": float(abs_delta.sum().item()),
                "delta_l2_norm": math.sqrt(delta_sq),
                "delta_linf_norm": float(abs_delta.max().item())
                if numel
                else float("nan"),
                "baseline_l2_norm": math.sqrt(baseline_sq),
                "tuned_l2_norm": math.sqrt(tuned_sq),
                "relative_l2": math.sqrt(delta_sq) / (math.sqrt(baseline_sq) + EPS),
                "cosine_similarity": dot
                / ((math.sqrt(tuned_sq) * math.sqrt(baseline_sq)) + EPS),
                "sign_flip_fraction": float(
                    ((tuned * baseline) < 0).float().mean().item()
                )
                if numel
                else float("nan"),
            }
            for suffix, value in abs_delta_quantiles.items():
                row[f"abs_delta_{suffix}"] = value
            for suffix, value in abs_z_quantiles.items():
                row[f"abs_delta_z_{suffix}"] = value
            for threshold in TAIL_THRESHOLDS:
                suffix = str(threshold).replace(".", "p")
                row[f"tail_abs_delta_z_gt_{suffix}_fraction"] = float(
                    (abs_z > threshold).float().mean().item()
                )
            # append the row as a single entry to the list.
            tensor_rows.append(row)

            # Accumulate running totals into three parallel aggregation buckets simultaneously, so
            # later statistic can be computed at different granularity levels without re-loading tensors
            # Bucket 1 - Component-level
            # groups all tensors belonging to the same component type (e.g., all conv tensors, all crf tensors)
            # later produces rows in weight_component_stats.csv
            component_key = (spec.name, "component", metadata["component"])
            component_agg = aggregate_map.setdefault(component_key, {})
            add_to_aggregate(
                component_agg, tuned=tuned, baseline=baseline, delta=delta, abs_z=abs_z
            )
            # Bucket 2 - Model-wide
            # group all tensor belonging to the same model
            all_key = (spec.name, "all", "all")
            all_agg = aggregate_map.setdefault(all_key, {})
            add_to_aggregate(
                all_agg, tuned=tuned, baseline=baseline, delta=delta, abs_z=abs_z
            )
            # Bucket 3 - Layer-wise
            # group all all tensor belonging to the same layer
            layer_key = (
                spec.name,
                "layer",
                metadata["component"],
                metadata["layer"],
                metadata["submodule_group"],
            )
            layer_agg = aggregate_map.setdefault(layer_key, {})
            add_to_aggregate(
                layer_agg, tuned=tuned, baseline=baseline, delta=delta, abs_z=abs_z
            )

            del tuned_raw, tuned, delta, abs_delta, abs_z
            gc.collect()

    component_rows: list[dict[str, object]] = []
    layer_rows: list[dict[str, object]] = []
    model_summary_rows: list[dict[str, object]] = []
    spec_lookup = {spec.name: spec for spec in specs}

    # compute total model change
    # it is needed later to compute delta_energy_share (what fraction of total change each component/layer contributes)
    model_delta_sq = {
        name: float(aggregate_map[(name, "all", "all")].get("delta_sq_sum", 0.0))
        for name in spec_lookup
        if (name, "all", "all") in aggregate_map
    }

    # iterate over every entry in aggregate_map and dispatches based on the key type
    for key, aggregate in aggregate_map.items():
        spec = spec_lookup[str(key[0])]
        # if the entry is for a component
        if key[1] == "component":
            # retrieve the component
            component = str(key[2])
            # convert accumulated sums into per-component statistics
            record = aggregate_record(
                {
                    **base_model_record(spec),
                    "component": component,
                    "component_label": COMPONENT_LABELS.get(component, component),
                    "component_order": COMPONENT_ORDER.get(component, 99),
                },
                aggregate,
                model_delta_sq_sum=model_delta_sq.get(spec.name),
            )
            component_rows.append(record)
            # generate histogram bin rows
            histogram_rows.extend(
                histogram_records(
                    spec,
                    {
                        "aggregate": "component",
                        "component": component,
                        "component_label": COMPONENT_LABELS.get(component, component),
                    },
                    aggregate,
                )
            )
        # aggregate record for models
        elif key[1] == "all":
            record = aggregate_record(
                {
                    **base_model_record(spec),
                    "component": "all",
                    "component_label": "All tensors",
                    "component_order": -1,
                },
                aggregate,
            )
            model_summary_rows.append(record)
            histogram_rows.extend(
                histogram_records(
                    spec,
                    {
                        "aggregate": "all",
                        "component": "all",
                        "component_label": "All tensors",
                    },
                    aggregate,
                )
            )
        # aggregate record for layers
        elif key[1] == "layer":
            component = str(key[2])
            layer = key[3]
            submodule_group = str(key[4])
            layer_rows.append(
                aggregate_record(
                    {
                        **base_model_record(spec),
                        "component": component,
                        "component_label": COMPONENT_LABELS.get(component, component),
                        "component_order": COMPONENT_ORDER.get(component, 99),
                        "layer": layer,
                        "submodule_group": submodule_group,
                    },
                    aggregate,
                    model_delta_sq_sum=model_delta_sq.get(spec.name),
                )
            )

    # convert row lists into dataframes
    component_df = pd.DataFrame.from_records(component_rows)
    model_summary_df = pd.DataFrame.from_records(model_summary_rows)
    tensor_df = pd.DataFrame.from_records(tensor_rows)
    layer_df = pd.DataFrame.from_records(layer_rows)
    histogram_df = pd.DataFrame.from_records(histogram_rows)

    # add per component energy share to model summary
    if not component_df.empty:
        pivot = component_df.pivot(
            index="name", columns="component", values="delta_energy_share"
        )
        for component in COMPONENT_ORDER:
            if component in pivot.columns:
                model_summary_df[f"{component}_delta_energy_share"] = model_summary_df[
                    "name"
                ].map(pivot[component])
        conv_share = model_summary_df.get("conv_delta_energy_share", 0)
        crf_share = model_summary_df.get("crf_delta_energy_share", 0)
        model_summary_df["conv_crf_delta_energy_share"] = conv_share.fillna(
            0
        ) + crf_share.fillna(0)

    # top-20 most changed tensors per model
    top_tensor_df = (
        # rank tensors by relative_l2 (noramlised delta magnitude) within each model
        tensor_df.sort_values(["name", "relative_l2"], ascending=[True, False])
        .groupby("name", sort=False)
        .head(20)
        .copy()
    )
    top_tensor_df["relative_l2_rank"] = top_tensor_df.groupby("name")[
        "relative_l2"
    ].rank(method="first", ascending=False)

    # sort all dataframe using sort_order and name
    # such that models appear as defined by MODEL_ORDER sequence
    # within each model, rows are ordered by component hierarchy -> layers -> tensor name
    sort_columns = ["sort_order", "name"]
    tensor_df = tensor_df.sort_values(
        sort_columns + ["component_order", "layer", "submodule", "tensor"],
        na_position="last",
    )
    component_df = component_df.sort_values(
        sort_columns + ["component_order", "component"]
    )
    layer_df = layer_df.sort_values(
        sort_columns + ["component_order", "layer", "submodule_group"],
        na_position="last",
    )
    model_summary_df = model_summary_df.sort_values(sort_columns)
    histogram_df = histogram_df.sort_values(
        sort_columns + ["aggregate", "component", "bin_index"]
    )
    top_tensor_df = top_tensor_df.sort_values(sort_columns + ["relative_l2_rank"])

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "tensor": out_dir / "weight_tensor_stats.csv",
        "component": out_dir / "weight_component_stats.csv",
        "layer": out_dir / "weight_layer_stats.csv",
        "histogram": out_dir / "weight_delta_histogram.csv",
        "model_summary": out_dir / "weight_model_summary.csv",
        "top_tensors": out_dir / "weight_top_tensors.csv",
    }
    # creates the output directory containing six CSV files at different granularity
    # Per-tensor comparison statistics
    tensor_df.to_csv(outputs["tensor"], index=False)
    # Component-level aggeegates with energy share
    component_df.to_csv(outputs["component"], index=False)
    # Per-layer aggregates with energy share
    layer_df.to_csv(outputs["layer"], index=False)
    # Z-score distribution histograms
    histogram_df.to_csv(outputs["histogram"], index=False)
    # One row per model + pivoted energy shares
    model_summary_df.to_csv(outputs["model_summary"], index=False)
    # Top 20 most-changed tensors per model.
    top_tensor_df.to_csv(outputs["top_tensors"], index=False)
    return outputs

# print out where each tables are being saved
def print_outputs(outputs: dict[str, Path]) -> None:
    for label, path in outputs.items():
        print(f"{label}: {path}")


# The argument parser defines the commandline input for this utility
def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    # The repo root containing the models and table directories.
    parser.add_argument("--repo-root", type=Path, default=None)
    # The model root containing subdirectories for each model, each with a dorado export subdirectory.
    parser.add_argument("--models-root", type=Path, default=None)
    # The baseline model directory containing .tensor files for the pretained sup baseline.
    parser.add_argument("--baseline-dorado-dir", type=Path, required=True)
    # The model metadata CSV file containing at least a "name" column, and optionally "family" and "variant" columns. Used to discover models and assign them to families/variants for analysis.
    parser.add_argument("--model-metadata", type=Path, default=None)
    # The output directory where the resulting CSV tables will be saved. Defaults to "tables/weights" under the repo root.
    parser.add_argument("--out-dir", type=Path, default=None)
    # Optional comma-separated list of model names to include in the analysis. If not provided, all models found in the metadata will be included.
    parser.add_argument(
        "--models", default=None, help="Comma-separated current model names to include"
    )
    # Epsilon value added to the baseline tensor standard deviation when calculating abs_delta_z to avoid division by zero. Default is 1e-8.
    parser.add_argument(
        "--z-eps",
        type=float,
        default=1e-8,
        help="Epsilon added to baseline tensor std for abs_delta_z",
    )
    # Optional arguments for including seeded models in the analysis. If --seeded-models-root is provided, it will look for subdirectories containing dorado exports of seeded runs.
    parser.add_argument("--seeded-models-root", type=Path, default=None)
    # Optional comma-separated list of seeded run names or base model names to include when --seeded-models-root is provided. If not provided, all seeded runs found under the seeded models root will be included.
    parser.add_argument("--seeded-out-dir", type=Path, default=None)
    parser.add_argument(
        "--seeded-models",
        default=None,
        help="Comma-separated seeded run names or base model names to include",
    )
    return parser.parse_args(argv)


# The driver function orchestrates the workflow:
# 1. Parse command-line arguments
# 2. Resolve paths for repo root, models root, metadata, output directory
# 3. Discover model specifications based on metadata and selected models
# 4. Write outputs for the discovered models
# 5. If seeded models are included, discover seeded specs and write outputs for them as well
def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    # repo_root is the root directory of the COMPX593-Thesis repository
    repo_root = args.repo_root.resolve()
    # model_root is the models directory containing subdirectories for each model
    models_root = args.models_root.resolve()
    # model_metadata is the model.csv file
    metadata_path = args.model_metadata.resolve()
    # output directory where computed metrics and deltas are saved.
    out_dir = args.out_dir.resolve()
    # pretaine sup model directory with .tensor weights.
    baseline_dorado_dir = args.baseline_dorado_dir.resolve()

    # a set of selected models, or empty set
    selected_models = split_csv_arg(args.models)

    # a list of model metadata retrieved from model.csv
    specs = discover_model_specs(
        models_root, metadata_path, selected_models=selected_models
    )

    # compute then print out where the output .csv files are being saved.
    print_outputs(write_outputs(specs, baseline_dorado_dir, out_dir, z_eps=args.z_eps))

    # if random seeded models are the target
    if args.seeded_models_root is not None:
        seeded_out_dir = (args.seeded_out_dir or out_dir / "seeded").resolve()
        # select all seeded models
        selected_seeded = split_csv_arg(args.seeded_models)
        seeded_specs = discover_seeded_specs(
            args.seeded_models_root.resolve(), selected_base_models=selected_seeded
        )
        # print a lager table with 10 seeded models per model type (i.e., Alpha_s01, Alpha_s02)
        print_outputs(
            write_outputs(
                seeded_specs, baseline_dorado_dir, seeded_out_dir, z_eps=args.z_eps
            )
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
