"""Deterministic rank-aware ICA and EOG proposal logic for DEMI EEG.

This module prepares the separate approximately 1--45 Hz ICA-fitting branch,
reuses the accepted reference/interpolation surface, estimates EEG rank, fits
extended Infomax with a fixed seed, scores HEO and VEO separately, and applies
the fixed ordinary-file routing rule. ID 86 is routed to its predeclared
component-review boundary without automatic exclusion or post-ICA output.

All operations are in memory. ICA objects and continuous derivatives are saved
by the storage layer only after routing is complete. This module does not build
epochs, run AutoReject, compute CSD, or make participant inclusion decisions.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import mne
import numpy as np

from .contracts import EEG_TARGET_CHANNELS, EOG_CHANNELS
from .stages import (
    ObjectiveStop,
    apply_fir_filter,
    apply_line_noise_branch,
    apply_reference,
    interpolate_global_bads,
    validate_finite_continuous,
)


def prepare_ica_fitting_copy(
    typed_montaged_source: mne.io.BaseRaw,
    *,
    global_bads: Sequence[str],
    reference: Mapping[str, Any],
    config: Mapping[str, Any],
) -> tuple[mne.io.BaseRaw, dict[str, Any]]:
    """Create the independent ICA branch with the accepted shared surface.

    Args:
        typed_montaged_source: Full-session source branch before filtering.
        global_bads: Accepted primary scalp bad channels.
        reference: Exact reference-source calculation used by the analysis branch.
        config: Validated complete production configuration.

    Returns:
        Prepared ICA-fitting ``Raw`` and complete preparation evidence.

    Side effects:
        Allocates and mutates an independent in-memory copy. Writes nothing.
    """

    ica_settings = config["ica"]
    branch = typed_montaged_source.copy()
    line_evidence = apply_line_noise_branch(branch, config)
    filter_evidence = apply_fir_filter(
        branch,
        high_pass_hz=float(ica_settings["high_pass_hz"]),
        low_pass_hz=float(ica_settings["low_pass_hz"]),
        settings=config["analysis_filter"],
    )
    reference_evidence = apply_reference(branch, reference)
    interpolation_evidence = interpolate_global_bads(branch, global_bads, config)
    validation_evidence = validate_finite_continuous(
        branch, float(config["validation"]["finite_scan_chunk_seconds"])
    )
    return branch, {
        "line_noise": line_evidence,
        "filter": filter_evidence,
        "reference_sources_match_analysis_branch": reference_evidence["reference_sources"]
        == reference["reference_sources"],
        "reference_targets": reference_evidence["reference_targets"],
        "interpolation_candidates_match_analysis_branch": interpolation_evidence[
            "interpolation_candidates"
        ]
        == sorted(global_bads),
        "interpolation": interpolation_evidence,
        "validation": validation_evidence,
        "full_session_sample_count": int(branch.n_times),
        "sampling_frequency_hz": float(branch.info["sfreq"]),
    }


def estimate_eeg_rank(raw: mne.io.BaseRaw, config: Mapping[str, Any]) -> dict[str, Any]:
    """Estimate and validate EEG rank after reference and interpolation.

    Args:
        raw: Prepared ICA-fitting branch.
        config: Validated complete production configuration. The fixed relative
            tolerance recognizes the exact average-reference dependency while
            remaining unchanged across recordings.

    Returns:
        Full MNE rank mapping plus the validated integer EEG rank.

    Raises:
        ObjectiveStop: If rank is missing, noninteger, nonpositive, or greater
            than the expected maximum after explicit average referencing.

    Side effects:
        Reads signal data for rank estimation. Writes nothing.
    """

    tolerance = float(config["ica"]["rank_tolerance"])
    tolerance_kind = config["ica"]["rank_tolerance_kind"]
    rank_map = mne.compute_rank(
        raw,
        rank=None,
        proj=False,
        tol=tolerance,
        tol_kind=tolerance_kind,
        on_rank_mismatch="raise",
        verbose="ERROR",
    )
    rank_value = rank_map.get("eeg")
    maximum_after_reference = len(EEG_TARGET_CHANNELS) - 1
    valid = (
        isinstance(rank_value, (int, np.integer))
        and int(rank_value) > 0
        and int(rank_value) <= maximum_after_reference
    )
    evidence = {
        "method": "mne.compute_rank",
        "tolerance": tolerance,
        "tolerance_kind": tolerance_kind,
        "tolerance_fixed_across_recordings": True,
        "rank_by_channel_type": {key: int(value) for key, value in rank_map.items()},
        "estimated_eeg_rank": int(rank_value) if isinstance(rank_value, (int, np.integer)) else None,
        "eeg_channel_count": len(EEG_TARGET_CHANNELS),
        "maximum_valid_rank_after_average_reference": maximum_after_reference,
        "valid": bool(valid),
    }
    if not valid:
        raise ObjectiveStop(
            "rank_estimation",
            "invalid_ica_rank",
            f"Invalid post-reference/interpolation EEG rank: {rank_value!r}.",
            evidence,
        )
    return evidence


def fit_ica(raw: mne.io.BaseRaw, rank: Mapping[str, Any], config: Mapping[str, Any]) -> tuple[mne.preprocessing.ICA, dict[str, Any]]:
    """Fit deterministic extended Infomax using the estimated EEG rank.

    Args:
        raw: Prepared ICA-fitting branch.
        rank: Validated output from :func:`estimate_eeg_rank`.
        config: Validated complete production configuration.

    Returns:
        Fitted MNE ICA object and method/fit evidence.

    Side effects:
        Performs CPU-intensive fitting in memory. Writes nothing.
    """

    settings = config["ica"]
    estimated_rank = int(rank["estimated_eeg_rank"])
    ica = mne.preprocessing.ICA(
        n_components=estimated_rank,
        method=settings["method"],
        fit_params={"extended": bool(settings["extended"])},
        random_state=int(settings["random_seed"]),
        max_iter="auto",
    )
    ica.fit(
        raw,
        picks=list(EEG_TARGET_CHANNELS),
        decim=int(settings["fit_decim"]),
        reject_by_annotation=True,
        verbose="ERROR",
    )
    if int(ica.n_components_) != estimated_rank:
        raise ObjectiveStop(
            "ica_fit",
            "ica_rank_component_conflict",
            f"ICA fit {ica.n_components_} components for estimated rank {estimated_rank}.",
        )
    return ica, {
        "method": settings["method"],
        "fit_params": {"extended": bool(settings["extended"])},
        "random_seed": int(settings["random_seed"]),
        "estimated_eeg_rank": estimated_rank,
        "n_components_fitted": int(ica.n_components_),
        "fit_decim": int(settings["fit_decim"]),
        "effective_fit_sampling_frequency_hz": float(raw.info["sfreq"])
        / int(settings["fit_decim"]),
        "fit_sample_count_before_decimation": int(raw.n_times),
        "fit_sample_count_approximate_after_decimation": int(
            np.ceil(raw.n_times / int(settings["fit_decim"]))
        ),
        "fit_channel_names": list(ica.ch_names),
        "n_iterations": int(ica.n_iter_),
        "current_fit": ica.current_fit,
    }


def score_eog_components(
    ica: mne.preprocessing.ICA,
    raw: mne.io.BaseRaw,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Score HEO and VEO independently and preserve one row per component.

    Args:
        ica: Fitted ICA solution.
        raw: The exact fitting branch.
        config: Validated complete production configuration.

    Returns:
        Per-EOG flagged indices and per-component score rows.

    Raises:
        ObjectiveStop: If score arrays are missing, mismatched, nonfinite, or
            produce an out-of-range component index.

    Side effects:
        Filters temporary signals internally through MNE's EOG scorer. Writes
        nothing and does not change the ICA exclusion list.
    """

    settings = config["ica"]
    expected_components = int(ica.n_components_)
    candidates_by_eog: dict[str, list[int]] = {}
    scores_by_eog: dict[str, list[float]] = {}
    for channel in EOG_CHANNELS:
        indices, scores = ica.find_bads_eog(
            raw,
            ch_name=channel,
            threshold=float(settings["eog_threshold"]),
            l_freq=float(settings["eog_score_low_hz"]),
            h_freq=float(settings["eog_score_high_hz"]),
            reject_by_annotation=True,
            measure=settings["eog_measure"],
            verbose="ERROR",
        )
        score_array = np.asarray(scores, dtype=float)
        invalid_indices = sorted(
            index for index in indices if not isinstance(index, (int, np.integer)) or not 0 <= int(index) < expected_components
        )
        if score_array.shape != (expected_components,) or not np.isfinite(score_array).all() or invalid_indices:
            evidence = {
                "channel": channel,
                "score_shape": list(score_array.shape),
                "expected_components": expected_components,
                "nonfinite_score_count": int(np.size(score_array) - np.isfinite(score_array).sum()),
                "invalid_indices": invalid_indices,
            }
            raise ObjectiveStop(
                "component_scoring",
                "ambiguous_automatic_ica_rule",
                f"Invalid EOG component evidence for {channel}: {evidence}",
                evidence,
            )
        candidates_by_eog[channel] = sorted(set(int(index) for index in indices))
        scores_by_eog[channel] = [float(value) for value in score_array]

    rows: list[dict[str, Any]] = []
    for component in range(expected_components):
        support = [
            channel for channel in EOG_CHANNELS if component in candidates_by_eog[channel]
        ]
        heo_score = scores_by_eog["HEO"][component]
        veo_score = scores_by_eog["VEO"][component]
        rows.append(
            {
                "component_index": component,
                "method": settings["method"],
                "random_seed": int(settings["random_seed"]),
                "estimated_rank": expected_components,
                "heo_score": heo_score,
                "veo_score": veo_score,
                "maximum_absolute_eog_score": max(abs(heo_score), abs(veo_score)),
                "eog_support_channels": support,
                "automatic_proposal": bool(support),
                "final_action": "pending_decision",
                "reason": "pending_exception_or_continuation_decision",
            }
        )
    return {
        "scoring_method": "mne.preprocessing.ICA.find_bads_eog",
        "measure": settings["eog_measure"],
        "threshold": float(settings["eog_threshold"]),
        "score_filter_hz": [
            float(settings["eog_score_low_hz"]),
            float(settings["eog_score_high_hz"]),
        ],
        "candidates_by_eog": candidates_by_eog,
        "component_rows": rows,
    }


def historical_eog_component_order(scores: Mapping[str, Any]) -> list[int]:
    """Reproduce historical MNE multi-EOG ranking and deduplication.

    The historical executable called ``find_bads_eog`` once without a channel
    override. MNE concatenated HEO and VEO findings, ranked occurrences by
    descending absolute score, and removed duplicate component indices while
    preserving that rank. The active scorer keeps channel evidence separate;
    this function reconstructs the historical combined result without a cap.
    """

    rows = {int(row["component_index"]): row for row in scores["component_rows"]}
    score_keys = {"HEO": "heo_score", "VEO": "veo_score"}
    occurrences: list[int] = []
    occurrence_scores: list[float] = []
    for channel in EOG_CHANNELS:
        for component in scores["candidates_by_eog"][channel]:
            component = int(component)
            occurrences.append(component)
            occurrence_scores.append(float(rows[component][score_keys[channel]]))
    if not occurrences:
        return []
    ranked = np.asarray(occurrences, dtype=int)[
        np.abs(np.asarray(occurrence_scores, dtype=float)).argsort()[::-1]
    ].tolist()
    selected: list[int] = []
    for component in ranked:
        if component not in selected:
            selected.append(component)
    return selected


def propose_eog_components(scores: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    """Recover the complete historical EOG-component selection.

    Args:
        scores: Output from :func:`score_eog_components`.
        config: Validated complete production configuration.

    Returns:
        All candidates and the score-ranked, deduplicated final list. This
        stage deliberately does not apply or route the proposal.

    Side effects:
        None.
    """

    del config
    proposed = historical_eog_component_order(scores)
    return {
        "proposal_rule": "historical_mne_find_bads_eog_all_ranked_deduplicated_components",
        "detected_candidates_by_eog": {
            channel: [int(index) for index in scores["candidates_by_eog"][channel]]
            for channel in EOG_CHANNELS
        },
        "proposed_components": proposed,
        "proposal_count": len(proposed),
        "all_returned_components_selected": True,
        "truncated_or_capped": False,
        "ordering": "descending_absolute_eog_score_across_HEO_then_VEO_with_stable_deduplication",
        "zero_component_path": not proposed,
    }


def decide_component_route(
    recording_id: int,
    scores: Mapping[str, Any],
    proposal: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Route an ordinary file or the predeclared ID-86 exception.

    Args:
        recording_id: EEG recording identifier parsed from the EDF name.
        scores: Component score evidence whose rows will receive final actions.
        proposal: Output of :func:`propose_eog_components`.
        config: Validated complete production configuration.

    Returns:
        Route, final exclusions, reason, and finalized component rows.

    Side effects:
        Mutates only copied row dictionaries returned in the result.
    """

    proposed = list(proposal["proposed_components"])
    rows = [dict(row) for row in scores["component_rows"]]
    exception_ids = set(config["ica"]["component_review_exception_ids"])
    if recording_id in exception_ids:
        route = "stop_for_component_review"
        exclusions: list[int] = []
        reason = "predeclared_id86_component_review_boundary"
        for row in rows:
            row["final_action"] = "review_required_no_automatic_application"
            row["reason"] = reason
    elif not proposed:
        route = "stop_for_exceptional_review"
        exclusions = []
        reason = "historical_no_eog_components_detected"
        for row in rows:
            row["final_action"] = "not_applied_exceptional_review"
            row["reason"] = reason
    else:
        route = "continue_automatic"
        exclusions = proposed
        reason = "historical_all_find_bads_eog_components"
        for row in rows:
            if row["component_index"] in exclusions:
                row["final_action"] = "exclude"
                row["reason"] = reason
            else:
                row["final_action"] = "retain"
                row["reason"] = "not_proposed_by_fixed_eog_rule"
    return {
        "recording_id": int(recording_id),
        "route": route,
        "automatic_application_authorized": route == "continue_automatic",
        "proposed_components": proposed,
        "final_exclusions": exclusions,
        "reason": reason,
        "component_rows": rows,
    }


def apply_ica_exclusions(
    ica: mne.preprocessing.ICA,
    analysis_raw: mne.io.BaseRaw,
    route: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply authorized exclusions to the analysis branch in place.

    Args:
        ica: Fitted solution.
        analysis_raw: Accepted 0.5--45 Hz referenced/interpolated branch.
        route: Output of :func:`decide_component_route`.

    Returns:
        Application status and exclusions.

    Raises:
        ObjectiveStop: If automatic application is not authorized.

    Side effects:
        Mutates ``analysis_raw`` only when at least one exclusion is authorized;
        always records the authorized list in ``ica.exclude``.
    """

    if not route["automatic_application_authorized"]:
        raise ObjectiveStop(
            "exception_or_continuation_decision",
            route["reason"],
            f"Recording stopped at component-review boundary: {route['reason']}",
            dict(route),
        )
    exclusions = list(route["final_exclusions"])
    ica.exclude = exclusions
    if exclusions:
        ica.apply(analysis_raw, exclude=exclusions, verbose="ERROR")
    return {
        "authorized": True,
        "final_exclusions": exclusions,
        "exclusion_count": len(exclusions),
        "ica_application_performed": bool(exclusions),
        "zero_component_branch_left_analysis_samples_unchanged": not exclusions,
    }
