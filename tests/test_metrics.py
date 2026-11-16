import numpy as np
import pandas as pd
import pytest

from scbatchflow.metrics import evaluate, kbet_acceptance, lisi_scores


def embedding(batch_shift, n=300, seed=0):
    """Three well-separated cell types in two batches; batch 2 is shifted by `batch_shift`."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(["a", "b", "c"], n // 3)
    batch = np.tile(["b1", "b2"], n // 2)
    centres = {"a": (0, 0), "b": (12, 0), "c": (0, 12)}
    emb = np.array([centres[lab] for lab in labels], dtype=float) + rng.normal(0, 1, (n, 2))
    emb[batch == "b2", 0] += batch_shift
    return emb, pd.Series(batch), pd.Series(labels)


def test_lisi_is_one_for_separated_and_two_for_mixed():
    mixed, batch, _ = embedding(0)
    split, _, _ = embedding(60)
    assert np.median(lisi_scores(mixed, batch)) > 1.7
    assert np.median(lisi_scores(split, batch)) == pytest.approx(1.0, abs=0.05)


def test_kbet_accepts_mixed_and_rejects_split():
    mixed, batch, _ = embedding(0)
    split, _, _ = embedding(60)
    assert kbet_acceptance(mixed, batch) > 0.8
    assert kbet_acceptance(split, batch) < 0.05


def test_evaluate_ranks_a_mixed_embedding_above_a_batch_split_one():
    mixed, batch, labels = embedding(0)
    split, _, _ = embedding(5)
    good, bad = evaluate(mixed, batch, labels), evaluate(split, batch, labels)
    assert good["batch_score"] > bad["batch_score"]
    assert good["overall_score"] > bad["overall_score"]
    assert good["ari"] == pytest.approx(1.0)
    for key in ("ilisi", "clisi", "label_asw", "batch_asw", "kbet_acceptance"):
        assert 0 <= good[key] <= 1


def test_overmixing_is_penalised_on_biology():
    """An embedding that erases cell types must lose on bio conservation."""
    mixed, batch, labels = embedding(0)
    rng = np.random.default_rng(1)
    collapsed = rng.normal(0, 1, mixed.shape)
    erased = evaluate(collapsed, batch, labels)["bio_score"]
    assert erased < evaluate(mixed, batch, labels)["bio_score"]


def test_unassigned_labels_do_not_crash():
    mixed, batch, _ = embedding(0)
    out = evaluate(mixed, batch, pd.Series(["Unassigned"] * len(batch)))
    assert np.isnan(out["bio_score"])
    assert out["overall_score"] == out["batch_score"]
