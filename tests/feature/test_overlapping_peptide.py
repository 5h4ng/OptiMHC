import networkx as nx
import pandas as pd

from optimhc.feature.overlapping_peptide import OverlappingPeptideFeatureGenerator
from optimhc.psm_container import PsmContainer


def _generator(peptides=(), min_overlap_length=3):
    return OverlappingPeptideFeatureGenerator(
        peptides=list(peptides),
        min_overlap_length=min_overlap_length,
        min_length=1,
        max_length=100,
        fill_missing="zero",
    )


def test_empty_input_returns_an_empty_feature_table():
    features = _generator().generate_features()

    assert features.empty
    assert features.columns.tolist() == [
        "Peptide",
        "contig_member_count",
        "contig_extension_ratio",
        "contig_member_rank",
        "contig_length",
    ]


def test_transitive_reduction_preserves_a_two_node_cycle():
    graph = nx.DiGraph()
    graph.add_edge("ABCDE", "DEABC", weight=2)
    graph.add_edge("DEABC", "ABCDE", weight=3)

    reduced = _generator()._remove_transitive_edges(graph)

    assert set(reduced.edges) == {("ABCDE", "DEABC"), ("DEABC", "ABCDE")}


def test_transitive_reduction_removes_an_offset_consistent_long_edge():
    graph = nx.DiGraph()
    chain = ["ABCDEFGH", "BCDEFGHI", "CDEFGHIJ", "DEFGHIJK", "EFGHIJKL"]
    for left, right in zip(chain, chain[1:]):
        graph.add_edge(left, right, weight=7)
    graph.add_edge(chain[0], chain[-1], weight=4)

    reduced = _generator()._remove_transitive_edges(graph)

    assert not reduced.has_edge(chain[0], chain[-1])
    assert all(reduced.has_edge(left, right) for left, right in zip(chain, chain[1:]))


def test_transitive_reduction_keeps_an_offset_inconsistent_reachable_edge():
    graph = nx.DiGraph()
    graph.add_edge("ABCDEFGH", "BCDEFGHI", weight=7)  # offset 1
    graph.add_edge("BCDEFGHI", "CDEFGHIJ", weight=7)  # offset 1
    graph.add_edge("ABCDEFGH", "CDEFGHIJ", weight=4)  # offset 4, not 2

    reduced = _generator()._remove_transitive_edges(graph)

    assert reduced.has_edge("ABCDEFGH", "CDEFGHIJ")


def test_cycle_is_emitted_as_a_canonical_contig():
    graph = nx.DiGraph()
    graph.add_edge("ABCDE", "DEABC", weight=2)
    graph.add_edge("DEABC", "ABCDE", weight=3)

    contigs = _generator()._simplify_graph_to_contigs(graph)

    assert contigs == [["ABCDE", "DEABC"]]


def test_linear_layout_is_independent_of_lexicographic_endpoint_order():
    source = "DAPSDQSDTSESDVDLG"
    sink = "APSDQSDTSESDVDLGDG"

    features = (
        _generator([source, sink], min_overlap_length=16)
        .generate_features()
        .set_index("Peptide")
    )

    assert features.loc[source, "contig_member_count"] == 2
    assert features.loc[sink, "contig_member_count"] == 2
    assert features.loc[source, "contig_length"] == 19
    assert features.loc[sink, "contig_length"] == 19


def test_branching_graph_does_not_drop_internal_nodes():
    graph = nx.DiGraph()
    graph.add_edges_from(
        [
            ("A", "B"),
            ("A", "C"),
            ("B", "D"),
        ]
    )

    contigs = _generator()._simplify_graph_to_contigs(graph)
    emitted_nodes = [node for contig in contigs for node in contig]

    assert sorted(emitted_nodes) == ["A", "B", "C", "D"]
    assert len(emitted_nodes) == len(set(emitted_nodes))


def test_raw_and_simplified_graphs_are_independent():
    generator = _generator(
        ["ABCDEFGH", "BCDEFGHI", "CDEFGHIJ"], min_overlap_length=6
    )

    generator.generate_features()

    assert generator.overlap_graph is not generator.simplified_graph
    assert generator.overlap_graph.has_edge("ABCDEFGH", "CDEFGHIJ")
    assert not generator.simplified_graph.has_edge("ABCDEFGH", "CDEFGHIJ")


def test_contig_rank_preserves_peptide_row_competition_ranking():
    generator = _generator(
        ["ABCDEF", "DEFGHI", "GHIKLM", "PQRSTU", "STUVWX"],
        min_overlap_length=3,
    )

    features = generator.generate_features().set_index("Peptide")

    assert features.loc["ABCDEF", "contig_member_rank"] == 1
    assert features.loc["PQRSTU", "contig_member_rank"] == 4


def test_containment_and_features_are_invariant_to_input_order():
    peptides = [
        "ABCDEFGH",
        "BCDEFGHI",
        "CDEFGHIJ",
        "PQRSTUVX",
        "QRSTUVXY",
        "RSTUVXYZ",
    ]

    expected = (
        _generator(peptides, min_overlap_length=6)
        .generate_features()
        .sort_values("Peptide")
        .reset_index(drop=True)
    )
    actual = (
        _generator(reversed(peptides), min_overlap_length=6)
        .generate_features()
        .sort_values("Peptide")
        .reset_index(drop=True)
    )

    pd.testing.assert_frame_equal(actual, expected)


def test_contained_peptide_uses_a_stable_largest_container():
    generator = _generator()

    accepted, contained_mapping = generator._remove_redundant_peptides(
        ["ABCDEYY", "ABCDE", "ABCDEXX"]
    )

    assert accepted == ["ABCDEXX", "ABCDEYY"]
    assert contained_mapping == {"ABCDE": "ABCDEXX"}


def test_from_config_deduplicates_peptides_in_stable_order():
    frame = pd.DataFrame(
        {
            "psm_id": [0, 1, 2],
            "run": ["run", "run", "run"],
            "scan": [1, 2, 3],
            "rank": [1, 1, 1],
            "sequence": ["PEPTIDE", "ABCDEFG", "PEPTIDE"],
            "mods": ["", "", ""],
            "mod_sites": ["", "", ""],
            "charge": [2, 2, 2],
            "proteins": ["P1", "P2", "P1"],
            "is_decoy": [False, False, False],
        }
    )
    psms = PsmContainer(frame)

    generator = OverlappingPeptideFeatureGenerator.from_config(
        psms,
        {"removePreNxtAA": False},
        {},
    )

    assert generator.original_peptides == ["ABCDEFG", "PEPTIDE"]
