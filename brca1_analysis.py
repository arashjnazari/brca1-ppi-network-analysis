"""Descriptive analysis of the supplied BRCA1-centered STRING association network.

The analyzed graph contains nonisolated endpoints of retained associations,
with the query protein excluded. It remains a selected protein neighborhood.
"""

from hashlib import sha256
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "string_interactions_short_501_nodes.tsv"
OUTPUT_DIR = BASE_DIR / "output"
SEED_PROTEIN = "BRCA1"
SCORE_THRESHOLD = 0.9
LAYOUT_SEED = 42


def load_graph(data_path=DATA_FILE):
    """Return the unweighted, undirected active-endpoint graph and its metadata.

    Retain combined STRING scores >= 0.9, omit self-loops, collapse duplicate
    unordered pairs, and exclude BRCA1. Nodes with no remaining association are
    omitted, matching the original analysis population. Their count is recorded.
    """
    data_path = Path(data_path)
    table = pd.read_csv(data_path, sep="\t")
    required = {"#node1", "node2", "combined_score"}
    if not required.issubset(table.columns):
        raise ValueError(f"Missing required columns: {sorted(required - set(table.columns))}")
    if table[["#node1", "node2", "combined_score"]].isna().any().any():
        raise ValueError("Endpoint labels and confidence scores must not be missing.")
    score = pd.to_numeric(table["combined_score"], errors="raise")
    if not np.isfinite(score).all() or not score.between(0, 1).all():
        raise ValueError("This export must use finite combined scores on the 0–1 scale.")

    raw_nodes = set(table["#node1"]) | set(table["node2"])
    if SEED_PROTEIN not in raw_nodes:
        raise ValueError(f"Query protein {SEED_PROTEIN} is absent from the supplied export.")
    raw_pairs = [tuple(sorted((a, b))) for a, b in zip(table["#node1"], table["node2"])]
    retained = table.loc[score >= SCORE_THRESHOLD]
    threshold_pairs = {
        tuple(sorted((a, b)))
        for a, b in zip(retained["#node1"], retained["node2"])
        if a != b
    }
    threshold_nodes = {node for pair in threshold_pairs for node in pair}
    analysis_pairs = sorted(pair for pair in threshold_pairs if SEED_PROTEIN not in pair)
    active_nodes = {node for pair in analysis_pairs for node in pair}

    graph = nx.Graph()
    # Stable insertion order also makes graph algorithms' tie handling repeatable.
    graph.add_nodes_from(sorted(active_nodes))
    graph.add_edges_from(analysis_pairs)
    if graph.number_of_edges() == 0:
        raise ValueError("No non-query associations remain after filtering.")
    graph.graph.update(
        data_file=data_path.name,
        input_sha256=sha256(data_path.read_bytes()).hexdigest(),
        query_protein=SEED_PROTEIN,
        score_threshold=SCORE_THRESHOLD,
        raw_rows=len(table),
        raw_nodes=len(raw_nodes),
        raw_self_loop_rows=sum(a == b for a, b in raw_pairs),
        raw_duplicate_undirected_rows=len(raw_pairs) - len(set(raw_pairs)),
        threshold_rows=len(retained),
        threshold_edges_before_query_removal=len(threshold_pairs),
        threshold_nodes_before_query_removal=len(threshold_nodes),
        query_incident_edges_removed=len(threshold_pairs) - len(analysis_pairs),
        isolates_created_by_query_removal=len(threshold_nodes - active_nodes - {SEED_PROTEIN}),
        nonquery_export_nodes_omitted=len(raw_nodes - active_nodes - {SEED_PROTEIN}),
        graph_type="simple_undirected_unweighted",
        node_population="nonisolated_endpoints_after_filtering_and_query_exclusion",
    )
    return graph


def degree_table(graph):
    """Return degree counts, probability masses, and the inclusive empirical CCDF."""
    degrees = np.array([degree for _, degree in graph.degree()], dtype=int)
    values, counts = np.unique(degrees, return_counts=True)
    return pd.DataFrame({
        "degree": values,
        "count": counts,
        "probability": counts / len(degrees),
        "ccdf": np.cumsum(counts[::-1])[::-1] / len(degrees),
    })


def analyze_graph(graph):
    """Compute descriptive metrics; communities are evaluated on the LCC only."""
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), tuple(sorted(c))))
    lcc = graph.subgraph(components[0]).copy()
    communities = sorted(
        greedy_modularity_communities(lcc),
        key=lambda community: (-len(community), tuple(sorted(community))),
    )
    component_of = {node: i for i, component in enumerate(components, 1) for node in component}
    community_of = {
        node: f"C{i:02d}" for i, community in enumerate(communities, 1) for node in community
    }
    degree_centrality = nx.degree_centrality(graph)
    betweenness = nx.betweenness_centrality(graph, normalized=True, weight=None)
    node_metrics = pd.DataFrame([
        {
            "protein": node,
            "degree": graph.degree(node),
            "degree_centrality": degree_centrality[node],
            "betweenness_centrality": betweenness[node],
            "component": component_of[node],
            "community_lcc": community_of.get(node, pd.NA),
        }
        for node in sorted(graph)
    ])
    membership = pd.DataFrame([
        {"protein": node, "community_lcc": community_of[node], "community_size": len(community)}
        for community in communities for node in sorted(community)
    ])
    density = nx.density(graph)
    clustering = nx.average_clustering(graph, weight=None)
    summary = dict(graph.graph)
    summary.update(
        nodes=graph.number_of_nodes(),
        edges=graph.number_of_edges(),
        mean_degree=2 * graph.number_of_edges() / graph.number_of_nodes(),
        average_clustering=clustering,
        density=density,
        # p is a theoretical density-matched ER clustering reference, not a
        # simulated null distribution or a test of small-world structure.
        er_clustering_reference_p=density,
        clustering_to_er_reference_ratio=clustering / density,
        connected_components=len(components),
        lcc_nodes=lcc.number_of_nodes(),
        lcc_edges=lcc.number_of_edges(),
        lcc_node_fraction=lcc.number_of_nodes() / graph.number_of_nodes(),
        lcc_communities=len(communities),
        lcc_modularity=modularity(lcc, communities, weight=None),
        layout_seed=LAYOUT_SEED,
    )
    return summary, node_metrics, membership, lcc, communities


def plot_degree_distribution(table, path):
    """Plot empirical PMF and CCDF without a log-log regression fit."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, column, title, ylabel, color in (
        (axes[0], "probability", "Empirical degree distribution", "P(K = k)", "#28639e"),
        (axes[1], "ccdf", "Empirical complementary cumulative distribution", "P(K >= k)", "#a93d55"),
    ):
        ax.loglog(table["degree"], table[column], "o", color=color, markersize=4.5, alpha=0.85)
        ax.set(xlabel="Degree k", ylabel=ylabel, title=title)
        ax.grid(True, which="major", color="#d8dde3", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("BRCA1-associated proteins: score >= 0.9, query protein excluded", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_communities(lcc, communities, node_metrics, modularity_score, path):
    """Draw LCC topology; colors identify algorithmic communities only."""
    fig, ax = plt.subplots(figsize=(12, 10.5))
    pos = nx.spring_layout(lcc, seed=LAYOUT_SEED, k=0.18, iterations=80)
    palette = plt.get_cmap("tab20", max(1, len(communities)))
    community_of = {node: i for i, community in enumerate(communities) for node in community}
    nx.draw_networkx_edges(lcc, pos, ax=ax, alpha=0.12, edge_color="#737d87", width=0.7)
    nx.draw_networkx_nodes(
        lcc, pos, ax=ax,
        node_color=[palette(community_of[node]) for node in lcc],
        node_size=[12 + lcc.degree(node) * 3.5 for node in lcc],
        linewidths=0.3, edgecolors="white", alpha=0.9,
    )
    ranked = node_metrics.sort_values(["degree", "protein"], ascending=[False, True])
    names = [node for node in ranked["protein"] if node in lcc][:8]
    # Spread hub labels around their shared center rather than overprinting the
    # tightly packed hub coordinates. Leader lines retain the node association.
    hub_center = np.mean([pos[node] for node in names], axis=0)
    ordered_names = sorted(
        names, key=lambda node: np.arctan2(*(pos[node] - hub_center)[::-1])
    )
    for i, node in enumerate(ordered_names):
        angle = -np.pi + 2 * np.pi * (i + 0.5) / len(ordered_names)
        label_pos = hub_center + np.array([0.18 * np.cos(angle), 0.14 * np.sin(angle)])
        ax.annotate(
            node, xy=pos[node], xytext=label_pos, fontsize=9, fontweight="bold",
            ha="center", va="center",
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none", "pad": 1.5},
            arrowprops={"arrowstyle": "-", "color": "#525c66", "lw": 0.7},
        )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=palette(i),
               markeredgecolor="none", markersize=7, label=f"C{i + 1:02d}: {len(community)} proteins")
        for i, community in enumerate(communities)
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.055),
              ncol=4, frameon=False, fontsize=8)
    ax.set_title(
        f"BRCA1 association network: topological communities in the LCC\n"
        f"{lcc.number_of_nodes()} proteins | {lcc.number_of_edges()} associations | Q = {modularity_score:.3f}",
        fontsize=13,
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    """Export reproducible descriptive tables and figures for the supplied data."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    graph = load_graph()
    summary, nodes, membership, lcc, communities = analyze_graph(graph)
    degrees = degree_table(graph)
    pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "network_summary.csv", index=False)
    nodes.to_csv(OUTPUT_DIR / "node_metrics.csv", index=False)
    membership.to_csv(OUTPUT_DIR / "community_membership.csv", index=False)
    degrees.to_csv(OUTPUT_DIR / "degree_counts.csv", index=False)
    plot_degree_distribution(degrees, BASE_DIR / "degree_distribution.png")
    plot_communities(lcc, communities, nodes, summary["lcc_modularity"], BASE_DIR / "network_communities.png")
    print(f"BRCA1: {summary['nodes']} active proteins, {summary['edges']} associations")
    print(f"Mean degree = {summary['mean_degree']:.3f}; average clustering = {summary['average_clustering']:.4f}")
    print(f"Density-matched theoretical ER reference p = {summary['er_clustering_reference_p']:.5f}")
    print(f"LCC: {summary['lcc_nodes']} proteins; {summary['lcc_communities']} topological communities; Q = {summary['lcc_modularity']:.4f}")
    print(f"Non-query export proteins with no retained association: {summary['nonquery_export_nodes_omitted']}")
    print("Top degree centralities (ties ordered by protein label):")
    print(nodes.sort_values(["degree", "protein"], ascending=[False, True]).head(10).to_string(index=False))
    print(f"Tables saved to {OUTPUT_DIR}")
    print("This selected association network does not establish scale-free, small-world, or disease-causal behavior.")


if __name__ == "__main__":
    main()
