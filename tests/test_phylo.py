import dendropy
import pytest

from lineage_charisma.phylo import (
    cophenetic_matrix,
    filter_tips,
    group_separation,
    is_ultrametric,
    load_tree,
    nearest_neighbors,
    parse_tip_label,
    prune_tree,
    relabel_tips,
    tip_table,
)

NEWICK = (
    "(((Panthera_leo_FELIDAE_CARNIVORA:2,Panthera_tigris_FELIDAE_CARNIVORA:2):3,"
    "Felis_catus_FELIDAE_CARNIVORA:5):5,"
    "((Canis_lupus_CANIDAE_CARNIVORA:1,Canis_aureus_CANIDAE_CARNIVORA:1):6,"
    "Vulpes_vulpes_CANIDAE_CARNIVORA:7):3,"
    "Homo_sapiens_HOMINIDAE_PRIMATES:10);"
)

NEXUS = """#NEXUS
Begin taxa;
    Dimensions ntax=3;
    Taxlabels
        _Anolis_carolinensis
        Panthera_leo_FELIDAE_CARNIVORA
        Canis_lupus_CANIDAE_CARNIVORA
        ;
End;
Begin trees;
    Translate
        1 _Anolis_carolinensis,
        2 Panthera_leo_FELIDAE_CARNIVORA,
        3 Canis_lupus_CANIDAE_CARNIVORA
        ;
tree TREE1 = [&R] (1:[&height=0]20,(2:[&height=0]8,3:[&height=0]8):12);
End;
"""


@pytest.fixture
def tree():
    return dendropy.Tree.get(data=NEWICK, schema="newick", preserve_underscores=True, rooting="default-rooted")


def test_parse_tip_label_standard():
    p = parse_tip_label("Panthera_leo_FELIDAE_CARNIVORA")
    assert p == {
        "tip_label": "Panthera_leo_FELIDAE_CARNIVORA",
        "binomial": "Panthera leo",
        "genus": "Panthera",
        "epithet": "leo",
        "family": "FELIDAE",
        "order": "CARNIVORA",
    }


def test_parse_tip_label_mixed_case_family_and_hyphen():
    p = parse_tip_label("Mammuthus_columbi_Mammutidae_PROBOSCIDEA")
    assert p["family"] == "MAMMUTIDAE"
    assert parse_tip_label("Ursus_arctos-horribilis_URSIDAE_CARNIVORA")["epithet"] == "arctos-horribilis"


def test_parse_tip_label_handles_quotes_and_spaces():
    assert parse_tip_label("'Panthera leo FELIDAE CARNIVORA'")["binomial"] == "Panthera leo"


def test_parse_tip_label_rejects_outgroup_and_junk():
    assert parse_tip_label("_Anolis_carolinensis") is None
    assert parse_tip_label("node123") is None


def test_parse_tip_label_bare_binomial():
    p = parse_tip_label("Panthera_leo")
    assert p["binomial"] == "Panthera leo" and p["family"] is None


def test_tip_table_and_filter(tree):
    tips = tip_table(tree)
    assert len(tips) == 7
    assert len(filter_tips(tips, orders=["Carnivora"])) == 6
    assert set(filter_tips(tips, families=["canidae"])["binomial"]) == {"Canis lupus", "Canis aureus", "Vulpes vulpes"}


def test_load_nexus_with_translate_and_annotations(tmp_path):
    p = tmp_path / "t.tre"
    p.write_text(NEXUS)
    t = load_tree(p)
    tips = tip_table(t)
    assert set(tips["binomial"]) == {"Panthera leo", "Canis lupus"}
    assert len(list(t.leaf_node_iter())) == 3


def test_load_newick_autodetect(tmp_path):
    p = tmp_path / "t.nwk"
    p.write_text(NEWICK)
    assert len(tip_table(load_tree(p))) == 7


def test_prune_keeps_only_requested(tree):
    keep = ["Panthera_leo_FELIDAE_CARNIVORA", "Panthera_tigris_FELIDAE_CARNIVORA", "Canis_lupus_CANIDAE_CARNIVORA"]
    pruned = prune_tree(tree, keep)
    assert sorted(l.taxon.label for l in pruned.leaf_node_iter()) == sorted(keep)
    assert len(pruned.taxon_namespace) == 3
    assert len(list(tree.leaf_node_iter())) == 7


def test_prune_missing_label_raises(tree):
    with pytest.raises(KeyError):
        prune_tree(tree, ["Nope_nope_X_Y"])


def test_prune_preserves_patristic_distances(tree):
    full = cophenetic_matrix(tree)
    keep = ["Panthera_leo_FELIDAE_CARNIVORA", "Canis_lupus_CANIDAE_CARNIVORA", "Felis_catus_FELIDAE_CARNIVORA"]
    pruned = cophenetic_matrix(prune_tree(tree, keep))
    for a in keep:
        for b in keep:
            assert pruned.loc[a, b] == pytest.approx(full.loc[a, b])


def test_cophenetic_values(tree):
    d = cophenetic_matrix(tree)
    assert d.loc["Panthera_leo_FELIDAE_CARNIVORA", "Panthera_tigris_FELIDAE_CARNIVORA"] == pytest.approx(4)
    assert d.loc["Panthera_leo_FELIDAE_CARNIVORA", "Felis_catus_FELIDAE_CARNIVORA"] == pytest.approx(10)
    assert d.loc["Panthera_leo_FELIDAE_CARNIVORA", "Canis_lupus_CANIDAE_CARNIVORA"] == pytest.approx(20)
    assert (d.values == d.values.T).all()
    assert (d.values.diagonal() == 0).all()


def test_group_separation_and_neighbors(tree):
    d = cophenetic_matrix(tree)
    cats = ["Panthera_leo_FELIDAE_CARNIVORA", "Panthera_tigris_FELIDAE_CARNIVORA", "Felis_catus_FELIDAE_CARNIVORA"]
    assert group_separation(d, cats)["clusters"] is True
    mixed = ["Panthera_leo_FELIDAE_CARNIVORA", "Canis_lupus_CANIDAE_CARNIVORA"]
    assert group_separation(d, mixed)["clusters"] is False
    assert nearest_neighbors(d, "Canis_lupus_CANIDAE_CARNIVORA", k=1)[0][0] == "Canis_aureus_CANIDAE_CARNIVORA"


def test_relabel_does_not_mutate_original(tree):
    out = relabel_tips(tree, {"Panthera_leo_FELIDAE_CARNIVORA": "Panthera_leo"})
    assert "Panthera_leo" in {l.taxon.label for l in out.leaf_node_iter()}
    assert "Panthera_leo" not in {l.taxon.label for l in tree.leaf_node_iter()}


def test_ultrametric(tree):
    assert is_ultrametric(tree)
