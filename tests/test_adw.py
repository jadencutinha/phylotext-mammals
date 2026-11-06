from pathlib import Path

from lineage_charisma.adw import SECTIONS, extract_text, parse_account, query_names

FIXTURES = Path(__file__).parent / "fixtures"


def page(name):
    return parse_account((FIXTURES / name).read_text())


def test_species_account_is_recognised_with_its_title():
    account = page("adw_species_excerpt.html")
    assert account.is_species_account
    assert account.title == "Panthera leo"


def test_genus_page_and_missing_page_are_not_species_accounts():
    assert not page("adw_genus.html").is_species_account
    missing = page("adw_not_found.html")
    assert missing.not_found and not missing.is_species_account


def test_text_rule_keeps_the_five_sections_in_order_and_nothing_else():
    text = extract_text(page("adw_species_excerpt.html"))
    order = ["sub-Saharan Africa", "plains or savanna", "Gir Forest", "short, tawny coat", "most social", "large ungulates"]
    positions = [text.index(s) for s in order]
    assert positions == sorted(positions)
    for excluded in ("Reproduction text", "Reference list", "photo caption", "table text", "Kingdom", "Ethiopian", "Geographic Range", "Erin Harrington"):
        assert excluded not in text


def test_links_keep_their_text_and_reference_markers_are_dropped():
    text = extract_text(page("adw_species_excerpt.html"))
    assert "African lions ( Panthera leo ) live" in text or "African lions (Panthera leo) live" in text
    assert "33da4d" not in text


def test_registered_section_list():
    assert SECTIONS == ("geographic_range", "habitat", "physical_description", "behavior", "food_habits")


def test_query_names_follow_the_registered_order_without_duplicates():
    names = query_names("Lontra canadensis", "Lutra canadensis", ["Lutra canadensis", "Lontra canadensis", "Lutra hudsonica"])
    assert names == ["Lontra canadensis", "Lutra canadensis", "Lutra hudsonica"]
