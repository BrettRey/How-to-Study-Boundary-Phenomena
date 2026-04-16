from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "matrix_clean.csv"
OUTPUT = ROOT / "data" / "inventory_annotations.csv"


PERSONAL_MASC = {"he", "him", "his", "himself"}
PERSONAL_FEM = {"she", "her_acc", "her_dep", "hers", "herself"}

PERSONAL_EPICENE = {
    "I",
    "me",
    "mine",
    "my",
    "myself",
    "we_pron",
    "us_pron",
    "our",
    "ours",
    "ourselves",
    "you_pron_sing",
    "you_pron_plur",
    "your_sing",
    "yours_sing",
    "your_plur",
    "yours_plur",
    "yourself",
    "yourselves",
    "they_sing",
    "them_sing",
    "their_sing",
    "theirs_sing",
    "themself",
    "one_pron",
    "one_s",
    "oneself",
    "who_int",
    "who_rel",
    "whom_int",
    "whom_rel",
    "whose_int",
    "whoever",
    "whoever_s",
    "whosever",
    "whosoever",
    "whosoever_s",
}

GENDER_NEUTRAL_PRON = {
    "they_plur",
    "them_plur",
    "their_plur",
    "theirs_plur",
    "themselves",
    "whose_rel",
    "each_other",
    "one_another",
}

NONPERSONAL_GENERAL = {
    "it_plain",
    "its",
    "itself",
    "whatever_pron",
    "much",
    "little",
    "a_little",
    "something",
    "anything",
    "everything",
    "nothing",
}

TEMPORAL_NOUNS = {
    "today",
    "tomorrow",
    "tonight",
    "yesterday",
}

NONPERSONAL_TEMPORAL = {
    "once",
    "twice",
    "thrice",
}

NONPERSONAL_LOCATIVE = {
    "somewhere",
    "anywhere",
    "everywhere",
    "nowhere",
}

PERSONAL_DETERMINATIVES = {
    "anybody",
    "anyone",
    "everybody",
    "everyone",
    "nobody",
    "no_one",
    "somebody",
    "someone",
    "we_det",
    "us_det",
    "you_det",
}

GENDER_NEUTRAL_DETERMINATIVES = {
    "all",
    "another",
    "any",
    "both",
    "certain",
    "each",
    "either",
    "enough",
    "every",
    "few",
    "fewer",
    "fewest",
    "many",
    "many_a",
    "most",
    "neither",
    "next",
    "last",
    "no",
    "none",
    "several",
    "some",
    "somewhat",
    "such",
    "sufficient",
    "the",
    "these",
    "those",
    "this",
    "that",
    "various_det",
    "zero_det",
    "one_det",
    "two_det",
    "three_det",
}

NOT_CORE_PROFORM = {"a", "a_certain", "a_few", "a_great_many"}
DUMMY_ONLY = {"it_dum"}
AMBIGUOUS_DUMMY_OR_LOCATIVE = set()
WHAT_FAMILY_DETERMINERS = {
    "what_det",
    "whatever_det",
}

WHAT_FAMILY_ITEMS = {
    "what_det",
    "what_pron",
    "whatever_det",
}

WHICH_FAMILY_ITEMS = {
    "which_rel_det",
    "which_int_det",
    "whichever_rel_det",
    "whichever_int_det",
    "which_rel_pron",
}

CLEAR_WHICH_DETERMINERS = {
    "which_rel_det",
    "which_int_det",
    "whichever_rel_det",
    "whichever_int_det",
}

CLEAR_PERSONAL_DETERMINERS = {
    "we_det",
    "us_det",
    "you_det",
}

CATEGORY_OVERRIDES = {
    "what_pron": {
        "revised_class": "determinative",
        "category_source": "user_instruction_wh_never_pronouns",
        "notes": "Seed row encoded pronoun use, but the lexical item is treated as determinative throughout. Keep this compatibility row out of clear anchors to avoid double counting.",
        "anchor_recommendation": "exclude",
    },
    "which_rel_pron": {
        "revised_class": "determinative",
        "category_source": "user_instruction_which_is_determinative",
        "notes": "Seed row encoded pronoun use, but the lexical item is treated as determinative throughout. Keep this compatibility row out of clear anchors to avoid double counting.",
        "anchor_recommendation": "exclude",
    },
    "more": {
        "revised_class": "determinative",
        "category_source": "reynolds2024_more_less_never_adverbs",
        "notes": "Retained as determinative under the unified analysis of more.",
    },
    "less": {
        "revised_class": "determinative",
        "category_source": "reynolds2024_more_less_never_adverbs",
        "notes": "Retained as determinative under the unified analysis of less.",
    },
    "today": {
        "revised_class": "noun",
        "category_source": "user_instruction_temporal_items_are_nouns",
        "notes": "Temporal noun; exclude from the pronoun/determinative inventory.",
        "anchor_recommendation": "exclude",
    },
    "tomorrow": {
        "revised_class": "noun",
        "category_source": "user_instruction_temporal_items_are_nouns",
        "notes": "Temporal noun; exclude from the pronoun/determinative inventory.",
        "anchor_recommendation": "exclude",
    },
    "tonight": {
        "revised_class": "noun",
        "category_source": "user_instruction_temporal_items_are_nouns",
        "notes": "Temporal noun; exclude from the pronoun/determinative inventory.",
        "anchor_recommendation": "exclude",
    },
    "yesterday": {
        "revised_class": "noun",
        "category_source": "user_instruction_temporal_items_are_nouns",
        "notes": "Temporal noun; exclude from the pronoun/determinative inventory.",
        "anchor_recommendation": "exclude",
    },
    "there": {
        "revised_class": "preposition",
        "category_source": "user_instruction_there_is_preposition_except_existential_and_tag_pronoun",
        "notes": "Seed row conflates existential/tag pronoun uses with the broader preposition use; exclude from the pronoun/determinative inventory until context-sensitive splitting is added.",
        "anchor_recommendation": "exclude",
    },
}

ANCHOR_REVIEWS = {
    "whatever_pron",
    "each_other",
    "one_another",
    "it_dum",
}


def classify_personhood(lemma: str, revised_class: str) -> str:
    if lemma in PERSONAL_MASC:
        return "personal_masculine"
    if lemma in PERSONAL_FEM:
        return "personal_feminine"
    if lemma in PERSONAL_EPICENE or lemma in PERSONAL_DETERMINATIVES:
        return "personal_epicene"
    if lemma in GENDER_NEUTRAL_PRON or lemma in GENDER_NEUTRAL_DETERMINATIVES:
        return "gender_neutral"
    if lemma in NONPERSONAL_GENERAL:
        return "non_personal_general"
    if lemma in TEMPORAL_NOUNS:
        return "temporal_noun_outside_system"
    if revised_class == "preposition":
        return "preposition_outside_system"
    if lemma in NONPERSONAL_TEMPORAL:
        return "non_personal_temporal"
    if lemma in NONPERSONAL_LOCATIVE:
        return "non_personal_locative"
    if lemma in DUMMY_ONLY:
        return "dummy_outside_system"
    if lemma in AMBIGUOUS_DUMMY_OR_LOCATIVE:
        return "dummy_or_locative_review"
    if lemma in NOT_CORE_PROFORM:
        return "not_core_pro_form"
    if lemma in WHAT_FAMILY_ITEMS or lemma in WHICH_FAMILY_ITEMS:
        return "non_personal_if_pro_form__neutral_with_overt_head"
    if revised_class == "pronoun":
        return "review_needed"
    return "gender_neutral_or_contextual"


def main() -> None:
    with INPUT.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    records = []
    for row in rows:
        lemma = row["lemma"]
        seed_class = row["class"]
        record = {
            "lemma": lemma,
            "seed_class": seed_class,
            "revised_class": seed_class,
            "category_source": "matrix_clean_seed",
            "personhood_profile": "",
            "anchor_recommendation": "retain",
            "notes": "",
        }

        if lemma in CATEGORY_OVERRIDES:
            override = CATEGORY_OVERRIDES[lemma]
            record["revised_class"] = override["revised_class"]
            record["category_source"] = override["category_source"]
            if "anchor_recommendation" in override:
                record["anchor_recommendation"] = override["anchor_recommendation"]
            if "notes" in override:
                record["notes"] = override["notes"]

        record["personhood_profile"] = classify_personhood(
            lemma, record["revised_class"]
        )

        if lemma in ANCHOR_REVIEWS and record["anchor_recommendation"] == "retain":
            record["anchor_recommendation"] = "review"

        if lemma in NOT_CORE_PROFORM:
            record["anchor_recommendation"] = "exclude"
            record["notes"] = (
                "Outside the core pro-form inventory in the personhood paper."
                if not record["notes"]
                else f"{record['notes']} Outside the core pro-form inventory in the personhood paper."
            )
        elif lemma in DUMMY_ONLY:
            record["anchor_recommendation"] = "exclude"
            record["notes"] = (
                "Dummy pronoun outside the personhood system."
                if not record["notes"]
                else f"{record['notes']} Dummy pronoun outside the personhood system."
            )
        elif lemma in AMBIGUOUS_DUMMY_OR_LOCATIVE:
            record["anchor_recommendation"] = "review"
            record["notes"] = (
                "Row may conflate locative and dummy uses."
                if not record["notes"]
                else f"{record['notes']} Row may conflate locative and dummy uses."
            )
        elif lemma in CLEAR_WHICH_DETERMINERS:
            record["notes"] = (
                "Treated as determinative throughout; retained as a clear anchor."
                if not record["notes"]
                else f"{record['notes']} Treated as determinative throughout; retained as a clear anchor."
            )
        elif lemma in WHAT_FAMILY_DETERMINERS:
            note = (
                "Treated as determinative throughout, including exclamative what a N; retained as a clear anchor."
                if lemma == "what_det"
                else "Treated as determinative throughout; retained as a clear anchor."
            )
            record["notes"] = (
                note if not record["notes"] else f"{record['notes']} {note}"
            )
        elif lemma in CLEAR_PERSONAL_DETERMINERS:
            record["notes"] = (
                "Treated as determinative throughout; retained as a clear anchor."
                if not record["notes"]
                else f"{record['notes']} Treated as determinative throughout; retained as a clear anchor."
            )
        elif lemma in {"each_other", "one_another"}:
            record["notes"] = (
                "Reciprocal boundary items; exclude from clear pronoun anchors."
                if not record["notes"]
                else f"{record['notes']} Reciprocal boundary items; exclude from clear pronoun anchors."
            )
            record["anchor_recommendation"] = "exclude"
        elif lemma in {"this", "that", "these", "those"}:
            record["notes"] = (
                "Generally gender-neutral with a non-personal default in fused-head use."
                if not record["notes"]
                else f"{record['notes']} Generally gender-neutral with a non-personal default in fused-head use."
            )

        records.append(record)

    with OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "lemma",
                "seed_class",
                "revised_class",
                "category_source",
                "personhood_profile",
                "anchor_recommendation",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
