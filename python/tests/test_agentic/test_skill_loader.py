"""Unit tests for mokume.agentic.skill_loader."""

import pytest

from mokume.agentic.skill_loader import (
    Skill,
    _parse_skill_md,
    extract_prompts,
    load_skill,
    load_skill_data,
    skill_body,
)


class TestParseSkillMd:
    """Tests for _parse_skill_md frontmatter parsing."""

    def test_parses_frontmatter_fields(self):
        text = "---\nname: test-skill\nversion: 1.2.3\ndescription: A test skill\n---\n\n# Body\n"
        skill = _parse_skill_md(text)
        assert skill.name == "test-skill"
        assert skill.version == "1.2.3"
        assert skill.description == "A test skill"

    def test_body_excludes_frontmatter(self):
        text = "---\nname: x\n---\n\n# Hello\nworld\n"
        skill = _parse_skill_md(text)
        assert "---" not in skill.body
        assert "# Hello" in skill.body
        assert "world" in skill.body

    def test_raises_on_missing_frontmatter(self):
        with pytest.raises(ValueError, match="missing YAML frontmatter"):
            _parse_skill_md("# No frontmatter here\n")

    def test_meta_contains_raw_dict(self):
        text = "---\nname: m\ncustom_key: 42\n---\n\nbody\n"
        skill = _parse_skill_md(text)
        assert skill.meta["custom_key"] == 42

    def test_defaults_for_missing_fields(self):
        text = "---\nname: minimal\n---\n\nbody\n"
        skill = _parse_skill_md(text)
        assert skill.version == "0.0.0"
        assert skill.description == ""


class TestExtractPrompts:
    """Tests for extract_prompts code-fence-aware parsing."""

    def test_extracts_system_and_user(self):
        body = (
            "# Title\n\n"
            "## System Prompt\n\n"
            "```\nYou are an expert.\n```\n\n"
            "## User Prompt\n\n"
            "```\nAnalyze {data}.\n```\n"
        )
        skill = Skill(name="t", version="0", description="", body=body, meta={})
        prompts = extract_prompts(skill)
        assert prompts["system"] == "You are an expert."
        assert prompts["user"] == "Analyze {data}."

    def test_ignores_headings_inside_code_blocks(self):
        body = (
            "## System Prompt\n\n"
            "```\n## This is NOT a heading\nContent here.\n```\n\n"
            "## User Prompt\n\n"
            "```\nUser content.\n```\n"
        )
        skill = Skill(name="t", version="0", description="", body=body, meta={})
        prompts = extract_prompts(skill)
        assert "## This is NOT a heading" in prompts["system"]
        assert prompts["user"] == "User content."

    def test_returns_empty_when_no_code_blocks(self):
        body = "## System Prompt\n\nNo code block here.\n"
        skill = Skill(name="t", version="0", description="", body=body, meta={})
        prompts = extract_prompts(skill)
        assert prompts == {}


class TestLoadSkill:
    """Tests for load_skill against real skill files."""

    def test_loads_proteomics_heuristics(self):
        skill = load_skill("proteomics-heuristics")
        assert skill.name == "proteomics-heuristics"
        assert skill.version == "0.1.0"
        assert "heuristics" in skill.description.lower()

    def test_loads_proteomics_proposal(self):
        skill = load_skill("proteomics-proposal")
        prompts = extract_prompts(skill)
        assert "system" in prompts
        assert "user" in prompts
        assert "{relevant_heuristics}" in prompts["system"]
        assert "{data_profile_json}" in prompts["user"]

    def test_raises_on_nonexistent_skill(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            load_skill("nonexistent-skill-xyz")


class TestLoadSkillData:
    """Tests for load_skill_data companion YAML loading."""

    def test_loads_heuristics_data(self):
        data = load_skill_data("proteomics-heuristics")
        assert "data_type_priors" in data
        assert "LFQ" in data["data_type_priors"]
        assert "condition_rules" in data

    def test_loads_benchmarks_data(self):
        data = load_skill_data("proteomics-benchmarks")
        assert "PXD001819" in data
        assert data["PXD001819"]["best_known"]["method"] == "deqms"

    def test_raises_on_missing_data_file(self):
        with pytest.raises(FileNotFoundError, match="No data.yaml"):
            load_skill_data("proteomics-proposal")


class TestSkillBody:
    """Tests for skill_body convenience function."""

    def test_returns_markdown_body(self):
        body = skill_body("proteomics-heuristics")
        assert "# Proteomics Heuristics" in body
        assert not body.startswith("---")
