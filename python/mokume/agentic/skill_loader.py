"""Skill loader: parse SKILL.md (YAML frontmatter + markdown body)."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import yaml

from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.skill_loader")

_SKILLS_DIR = Path(__file__).parent / "skills"

_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)
_CODEBLOCK_RE = re.compile(r"```[^\n]*\n(.*?)```", re.DOTALL)


@dataclass(frozen=True)
class Skill:
    """Parsed skill from a SKILL.md file."""

    name: str
    version: str
    description: str
    body: str
    meta: dict = field(default_factory=dict)

    def code_blocks(self) -> list[str]:
        """Extract all fenced code blocks from the markdown body."""
        return _CODEBLOCK_RE.findall(self.body)


def _parse_skill_md(text: str) -> Skill:
    """Parse a SKILL.md file into a Skill dataclass."""
    match = _FRONTMATTER_RE.match(text)
    if not match:
        raise ValueError("SKILL.md missing YAML frontmatter")
    meta = yaml.safe_load(match.group(1)) or {}
    body = text[match.end() :]
    return Skill(
        name=meta.get("name", ""),
        version=meta.get("version", "0.0.0"),
        description=meta.get("description", ""),
        body=body,
        meta=meta,
    )


@lru_cache(maxsize=16)
def load_skill(name: str) -> Skill:
    """Load a skill by directory name under skills/."""
    skill_dir = _SKILLS_DIR / name
    skill_file = skill_dir / "SKILL.md"
    if not skill_file.exists():
        raise FileNotFoundError(f"Skill '{name}' not found at {skill_file}")
    text = skill_file.read_text(encoding="utf-8")
    skill = _parse_skill_md(text)
    logger.debug("Loaded skill '%s' v%s", skill.name, skill.version)
    return skill


@lru_cache(maxsize=16)
def load_skill_data(name: str) -> dict:
    """Load companion data.yaml for a skill (structured data).

    The returned dict is cached; callers MUST NOT mutate it.
    """
    data_file = _SKILLS_DIR / name / "data.yaml"
    if not data_file.exists():
        raise FileNotFoundError(f"No data.yaml for skill '{name}' at {data_file}")
    with open(data_file, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def extract_prompts(skill: Skill) -> dict[str, str]:
    """Extract prompt templates from a prompt skill's code blocks.

    Expects sections named '## System Prompt' and '## User Prompt',
    each containing exactly one fenced code block.  Code-fence aware:
    headings inside fenced blocks are ignored.
    """
    lines = skill.body.splitlines()
    sections: dict[str, str] = {}
    current_heading: str | None = None
    current_lines: list[str] = []
    in_code_block = False

    for line in lines:
        if line.startswith("```"):
            in_code_block = not in_code_block
        if not in_code_block and line.startswith("## "):
            if current_heading is not None:
                sections[current_heading] = "\n".join(current_lines)
            current_heading = line[3:].strip().lower()
            current_lines = []
        else:
            current_lines.append(line)

    if current_heading is not None:
        sections[current_heading] = "\n".join(current_lines)

    prompts: dict[str, str] = {}
    for heading, content in sections.items():
        blocks = _CODEBLOCK_RE.findall(content)
        if not blocks:
            continue
        if "system" in heading:
            prompts["system"] = blocks[0].strip()
        elif "user" in heading:
            prompts["user"] = blocks[0].strip()
    return prompts


def skill_body(name: str) -> str:
    """Convenience: return just the markdown body of a skill."""
    return load_skill(name).body
