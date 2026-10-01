"""Sphinx configuration for the APOST-3D documentation site."""

project = "APOST-3D"
copyright = "P. Salvador (University of Girona), M. Gimferrer (University of Göttingen) and collaborators"
author = "P. Salvador, M. Gimferrer and collaborators"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

source_suffix = {
    ".md": "markdown",
}

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
]

# Anchors for headings up to level 3, so pages can link to sections
# ("page.md#section-title").
myst_heading_anchors = 3

# Copy buttons skip shell prompts and output lines.
copybutton_exclude = ".linenos, .gp, .go"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- HTML output -------------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/logo-apost.png"
html_title = "APOST-3D documentation"

html_theme_options = {
    "sidebar_hide_name": True,
    "source_repository": "https://github.com/mgimferrer/APOST3D/",
    "source_branch": "MAJOR-UPDATE",
    "source_directory": "docs/source/",
}
