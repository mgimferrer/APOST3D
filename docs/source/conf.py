"""Sphinx configuration for the APOST-3D documentation site."""

project = "APOST-3D"
copyright = "P. Salvador and collaborators, Universitat de Girona"
author = "P. Salvador and collaborators"

extensions = [
    "myst_parser",
]

source_suffix = {
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

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
