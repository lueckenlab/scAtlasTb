# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import sys
import os
from pdf2image import convert_from_path
import shutil
from datetime import datetime
from importlib.metadata import metadata
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))


# -- Project information -----------------------------------------------------

project_name = 'scAtlasTb'
author = 'Michaela Mueller'
copyright = f"{datetime.now():%Y}, {author}."
version = 'dev'
repository_url = f"https://github.com/HCA-integration/scAtlasTb.git"

# The full version, including alpha/beta/rc tags
release = version

bibtex_bibfiles = ["references.bib"]
templates_path = ["_templates"]
nitpicky = True  # Warn about broken links
needs_sphinx = "4.0"

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "mumichae",  # Username
    "github_repo": project_name,  # Repo name
    "github_version": "main",  # Version
    "conf_py_path": "/docs/",  # Path in the checkout to the docs root
}

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_nb",
    "sphinxemoji.sphinxemoji",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx_autodoc_typehints",
    "sphinx.ext.mathjax",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_design",
    "sphinxcontrib.mermaid",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
myst_heading_anchors = 6  # create anchors for h1-h6
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
myst_fence_as_directive = ["mermaid"]  # render ```mermaid fences (as on GitHub)
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True
typehints_defaults = "braces"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = ["js/table-wrap.js"]
html_title = project_name

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "path_to_docs": "docs/",
}

pygments_style = "default"

nitpick_ignore = [
    # If building the documentation fails because of a missing link that is outside your control,
    # you can add an exception to this list.
    #     ("py:class", "igraph.Graph"),
]


def convert_pdfs(app):
    base_dir = os.path.join(app.confdir, '_static', 'images')
    for root, _, files in os.walk(base_dir):
        for f in files:
            if f.lower().endswith('.pdf'):
                pdf_file = os.path.join(root, f)
                png_file = os.path.splitext(pdf_file)[0] + '.png'
                if not os.path.exists(png_file):
                    images = convert_from_path(pdf_file, first_page=1, last_page=1)
                    images[0].save(png_file, 'PNG')


def setup(app):
    """App setup hook."""
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
    )
    app.connect('builder-inited', convert_pdfs)
