# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os, sys
import json
import urllib.request
    
sys.path.insert(0, os.path.abspath('.'))

project = 'HOPPET'
copyright = '2025, Alexander Karlberg, Paolo Nason, Gavin Salam, Giulia Zanderighi, Frederic Dreyer, Juan Rojo'
author = 'Alexander Karlberg, Paolo Nason, Gavin Salam, Giulia Zanderighi, Frederic Dreyer, Juan Rojo'
release = "2.1.4"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',   # for Google-style or NumPy-style docstrings
    'sphinx.ext.viewcode',   # adds “View source” links
    'sphinx.ext.mathjax',
    ]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']


# --- Get latest HOPPET release PDF URL from GitHub -----------------
def get_latest_hoppet_doc_url():
    """
    Returns the URL of the latest HOPPET-doc.pdf from the GitHub releases API.
    Falls back to a fixed URL if the request fails.
    """
    api_url = "https://api.github.com/repos/hoppet-code/hoppet/releases/latest"
    try:
        with urllib.request.urlopen(api_url, timeout=10) as resp:
            data = json.load(resp)
        tag = data.get("tag_name")
        if not tag:
            return "https://github.com/hoppet-code/hoppet"
        # Construct the release asset URL (adjust if your naming changes)
        return f"https://github.com/hoppet-code/hoppet/releases/download/{tag}/HOPPET-doc.pdf"
    except Exception as e:
        print(f"[WARN] Could not fetch latest HOPPET release: {e}")
        return "https://github.com/hoppet-code/hoppet"

def _suppress_hoppet_data_docstrings(app, what, name, obj, options, lines):
    """
    Sphinx autodoc hook to suppress docstrings for data (autodata) entries
    from the 'hoppet' package — prevents printing the builtin int doc.
    - what: 'module', 'class', 'exception', 'function', 'method', 'attribute', 'data'
    - name: full dotted name of the object
    - obj: the object itself
    - lines: list of doc lines (modifiable in-place)
    """
    # Only act on module-level data items (autodata produces what == 'data')
    if what != "data":
        return

    # Only target items in the hoppet package (change prefix if needed)
    if not name.startswith("hoppet."):
        return

    # Optionally restrict to ints (or remove this check to target any data)
    try:
        if not isinstance(obj, int):
            return
    except Exception:
        # some extension objects might raise on isinstance; be conservative
        return

    # Clear the docstring lines so nothing is printed
    lines[:] = []

# conf.py — robust generator that works across Sphinx versions

def _generate_hoppet_constants(app):
    import importlib, inspect, os
    MODULE = "hoppet"
    OUT_DIR = os.path.join(os.path.dirname(__file__), "_generated")
    OUT_FILE = os.path.join(OUT_DIR, "hoppet_constants.rst")

    # simple logger helper that is safe across Sphinx versions
    def log(level, msg):
        # try sphinx.util.logging first (newer Sphinx)
        try:
            from sphinx.util import logging as _slogging
            logger = _slogging.getLogger(__name__)
            if level == "info":
                logger.info(msg)
            elif level == "warning":
                logger.warning(msg)
            elif level == "error":
                logger.error(msg)
            else:
                logger.info(msg)
            return
        except Exception:
            pass
        # fallback to legacy app.warn/app.info if available
        try:
            if hasattr(app, "warn") and level in ("warning", "error"):
                app.warn(msg)
                return
            if hasattr(app, "info") and level == "info":
                app.info(msg)
                return
        except Exception:
            pass
        # final fallback to print
        print(f"[{level.upper()}] {msg}")

    # Configure prefix groups here
    PREFIX_GROUPS = [
        ("if", "PDF indices"),
        ("factscheme", "Factorisation schemes"),
        ("nnlo_splitting", "NNLO splitting function choices"),
        ("nnlo_nfthreshold", "NNLO mass threshold choices"),
        ("n3lo_splitting_approx", "N3LO splitting function approximations"),
        ("n3lo_splitting", "N3LO splitting function choices"),
        ("n3lo_nfthreshold", "N3LO mass threshold choices"),
        ("cc_", "Convolution communicator pieces"),
        ("iF", "Structure function indices"),
        ("scale_choice", "Structure function scale choices"),
        # add other groups as needed
    ]
    OTHER_GROUP = "Other constants"

    # Optional extra text to include under specific section titles
    GROUP_DESCRIPTIONS = {
        "PDF indices": """ These are the PDF indices of lists that are returned by e.g. :func:`Eval` or that are to be passed to :func:`EvalIFlv`.
        Note that index 13, not listed below, is used to encode information about the PDF flavour representation. 
        """,        
        "Structure function indices": """These constants label the various structure functions  used in HOPPET

            * F1 W+ : :math:`d+\\bar{u}`                                                       
            * F2 W+ : :math:`d + \\bar{u}`                                                      
            * F3 W+ : :math:`d + \\bar{u}`                                                      
            * F1 W- : :math:`\\bar{d} + u`                                                      
            * F2 W- : :math:`\\bar{d} + u`                                                      
            * F3 W- : :math:`\\bar{d} + u`                                                      
            * F1 Z  : :math:`(d + \\bar{d})  v_d^2 a_d^2 + (u + \\bar{u})  v_u^2 a_u^2`     
            * F2 Z  : :math:`(d + \\bar{d})  v_d^2 a_d^2 + (u + \\bar{u})  v_u^2a_u^2`     
            * F3 Z  : :math:`(d + \\bar{d})  2 v_d a_d + (u + \\bar{u})  2 v_u a_u`           
            * F1 γ  : :math:`(d + \\bar{d})  e^2_d + (u + \\bar{u})  e^2_u`                     
            * F2 γ  : :math:`(d + \\bar{d})  e^2_d + (u + \\bar{u})  e^2_u`                     
            * F1 γZ : :math:`(d + \\bar{d})  e_d  2v_d + (u + \\bar{u})  e_u  2v_u` 
            * F2 γZ : :math:`(d + \\bar{d})  e_d  2v_d + (u + \\bar{u})  e_u  2v_u`
            * F3 γZ : :math:`(d + \\bar{d})  e_d  2a_d + (u + \\bar{u})  e_u  2a_u` 

        where the :math:`u` and :math:`d` refer to the sum over up- and down-type quarks, and :math:`v_i`, :math:`a_i`, and :math:`e_i` are the vector, axial, and electric charges of an :math:`i`-type quark. 
""",
        "Geometry indices":
            "Geometry-related integer constants used for grid layout and transforms.",
        # Add more descriptions as you wish
    }

    def is_constant_name(name):
        return not name.startswith("_")

    def is_constant_obj(obj):
        return isinstance(obj, (int, float, str, bool, tuple))

    try:
        mod = importlib.import_module(MODULE)
    except Exception as e:
        log("warning", f"Could not import {MODULE} to generate constants: {e}")
        return

    names = sorted(name for name in dir(mod) if is_constant_name(name))
    constants = []
    for n in names:
        try:
            obj = getattr(mod, n)
        except Exception:
            continue
        if inspect.ismodule(obj) or inspect.isroutine(obj) or inspect.isclass(obj):
            continue
        if is_constant_obj(obj):
            constants.append((n, obj))

    groups = {title: [] for _, title in PREFIX_GROUPS}
    groups[OTHER_GROUP] = []
    for name, obj in constants:
        matched = False
        for prefix, title in PREFIX_GROUPS:
            if name.startswith(prefix):
                groups[title].append((name, obj))
                matched = True
                break
        if not matched:
            groups[OTHER_GROUP].append((name, obj))

    # prune empty groups
    groups = {k: v for k, v in groups.items() if v}

    os.makedirs(OUT_DIR, exist_ok=True)
    try:
        with open(OUT_FILE, "w", encoding="utf8") as f:
            TITLE = "HOPPET constants"
            f.write(TITLE + "\n")
            f.write("=" * len(TITLE) + "\n\n")
            f.write("Auto-generated constants reference for the `hoppet` module.\n\n")

            for group_title, items in groups.items():
                f.write(group_title + "\n")
                f.write("-" * len(group_title) + "\n\n")

                # Add the custom text if available
                desc = GROUP_DESCRIPTIONS.get(group_title)
                if desc:
                    f.write(desc.strip() + "\n\n")
                
                for name, obj in sorted(items, key=lambda x: x[0]):
                    f.write(f".. autodata:: {MODULE}.{name}\n")
                    # show value and type by default - modify as you prefer:
                    #f.write("   :annotation: " + (type(obj).__name__) + "\n")
                    f.write("\n")
    except Exception as e:
        log("error", f"Failed to write {OUT_FILE}: {e}")
        return

    total = sum(len(v) for v in groups.values())
    log("info", f"Generated {OUT_FILE} with {total} constants.")

def _generate_latest_doc_link(app):
    import os
    OUT_DIR = os.path.join(os.path.dirname(__file__), "_generated")
    os.makedirs(OUT_DIR, exist_ok=True)
    outfile = os.path.join(OUT_DIR, "github_links.rst")

    url = get_latest_hoppet_doc_url()
    with open(outfile, "w", encoding="utf8") as f:
        f.write("Manual in PDF format\n")
        f.write("=====================================\n\n")
        f.write(f"`Download latest PDF manual <{url}>`_\n")

    print(f"[INFO] Generated {outfile} with latest doc URL: {url}")
    
def setup(app):
    app.connect("autodoc-process-docstring", _suppress_hoppet_data_docstrings)
    app.connect("builder-inited", _generate_hoppet_constants)
    app.connect("builder-inited", _generate_latest_doc_link)



