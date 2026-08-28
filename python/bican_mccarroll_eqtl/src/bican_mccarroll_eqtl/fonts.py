"""Manuscript font handling for matplotlib-based plots.

Mirrors the R-side ``fonts.R`` in ``bican.mccarroll.eqtl``: a single place
that pins the font used for manuscript figures, and a loud check for when
that font is not actually resolvable on the current system.
"""

import warnings
from contextlib import contextmanager

import matplotlib as mpl
from matplotlib import font_manager

#: Font used for all manuscript figures.
MANUSCRIPT_FONT = "Arial"


@contextmanager
def manuscript_font():
    """Context manager that renders matplotlib text in the manuscript font.

    Sets ``font.family``/``font.sans-serif`` to :data:`MANUSCRIPT_FONT` and
    ``svg.fonttype`` to ``"none"`` so that SVG output keeps real ``<text>``
    elements (editable, with the font name embedded) rather than converting
    glyphs to vector outlines. Settings are restored on exit via
    ``matplotlib.rc_context``, even if an exception is raised.
    """
    with mpl.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [MANUSCRIPT_FONT],
            "svg.fonttype": "none",
        }
    ):
        yield


def check_manuscript_font():
    """Warn if the manuscript font is not resolvable on this system.

    matplotlib silently substitutes a fallback font if
    :data:`MANUSCRIPT_FONT` is not installed, so this check exists to
    surface that condition loudly instead of letting figures regenerate in
    the wrong font unnoticed.
    """
    try:
        found_path = font_manager.findfont(
            MANUSCRIPT_FONT, fallback_to_default=False
        )
    except Exception:
        found_path = None

    if found_path is None or MANUSCRIPT_FONT.lower() not in found_path.lower():
        fallback = found_path or font_manager.findfont(MANUSCRIPT_FONT)
        warnings.warn(
            f"Manuscript font '{MANUSCRIPT_FONT}' not found on this system; "
            f"figures will fall back to '{fallback}'.",
            stacklevel=2,
        )
