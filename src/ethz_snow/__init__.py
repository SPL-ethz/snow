import sys
import os

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover


# Set the bibtex entry to the article referenced in CITATION.
def _get_bibtex():
    citation_file = os.path.join(os.path.dirname(__file__), "CITATION")

    with open(citation_file, "r") as citation:
        refs = citation.read().split("@article")[1:]
        if len(refs) == 0:
            return ""
        bibtexreference = f"@article{refs[0]}"
    return bibtexreference


__citation__ = __bibtex__ = _get_bibtex()

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "ethz_snow"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
