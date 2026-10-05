"""Minimal stand-in for the stdlib cgi module, removed in Python 3.13
(PEP 594). ete3's webplugin submodule (imported unconditionally by
ete3/__init__.py, even though nothing in this repo uses ete3's web
features) only calls cgi.parse_qs, which was itself just a thin wrapper
around urllib.parse.parse_qs in the versions of cgi that still had it.
"""
from urllib.parse import parse_qs
