"""Feature-grouped, surface-aware CLI command registry (issues #144/#660)."""

from __future__ import annotations

from argparse import _SubParsersAction
from functools import lru_cache
from importlib import import_module
from types import ModuleType
from typing import Any, Callable, Iterable, Iterator
from ..command_surface import COMMAND_SURFACE_BY_NAME, records_for

CORE_MODULE_NAMES = (
    "assembly", "assays", "blast_panel", "budget", "budget_llm", "comparative",
    "diagnostics", "discovery", "evidence", "extraction", "family",
    "homolog_linking", "llm", "mapping", "orchestration", "orthology", "package",
    "ranking", "scenarios", "tnt1",
    "variants", "wet", "production_research", "production_direct_c",
)
# Public release: truth-scoring benchmark command modules are not shipped.
BENCHMARK_MODULE_NAMES = tuple(dict.fromkeys(CORE_MODULE_NAMES + ("research", "retrieval")))
SURFACES = frozenset({"core", "all", "benchmark"})


def _surface_module_names(surface: str) -> tuple[str, ...]:
    if surface not in SURFACES:
        raise ValueError(f"unknown CLI surface: {surface!r}")
    if surface == "core":
        return CORE_MODULE_NAMES
    if surface == "benchmark":
        return BENCHMARK_MODULE_NAMES
    return tuple(dict.fromkeys(CORE_MODULE_NAMES + BENCHMARK_MODULE_NAMES))


@lru_cache(maxsize=None)
def _load_modules(surface: str) -> tuple[ModuleType, ...]:
    return tuple(import_module(f"{__name__}.{name}") for name in _surface_module_names(surface))


def _collect_handlers(modules: Iterable[Any]) -> dict[str, Callable[[Any], dict]]:
    handlers: dict[str, Callable[[Any], dict]] = {}
    for module in modules:
        module_handlers = dict(module.HANDLERS)
        duplicate = handlers.keys() & module_handlers.keys()
        # Production-owned lazy registries intentionally shadow the historical
        # monolithic research/retrieval registrations on the bench superset.
        retained = {
            name for name in duplicate
            if COMMAND_SURFACE_BY_NAME.get(name, None)
            and COMMAND_SURFACE_BY_NAME[name].feature_module.startswith("production_")
        }
        for name in retained:
            module_handlers.pop(name, None)
        duplicate -= retained
        if duplicate:
            raise RuntimeError(f"duplicate CLI command handler(s): {sorted(duplicate)}")
        handlers.update(module_handlers)
    return handlers


@lru_cache(maxsize=None)
def handlers_for_surface(surface: str = "all") -> dict[str, Callable[[Any], dict]]:
    handlers = _collect_handlers(_load_modules(surface))
    allowed = {record.command for record in records_for(surface)}
    unknown = set(handlers) - set(COMMAND_SURFACE_BY_NAME)
    if unknown:
        raise RuntimeError(f"unclassified CLI command(s): {sorted(unknown)}")
    return {name: handler for name, handler in handlers.items() if name in allowed}


def command_names(surface: str = "all") -> frozenset[str]:
    return frozenset(handlers_for_surface(surface))


class _LazyCommandNames:
    """Set-like compatibility view without importing benchmark modules eagerly."""

    def _value(self) -> frozenset[str]:
        return command_names("all")

    def __iter__(self) -> Iterator[str]:
        return iter(self._value())

    def __contains__(self, item: object) -> bool:
        return item in self._value()

    def __len__(self) -> int:
        return len(self._value())

    def difference(self, *others: object) -> set[str]:
        return self._value().difference(*others)

    def __repr__(self) -> str:
        return repr(self._value())

class _NoopParser:
    """Sink parser options when a legacy module shadows a production command."""
    def add_argument(self, *args: Any, **kwargs: Any) -> None: return None
    def set_defaults(self, **kwargs: Any) -> None: return None

class _DuplicateSafeSubparsers:
    def __init__(self, subparsers: _SubParsersAction) -> None: self._subparsers = subparsers
    @property
    def choices(self) -> dict[str, Any]: return self._subparsers.choices
    def add_parser(self, name: str, *args: Any, **kwargs: Any) -> Any:
        if name in self._subparsers.choices: return _NoopParser()
        return self._subparsers.add_parser(name, *args, **kwargs)


COMMAND_NAMES = _LazyCommandNames()


def register_all(subparsers: _SubParsersAction, surface: str = "all") -> None:
    """Register the selected command surface exactly once."""
    modules = _load_modules(surface)
    handlers = handlers_for_surface(surface)
    safe_subparsers = _DuplicateSafeSubparsers(subparsers)
    for module in modules:
        module.register(safe_subparsers)
    # Mixed feature modules (for example ranking and assays) contain both
    # production and benchmark commands.  They share argparse setup but the
    # command-surface registry is authoritative about entrypoint ownership.
    for command in set(subparsers.choices) - set(handlers):
        subparsers.choices.pop(command, None)
    registered = set(subparsers.choices)
    expected = set(handlers)
    if registered != expected:
        missing = sorted(expected - registered)
        extra = sorted(registered - expected)
        raise RuntimeError(f"CLI registry mismatch: missing={missing}, extra={extra}")


def handler_for(command: str, surface: str = "all") -> Callable[[Any], dict] | None:
    return handlers_for_surface(surface).get(command)
