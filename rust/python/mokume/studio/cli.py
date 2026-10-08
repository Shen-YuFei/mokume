"""Command-line startup for Mokume Studio."""

from __future__ import annotations

import argparse
import importlib
import secrets
import socket
import sys
import threading
import time
import webbrowser
from collections.abc import Callable, Sequence
from typing import Any

HOST = "127.0.0.1"
DEFAULT_PORT = 8765
PORT_ATTEMPTS = 50
_MISSING_EXTRA = 'Mokume Studio dependencies are missing; install "mokume[studio]".'


def _port_number(value: str) -> int:
    try:
        port = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("port must be an integer") from exc
    if not 1 <= port <= 65535:
        raise argparse.ArgumentTypeError("port must be between 1 and 65535")
    return port


def build_parser() -> argparse.ArgumentParser:
    """Build the public ``mokume studio`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="mokume studio",
        description="Launch the local Mokume Studio web application.",
    )
    parser.add_argument(
        "--port",
        type=_port_number,
        metavar="PORT",
        help=f"bind strictly to PORT (default: first free port from {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="do not open Studio in the default browser",
    )
    return parser


def _bind_port(port: int) -> socket.socket:
    listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        listener.bind((HOST, port))
        listener.listen()
    except OSError:
        listener.close()
        raise
    return listener


def _select_socket(requested_port: int | None) -> tuple[socket.socket, int]:
    if requested_port is not None:
        try:
            return _bind_port(requested_port), requested_port
        except OSError as exc:
            raise RuntimeError(
                f"cannot bind {HOST}:{requested_port}: {exc.strerror or exc}"
            ) from exc

    for port in range(DEFAULT_PORT, DEFAULT_PORT + PORT_ATTEMPTS):
        try:
            return _bind_port(port), port
        except OSError:
            continue
    last_port = DEFAULT_PORT + PORT_ATTEMPTS - 1
    raise RuntimeError(
        f"no free Studio port found on {HOST} in range {DEFAULT_PORT}-{last_port}"
    )


def _load_runtime() -> tuple[Any, Callable[[], Any]]:
    try:
        uvicorn = importlib.import_module("uvicorn")
        create_app = importlib.import_module("mokume.studio.app").create_app
    except ModuleNotFoundError as exc:
        raise RuntimeError(_MISSING_EXTRA) from exc
    return uvicorn, create_app


def _open_when_ready(server: Any, url: str) -> None:
    while not server.started and not server.should_exit:
        time.sleep(0.01)
    if server.started:
        try:
            webbrowser.open(url)
        except (OSError, webbrowser.Error) as exc:
            print(f"warning: could not open browser: {exc}", file=sys.stderr)


def _serve(
    uvicorn: Any,
    create_app: Callable[[], Any],
    listener: socket.socket,
    port: int,
    *,
    open_browser: bool,
) -> None:
    startup_token = secrets.token_urlsafe(32)
    server_holder: dict[str, Any] = {}

    def request_shutdown() -> None:
        server_holder["server"].should_exit = True

    app = create_app(
        startup_token=startup_token,
        shutdown_callback=request_shutdown,
    )
    config = uvicorn.Config(app, host=HOST, port=port)
    server = uvicorn.Server(config)
    server_holder["server"] = server
    url = f"http://{HOST}:{port}/?token={startup_token}"
    print(f"Mokume Studio: {url}")
    if open_browser:
        threading.Thread(
            target=_open_when_ready,
            args=(server, url),
            daemon=True,
        ).start()
    try:
        server.run(sockets=[listener])
    except KeyboardInterrupt:
        # Python 3.13's asyncio runner can re-raise after Uvicorn shuts down.
        pass


def main(argv: Sequence[str] | None = None) -> int:
    """Parse Studio options and run its single-process local ASGI server."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        uvicorn, create_app = _load_runtime()
        listener, port = _select_socket(args.port)
    except RuntimeError as exc:
        parser.exit(1, f"mokume studio: error: {exc}\n")

    try:
        _serve(
            uvicorn,
            create_app,
            listener,
            port,
            open_browser=not args.no_browser,
        )
    finally:
        listener.close()
    return 0
